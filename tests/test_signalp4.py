"""SignalP 4 / MHCscan Phobius adapters; fabricated predictor fixtures are labeled."""
import json
from pathlib import Path
import shutil
import sys
from argparse import Namespace
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import localization as loc
import genome_pipeline as gp
from predictor_setup import configure
from test_workflow import make_hmm, read


SIGNALP4 = '''# SignalP-4.1 gram+ predictions
# name Cmax pos Ymax pos Smax pos Smean D ? Dmaxcut Networks-used
positive 0.600 21 0.700 21 0.900 10 0.800 0.800 Y 0.570 SignalP-TM
negative 0.100 20 0.100 20 0.200 8 0.100 0.100 N 0.570 SignalP-TM
'''


def test_signalp4_native_scores_cleavage_and_limitations(tmp_path):
    path=tmp_path/'signalp.txt';path.write_text(SIGNALP4)
    rows=loc.parse_native('signalp',path)
    assert rows['positive']['generation']=='4' and rows['positive']['prediction']=='SP'
    assert rows['positive']['cleavage']=='20-21'  # Ymax=21 is first mature residue.
    assert rows['positive']['d_score']=='0.800' and rows['positive']['d_cutoff']=='0.570'
    assert rows['positive']['networks']=='SignalP-TM'
    assert rows['negative']['prediction']=='OTHER' and rows['negative']['cleavage']==''
    r=loc.infer(dict(signalp=rows['positive']),'positive')
    assert r['lipoprotein_assessment']=='unavailable' and r['membrane_attachment_evidence']=='unavailable'
    assert r['extracellular_exposure_evidence']=='unresolved'
    assert r['localization_interpretation']=='export_signal_detected_destination_unresolved'
    assert r['signalp_d_score']=='0.800' and r['signalp_generation']=='4'
    negative=loc.infer(dict(signalp=rows['negative'],phobius=dict(tm_count=0,signal_peptide=False,topology='i')),'positive')
    assert negative['localization_interpretation']=='unresolved_no_positive_targeting_evidence'
    assert negative['phobius_status']=='assessed' and negative['lipoprotein_assessment']=='unavailable'


@pytest.mark.parametrize('bad', [
    SIGNALP4.replace('0.800 Y', 'nan Y'),
    SIGNALP4.replace(' Y ', ' MAYBE '),
    SIGNALP4+'positive 0.600 21 0.700 21 0.900 10 0.800 0.800 Y 0.570 SignalP-TM\n',
])
def test_bad_signalp4_output_fails(tmp_path,bad):
    path=tmp_path/'bad';path.write_text(bad)
    with pytest.raises(ValueError):loc.parse_native('signalp',path)


def test_signalp4_gram_arguments_and_version(tmp_path):
    for gram,flag in [('positive','gram+'),('negative','gram-')]:
        argv,fmt=loc.default_command('signalp',tmp_path/'input.faa',tmp_path/'out',gram,4,'fast','signalp4')
        assert argv==['signalp','-t',flag,'-f','short','-s','best','-c','70',str(tmp_path/'input.faa')]
        assert fmt=='stdout'


def fake_conda(tmp_path):
    """An executable test double that forwards conda run to a synthetic predictor."""
    prefix=tmp_path/'environments/signalP';(prefix/'bin').mkdir(parents=True)
    predictor=prefix/'bin/signalp'
    predictor.write_text('#!'+sys.executable+'\n'+'''import sys
from Bio import SeqIO
if '-V' in sys.argv:
    print('SignalP-4.1 (TEST DOUBLE)');sys.exit(0)
gram=sys.argv[sys.argv.index('-t')+1]
assert gram in ['gram+','gram-']
assert '-f' in sys.argv and sys.argv[sys.argv.index('-f')+1]=='short'
print('# SignalP-4.1 '+gram+' predictions')
print('# name Cmax pos Ymax pos Smax pos Smean D ? Dmaxcut Networks-used')
for r in SeqIO.parse(sys.argv[-1],'fasta'):
    print(r.id+' 0.600 21 0.700 21 0.900 10 0.800 0.800 Y 0.570 SignalP-TM')
''');predictor.chmod(0o755)
    conda=tmp_path/'fake-conda'
    conda.write_text('#!'+sys.executable+'\n'+f'''import sys,json,subprocess
from pathlib import Path
if sys.argv[1:]==['env','list','--json']:
    print(json.dumps({{'envs':[{str(prefix)!r}]}}));sys.exit(0)
assert sys.argv[1:5]==['run','--no-capture-output','-p',{str(prefix)!r}]
sys.exit(subprocess.run([str(Path(sys.argv[4])/'bin'/sys.argv[5])]+sys.argv[6:]).returncode)
''');conda.chmod(0o755)
    return conda,prefix


def test_configure_existing_conda_environment(tmp_path):
    conda,prefix=fake_conda(tmp_path)
    args=Namespace(signalp='run',signalp_version='4',signalp4_env='signalP',conda_executable=str(conda),phobius='skip',phobius_dir=None)
    config=configure(args,{})['signalp']
    assert config['argv'][:6]==[str(conda),'run','--no-capture-output','-p',str(prefix),'signalp']
    assert config['version_argv'][-1]=='-V' and config['format']=='signalp4'
    assert str(prefix/'bin/signalp') in config['model_files']
    args.signalp_version='6'
    with pytest.raises(ValueError,match='requires'):configure(args,{})
    args.signalp='skip'
    assert configure(args,{})=={}  # skipped tools never access their installation


def test_separate_phobius_stderr_and_use_stdin_with_space_paths(tmp_path):
    source=tmp_path/'input proteins.faa';source.write_text('>safe\nMCAACH\n')
    script=tmp_path/'fake_phobius.py'
    script.write_text('import sys\nassert ">safe" in sys.stdin.read()\nprint("Phobius ver 1.01\\nwarning from decoder",file=sys.stderr)\nprint("SEQENCE ID TM SP PREDICTION\\nsafe 0 Y n4-15c20/21o")\n')
    out=tmp_path/'predictions';log=tmp_path/'stderr.log';manifest=dict(commands=[])
    loc.run_prediction([sys.executable,str(script)],source,out,log,manifest)
    assert loc.parse_native('phobius',out)['safe']['signal_peptide']
    assert 'Phobius ver' in log.read_text() and 'Phobius ver' not in out.read_text()


def test_phobius_bundle_validation_and_model_tracking(tmp_path):
    folder=tmp_path/'phobius';folder.mkdir()
    for name in ['phobius.pl','decodeanhmm.64bit','phobius.model','phobius.options']:
        (folder/name).write_text('# phobius.pl 1.01\n')
    args=Namespace(signalp='skip',phobius='run',phobius_dir=folder)
    with pytest.raises(ValueError,match='chmod'):configure(args,{})
    (folder/'decodeanhmm.64bit').chmod(0o755)
    config=configure(args,{})['phobius']
    assert config['argv'][-1]=='-short' and config['stdin'] is True
    assert config['version']=='Phobius 1.01' and len(config['model_files'])==4


@pytest.mark.skipif(not all(shutil.which(t) for t in ['diamond','mcl','mcxload']),reason='external clustering binaries not installed')
def test_wrapper_signalp4_env_and_no_psortb_or_lipop(tmp_path):
    conda,prefix=fake_conda(tmp_path)
    hmm=tmp_path/'models.hmm';seqs=make_hmm(hmm)
    inp=tmp_path/'proteomes';inp.mkdir();(inp/'one.faa').write_text('>original\n'+seqs['A']+'CAACHCAACH\n')
    imports=tmp_path/'imports';imports.mkdir();(imports/'one.phobius.txt').write_text('SEQENCE ID TM SP PREDICTION\noriginal 0 Y n4-15c20/21o\n')
    out=tmp_path/'analysis with spaces'
    argv=['-i',str(inp),'--input-type','proteins','-o',str(out),'--hmm',str(hmm),'--gram','positive','-t','1','--no-plots',
          '--signalp','run','--signalp4-env','signalP','--conda-executable',str(conda),'--localization-import',str(imports)]
    assert gp.main(argv)==0
    rows=read(out/'protein_annotations.tsv');r=rows[0]
    assert r['signalp_generation']=='4' and r['signalp_d_score']=='0.800' and r['signalp_cleavage']=='20-21'
    assert r['phobius_status']=='assessed' and r['lipoprotein_assessment']=='unavailable'
    assert r['psortb_status']==r['lipop_status']=='unavailable'
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['parameters']['psortb']==manifest['parameters']['lipop']=='skip'
    native=next(e for e in manifest['localization'] if e['tool']=='signalp')
    assert native['version']=='SignalP-4.1 (TEST DOUBLE)' and native['detected_generations']==['4']
    assert '-V' not in native['argv'] and 'gram+' in native['argv']
    assert not any(e['status']=='ran' for e in manifest['localization'] if e['tool'] in ['psortb','lipop'])
    assert gp.main(argv+['--resume'])==0
    manifest=json.loads((out/'manifest.json').read_text());assert next(e for e in manifest['localization'] if e['tool']=='signalp')['cache_reused']
    assert 'signalp_d_score' in (out/'index.html').read_text()
