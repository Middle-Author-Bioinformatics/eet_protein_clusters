"""Consequential genome, motif, adapter and end-to-end workflow checks."""
import csv
import json
from pathlib import Path
import random
import shutil
import sys

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, BeforePosition
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import genome_inputs as gi
import genome_pipeline as gp
import localization as loc
from protein_motifs import characterize
from test_workflow import make_hmm, read


def test_source_motif_patterns_positions_and_nonoverlap():
    seq = ('G'*40).join(['CAACH', 'CAAACH', 'CAAAACH', 'C'+'A'*14+'CH', 'C'+'A'*15+'CH', 'GAHRCPG', 'SASPT', 'SASPT'])
    counts, hits = characterize(seq.lower())
    assert {k: counts[k] for k in ['cxxch_count', 'cxxxch_count', 'cxxxxch_count', 'cx14ch_count', 'cx15ch_count', 'heme_b_count', 'hematite_count']} == dict(cxxch_count=1, cxxxch_count=1, cxxxxch_count=1, cx14ch_count=1, cx15ch_count=1, heme_b_count=1, hematite_count=2)
    assert counts['mhcscan_hemes']==2 and counts['fegenie_heme_c_total']==5
    assert all(seq[h['start']-1:h['end']]==h['match'] for h in hits)
    assert characterize('M'+'P'*100)[0]['cxxch_count']==0
    # Shared cysteine overlaps are not counted twice by source re.findall semantics.
    assert characterize('CAACHCAACH')[0]['cxxch_count']==2
    assert characterize('CAACHACH')[0]['cxxch_count']==1


def gbk(path, proteins):
    from Bio.Data import CodonTable
    codons = CodonTable.unambiguous_dna_by_id[11].back_table
    chunks = []; features = []; offset = 0
    for i, (pid, seq) in enumerate(proteins):
        coding = ''.join(codons[aa] for aa in seq)+'TAA'
        genomic = coding if i%2==0 else str(Seq(coding).reverse_complement())
        features.append(SeqFeature(SimpleLocation(offset,offset+len(coding),strand=(-1 if i%2 else 1)),type='CDS',
                                   qualifiers={'protein_id':[pid], 'translation':[seq], 'transl_table':['11']}))
        chunks.append(genomic+'A'*30); offset += len(coding)+30
    rec = SeqRecord(Seq(''.join(chunks)), id='ctg', description='synthetic'); rec.annotations['molecule_type']='DNA'
    rec.features = features
    SeqIO.write([rec], path, 'genbank')


def test_genbank_preserves_annotation_ids_reverse_and_partial(tmp_path):
    inp=tmp_path/'genomes'; inp.mkdir(); gbk(inp/'one.gbk', [('same','MCAACH'),('same','MSASPT')]); gbk(inp/'two.gbk', [('same','MCAACH')])
    files=gi.discover(inp); grams=gi.gram_assignments(files, 'positive')
    result=gi.prepare(files, grams, tmp_path/'prepared')
    rows=read(tmp_path/'prepared/protein_mapping.tsv')
    assert len(rows)==3 and len({r['protein_id'] for r in rows})==3
    assert rows[1]['start']=='52' and rows[1]['end']=='72' and rows[1]['strand']=='-'
    assert all(r['annotation_source']=='genbank_translation' for r in rows)
    assert sum(g['proteins'] for g in result)==3
    with pytest.raises(ValueError,match='ambiguous'):loc.map_predictions({'same':{'prediction':'OTHER'}},rows[:2])


def test_genbank_translation_join_and_pseudogenes(tmp_path):
    rec=SeqRecord(Seq('GTGAAATAG'+'A'*10+'ATGGCTTAA'),id='ctg');rec.annotations['molecule_type']='DNA'
    rec.features=[SeqFeature(SimpleLocation(0,9),type='CDS',qualifiers={'locus_tag':['first']}),
                  SeqFeature(CompoundLocation([SimpleLocation(BeforePosition(0),3),SimpleLocation(19,25)]),type='CDS',qualifiers={'locus_tag':['joined'],'translation':['MMA']}),
                  SeqFeature(SimpleLocation(0,9),type='CDS',qualifiers={'pseudo':['']})]
    path=tmp_path/'g.gbk';SeqIO.write(rec,path,'genbank')
    rows=list(gi.records(path,'genbank'))
    assert rows[0][1]=='MK*' and rows[0][0]['annotation_source']=='genbank_CDS_translated'
    assert rows[1][0]['partial']=='true' and 'join' in rows[1][0]['location'] and len(rows)==2


def test_metadata_and_ambiguous_inputs(tmp_path):
    (tmp_path/'x.fna').write_text('>ctg\nATGAAA\n')
    files=gi.discover(tmp_path)
    with pytest.raises(ValueError,match='never guessed'):gi.gram_assignments(files)
    meta=tmp_path/'meta.tsv';meta.write_text('genome\tgram\nx\tnegative\n')
    assert gi.gram_assignments(files,metadata=meta)=={'x':'negative'}
    (tmp_path/'x.faa').write_text('>p\nMKK\n')
    with pytest.raises(ValueError,match='ambiguous'):gi.discover(tmp_path)
    assert gi.discover(tmp_path,'assembly')[0][2]=='assembly'
    (tmp_path/'x.gbk').write_text('')
    with pytest.raises(ValueError,match='Duplicate'):gi.discover(tmp_path,'auto')


def test_real_pyrodigal_prediction_reproducible(tmp_path):
    pytest.importorskip('pyrodigal')
    rng=random.Random(42)
    codons=['GCT','GAA','AAA','GGT','GTT','GAC','TCT','GCC','CTG','ACG']
    dna=''.join('TAA'+'A'*15+'AGGAGG'+'A'*7+'ATG'+''.join(rng.choices(codons,k=180))+'TAA'+'C'*30 for _ in range(12))
    fasta=tmp_path/'b.fna';fasta.write_text('>contig\n'+dna+'\n')
    a=list(gi.records(fasta,'assembly'));b=list(gi.records(fasta,'assembly'))
    assert a==b and len(a)>5
    assert all(m['translation_table']==11 and m['annotation_source']=='pyrodigal_meta' and int(m['start'])<=int(m['end']) for m,s in a)
    assert any(len(s)>170 for m,s in a)


def test_native_parsers_and_id_validation(tmp_path):
    ps=tmp_path/'p';ps.write_text('SeqID\tFinal_Localization\tFinal_Score\tSecondary_Localization\na desc\tUnknown\t0.0\t\nb\tCytoplasmic\t9.8\t\n')
    assert loc.parse_native('psortb',ps)['a']['localization']=='Unknown'
    sp=tmp_path/'s';sp.write_text('# SignalP-6.0\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\na\tLIPO\t0.01\t0.01\tCS pos: 20-21\n')
    assert loc.parse_native('signalp',sp)['a']['prediction']=='LIPO'
    ph=tmp_path/'t';ph.write_text('SEQENCE ID TM SP PREDICTION\na 0 Y n4-16c20/21o\nb 2 0 i5-25o60-80i\n')
    assert loc.parse_native('phobius',ph)['b']['tm_count']==2
    lp=tmp_path/'l';lp.write_text('# a SpII score=9.2 margin=4.0 cleavage=20-21 Pos+2=G\n# Cut-off=-3\n')
    assert loc.parse_native('lipop',lp)['a']['prediction']=='SpII'
    mapped=loc.map_predictions({'a':{'tm_count':2}},[dict(protein_id='safe',original_id='a')]); assert mapped=={'safe':{'tm_count':2}}
    with pytest.raises(ValueError,match='Unknown'):loc.map_predictions({'missing':{}},[dict(protein_id='safe',original_id='a')])
    ps.write_text(ps.read_text()+'a\tUnknown\t0\t\n')
    with pytest.raises(ValueError,match='duplicate'):loc.parse_native('psortb',ps)
    sp.write_text('ID\tPrediction\na\tSP5\n')
    with pytest.raises(ValueError,match='unsupported'):loc.parse_native('signalp',sp)


def test_localization_unknown_refinement_separate_from_exposure():
    unknown={'psortb':dict(localization='Unknown',score='0')}
    r=loc.infer(dict(unknown,signalp=dict(prediction='LIPO')), 'positive')
    assert r['psortb_localization']=='Unknown' and r['localization_interpretation']=='putative_membrane_attached_lipoprotein'
    assert r['extracellular_exposure_evidence']=='unresolved' and r['phobius_tm_count']=='unavailable'
    r=loc.infer(dict(unknown,signalp=dict(prediction='SP'),phobius=dict(tm_count=0,signal_peptide=True,topology='n5-15c20/21o')), 'negative')
    assert r['localization_interpretation']=='export_signal_detected_destination_unresolved'
    assert r['cell_wall_evidence']=='unresolved'
    r=loc.infer(dict(psortb=dict(localization='Cytoplasmic'),signalp=dict(prediction='LIPO')), 'positive')
    assert r['localization_interpretation']=='Cytoplasmic' and r['localization_conflict']
    r=loc.infer(unknown, 'positive'); assert r['localization_interpretation']=='Unknown'
    r=loc.infer({}, 'negative'); assert r['psortb_status']=='unavailable' and r['signalp_status']=='unavailable'
    assert loc.infer(dict(unknown,phobius=dict(tm_count=2,signal_peptide=False,topology='i5-20o40-60i')), 'positive')['localization_interpretation']=='putative_membrane_protein'


def predictor_stub(tmp_path):
    path=tmp_path/'fixture_predictor.py'
    path.write_text('import sys\nfrom Bio import SeqIO\nprint("SeqID\\tFinal_Localization\\tFinal_Score")\nfor r in SeqIO.parse(sys.argv[1],"fasta"): print(r.id+"\\tUnknown\\t0.0")\n')
    return dict(psortb=dict(argv=[sys.executable,str(path),'{input}'], output='stdout', version='TEST DOUBLE; not PSORTb'))


def test_predictor_execution_cache_and_partial_import(tmp_path):
    from argparse import Namespace
    fasta=tmp_path/'p.faa';fasta.write_text('>safeA\nMKK\n>safeB\nMKK\n')
    genome=dict(genome='g',fasta=str(fasta),gram='positive')
    mapping=[dict(protein_id='safeA',original_id='a'),dict(protein_id='safeB',original_id='b')]
    args=Namespace(psortb='run',signalp='skip',phobius='skip',lipop='auto', localization_import=None,threads=1,signalp_mode='fast',resume=False,force=False,predictor_cache_tag='')
    manifest=dict(localization=[],skipped_analyses=[],warnings=[],commands=[])
    config=predictor_stub(tmp_path)
    result=loc.collect(genome,mapping,args,config,tmp_path/'out',manifest)
    assert len(result['psortb'])==2 and not result['phobius'] and len(manifest['commands'])==1
    args.resume=True
    loc.collect(genome,mapping,args,config,tmp_path/'out',manifest)
    assert manifest['localization'][-4]['cache_reused'] and len(manifest['commands'])==1
    (tmp_path/'out/psortb.txt').write_text('corrupt')
    loc.collect(genome,mapping,args,config,tmp_path/'out',manifest)
    assert len(manifest['commands'])==2
    imports=tmp_path/'imports';imports.mkdir();(imports/'g.psortb.tsv').write_text('SeqID\tLocalization\tScore\na\tUnknown\t0\n')
    args.localization_import=imports
    result=loc.collect(genome,mapping,args,config,tmp_path/'out',manifest)
    assert len(result['psortb'])==1 and 'safeB' not in result['psortb'] and manifest['localization'][-4]['status']=='imported'


def test_missing_requested_predictor_has_actionable_error(tmp_path,monkeypatch):
    from argparse import Namespace
    monkeypatch.setattr(shutil,'which',lambda _:None)
    args=Namespace(psortb='run',signalp='skip',phobius='skip',lipop='skip', localization_import=None,threads=1,signalp_mode='fast',resume=False,force=False,predictor_cache_tag='')
    with pytest.raises(ValueError,match='--psortb skip'):
        loc.collect(dict(genome='g',fasta='unused',gram='positive'),[],args,{},tmp_path,dict(localization=[],skipped_analyses=[],commands=[],warnings=[]))


@pytest.mark.skipif(not all(shutil.which(t) for t in ['diamond','mcl','mcxload']),reason='external clustering binaries not installed')
def test_real_genome_to_domains_motifs_localization_wrapper(tmp_path):
    hmm=tmp_path/'models.hmm'; seqs=make_hmm(hmm); a=seqs['A'];b=seqs['B']
    proteins=[('original',a+'GGGGGG'+a+'GGGGGG'+b+'GGCAACHGGCAACHGGSASPT'),('other',a+'GGGGGG'+a+'GGGGGG'+b+'GGCAACHGGCAACHGGSASPT'),('negative','M'+'P'*200)]
    inp=tmp_path/'genomes';inp.mkdir();gbk(inp/'one.gbk',proteins);gbk(inp/'two.gbk',[('original',b)])
    meta=tmp_path/'gram.tsv';meta.write_text('genome\tgram\none\tpositive\ntwo\tnegative\n')
    imports=tmp_path/'native';imports.mkdir()
    # Partial native imports deliberately leave the last protein unavailable.
    (imports/'one.signalp.tsv').write_text('# ID\tPrediction\tCS Position\noriginal\tLIPO\tCS pos: 20-21\nother\tOTHER\t\n')
    (imports/'one.phobius.txt').write_text('SEQENCE ID TM SP PREDICTION\noriginal 0 Y n4-15c20/21o\nother 0 0 o\nnegative 0 0 o\n')
    config=tmp_path/'config.json';config.write_text(json.dumps(predictor_stub(tmp_path)))
    out=tmp_path/'analysis'
    argv=['-i',str(inp),'-o',str(out),'--hmm',str(hmm),'--genome-metadata',str(meta),'-t','1','--no-plots','--localization-config',str(config),'--localization-import',str(imports),'--signalp','auto','--phobius','auto','--psortb','run']
    assert gp.main(argv)==0
    rows=read(out/'protein_annotations.tsv'); assert len(rows)==4
    first=rows[0]; assert first['original_id']=='original' and first['cxxch_count']=='2' and first['hematite_count']=='1'
    assert first['psortb_localization']=='Unknown' and first['localization_interpretation']=='putative_membrane_attached_lipoprotein'
    assert first['extracellular_exposure_evidence']=='unresolved'
    assert int(first['n_domain_hits'])==3 and first['architecture'].count('TEST_A.1')==2
    assert rows[2]['n_domain_hits']=='0' and rows[2]['signalp_status']=='unavailable'
    assert rows[0]['cluster']==rows[1]['cluster']
    summaries=read(out/'genome_summary.tsv'); assert summaries[0]['proteins']=='3' and summaries[0]['psortb_unknown']=='3'
    assert summaries[1]['gram']=='negative'
    dom=read(out/'domains/domain_prevalence.tsv');assert dom
    manifest=json.loads((out/'manifest.json').read_text());assert manifest['counts']['proteins']==4 and manifest['motif_sources']['FeGenie']['commit']
    assert any(s['tool']=='lipop' for s in manifest['skipped_analyses'])
    assert manifest['localization'][0]['argv'][1].endswith('fixture_predictor.py')
    raw=out/'domains/work/domains.domtblout';before=raw.stat().st_mtime_ns
    assert gp.main(argv+['--resume'])==0
    manifest=json.loads((out/'manifest.json').read_text()); assert manifest['stages']['gene_extraction']['cache_reused']
    assert all(e['cache_reused'] for e in manifest['localization'] if e['status']=='ran')
    assert raw.stat().st_mtime_ns==before and manifest['stages']['cluster_domains.py']['provenance']['search_definition']['cache_reused']
    # Genome changes invalidate preparation, alignment, HMM and localization caches.
    gbk(inp/'two.gbk',[('original',b),('new','MCAACH')])
    assert gp.main(argv+['--resume'])==0
    manifest=json.loads((out/'manifest.json').read_text());assert not manifest['stages']['gene_extraction']['cache_reused']
    assert len(read(out/'protein_annotations.tsv'))==5
    assert not manifest['stages']['cluster_domains.py']['provenance']['search_definition']['cache_reused']
    assert any(e['status']=='ran' and not e['cache_reused'] for e in manifest['localization'])
    assert 'Per-cluster domains and networks' in (out/'index.html').read_text()


def test_protein_folder_no_coordinates_and_same_ids_across_genomes(tmp_path):
    inp=tmp_path/'in';inp.mkdir();(inp/'a.faa').write_text('>same description\nMCAACH\n');(inp/'b.faa').write_text('>same description\nMCAACH\n')
    files=gi.discover(inp, 'proteins');gi.prepare(files, gi.gram_assignments(files, 'negative'), tmp_path/'out')
    rows=read(tmp_path/'out/protein_mapping.tsv')
    assert rows[0]['original_id']==rows[1]['original_id']=='same'
    assert rows[0]['protein_id']!=rows[1]['protein_id']
    assert all(r['start']==r['end']==r['contig']=='' and r['partial']=='unknown' for r in rows)


def test_default_localization_gram_selection_and_signalp_threads(tmp_path):
    positive, format_p=loc.default_command('psortb', tmp_path/'p.faa', tmp_path/'out', 'positive', 4, 'fast')
    negative, _=loc.default_command('psortb', tmp_path/'p.faa', tmp_path/'out', 'negative', 4, 'fast')
    assert positive[1]=='-p' and negative[1]=='-n' and positive[2:4]==['-o','long'] and format_p=='stdout'
    argv, name=loc.default_command('signalp', tmp_path/'p.faa', tmp_path/'out', 'positive', 4, 'slow-sequential')
    assert argv[argv.index('--organism')+1]=='other' and argv[argv.index('--torch_num_threads')+1]=='4'
    assert argv[argv.index('--mode')+1]=='slow-sequential' and name=='prediction_results.txt'


def test_config_rejects_invalid_argv_and_hashes_relative_models(tmp_path):
    config=tmp_path/'c.json';config.write_text('{"psortb":{"argv":[]}}')
    with pytest.raises(ValueError,match='nonempty'):gp.load_config(config)
    model=tmp_path/'model';model.write_text('x')
    config.write_text('{"psortb":{"version":"3.0.6","model_files":["model"]}}')
    assert gp.load_config(config)['psortb']['model_files']==[str(model.resolve())]
