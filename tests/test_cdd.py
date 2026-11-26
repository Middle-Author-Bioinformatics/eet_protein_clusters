"""Controlled CDD fixtures, plus optional execution against an installed snapshot."""
import csv
import io
import json
import os
from pathlib import Path
import tarfile
import pytest
import cluster_cdd as cdd
import download_cdd as setup
from common import dump,sha256
from domain_members import DomainMembers
import network_reports


def rows(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def fixture(tmp_path,single=False):
    database=tmp_path/'CDD';data=database/'snapshot/data';data.mkdir(parents=True)
    (data/'cddid.tbl').write_text('100\tcd00001\tDomain_A\tTest domain A\t40\n101\tcd00002\tDomain_B\tTest domain B\t40\n200\tcl00001\tSuper_A\tSuperfamily A\t0\n300\tPRK00001\tWhole_AB\tWhole protein\t200\n')
    (data/'bitscore_specific.txt').write_text('100\tcd00001\t25\n101\tcd00002\t30\n')
    (data/'cdtrack.txt').write_text('')
    db=database/'snapshot/db';db.mkdir();(db/'Cdd.rps').write_text('fixture database, not executable')
    dump(database/'cdd_manifest.json',dict(schema=1,subset='full',release_info='synthetic',db_prefix='snapshot/db/Cdd',data_dir='snapshot/data',files={str(p.relative_to(database)):sha256(p) for p in database.rglob('*') if p.is_file()}))
    proteins=tmp_path/'clusters';proteins.mkdir()
    # A-A-B, A only, B only, negative. Separate c2 carriers must not form a pair.
    (proteins/'c1.faa').write_text('>pair annotated protein\n'+'A'*200+'\n')
    (proteins/'c2.faa').write_text('>only_a\n'+'A'*200+'\n>only_b\n'+'A'*200+'\n>negative\n'+'A'*200+'\n')
    if single:
        path=tmp_path/'all.faa';path.write_text(''.join(p.read_text() for p in sorted(proteins.glob('*.faa'))))
    else:path=proteins
    align=[];blocks=[]
    specs=[[(100,1,40),(100,60,99),(101,130,169)],[(100,1,40)],[(101,1,40)],[]]
    names={100:('cd00001','Domain_A'),101:('cd00002','Domain_B')}
    for n,hits in enumerate(specs,1):
        q=f'p{n:010d}';ql=f'Query_{n}'
        blocks+=['QUERY\t'+ql+'\tPeptide\t200\t'+q,'DOMAINS']
        for uid,a,b in hits:
            acc,name=names[uid]
            align.append(f'{q}\tgnl|CDD|{uid}\t{a}\t{b}\t1\t40\t1e-10\t50\t200\t40')
            blocks.append(f'1\t{ql}\tSpecific\t{uid}\t{a}\t{b}\t1e-10\t50\t{acc}\t{name}\t-\t'+('200' if uid==100 else '-'))
        if n==1:
            blocks+=['1\t'+ql+'\tSuperfamily\t200\t1\t40\t1e-9\t49\tcl00001\tSuper_A\t-\t-',
                     '1\t'+ql+'\tMultidom\t300\t1\t200\t1e-20\t100\tPRK00001\tWhole_AB\t-\t-']
        blocks+=['ENDDOMAINS','ENDQUERY']
    native=tmp_path/'native.out';native.write_text('\n'.join(['#Native representative fixture','DATA','SESSION\t1\tblastp\t2.12.0+\tCdd\tBLOSUM62\t0.01',*blocks,'ENDSESSION\t1','ENDDATA'])+'\n')
    raw=tmp_path/'raw.tsv';raw.write_text('\n'.join(align)+'\n')
    return database,path,native,raw


def args(tmp_path,single=False):
    database,proteins,native,raw=fixture(tmp_path,single)
    return ['--input',str(proteins),'--cdd',str(database),'--rpsbproc-output',str(native),'--raw-hits',str(raw),'--output',str(tmp_path/'out'),'--no-plots'],native,raw


def test_cdd_repeats_hierarchy_negative_counts_and_network_refresh(tmp_path):
    argv,native,raw=args(tmp_path);assert cdd.main(argv)==0
    out=tmp_path/'out';hits=rows(out/'domain_hits.tsv');selected=[h for h in hits if h['selected']=='True']
    assert len(selected)==5 and sum(h['accession']=='cd00001' for h in selected)==3
    assert {h['selection_reason'] for h in hits if h['selected']=='False'}=={'superfamily_redundant_with_selected_domain','whole_protein_multidomain_model_not_a_domain'}
    assert all(h['model_coverage']=='1.0' for h in selected)
    assert all(h['specific_bitscore_cutoff'] in ['25.0','30.0'] for h in selected)
    assert 'hmm_from' not in hits[0] and 'independent_evalue' not in hits[0]
    arch=rows(out/'protein_architectures.tsv');assert len(arch)==4
    assert next(x for x in arch if x['protein_id']=='negative')['domain_set']==''
    assert not rows(out/'clusters/c2/cooccurrence_edges.tsv')
    assert rows(out/'clusters/c1/cooccurrence_edges.tsv')[0]['proteins_with_both']=='1'
    import networkx as nx
    graph=nx.read_graphml(out/'all_domains_network.graphml')
    assert graph.graph['domain_source']=='CDD' and graph.number_of_nodes()==2 and graph.number_of_edges()==1
    assert graph.nodes['cd00001']['cdd_pssm_id']=='100'
    assert graph['cd00001']['cd00002']['jaccard']==pytest.approx(1/3)
    assert network_reports.main(['--domains',str(out)])==0
    assert DomainMembers(out).export_pair(['cd00001','cd00002'],tmp_path/'pair.tsv','c2')==0
    assert cdd.main(argv+['--resume'])==0
    assert json.loads((out/'manifest.json').read_text())['search_definition']['cache_reused']
    native.write_text(native.read_text()+'# changed\n')
    assert cdd.main(argv+['--resume'])==0
    assert not json.loads((out/'manifest.json').read_text())['search_definition']['cache_reused']


def test_single_fasta_native_schema_and_existing_membership(tmp_path):
    argv,native,raw=args(tmp_path,True);assert cdd.main(argv)==0
    out=tmp_path/'out';summary=rows(out/'collection_summary.tsv')[0]
    assert summary['total_proteins']=='4' and summary['proteins_without_domains']=='1'
    assert not (out/'clusters').exists() and not (out/'cluster_summary.tsv').exists()
    assert 'cluster' not in rows(out/'protein_architectures.tsv')[0]
    a=next(x for x in rows(out/'domain_prevalence.tsv') if x['accession']=='cd00001')
    assert a['proteins_with_domain']=='2' and a['domain_copies']=='3' and a['fraction_proteins']=='0.5'
    mapping=tmp_path/'membership.tsv';mapping.write_text('protein_id\tcluster\npair\tc1\nonly_a\tc2\nonly_b\tc2\nnegative\tc3\n')
    assert network_reports.main(['--domains',str(out),'--membership',str(mapping)])==0


@pytest.mark.parametrize('change',[
    lambda s:s.replace('p0000000001','renamed'),
    lambda s:s.replace('Peptide\t200','Peptide\t201'),
    lambda s:s.replace('\t1\t40\t1e-10','\t1\t201\t1e-10'),
    lambda s:s.replace('\tcd00001\t','\tcd99999\t'),
    lambda s:s.replace('ENDDATA',''),
    lambda s:s.replace('\tblastp\t','\tRPSTBLASTN\t'),
    lambda s:s[:s.index('QUERY\tQuery_4')]+ 'ENDSESSION\t1\nENDDATA\n',
])
def test_bad_native_data_rejected_without_successful_cache(tmp_path,change):
    argv,native,raw=args(tmp_path);native.write_text(change(native.read_text()))
    assert cdd.main(argv)==2
    assert not (tmp_path/'out/manifest.json').exists() and not (tmp_path/'out/work/cdd.stage.json').exists()


def test_missing_coverage_remains_explicit(tmp_path):
    argv,native,raw=args(tmp_path);raw.write_text('')
    assert cdd.main(argv+['--min-model-coverage','0.5'])==0
    hits=rows(tmp_path/'out/domain_hits.tsv')
    assert all(h['model_coverage']=='' and h['selected']=='False' for h in hits)
    assert rows(tmp_path/'out/cluster_summary.tsv')[0]['distinct_domains']=='0'


def test_duplicate_alignment_regions_not_repeats(tmp_path):
    argv,native,raw=args(tmp_path)
    first=next(x for x in native.read_text().splitlines() if '\tSpecific\t100\t1\t40\t' in x)
    native.write_text(native.read_text().replace(first,first+'\n'+first,1))
    assert cdd.main(argv)==0
    assert len([x for x in rows(tmp_path/'out/domain_hits.tsv') if x['selected']=='True'])==5


def test_derived_superfamily_and_legacy_session_format(tmp_path):
    argv,native,raw=args(tmp_path)
    s=native.read_text().replace('blastp\t2.12.0+','RPSBLAST 2.12.0+')
    s=s.replace('\tSuperfamily\t200\t','\tSuperfamily\t0\t')
    native.write_text(s)
    assert cdd.main(argv)==0
    h=next(x for x in rows(tmp_path/'out/domain_hits.tsv') if x['hit_type']=='Superfamily')
    assert h['reported_pssm_id']=='0' and h['pssm_id']=='200' and h['model_length']=='0'
    assert h['model_coverage']=='' and h['selected']=='False'


@pytest.mark.parametrize('nonspecific',['Non-Specific','Non-specific'])
def test_genuine_superfamily_fallback_and_nonspecific_domain(tmp_path,nonspecific):
    argv,native,raw=args(tmp_path)
    lines=native.read_text().splitlines()
    # Generic region in another part of the sequence is retained as a generic hit.
    for i,line in enumerate(lines):
        if '\tSuperfamily\t200\t' in line:
            lines[i]=line.replace('\t1\t40\t','\t175\t200\t')
    native.write_text('\n'.join(lines).replace('\tSpecific\t101\t','\t'+nonspecific+'\t101\t')+'\n')
    assert cdd.main(argv)==0
    hits=rows(tmp_path/'out/domain_hits.tsv')
    h=next(x for x in hits if x['hit_type']=='Superfamily')
    assert h['selected']=='True' and h['model_coverage']==''
    assert len([x for x in hits if x['hit_type']==nonspecific and x['selected']=='True'])==2


@pytest.mark.parametrize('kind',['path','symlink'])
def test_download_archive_validation(tmp_path,kind):
    archive=tmp_path/'bad.tgz'
    with tarfile.open(archive,'w:gz') as t:
        m=tarfile.TarInfo('../outside' if kind=='path' else 'link')
        if kind=='symlink':m.type=tarfile.SYMTYPE;m.linkname='/etc/passwd'
        t.addfile(m)
    with pytest.raises(ValueError):setup.extract(archive,tmp_path/'db')


def test_snapshot_changed_model_invalidates(tmp_path):
    database,_,_,_=fixture(tmp_path);assert setup.checked_snapshot(database)['subset']=='full'
    (database/'snapshot/db/Cdd.rps').write_text('changed')
    with pytest.raises(ValueError,match='changed'):setup.checked_snapshot(database)


@pytest.mark.parametrize('change',['missing','edited'])
def test_readme_drift_is_recorded_without_altering_snapshot(tmp_path,caplog,change):
    argv,native,raw=args(tmp_path)
    database=tmp_path/'CDD';manifest_path=database/'cdd_manifest.json'
    original=json.loads(manifest_path.read_text())
    for name in ['CDD_README.txt','rpsbproc_README.txt']:
        p=database/'snapshot'/name;p.write_text('original documentation')
        original['files'][str(p.relative_to(database))]=sha256(p)
    dump(manifest_path,original);saved=manifest_path.read_bytes()
    for name in ['CDD_README.txt','rpsbproc_README.txt']:
        p=database/'snapshot'/name
        if change=='missing':p.unlink()
        else:p.write_text('edited documentation')
    # Exercise the actual analysis path: all other snapshot checks still run.
    assert cdd.main(argv)==0
    output=json.loads((tmp_path/'out/manifest.json').read_text())
    verification=output['database']['verification']
    assert verification['search_annotation_files_verified']==len(original['files'])-2
    assert len(verification['documentation'])==2
    expected='missing' if change=='missing' else 'checksum_mismatch'
    assert {r['status'] for r in verification['documentation']}=={expected}
    assert all(r['expected_sha256'] for r in verification['documentation'])
    assert 'documentation' in caplog.text and manifest_path.read_bytes()==saved


def test_readme_warning_never_hides_corrupted_scientific_data(tmp_path):
    database,_,_,_=fixture(tmp_path)
    manifest_path=database/'cdd_manifest.json';m=json.loads(manifest_path.read_text())
    m['files']['snapshot/CDD_README.txt']='0'*64;dump(manifest_path,m)
    (database/'snapshot/db/Cdd.rps').write_text('corrupted')
    (database/'snapshot/data/cddid.tbl').unlink()
    with pytest.raises(ValueError) as error:setup.checked_snapshot(database)
    assert 'Cdd.rps: checksum_mismatch' in str(error.value)
    assert 'cddid.tbl: missing' in str(error.value)


def test_only_known_snapshot_root_readmes_are_noncritical(tmp_path):
    database,_,_,_=fixture(tmp_path)
    manifest_path=database/'cdd_manifest.json';m=json.loads(manifest_path.read_text())
    m['files']['snapshot/db/CDD_README.txt']='0'*64;dump(manifest_path,m)
    with pytest.raises(ValueError,match='snapshot/db/CDD_README.txt: missing'):
        setup.checked_snapshot(database)


def test_download_complete_snapshot_resume_and_failed_update(tmp_path,monkeypatch):
    calls=[]
    def fetch(url,path):
        calls.append(url)
        if url.endswith('Cdd_LE.tar.gz'):
            with tarfile.open(path,'w:gz') as t:
                for name in ['Cdd.rps','Cdd.loo']:
                    b=b'database';m=tarfile.TarInfo(name);m.size=len(b);t.addfile(m,io.BytesIO(b))
        elif url.endswith('cdd.info'):path.write_text('cdd version fixture\n')
        elif url.endswith('.gz'):
            import gzip
            with gzip.open(path,'wb') as f:f.write(b'fixture\n')
        else:path.write_text('metadata\n')
        return dict(url=url,sha256=sha256(path),bytes=path.stat().st_size)
    monkeypatch.setattr(setup,'download',fetch)
    out=tmp_path/'CDD';assert setup.main(['-o',str(out)])==0
    manifest=setup.checked_snapshot(out);assert (out/manifest['db_prefix']).with_suffix('.rps').exists()
    n=len(calls);assert setup.main(['-o',str(out),'--resume'])==0 and len(calls)==n
    assert setup.main(['-o',str(out),'--subset','ncbi','--resume'])==2
    assert setup.main(['-o',str(out),'--force'])==0
    assert len(list((out/'snapshots').iterdir()))==2
    old=(out/'cdd_manifest.json').read_bytes()
    monkeypatch.setattr(setup,'download',lambda *a:(_ for _ in ()).throw(OSError('injected disk full')))
    assert setup.main(['-o',str(out),'--force'])==2
    assert (out/'cdd_manifest.json').read_bytes()==old and not (out/'.download.lock').exists()


def test_pfams_cannot_be_overwritten(tmp_path):
    argv,_,_=args(tmp_path);out=tmp_path/'out';out.mkdir();dump(out/'manifest.json',{'hmm_sha256':'pfam'})
    assert cdd.main(argv+['--resume'])==2
    assert json.loads((out/'manifest.json').read_text())=={'hmm_sha256':'pfam'}


def test_real_rpsblast_rpsbproc_snapshot(tmp_path):
    snapshot=os.environ.get('EET_CDD_TEST_DATABASE');proteins=os.environ.get('EET_CDD_TEST_PROTEINS')
    if not snapshot or not proteins:pytest.skip('Set EET_CDD_TEST_DATABASE and EET_CDD_TEST_PROTEINS for the downloaded SMART control snapshot')
    argv=['--cdd',snapshot,'--input',proteins,'--output',str(tmp_path/'out'),'--threads','2','--no-plots']
    assert cdd.main(argv)==0
    out=tmp_path/'out';summary=rows(out/'cluster_summary.tsv')
    assert sum(int(r['n_proteins']) for r in summary)==4
    hits=[r for r in rows(out/'domain_hits.tsv') if r['selected']=='True']
    assert len(hits)==5 and {r['accession'] for r in hits}=={'smart00257','smart00287'}
    assert len([r for r in hits if r['protein_id']=='repeat_pair'])==3
    assert all(r['model_coverage']=='1.0' for r in hits)
    assert rows(out/'clusters/c1/cooccurrence_edges.tsv')[0]['proteins_with_both']=='1'
    assert not rows(out/'clusters/c2/cooccurrence_edges.tsv')
    stamp=(out/'work/rpsblast.asn').stat().st_mtime_ns
    assert cdd.main(argv+['--resume'])==0
    assert (out/'work/rpsblast.asn').stat().st_mtime_ns==stamp
