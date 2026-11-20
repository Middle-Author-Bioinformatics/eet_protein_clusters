"""Single-FASTA analysis must match the equivalent one-cluster analysis."""
import gzip
import json
from pathlib import Path
import subprocess
import sys

import networkx as nx
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import cluster_domains as cd
from test_workflow import domain_fixture, read


def test_standalone_matches_folder_statistics_and_exports(tmp_path):
    folder,hmm=domain_fixture(tmp_path)
    fasta=folder/'family_1.faa'
    # Remove the second cluster so the target database and denominator are identical.
    (folder/'family_2.faa').unlink()
    single=tmp_path/'single';clustered=tmp_path/'clustered'
    assert cd.main(['--input',str(fasta),'--hmm',str(hmm),'-o',str(single),'-t','2'])==0
    assert cd.main(['--clusters',str(folder),'--hmm',str(hmm),'-o',str(clustered),'--no-plots','-t','2'])==0
    for name in ['domain_hits.tsv','protein_architectures.tsv']:
        assert read(single/name)==[{k:v for k,v in row.items() if k!='cluster'} for row in read(clustered/name)]
    folder_domains={r['accession']:r for r in read(clustered/'domain_prevalence.tsv')}
    for r in read(single/'domain_prevalence.tsv'):
        original=folder_domains[r['accession']]
        assert r['fraction_proteins']==original['fraction_cluster_proteins']
        assert r['domain_copies']==original['domain_copies']
    domains={r['accession']:r for r in read(single/'domain_prevalence.tsv')}
    assert float(domains['TEST_A.1']['fraction_proteins'])==pytest.approx(2/3)
    assert domains['TEST_A.1']['domain_copies']=='3'
    assert all(not any('cluster' in k or 'outside' in k for k in r) for r in domains.values())
    assert not (single/'cluster_summary.tsv').exists() and not (single/'clusters').exists()
    summary=read(single/'collection_summary.tsv')[0]
    assert summary['dataset']=='family_1' and summary['total_proteins']=='3'
    assert summary['proteins_without_domains']=='1'
    for name in ['domain_hits.tsv','protein_architectures.tsv','id_mapping.tsv','cooccurrence_edges.tsv']:
        assert not any('cluster' in k for k in read(single/name)[0])
    assert read(single/'cooccurrence_edges.tsv')[0]['proteins_with_both']=='1'
    graph=nx.read_graphml(single/'domain_network.graphml')
    assert len(graph)==2 and graph.number_of_edges()==1
    assert graph.graph['input_mode']=='proteins'
    assert (single/'domain_network.svg').read_text().startswith('<?xml')
    assert graph.graph['dataset']=='family_1' and graph.graph['total_proteins']==3
    assert not any('cluster' in k for k in graph.graph)
    assert all(not any('cluster' in k for k in attrs) for _,attrs in graph.nodes(data=True))
    assert all(not any('cluster' in k for k in attrs) for _,_,attrs in graph.edges(data=True))
    report=(single/'index.html').read_text()
    assert 'Protein collection domain profile' in report and 'TEST_A.1' in report
    assert 'cluster' not in report.lower()
    assert (single/'domain_prevalence.svg').exists()
    manifest=json.loads((single/'manifest.json').read_text())
    assert manifest['input_mode']=='proteins' and manifest['protein_collections']==['family_1']
    assert manifest['counts']['proteins']==3
    # Folder mode keeps its existing per-cluster output layout.
    assert not (clustered/'domain_network.graphml').exists()
    assert (clustered/'clusters/family_1/domain_network.graphml').exists()


def test_gzipped_fasta_subprocess_cli_and_resume(tmp_path):
    folder,hmm=domain_fixture(tmp_path)
    fasta=tmp_path/'candidate proteins.faa.gz'
    with gzip.open(fasta,'wt') as f:f.write((folder/'family_1.faa').read_text())
    out=tmp_path/'profile'
    command=[sys.executable,str(Path(cd.__file__)), '-i',str(fasta),'--hmm',str(hmm),'-o',str(out),'--no-plots']
    result=subprocess.run(command,text=True,capture_output=True)
    assert result.returncode==0,result.stderr
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['protein_collections']==['candidate proteins']
    assert manifest['counts']['proteins']==3
    assert not (out/'domain_network.svg').exists()
    assert (out/'domain_network.graphml').exists()
    assert cd.main(command[2:]+['--resume'])==0
    assert json.loads((out/'manifest.json').read_text())['search_definition']['cache_reused']
    assert read(out/'cooccurrence_edges.tsv')[0]['proteins_with_both']=='1'


def test_standalone_resume_cleans_old_plot_and_folder_mode_exports(tmp_path):
    folder,hmm=domain_fixture(tmp_path);fasta=folder/'family_1.faa';out=tmp_path/'profile'
    args=['-i',str(fasta),'--hmm',str(hmm),'-o',str(out)]
    assert cd.main(args)==0 and (out/'domain_network.svg').exists()
    assert cd.main(args+['--resume','--no-plots'])==0
    assert not (out/'domain_network.svg').exists()
    assert (out/'domain_network.graphml').exists()
    assert not (out/'domain_prevalence.svg').exists()
    assert cd.main(['--input',str(folder),'--hmm',str(hmm),'-o',str(out),'--resume','--no-plots'])==0
    for name in ['domain_network.svg','domain_network.graphml','cooccurrence_edges.tsv']:
        assert not (out/name).exists()
    assert len(read(out/'cluster_summary.tsv'))==2
    assert not (out/'collection_summary.tsv').exists()
    # Switching back must remove generated cluster folders and old schemas.
    assert cd.main(args+['--resume','--no-plots'])==0
    assert not (out/'clusters').exists() and not (out/'cluster_summary.tsv').exists()
    assert (out/'collection_summary.tsv').exists()
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['input_mode']=='proteins' and not manifest['search_definition']['cache_reused']


def test_domain_negative_single_fasta(tmp_path):
    folder,hmm=domain_fixture(tmp_path)
    fasta=tmp_path/'negative.faa';fasta.write_text('>unannotated\nM'+'P'*90+'\n')
    out=tmp_path/'profile'
    assert cd.main(['--input',str(fasta),'--hmm',str(hmm),'-o',str(out)])==0
    summary=read(out/'collection_summary.tsv')[0]
    assert summary['total_proteins']=='1' and summary['proteins_with_domains']=='0'
    assert not read(out/'domain_prevalence.tsv')
    assert not read(out/'cooccurrence_edges.tsv')
    assert len(nx.read_graphml(out/'domain_network.graphml'))==0
    assert (out/'domain_network.svg').exists()
    assert read(out/'protein_architectures.tsv')[0]['protein_id']=='unannotated'


@pytest.mark.parametrize('content', ['', '>duplicate\nMACDE\n>duplicate\nMFGHI\n'])
def test_empty_or_duplicate_single_fasta_fails(tmp_path,content):
    folder,hmm=domain_fixture(tmp_path)
    fasta=tmp_path/'invalid.faa';fasta.write_text(content);out=tmp_path/'profile'
    assert cd.main(['--input',str(fasta),'--hmm',str(hmm),'-o',str(out),'--no-plots'])==2
    assert not (out/'manifest.json').exists()


def test_missing_input_and_output_containing_input_fail_safely(tmp_path):
    folder,hmm=domain_fixture(tmp_path)
    missing=tmp_path/'absent.faa';out=tmp_path/'missing_profile'
    assert cd.main(['--input',str(missing),'--hmm',str(hmm),'-o',str(out)])==2
    assert not out.exists()
    fasta=folder/'family_1.faa';before=fasta.read_bytes()
    assert cd.main(['--input',str(fasta),'--hmm',str(hmm),'-o',str(folder),'--force'])==2
    assert fasta.read_bytes()==before


def test_published_cache_migration_and_report_only_resume(tmp_path,monkeypatch):
    import pyhmmer
    from common import dump, sha256, signature
    folder,hmm=domain_fixture(tmp_path);fasta=folder/'family_1.faa';out=tmp_path/'profile'
    args=['-i',str(fasta),'--hmm',str(hmm),'-o',str(out),'--no-plots']
    assert cd.main(args)==0
    work=out/'work';raw=work/'domains.domtblout';raw_time=raw.stat().st_mtime_ns
    stamp=json.loads((work/'search.stage.json').read_text())
    # Recreate the exact key layout used by the published v1.1.0 script.
    stamp['key']=signature(dict(script='dd77ecf24746035c42de6036932c397c67eb498f3b7df6eab8eb9255770b5d3f',common='ddce5e4ebc029d185dd9e45046d13d8daa1bbd1810890bd7b638c1b40403f681',memberships=[('family_1',sha256(fasta))],proteins=sha256(work/'proteins.faa'),hmm=sha256(hmm),pyhmmer=pyhmmer.__version__,fallback=None))
    dump(work/'search.stage.json',stamp)
    (out/'cluster_summary.tsv').write_text('legacy output')
    legacy=out/'clusters/family_1';legacy.mkdir(parents=True)
    (legacy/'domain_prevalence.tsv').write_text('legacy output')
    def cannot_search(*args,**kwargs):raise AssertionError('Presentation refresh must not rerun HMMER')
    monkeypatch.setattr(pyhmmer.hmmer,'hmmsearch',cannot_search)
    assert cd.main(args+['--resume'])==0
    assert raw.stat().st_mtime_ns==raw_time
    result=json.loads((out/'manifest.json').read_text())['search_definition']
    assert result['cache_reused'] and result['legacy_cache_migrated']
    assert not (out/'cluster_summary.tsv').exists() and not (out/'clusters').exists()
    assert (out/'collection_summary.tsv').exists()
    assert cd.main(args+['--resume'])==0
    result=json.loads((out/'manifest.json').read_text())['search_definition']
    assert result['cache_reused'] and not result['legacy_cache_migrated']
