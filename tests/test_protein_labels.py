"""Annotated headers are display metadata; scientific IDs/graphs remain unchanged."""
from argparse import Namespace
import csv
import gzip
import json
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from protein_labels import ProteinLabels, header_label
from domain_members import DomainMembers
from network_reports import main
from test_domain_members import fixture, graph
from test_networks import opts, payload
from test_workflow import read, domain_fixture, normalized
import cluster_domains as cd
import cluster_proteins as cp


def labels(folder):
    folder.mkdir()
    (folder/'c1.faa').write_text('>pair | nr_hit=WP_A | nr_annotation=Putative cytochrome c protein [Species]\nMACDE\n'
                               '>only_a | nr_status=no_qualifying_hit\nMACDE\n')
    with gzip.open(folder/'c2.faa.gz','wt') as f:
        f.write('>pair_other | nr_hit=WP_B | nr_annotation=Cell wall protein [Other species]\nMACDE\n')
    return folder


def test_header_import_keeps_ids_annotations_and_missing_states():
    pid,row=header_label('ID|one old note | nr_hit=WP_1 | nr_annotation=SLH protein [Species]')
    assert pid=='ID|one' and row['annotation']=='SLH protein [Species]' and row['nr_hit']=='WP_1'
    assert header_label('id | nr_status=no_qualifying_hit')[1]['annotation']==''
    assert header_label('id | nr_hit=WP | nr_annotation=[title unavailable in DIAMOND database]')[1]['annotation_status']=='title_unavailable'
    assert header_label('id plain descriptive name')[1]['annotation']=='plain descriptive name'
    assert header_label('id')[1]['annotation']==''


def test_domain_carriers_and_pair_export_include_annotations_without_changing_counts(tmp_path):
    domains=tmp_path/'domains';domains.mkdir();fixture(domains)
    source=labels(tmp_path/'annotated')
    imported=ProteinLabels(source,tmp_path/'cache')
    members=DomainMembers(domains,imported)
    data=members.payload(graph(),list(graph()),20000)
    proteins={p['protein_id']:p for p in data['proteins']}
    assert len(proteins)==4 and data['matching_proteins']==4
    assert proteins['pair']['annotation']=='Putative cytochrome c protein [Species]'
    assert proteins['pair']['domain_names']=='SLH;LysM'
    assert 'annotation' not in proteins['only_b']  # partial imports retain unannotated IDs
    assert not proteins['only_a']['annotation']
    members.export_pair(['PF00395','PF01476'],tmp_path/'pair.tsv')
    exported=read(tmp_path/'pair.tsv')
    assert len(exported)==2 and {r['nr_hit'] for r in exported}=={'WP_A','WP_B'}
    assert imported.provenance['matches']['domain_architectures']['matched']==3
    assert imported.provenance['matches']['domain_architectures']['without_imported_header']==2


def test_header_cache_refreshes_and_rejects_duplicate_or_renamed_ids(tmp_path):
    source=labels(tmp_path/'annotated');work=tmp_path/'cache';imported=ProteinLabels(source,work)
    before=imported.path.stat().st_mtime_ns
    assert ProteinLabels(source,work).path.stat().st_mtime_ns==before
    file=source/'c1.faa';file.write_text(file.read_text().replace('Putative cytochrome c protein','Updated description'))
    newer=ProteinLabels(source,work)
    rows=[dict(protein_id='pair'),dict(protein_id='pair_not_real')];newer.decorate(rows)
    assert rows[0]['annotation']=='Updated description [Species]' and 'annotation' not in rows[1]
    domains=tmp_path/'domains';domains.mkdir();fixture(domains)
    wrong=tmp_path/'renamed.faa';wrong.write_text('>changed_ID new name\nMACDE\n')
    with pytest.raises(ValueError,match='first whitespace-delimited'):
        DomainMembers(domains,ProteinLabels(wrong,work))
    (source/'duplicate.faa').write_text('>pair duplicated ID\nMACDE\n')
    with pytest.raises(ValueError,match='globally unique'):ProteinLabels(source,work)


def test_report_regeneration_reads_annotated_fasta_without_rescanning(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'domains'
    assert cd.main(['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'--no-plots'])==0
    arch_before=(out/'protein_architectures.tsv').read_bytes()
    raw=out/'work/domains.domtblout';raw_time=raw.stat().st_mtime_ns
    annotated=tmp_path/'annotated';annotated.mkdir()
    for source in inp.glob('*.faa'):
        lines=source.read_text().splitlines()
        (annotated/source.name).write_text('\n'.join(line+' | nr_hit=WP_TEST | nr_annotation=Protein display name' if line.startswith('>') else line for line in lines)+'\n')
    assert main(['--domains',str(out),'--protein-fastas',str(annotated)])==0
    data=payload(out/'all_domains_network.html')['members']
    assert next(p for p in data['proteins'] if p['protein_id']=='one')['annotation']=='Protein display name'
    assert raw.stat().st_mtime_ns==raw_time and (out/'protein_architectures.tsv').read_bytes()==arch_before
    manifest=json.loads((out/'network_report_manifest.json').read_text())
    assert manifest['protein_header_import']['matches']['domain_architectures']['matched']==4
    assert main(['--domains',str(out),'--protein-fastas',str(annotated),'--pair','TEST_A','TEST_B','--proteins-out',str(tmp_path/'pair.tsv')])==0
    assert read(tmp_path/'pair.tsv')[0]['annotation']=='Protein display name'


def test_existing_clustering_folder_is_automatically_read_and_graph_retained(tmp_path):
    clustering=tmp_path/'clustering';clustering.mkdir()
    db,_=normalized(clustering,{'a':'M'*100,'b':'M'*100})
    hits=clustering/'hits';hits.write_text('a b 75 90 1 90 1 90 1e-20 100 100 100\n')
    args=Namespace(evalue=1e-5,min_identity=30,min_coverage=.5,weight='neglog10-evalue',evalue_cap=200)
    cp.build_graph(db,hits,clustering/'edges',args);cp.write_clusters(db,clustering/'norm.faa',[[1,2]],clustering);db.close()
    work=clustering/'work';work.mkdir();(clustering/'index.db').rename(work/'proteins.sqlite')
    fasta=clustering/'clusters/cluster_000001.faa'
    fasta.write_text(fasta.read_text().replace('>a original header','>a | nr_hit=WP_A | nr_annotation=Cytochrome c protein'))
    assert main(['--clustering',str(clustering)])==0
    data=payload(clustering/'protein_network.html')
    assert data['nodes'][0]['protein_id']=='a' and data['nodes'][0]['label']=='a'
    assert data['nodes'][0]['annotation']=='Cytochrome c protein'
    assert data['edges'][0]['identity_percent']==75
    assert read(clustering/'protein_network_nodes.tsv')[0]['annotation']=='Cytochrome c protein'


def test_unsafe_annotation_is_only_embedded_as_escaped_data(tmp_path):
    domains=tmp_path/'domains';domains.mkdir();fixture(domains)
    fasta=tmp_path/'unsafe.faa';fasta.write_text('>pair | nr_hit=WP_X | nr_annotation=</script><img src=x onerror=alert(1)>\nMACDE\n')
    from interactive_network import domain_network
    domain_network(graph(),domains,opts(),members=DomainMembers(domains,ProteinLabels(fasta,tmp_path/'cache')))
    html=(domains/'domain_network.html').read_text()
    assert '<img src=x' not in html
    assert '</script><img' in next(p for p in payload(domains/'domain_network.html')['members']['proteins'] if p['protein_id']=='pair')['annotation']
