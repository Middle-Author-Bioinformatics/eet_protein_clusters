"""Cluster highlights use complete architectures, not the sampled HTML proteins."""
import csv
import json
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from domain_members import DomainMembers
from network_reports import main
from test_domain_members import fixture,graph
from test_networks import payload
from test_workflow import domain_fixture,read
import cluster_domains as cd
import cluster_proteins as cp


def membership(path,rows):
    with open(path,'w') as stream:
        w=csv.writer(stream,delimiter='\t');w.writerow(['protein_id','cluster']);w.writerows(rows)
    return path


def clusters(data):
    return {row['id']:row for row in data['cluster_index']['clusters']}


@pytest.mark.parametrize('limit',[0,1,20000])
def test_cluster_presence_and_counts_are_independent_of_embedding_cap(tmp_path,limit):
    index=fixture(tmp_path)
    data=index.payload(graph(),list(graph()),limit)
    assert data['embedded_count']<=limit
    c=clusters(data)
    assert c['c1']['proteins']==3 and c['c2']['proteins']==2
    assert c['c1']['domains']=={'PF00395.29':2,'PF01476.21':2}
    assert c['c2']['domains']=={'PF00395.29':1,'PF01476.21':1}
    assert all(r['cluster'] in c for r in data['proteins'])


def test_node_caps_negative_clusters_and_per_cluster_scope(tmp_path):
    fixture(tmp_path)
    with open(tmp_path/'protein_architectures.tsv','a') as stream:stream.write('negative_new\tc3\t\t\t\n')
    index=DomainMembers(tmp_path)
    data=index.payload(graph(),['PF00395.29'],0)
    c=clusters(data)
    assert c['c1']['total_domains']==2 and c['c1']['domains']=={'PF00395.29':2}
    assert c['c3']['proteins']==1 and c['c3']['domains']=={} and c['c3']['total_domains']==0
    focused=index.payload(graph('c1'),list(graph()),20)
    assert set(clusters(focused))=={'c1'}


def test_native_architectures_get_known_cluster_ids_from_membership(tmp_path):
    fixture(tmp_path,False)
    mapping=membership(tmp_path/'membership.tsv',[('pair','mcl_1'),('only_a','mcl_2'),('only_b','mcl_2'),('pair_other','mcl_3'),('negative','mcl_4')])
    index=DomainMembers(tmp_path,membership=mapping)
    data=index.payload(graph(),list(graph()),20000)
    assert {p['protein_id']:p['cluster'] for p in data['proteins']}=={'pair':'mcl_1','only_a':'mcl_2','only_b':'mcl_2','pair_other':'mcl_3'}
    assert clusters(data)['mcl_2']['domains']=={'PF00395.29':1,'PF01476.21':1}
    # Both domains are present in mcl_2, but they occur in different proteins.
    # A highlight does not create a within-cluster co-occurrence edge.
    assert index.export_pair(['PF00395','PF01476'],tmp_path/'pair.tsv','mcl_2')==0
    assert not read(tmp_path/'pair.tsv')
    assert index.export_pair(['PF00395','PF01476'],tmp_path/'pair.tsv','mcl_1')==1
    assert read(tmp_path/'pair.tsv')[0]['cluster']=='mcl_1'
    assert clusters(data)['mcl_4']['total_domains']==0
    assert index.cluster_provenance['membership']['sha256']


def test_missing_partial_updated_and_conflicting_membership(tmp_path):
    fixture(tmp_path,False)
    index=DomainMembers(tmp_path)
    assert not index.payload(graph(),list(graph()),100)['cluster_index']['available']
    mapping=membership(tmp_path/'membership.tsv',[('pair','mcl_1')])
    index=DomainMembers(tmp_path,membership=mapping)
    data=index.payload(graph(),list(graph()),100)
    assert data['cluster_index']['unassigned_proteins']==4
    assert next(p for p in data['proteins'] if p['protein_id']=='only_a')['cluster']==''
    mapping=membership(mapping,[('pair','mcl_9')])
    changed=DomainMembers(tmp_path,membership=mapping)
    assert set(clusters(changed.payload(graph(),list(graph()),0)))=={'mcl_9'}
    fixture(tmp_path,True)
    with pytest.raises(ValueError,match='conflicts'):DomainMembers(tmp_path,membership=mapping)


@pytest.mark.parametrize('rows,match',[
    ([('pair','c1'),('pair','c1')],'Duplicate'),
    ([('pair','')],'Empty'),
    ([('not_an_id','c1')],'do not match'),
])
def test_invalid_mapping_never_assigns_guessed_cluster(tmp_path,rows,match):
    fixture(tmp_path,False)
    path=membership(tmp_path/'membership.tsv',rows)
    with pytest.raises(ValueError,match=match):DomainMembers(tmp_path,membership=path)


def test_report_only_single_fasta_mapping_preserves_domain_search(tmp_path):
    folder,hmm=domain_fixture(tmp_path);out=tmp_path/'domains'
    assert cd.main(['-i',str(folder/'family_1.faa'),'--hmm',str(hmm),'-o',str(out),'--no-plots'])==0
    raw=out/'work/domains.domtblout';stamp=raw.stat().st_mtime_ns
    arch=(out/'protein_architectures.tsv').read_bytes()
    mapping=membership(tmp_path/'membership.tsv',[('one','cluster_1'),('two','cluster_2'),('no_hits','cluster_3')])
    assert main(['--domains',str(out),'--membership',str(mapping),'--network-max-proteins','0'])==0
    data=payload(out/'domain_network.html')['members']
    assert not data['proteins'] and set(clusters(data))=={'cluster_1','cluster_2','cluster_3'}
    assert raw.stat().st_mtime_ns==stamp and (out/'protein_architectures.tsv').read_bytes()==arch
    manifest=json.loads((out/'network_report_manifest.json').read_text())
    assert manifest['domain_cluster_mapping']['matched']==3
    assert manifest['domain_cluster_mapping']['source']=='membership.tsv'
    assert main(['--domains',str(out),'--membership',str(mapping),'--pair','TEST_A','TEST_B','--cluster','cluster_1','--proteins-out',str(tmp_path/'pair.tsv')])==0
    assert read(tmp_path/'pair.tsv')[0]['protein_id']=='one'


def test_clustering_argument_automatically_imports_membership_for_native_domains(tmp_path):
    from Bio import SeqIO
    from argparse import Namespace
    from test_workflow import normalized
    folder,hmm=domain_fixture(tmp_path);out=tmp_path/'domains';source=folder/'family_1.faa'
    assert cd.main(['-i',str(source),'--hmm',str(hmm),'-o',str(out),'--no-plots'])==0
    root=tmp_path/'clustering';root.mkdir()
    sequences={r.id:str(r.seq) for r in SeqIO.parse(source,'fasta')}
    db,_=normalized(root,sequences)
    hits=root/'empty_hits';hits.write_text('')
    cp.build_graph(db,hits,root/'edges',Namespace(evalue=1e-5,min_identity=30,min_coverage=.5,weight='neglog10-evalue',evalue_cap=200))
    cp.write_clusters(db,root/'norm.faa',[[1,2],[3]],root);db.close()
    work=root/'work';work.mkdir();(root/'index.db').rename(work/'proteins.sqlite')
    assert main(['--domains',str(out),'--clustering',str(root)])==0
    data=payload(out/'domain_network.html')['members']
    assert set(clusters(data))=={'cluster_000001','cluster_000002'}
    assert clusters(data)['cluster_000002']['total_domains']==0
    assert next(p for p in data['proteins'] if p['protein_id']=='one')['cluster']=='cluster_000001'
