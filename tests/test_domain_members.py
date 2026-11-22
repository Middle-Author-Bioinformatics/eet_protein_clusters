"""Protein identity in domain views must come from same-protein architectures."""
import csv
from pathlib import Path
import sqlite3
import sys
import networkx as nx
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from domain_members import DomainMembers
from interactive_network import domain_network, validate
from network_reports import main
from test_networks import opts, payload
from test_workflow import read, domain_fixture
import cluster_domains as cd


def fixture(folder, clustered=True):
    rows=[('pair', 'c1', 'PF00395.29;PF01476.21;PF00395.29'),
          ('only_a', 'c1', 'PF00395.29'), ('only_b', 'c1', 'PF01476.21'),
          ('pair_other', 'c2', 'PF00395.29;PF01476.21'), ('negative', 'c2', '')]
    with open(folder/'protein_architectures.tsv','w') as f:
        fields=['protein_id']+(['cluster'] if clustered else [])+['domain_set','domain_names','architecture_names']
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
        for pid, group, domains in rows:
            row=dict(protein_id=pid,domain_set=domains,domain_names='SLH;LysM' if domains else '',architecture_names='Named architecture')
            if clustered:row['cluster']=group
            w.writerow(row)
    return DomainMembers(folder)


def graph(cluster=None):
    g=nx.Graph(input_mode='clusters' if cluster else 'proteins')
    if cluster:g.graph['cluster']=cluster
    for a,name in [('PF00395.29','SLH'),('PF01476.21','LysM')]:g.add_node(a,name=name,proteins_with_domain=2 if cluster else 3,domain_copies=4)
    g.add_edge('PF00395.29','PF01476.21',jaccard=.4,proteins_with_both=1 if cluster else 2,lift=1)
    return g


def pair_ids(data):
    a,b=data['by_domain'].values()
    return {data['proteins'][i]['protein_id'] for i in set(a)&set(b)}


def test_carriers_repeat_dedup_scope_and_no_cross_protein_pairs(tmp_path):
    index=fixture(tmp_path)
    pooled=index.payload(graph(),list(graph()),20000)
    assert pooled['matching_proteins']==4 and not pooled['truncated']
    assert pair_ids(pooled)=={'pair','pair_other'}
    single=index.payload(graph('c1'),list(graph()),20000)
    assert pair_ids(single)=={'pair'} and single['matching_proteins']==3
    assert single['source_sha256']==index.checksum
    index.export_pair(['PF00395','PF01476'],tmp_path/'pair.tsv','c1')
    rows=read(tmp_path/'pair.tsv');assert [r['protein_id'] for r in rows]==['pair']
    assert rows[0]['domain_a']=='PF00395.29' and rows[0]['cluster']=='c1'


def test_native_fasta_caps_zero_full_export_and_cache_invalidation(tmp_path):
    index=fixture(tmp_path,False)
    before=index.path.stat().st_mtime_ns
    assert DomainMembers(tmp_path).path.stat().st_mtime_ns==before
    capped=index.payload(graph(),list(graph()),1)
    assert capped['truncated'] and capped['embedded_count']==1 and capped['matching_proteins']==4
    disabled=index.payload(graph(),list(graph()),0)
    assert disabled['embedded_count']==0 and disabled['truncated']
    assert main(['--domains',str(tmp_path),'--pair','PF00395','PF01476','--proteins-out',str(tmp_path/'pair.tsv')])==0
    rows=read(tmp_path/'pair.tsv');assert len(rows)==2 and 'cluster' not in rows[0]
    assert not (tmp_path/'network_report_manifest.json').exists()  # export does not refresh HTML
    with pytest.raises(ValueError,match='single-FASTA'):index.export_pair(['PF00395','PF01476'],tmp_path/'x.tsv','c1')
    with open(index.source,'a') as f:f.write('extra\tPF00395.29\tSLH\tSLH\n')
    newer=DomainMembers(tmp_path);assert newer.checksum!=index.checksum
    assert not newer.payload(graph(),list(graph()),100)['available']


def test_missing_or_mismatched_architectures_never_invent_protein_ids(tmp_path):
    index=DomainMembers(tmp_path)
    domain_network(graph(),tmp_path,opts(),members=index)
    assert not payload(tmp_path/'domain_network.html')['members']['available']
    with pytest.raises(ValueError,match='GraphML'):index.export_pair(['A','B'],tmp_path/'x')
    index=fixture(tmp_path)
    g=graph();g.edges['PF00395.29','PF01476.21']['proteins_with_both']=1
    result=index.payload(g,list(g),20)
    assert not result['available'] and 'pair counts' in result['reason']


def test_duplicate_proteins_and_bad_pair_requests_fail(tmp_path):
    index=fixture(tmp_path)
    with pytest.raises(ValueError,match='Unknown'):index.export_pair(['A','B'],tmp_path/'x')
    with pytest.raises(ValueError,match='different'):index.export_pair(['PF00395','PF00395'],tmp_path/'x')
    with pytest.raises(ValueError,match='Unknown cluster'):index.export_pair(['PF00395','PF01476'],tmp_path/'x','missing')
    with open(index.source,'a') as f:f.write('pair\tc1\tPF00395.29\tSLH\tSLH\n')
    with pytest.raises(sqlite3.IntegrityError):DomainMembers(tmp_path)
    assert not list((tmp_path/'work').glob('network_members-*.sqlite'))


def test_normal_search_and_report_only_rebuild_embed_actual_carriers(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'out'
    args=['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'--no-plots']
    assert cd.main(args)==0
    pooled=payload(out/'all_domains_network.html')['members']
    assert pair_ids(pooled)=={'one'}
    single=payload(out/'clusters/family_1/domain_network.html')['members']
    assert pair_ids(single)=={'one'}
    raw=out/'work/domains.domtblout';stamp=raw.stat().st_mtime_ns
    assert main(['--domains',str(out),'--network-max-proteins','0'])==0
    assert payload(out/'all_domains_network.html')['members']['embedded_count']==0
    assert raw.stat().st_mtime_ns==stamp
    assert main(['--domains',str(out),'--pair','TEST_A','TEST_B','--proteins-out',str(out/'pair.tsv')])==0
    assert [r['protein_id'] for r in read(out/'pair.tsv')]==['one']


def test_cross_view_links_and_limits(tmp_path):
    domains=tmp_path/'domains';domains.mkdir()
    index=fixture(domains)
    args=opts();args.clustering=tmp_path/'clustering';args.network_max_proteins=20000
    domain_network(graph(),domains,args,members=index)
    text=(domains/'domain_network.html').read_text()
    assert 'href="../clustering/protein_network.html"' in text and 'Complete protein architectures' in text
    args.network_max_proteins=-1
    with pytest.raises(ValueError,match='network-max-proteins'):validate(args)


def test_single_fasta_search_preserves_native_carrier_output(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'native'
    assert cd.main(['-i',str(inp/'family_1.faa'),'--hmm',str(hmm),'-o',str(out),'--no-plots'])==0
    data=payload(out/'domain_network.html')['members']
    assert pair_ids(data)=={'one'} and not data['truncated']
    assert all(not row['cluster'] for row in data['proteins'])
    assert 'cluster' not in read(out/'protein_architectures.tsv')[0]
    assert main(['--domains',str(out),'--pair','TEST_A','TEST_B','--proteins-out',str(out/'pair.tsv')])==0
    assert 'cluster' not in read(out/'pair.tsv')[0]
