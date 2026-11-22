"""Scientific edge metadata, named domains, bounded exports and report-only rebuilds."""
from argparse import Namespace
import json
from pathlib import Path
import re
import sqlite3
import sys
import pytest
import networkx as nx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import cluster_proteins as cp
import cluster_domains as cd
from interactive_network import write_html, domain_network, domain_communities
from network_reports import refresh_proteins
from test_workflow import normalized, domain_fixture, read


def payload(path):
    text=path.read_text();return json.loads(re.search(r'<script id="network-data" type="application/json">(.*?)</script>',text,re.S)[1])


def opts(**values):
    return Namespace(network_max_nodes=1000,network_max_edges=5000,no_interactive=False,network_cluster=None,**values)


def test_graph_metadata_follows_selected_weight_and_swaps_coverage(tmp_path):
    db,_=normalized(tmp_path,{'a':'M'*100,'b':'M'*200,'iso':'M'*50})
    hits=tmp_path/'hits';hits.write_text('a b 80 80 1 80 1 100 1e-10 80 100 200\nb a 70 90 1 180 1 90 1e-30 140 200 100\n')
    args=Namespace(evalue=1e-5,min_identity=30,min_coverage=.5,weight='neglog10-evalue',evalue_cap=200)
    cp.build_graph(db,hits,tmp_path/'edges',args)
    row=db.execute('SELECT weight,identity_percent,coverage_a,coverage_b,evalue,bitscore FROM edges').fetchone()
    assert row==(30,70,.9,.9,1e-30,140)
    # Metadata belongs to the max-weight qualifying HSP, not the highest identity.
    db.close()


def test_protein_colors_keep_mcl_membership_and_isolates(tmp_path):
    db,_=normalized(tmp_path,{'a':'M'*100,'b':'M'*100,'iso':'M'*50})
    hits=tmp_path/'hits';hits.write_text('a b 75 90 1 90 1 90 1e-20 100 100 100\n')
    args=Namespace(evalue=1e-5,min_identity=30,min_coverage=.5,weight='neglog10-evalue',evalue_cap=200)
    cp.build_graph(db,hits,tmp_path/'edges',args);cp.write_clusters(db,tmp_path/'norm.faa',[[1,2],[3]],tmp_path)
    from interactive_network import protein_network
    protein_network(db,tmp_path,opts());p=payload(tmp_path/'protein_network.html')
    assert len(p['nodes'])==3 and len(p['edges'])==1 and p['edges'][0]['identity_percent']==75
    assert [n['group'] for n in p['nodes']]==['cluster_000001','cluster_000001','cluster_000002']
    assert 'not genome-wide' in p['definition']
    assert len(read(tmp_path/'protein_network_nodes.tsv'))==3
    db.close()


def test_domain_communities_and_empty_graph(tmp_path):
    g=nx.Graph(input_mode='proteins')
    for a,name in [('A','Cell wall binding'),('B','Heme binding'),('C','Unpaired')]:g.add_node(a,name=name,proteins_with_domain=2,domain_copies=3,fraction_proteins=.5)
    g.add_edge('A','B',jaccard=.5,proteins_with_both=1,lift=2)
    domain_network(g,tmp_path,opts())
    p=payload(tmp_path/'domain_network.html');assert p['nodes'][0]['label']=='Cell wall binding (A)'
    assert g.nodes['A']['community']==g.nodes['B']['community']!=g.nodes['C']['community']
    assert p['edges'][0]['domain_a_name']=='Cell wall binding'
    empty=nx.Graph();domain_network(empty,tmp_path,opts());assert payload(tmp_path/'domain_network.html')['nodes']==[]
    assert not read(tmp_path/'domain_network_communities.tsv')


def test_html_embeds_safe_json_without_external_assets(tmp_path):
    label='LysM </script><img src=x onerror=alert(1)> @@SCRIPT@@ & binding'
    file=tmp_path/'network.html';write_html(file,[dict(id='A',label=label,group='C',name=label)],[],'<Title>','domains','Definition')
    assert payload(file)['nodes'][0]['label']==label
    content=file.read_text();assert '<img src=x' not in content and '<script src=' not in content
    assert '<title>&lt;Title&gt;</title>' in content and 'fetch(' not in content


def test_named_domain_outputs_pooled_counts_and_limits(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'out'
    args=['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'--no-plots','--network-max-nodes','1']
    assert cd.main(args)==0
    rows=read(out/'protein_architectures.tsv');r=next(r for r in rows if r['protein_id']=='one')
    assert r['domain_names']=='A;B' and r['architecture_names'].count('A (TEST_A.1)')==2
    edge=read(out/'clusters/family_1/cooccurrence_edges.tsv')[0]
    assert edge['domain_a_name']=='A' and edge['domain_b_name']=='B'
    pooled=nx.read_graphml(out/'all_domains_network.graphml')
    assert pooled.edges['TEST_A.1','TEST_B.1']['proteins_with_both']==1
    assert pooled.nodes['TEST_A.1']['proteins_with_domain']==2 and pooled.nodes['TEST_B.1']['proteins_with_domain']==2
    p=payload(out/'all_domains_network.html');assert p['total_nodes']==2 and len(p['nodes'])==1
    before=(out/'work/domains.domtblout').stat().st_mtime_ns
    assert cd.main(args+['--resume','--no-interactive'])==0
    assert not (out/'all_domains_network.html').exists() and not (out/'clusters/family_1/domain_network.html').exists()
    assert (out/'all_domains_network.graphml').exists() and (out/'work/domains.domtblout').stat().st_mtime_ns==before


def test_report_only_migrates_old_graph_without_diamond(tmp_path):
    db,_=normalized(tmp_path,{'a':'M'*100,'b':'M'*100,'iso':'M'*100})
    db.execute('CREATE TABLE edges (a INTEGER,b INTEGER,weight REAL)');db.execute('INSERT INTO edges VALUES (1,2,20)');db.commit()
    cp.write_clusters(db,tmp_path/'norm.faa',[[1,2],[3]],tmp_path);db.close()
    work=tmp_path/'work';work.mkdir();(tmp_path/'index.db').rename(work/'proteins.sqlite')
    hits=work/'diamond_hits.tsv';hits.write_text('p0000000001 p0000000002 85 90 1 90 1 90 1e-20 100 100 100\n')
    from common import dump,mark_stage
    params=dict(evalue=1e-5,min_identity=30,min_coverage=.5,weight='neglog10-evalue',evalue_cap=200,hits=None)
    dump(tmp_path/'manifest.json',dict(parameters=params));mark_stage(work/'diamond.stage.json','legacy',[hits])
    refresh_proteins(tmp_path,opts());assert payload(tmp_path/'protein_network.html')['edges'][0]['identity_percent']==85
    focused=opts();focused.network_cluster='cluster_000002';refresh_proteins(tmp_path,focused)
    p=payload(tmp_path/'protein_network.html');assert len(p['nodes'])==1 and not p['edges']
    focused.network_cluster='missing'
    with pytest.raises(ValueError,match='Unknown'):refresh_proteins(tmp_path,focused)
