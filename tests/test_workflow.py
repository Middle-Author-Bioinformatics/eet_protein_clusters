import csv
import json
import math
import random
import shutil
import sqlite3
import sys
from argparse import Namespace
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import common
import cluster_proteins as cp
import cluster_domains as cd


def read(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


def normalized(tmp_path,seqs):
    path=tmp_path/'p.faa';path.write_text(''.join(f'>{p} original header\n{s}\n' for p,s in seqs.items()))
    db,n=common.normalize([('',path)],tmp_path/'norm.faa',tmp_path/'index.db',tmp_path/'map.tsv')
    return db,n


def test_normalize_preserves_headers_and_duplicates_fail(tmp_path):
    db,n=normalized(tmp_path,{'original|A':'MACDE*','p0000000001':'MFGH'})
    assert n==2
    assert db.execute('SELECT header FROM proteins WHERE n=1').fetchone()[0]=='original|A original header'
    assert common.resolve(db,'p0000000001',False)[0]==2
    assert common.resolve(db,'p0000000001',True)[0]==1
    db.close()
    path=tmp_path/'duplicate.faa';path.write_text('>x\nMAC\n>x\nMAC\n')
    with pytest.raises(ValueError,match='Duplicate'):common.normalize([('',path)],tmp_path/'d.faa',tmp_path/'d.db',tmp_path/'d.tsv')


def test_graph_filters_reciprocal_dedup_and_zero_evalue(tmp_path):
    db,n=normalized(tmp_path,{'a':'M'*100,'b':'M'*100,'c':'M'*200,'iso':'M'*100})
    hits=tmp_path/'hits.tsv'
    hits.write_text('a b 80 90 1 90 1 90 1e-20 100 100 100\nb a 80 90 1 90 1 90 1e-30 110 100 100\na c 80 90 1 90 1 90 1e-20 100 100 200\na a 100 100 1 100 1 100 0 200 100 100\nb c 80 100 1 100 1 100 0 200 100 200\n')
    args=Namespace(evalue=1e-5,min_identity=30,min_coverage=.5,weight='neglog10-evalue',evalue_cap=200)
    stats=cp.build_graph(db,hits,tmp_path/'edges',args)
    assert stats['undirected_edges']==2
    assert db.execute('SELECT weight FROM edges WHERE a=1 AND b=2').fetchone()[0]==30
    assert db.execute('SELECT weight FROM edges WHERE a=2 AND b=3').fetchone()[0]==200
    assert stats['filtered_rows']==1 and stats['self_hits']==1
    mcl=tmp_path/'mcl';mcl.write_text('p0000000001\tp0000000002\tp0000000003\n')
    assert cp.memberships(db,mcl)==[[1,2,3],[4]]
    cp.write_clusters(db,tmp_path/'norm.faa',cp.memberships(db,mcl),tmp_path)
    assert len(list((tmp_path/'clusters').glob('*.faa')))==2
    assert len(read(tmp_path/'membership.tsv'))==4
    assert '>a original header' in (tmp_path/'clusters/cluster_000001.faa').read_text()
    db.close()


@pytest.mark.parametrize('line',['a b nan 20 1 20 1 20 1e-5 50 100 100','a b 40 20 1 200 1 20 1e-5 50 100 100','unknown b 40 20 1 20 1 20 1e-5 50 100 100'])
def test_invalid_hits_rejected(tmp_path,line):
    db,_=normalized(tmp_path,{'a':'M'*100,'b':'M'*100});p=tmp_path/'hits';p.write_text(line+'\n')
    args=Namespace(evalue=1e-5,min_identity=30,min_coverage=.5,weight='bitscore',evalue_cap=200)
    with pytest.raises(ValueError):cp.build_graph(db,p,tmp_path/'edges',args)
    db.close()


def test_mcl_duplicate_partition_rejected(tmp_path):
    db,_=normalized(tmp_path,{'a':'MA','b':'MA'});p=tmp_path/'mcl';p.write_text('p0000000001\np0000000001\n')
    with pytest.raises(ValueError,match='repeats'):cp.memberships(db,p)
    db.close()


def domain(a,start,end,evalue=1e-10,score=50,cov=1):
    return dict(accession=a,start=start,end=end,independent_evalue=evalue,domain_bitscore=score,model_coverage=cov)


def test_overlap_resolution_retains_repeats():
    a=domain('A.1',1,50);ar=domain('A.1',60,110);alt=domain('B.1',2,48,1e-5)
    kept,decisions=cd.select_domains([a,ar,alt])
    assert kept==[a,ar]
    assert any('overlaps_stronger' in r['selection_reason'] for r in decisions)
    assert len(cd.select_domains([a,ar,alt],'all')[0])==3
    assert not cd.select_domains([domain('A',1,10,cov=.2)],min_coverage=.5)[0]


def make_hmm(path,no_ga=False):
    import pyhmmer as ph
    rng=random.Random(19);aa='ACDEFGHIKLMNPQRSTVWY'
    seqs={a:''.join(rng.choices(aa,k=70)) for a in ['A','B']}
    alphabet=ph.easel.Alphabet.amino()
    with open(path,'wb') as f:
        for a,s in seqs.items():
            msa=ph.easel.TextMSA(name=a,sequences=[ph.easel.TextSequence(name='a',sequence=s),ph.easel.TextSequence(name='b',sequence=s)])
            hmm,_,_=ph.plan7.Builder(alphabet).build_msa(msa.digitize(alphabet),ph.plan7.Background(alphabet))
            hmm.accession=f'TEST_{a}.1'
            if not no_ga:hmm.cutoffs.gathering=(10.,10.)
            hmm.write(f)
    return seqs


def domain_fixture(tmp_path,no_ga=False):
    hmm=tmp_path/'models.hmm';seqs=make_hmm(hmm,no_ga)
    inp=tmp_path/'clusters';inp.mkdir()
    a,b=seqs['A'],seqs['B'];gap='GGGGGGGGGG'
    (inp/'family_1.faa').write_text(f'>one\n{a}{gap}{a}{gap}{b}\n>two\n{a}\n>no_hits\n'+'M'+'P'*90+'\n')
    (inp/'family_2.faa').write_text(f'>three\n{b}\n')
    return inp,hmm


def test_real_hmmer_prevalence_repeats_network_and_zero_hits(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'out'
    args=['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'-t','2','--model-batch-size','1']
    assert cd.main(args)==0
    rows=read(out/'clusters/family_1/domain_prevalence.tsv')
    byid={r['accession']:r for r in rows}
    assert int(byid['TEST_A.1']['proteins_with_domain'])==2
    assert int(byid['TEST_A.1']['domain_copies'])==3
    assert float(byid['TEST_A.1']['fraction_cluster_proteins'])==pytest.approx(2/3)
    assert float(byid['TEST_A.1']['mean_copies_when_present'])==1.5
    edge=read(out/'clusters/family_1/cooccurrence_edges.tsv')[0]
    assert int(edge['proteins_with_both'])==1
    assert float(edge['jaccard'])==.5
    assert float(edge['lift'])==1.5
    assert read(out/'cluster_summary.tsv')[0]['n_proteins']=='3'
    assert any(r['protein_id']=='no_hits' and r['n_domain_hits']=='0' for r in read(out/'protein_architectures.tsv'))
    assert (out/'clusters/family_1/domain_network.svg').exists()
    import networkx as nx
    graph=nx.read_graphml(out/'clusters/family_1/domain_network.graphml')
    assert graph.number_of_edges()==1 and len(graph)==2
    assert cd.main(args+['--resume'])==0
    assert json.load(open(out/'manifest.json'))['search_definition']['cache_reused']
    # A changed cluster assignment/header invalidates the search DB backup.
    p=inp/'family_2.faa';p.rename(inp/'family_renamed.faa')
    assert cd.main(args+['--resume','--no-plots'])==0
    assert not json.load(open(out/'manifest.json'))['search_definition']['cache_reused']
    assert not (out/'clusters/family_2').exists()
    assert (out/'clusters/family_renamed').exists()


def test_missing_ga_requires_explicit_fallback(tmp_path):
    inp,hmm=domain_fixture(tmp_path,True);out=tmp_path/'out'
    args=['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'--no-plots']
    assert cd.main(args)==2
    assert not (out/'manifest.json').exists()
    assert cd.main(args+['--force','--domain-evalue','1e-4'])==0
    assert all(r['mode']=='Evalue' for r in read(out/'model_inventory.tsv'))


def test_low_memory_search(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'out'
    assert cd.main(['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'--low-memory','--no-plots','-t','2'])==0
    assert len(read(out/'domain_hits.tsv'))>=4


def test_graphml_empty_edges_single_domain(tmp_path):
    inp,hmm=domain_fixture(tmp_path);out=tmp_path/'out'
    assert cd.main(['-i',str(inp),'--hmm',str(hmm),'-o',str(out),'--min-pair-count','2','--no-plots'])==0
    import networkx as nx
    assert len(read(out/'clusters/family_1/cooccurrence_edges.tsv'))==1
    assert nx.read_graphml(out/'clusters/family_1/domain_network.graphml').number_of_edges()==0
    assert len(read(out/'clusters/family_2/domain_prevalence.tsv'))==1


def test_external_diamond_mcl_end_to_end(tmp_path):
    if not all(shutil.which(t) for t in ['diamond','mcl','mcxload']):pytest.skip('DIAMOND/MCL binaries not available')
    rng=random.Random(34);aa='ACDEFGHIKLMNPQRSTVWY'
    a=''.join(rng.choices(aa,k=200));b=''.join(rng.choices(aa,k=200));c=''.join(rng.choices(aa,k=200))
    inp=tmp_path/'all.faa';inp.write_text(f'>a1 descriptive\n{a}\n>a2\n{a}\n>b1\n{b}\n>b2\n{b}\n>isolate\n{c}\n')
    out=tmp_path/'out'
    assert cp.main(['-i',str(inp),'-o',str(out),'-t','2'])==0
    summary=read(out/'cluster_summary.tsv')
    assert [int(r['n_proteins']) for r in summary]==[2,2,1]
    assert len(read(out/'membership.tsv'))==5
    stamp=(out/'work/diamond_hits.tsv').stat().st_mtime_ns
    assert cp.main(['-i',str(inp),'-o',str(out),'-t','2','--resume','--inflation','3.0'])==0
    assert (out/'work/diamond_hits.tsv').stat().st_mtime_ns==stamp


def test_omitted_nonisolated_mcl_vertex_fails(tmp_path):
    db,_=normalized(tmp_path,{'a':'MA','b':'MA'})
    db.execute('CREATE TABLE edges (a INTEGER,b INTEGER,weight REAL)');db.execute('INSERT INTO edges VALUES(1,2,20)')
    path=tmp_path/'partial';path.write_text('p0000000001\n')
    with pytest.raises(ValueError,match='omitted'):cp.memberships(db,path)
    db.close()


def test_all_isolates_cli(tmp_path):
    if not shutil.which('mcl'):pytest.skip('MCL binary not available')
    inp=tmp_path/'all.faa';inp.write_text('>a\nMACDEFG\n>b\nMWYVTSR\n')
    hits=tmp_path/'empty_hits.tsv';hits.write_text('')
    out=tmp_path/'out'
    assert cp.main(['-i',str(inp),'-o',str(out),'--hits',str(hits)])==0
    assert [r['n_proteins'] for r in read(out/'cluster_summary.tsv')]==['1','1']
    assert len(list((out/'clusters').glob('*.faa')))==2


def test_domain_negative_cluster_and_no_false_cross_protein_edges(tmp_path):
    inp,hmm=domain_fixture(tmp_path)
    # A and B are present in the cluster, but never in the same protein.
    seqs={r.id:str(r.seq) for r in __import__('Bio.SeqIO',fromlist=['SeqIO']).parse(inp/'family_1.faa','fasta')}
    b=str(next(__import__('Bio.SeqIO',fromlist=['SeqIO']).parse(inp/'family_2.faa','fasta')).seq)
    (inp/'family_1.faa').write_text(f'>two\n{seqs["two"]}\n>separateB\n{b}\n')
    (inp/'negative.faa').write_text('>negative_only\nM'+'P'*90+'\n')
    out=tmp_path/'out'
    assert cd.main(['-i',str(inp),'--hmm',str(hmm),'-o',str(out)])==0
    assert not read(out/'clusters/family_1/cooccurrence_edges.tsv')
    assert not read(out/'clusters/negative/domain_prevalence.tsv')
    assert (out/'clusters/negative/domain_network.svg').exists()
    negative=next(r for r in read(out/'cluster_summary.tsv') if r['cluster']=='negative')
    assert negative['n_proteins']=='1' and negative['fraction_with_domains']=='0.0'
