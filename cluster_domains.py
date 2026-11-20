#!/usr/bin/env python3
"""Profile domains and co-occurrence in one protein FASTA or a folder of clusters."""
import argparse
from collections import Counter,defaultdict
import csv
import gzip
import html
import itertools
import inspect
import json
import logging
import math
from pathlib import Path
import sqlite3
import sys
from common import (LOG,normalize,sha256,signature,writer,dump,versions,
                    prepare_output,base_manifest,finish,stage_reusable,mark_stage)

HIT_FIELDS=['protein_id','internal_id','cluster','accession','model_name','model_length','hmm_from','hmm_to','start','end',
            'envelope_start','envelope_end','sequence_bitscore','sequence_evalue','domain_bitscore',
            'independent_evalue','conditional_evalue','model_coverage','protein_coverage','threshold']


def label(value):
    return value.decode() if isinstance(value,bytes) else str(value or '')


def fasta_label(path):
    name=path.name.removesuffix('.gz').rsplit('.',1)[0]
    if name in ['', '.', '..']:raise ValueError(f'Unsafe/empty protein collection name: {path.name}')
    return name


def find_clusters(directory):
    paths=sorted(p for p in directory.iterdir() if p.is_file() and any(p.name.endswith(s) for s in ['.faa','.fa','.fasta','.faa.gz','.fa.gz','.fasta.gz']))
    if not paths:raise ValueError(f'{directory}: no cluster FASTA files found')
    files=[];seen=set()
    for path in paths:
        name=fasta_label(path)
        if name in seen:raise ValueError(f'Duplicate cluster stem: {name}')
        seen.add(name);files.append((name,path))
    return files


def resolve_input(path):
    """A single FASTA is one collection, not a request to cluster its proteins."""
    if path.is_file():return 'proteins',[(fasta_label(path),path)]
    if path.is_dir():return 'clusters',find_clusters(path)
    raise ValueError(f'{path}: input must be an existing protein FASTA file or cluster folder')


def model_inventory(hmm_path,db,evalue):
    import pyhmmer
    db.execute('DROP TABLE IF EXISTS models')
    db.execute('CREATE TABLE models (accession TEXT PRIMARY KEY, name TEXT, length INTEGER, threshold TEXT, mode TEXT)')
    seen=set()
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as f:
        for hmm in f:
            accession=label(hmm.accession) or label(hmm.name)
            base=accession.split('.')[0]
            if base in seen:raise ValueError(f'Duplicate HMM family/version: {accession}')
            seen.add(base)
            ga=hmm.cutoffs.gathering
            if ga is None and evalue is None:
                raise ValueError(f'{accession} lacks gathering thresholds. Supply --domain-evalue explicitly for models lacking GA.')
            mode='GA' if ga is not None else 'Evalue'
            threshold=f'GA sequence={ga[0]} domain={ga[1]}' if ga is not None else f'E-value sequence+independent-domain <= {evalue}'
            try:db.execute('INSERT INTO models VALUES (?,?,?,?,?)',(accession,label(hmm.name),hmm.M,threshold,mode))
            except sqlite3.IntegrityError as e:raise ValueError(f'Duplicate HMM accession/name: {accession}') from e
    count=db.execute('SELECT COUNT(*) FROM models').fetchone()[0]
    if not count:raise ValueError('Empty HMM database')
    db.commit();return count


def scan(db,hmm_path,protein_path,raw_path,args,manifest):
    import pyhmmer
    from contextlib import ExitStack
    db.execute('DROP TABLE IF EXISTS hits')
    db.execute('CREATE TABLE hits (protein_n INTEGER, accession TEXT, name TEXT, model_length INTEGER,hmm_from INTEGER,hmm_to INTEGER,start INTEGER,end INTEGER,env_from INTEGER,env_to INTEGER,seq_score REAL,seq_evalue REAL,dom_score REAL,i_evalue REAL,c_evalue REAL,model_coverage REAL,protein_coverage REAL,threshold TEXT)')
    db.execute('CREATE INDEX hits_protein_idx ON hits(protein_n)')
    alphabet=pyhmmer.easel.Alphabet.amino()
    completed=0;hit_count=0
    with ExitStack() as stack:
        sf=stack.enter_context(pyhmmer.easel.SequenceFile(str(protein_path),digital=True,alphabet=alphabet))
        targets=sf if args.low_memory else sf.read_block()
        hmmfile=stack.enter_context(pyhmmer.plan7.HMMFile(str(hmm_path)))
        raw=stack.enter_context(open(raw_path,'wb'))
        iterator=iter(hmmfile)
        while True:
            batch=list(itertools.islice(iterator,args.model_batch_size))
            if not batch:break
            for mode in ['GA','Evalue']:
                subset=[m for m in batch if (m.cutoffs.gathering is not None)==(mode=='GA')]
                if not subset:continue
                options=dict(bit_cutoffs='gathering') if mode=='GA' else dict(E=args.domain_evalue,domE=args.domain_evalue,incE=args.domain_evalue,incdomE=args.domain_evalue)
                manifest['search_options'][mode]=options
                for model,top in zip(subset,pyhmmer.hmmer.hmmsearch(subset,targets,cpus=args.threads,**options)):
                    accession=label(model.accession) or label(model.name)
                    threshold=db.execute('SELECT threshold FROM models WHERE accession=?',(accession,)).fetchone()[0]
                    top.write(raw,format='domains',header=True)
                    for hit in top:
                        if not hit.included:continue
                        internal=label(hit.name)
                        n=int(internal[1:])
                        for domain in hit.domains:
                            if not domain.included:continue
                            aln=domain.alignment
                            cov=(aln.hmm_to-aln.hmm_from+1)/model.M
                            pcov=(aln.target_to-aln.target_from+1)/hit.length
                            db.execute('INSERT INTO hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                                       (n,accession,label(model.name),model.M,aln.hmm_from,aln.hmm_to,
                                        aln.target_from,aln.target_to,domain.env_from,domain.env_to,
                                        hit.score,hit.evalue,domain.score,domain.i_evalue,domain.c_evalue,cov,pcov,threshold))
                            hit_count+=1
                    completed+=1
            db.commit();LOG.info('Searched %d HMMs; %d threshold-passing domain hits',completed,hit_count)
    manifest['models_searched']=completed
    return hit_count


def hit_rows(db):
    sql='''SELECT p.n,p.internal,p.original,p.cluster,h.accession,h.name,h.model_length,h.hmm_from,h.hmm_to,h.start,h.end,h.env_from,h.env_to,h.seq_score,h.seq_evalue,h.dom_score,h.i_evalue,h.c_evalue,h.model_coverage,h.protein_coverage,h.threshold
           FROM proteins p JOIN hits h ON p.n=h.protein_n ORDER BY p.n,h.start,h.end,h.accession'''
    for row in db.execute(sql):
        n,internal,original,cluster,*values=row
        keys=['accession','model_name','model_length','hmm_from','hmm_to','start','end','envelope_start','envelope_end',
              'sequence_bitscore','sequence_evalue','domain_bitscore','independent_evalue','conditional_evalue','model_coverage','protein_coverage','threshold']
        yield n,dict(protein_id=original,internal_id=internal,cluster=cluster,**dict(zip(keys,values)))


def select_domains(hits,policy='best',overlap_fraction=.5,min_coverage=0.):
    """Greedy overlap heuristic, not the PfamScan clan/nesting algorithm."""
    kept=[];decisions=[]
    for hit in sorted(hits,key=lambda h:(h['independent_evalue'],-h['domain_bitscore'],h['accession'],h['start'],h['end'])):
        reason='accepted'
        if hit['model_coverage']<min_coverage:reason='below_min_model_coverage'
        elif policy=='best':
            for old in kept:
                if hit['accession'].split('.')[0]==old['accession'].split('.')[0]:continue
                overlap=max(0,min(hit['end'],old['end'])-max(hit['start'],old['start'])+1)
                if overlap/min(hit['end']-hit['start']+1,old['end']-old['start']+1)>=overlap_fraction:
                    reason='overlaps_stronger_other_model:'+old['accession'];break
        if reason=='accepted':kept.append(hit)
        decisions.append(dict(**hit,selected=reason=='accepted',selection_reason=reason))
    return sorted(kept,key=lambda h:(h['start'],h['end'],h['accession'])),decisions


def search_cache_keys(files,combined,hmm_sha,pyhmmer_version,fallback):
    """Presentation changes must not invalidate an otherwise identical HMM search."""
    data=dict(memberships=[(c,sha256(p)) for c,p in files],proteins=sha256(combined),
              hmm=hmm_sha,pyhmmer=pyhmmer_version,fallback=fallback)
    scientific_code={f.__name__:inspect.getsource(f) for f in [normalize,model_inventory,scan]}
    key=signature(dict(search_schema=1,scientific_code=scientific_code,**data))
    # Exact source hashes of the published v1.1.0 implementation, whose search
    # and normalization semantics are unchanged. Raw/cache checksums still apply.
    legacy=signature(dict(script='dd77ecf24746035c42de6036932c397c67eb498f3b7df6eab8eb9255770b5d3f',
                          common='ddce5e4ebc029d185dd9e45046d13d8daa1bbd1810890bd7b638c1b40403f681',**data))
    # Do not retain migration eligibility if a future edit changes the search core.
    if signature(scientific_code)!='8280ecc153383d41c0ab05c4ac09baaa658132d78fd45a65a7cdfd04fadcd8a5':
        legacy=None
    return key,legacy


def draw_network(graph,path,max_nodes,max_edges,collection,total):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import networkx as nx
    keep=sorted(graph,key=lambda n:(-graph.nodes[n]['proteins_with_domain'],n))[:max_nodes]
    view=graph.subgraph(keep).copy()
    if view.number_of_edges()>max_edges:
        ranked=sorted(view.edges(data=True),key=lambda e:(-e[2]['proteins_with_both'],e[0],e[1]))
        view.remove_edges_from([(u,v) for u,v,_ in ranked[max_edges:]])
    fig,ax=plt.subplots(figsize=(10,8))
    if not len(view):
        ax.text(.5,.5,'No domains passed the selected thresholds',ha='center',va='center',transform=ax.transAxes)
    else:
        pos=nx.spring_layout(view,seed=17,weight='proteins_with_both') if view.number_of_edges() else nx.circular_layout(view)
        fraction_key='fraction_proteins' if graph.graph.get('input_mode')=='proteins' else 'fraction_cluster_proteins'
        sizes=[250+1600*view.nodes[n][fraction_key] for n in view]
        widths=[.5+4*e[2]['jaccard'] for e in view.edges(data=True)]
        nx.draw_networkx(view,pos,ax=ax,node_size=sizes,width=widths,node_color='#8ac6be',edge_color='#7b8e9a',
                         labels={n:n+'\n'+view.nodes[n]['name'] for n in view},font_size=8)
    ax.set_title(f'{collection}: within-protein domain co-occurrence\n{total} proteins; shown {len(view)}/{len(graph)} domains, {view.number_of_edges()}/{graph.number_of_edges()} edges',fontsize=11)
    ax.axis('off');fig.tight_layout();fig.savefig(path,format='svg',bbox_inches='tight');plt.close(fig)


def summarize(db,out,args,input_mode='clusters'):
    import networkx as nx
    standalone=input_mode=='proteins'
    # Clear convenience exports when switching modes or disabling plots on resume.
    for name in ['cooccurrence_edges.tsv','domain_network.graphml','domain_network.svg','domain_prevalence.svg','cluster_summary.tsv','collection_summary.tsv']:
        (out/name).unlink(missing_ok=True)
    if (out/'clusters').exists():
        for old in (out/'clusters').iterdir():
            if old.is_dir():
                for name in ['domain_prevalence.tsv','cooccurrence_edges.tsv','domain_network.graphml','domain_network.svg']:
                    (old/name).unlink(missing_ok=True)
                try:old.rmdir()
                except OSError:pass
        try:(out/'clusters').rmdir()
        except OSError:pass
    clusters={c:n for c,n in db.execute('SELECT cluster,COUNT(*) FROM proteins GROUP BY cluster ORDER BY cluster')}
    if standalone and len(clusters)!=1:raise ValueError('Single-FASTA mode requires exactly one input collection')
    total=sum(clusters.values())
    # Domain summaries and pair counts are aggregated once, then written per cluster.
    counts=defaultdict(Counter);copies=defaultdict(Counter);pairs=defaultdict(Counter)
    names={a:n for a,n in db.execute('SELECT accession,name FROM models')}
    with_domain=Counter()
    hit_fields=[f for f in HIT_FIELDS if not standalone or f!='cluster']
    af,aw=writer(out/'domain_hits.tsv',hit_fields+['selected','selection_reason'])
    pf,pw=writer(out/'protein_architectures.tsv',([] if standalone else ['cluster'])+['protein_id','n_domain_hits','n_distinct_domains','architecture','domain_set'])
    db.execute('DROP TABLE IF EXISTS architectures')
    db.execute('CREATE TABLE architectures (n INTEGER PRIMARY KEY, hits INTEGER, distinct_n INTEGER, architecture TEXT, domain_set TEXT)')
    try:
        for n,group in itertools.groupby(hit_rows(db),key=lambda row:row[0]):
            hits=[h for _,h in group]
            kept,decisions=select_domains(hits,args.overlap_policy,args.overlap_fraction,args.min_model_coverage)
            aw.writerows(decisions)
            if not kept:continue
            c=kept[0]['cluster'];distinct=sorted({h['accession'] for h in kept})
            with_domain[c]+=1
            counts[c].update(distinct);copies[c].update(h['accession'] for h in kept)
            pairs[c].update(itertools.combinations(distinct,2))
            db.execute('INSERT INTO architectures VALUES (?,?,?,?,?)',(n,len(kept),len(distinct),';'.join(f'{h["accession"]}@{h["start"]}-{h["end"]}' for h in kept),';'.join(distinct)))
        # Include domain-negative proteins in both the denominators and architecture table.
        for original,c,hits_n,distinct_n,arch,dset in db.execute('SELECT p.original,p.cluster,COALESCE(a.hits,0),COALESCE(a.distinct_n,0),COALESCE(a.architecture,""),COALESCE(a.domain_set,"") FROM proteins p LEFT JOIN architectures a ON p.n=a.n ORDER BY p.n'):
            pw.writerow(dict(cluster=c,protein_id=original,n_domain_hits=hits_n,n_distinct_domains=distinct_n,architecture=arch,domain_set=dset))
    finally:af.close();pf.close()
    db.commit()
    global_counts=Counter()
    for counter in counts.values():global_counts.update(counter)
    summary_fields=['dataset','total_proteins','proteins_with_domains','proteins_without_domains','fraction_with_domains','distinct_domains','domain_copies','cooccurrence_pairs','top_domain','top_domain_fraction'] if standalone else ['cluster','n_proteins','proteins_with_domains','fraction_with_domains','distinct_domains','domain_copies','cooccurrence_pairs','top_domain','top_domain_fraction']
    sf,sw=writer(out/('collection_summary.tsv' if standalone else 'cluster_summary.tsv'),summary_fields)
    prevalence_fields=['accession','name','total_proteins','proteins_with_domain','fraction_proteins','percent_proteins','domain_copies','mean_copies_per_protein','mean_copies_when_present'] if standalone else ['cluster','accession','name','cluster_proteins','proteins_with_domain','fraction_cluster_proteins','percent_cluster_proteins','domain_copies','mean_copies_per_member','mean_copies_when_present','outside_cluster_proteins','outside_proteins_with_domain','outside_prevalence','prevalence_ratio_vs_outside']
    fraction_key='fraction_proteins' if standalone else 'fraction_cluster_proteins'
    lf,lw=writer(out/'domain_prevalence.tsv',prevalence_fields)
    collection_rows=[]
    rows_for_html=[]
    try:
        for c,n in clusters.items():
            folder=out if standalone else out/'clusters'/c
            folder.mkdir(parents=True,exist_ok=True)
            graph_metadata=dict(dataset=c,total_proteins=n) if standalone else dict(cluster=c,cluster_proteins=n)
            graph=nx.Graph(**graph_metadata,input_mode=input_mode,definition='Domains co-occur in the same protein; repeats count once per protein for edges')
            rows=[]
            for a,k in sorted(counts[c].items(),key=lambda item:(-item[1],item[0])):
                outside_n=total-n;outside_k=global_counts[a]-k
                outside_p=outside_k/outside_n if outside_n else None
                prevalence=k/n
                row=dict(cluster=c,accession=a,name=names[a],cluster_proteins=n,proteins_with_domain=k,
                         fraction_cluster_proteins=prevalence,percent_cluster_proteins=100*prevalence,domain_copies=copies[c][a],
                         mean_copies_per_member=copies[c][a]/n,mean_copies_when_present=copies[c][a]/k,
                         outside_cluster_proteins=outside_n,outside_proteins_with_domain=outside_k,
                         outside_prevalence=outside_p if outside_p is not None else 'NA',
                         prevalence_ratio_vs_outside=prevalence/outside_p if outside_p else ('Inf' if outside_n else 'NA'))
                if standalone:
                    row=dict(accession=a,name=names[a],total_proteins=n,proteins_with_domain=k,fraction_proteins=prevalence,percent_proteins=100*prevalence,domain_copies=copies[c][a],mean_copies_per_protein=copies[c][a]/n,mean_copies_when_present=copies[c][a]/k)
                rows.append(row);graph.add_node(a,**{key:row[key] for key in ['name','proteins_with_domain',fraction_key,'domain_copies']})
            if not standalone:
                df,dw=writer(folder/'domain_prevalence.tsv',prevalence_fields);dw.writerows(rows);df.close()
            lw.writerows(rows)
            if standalone:collection_rows=rows
            edge_fields=['domain_a','domain_b','proteins_with_both',fraction_key,'jaccard','lift','p_b_given_a','p_a_given_b']
            ef,ew=writer(folder/'cooccurrence_edges.tsv',edge_fields)
            for (a,b),k in sorted(pairs[c].items()):
                na,nb=counts[c][a],counts[c][b]
                row=dict(domain_a=a,domain_b=b,proteins_with_both=k,fraction_cluster_proteins=k/n,jaccard=k/(na+nb-k),
                         lift=k*n/(na*nb),p_b_given_a=k/na,p_a_given_b=k/nb)
                if standalone:row['fraction_proteins']=row.pop('fraction_cluster_proteins')
                ew.writerow(row)
                if k>=args.min_pair_count:graph.add_edge(a,b,**{key:val for key,val in row.items() if key not in ['domain_a','domain_b']})
            ef.close()
            nx.write_graphml(graph,folder/'domain_network.graphml')
            if not args.no_plots:draw_network(graph,folder/'domain_network.svg',args.max_plot_nodes,args.max_plot_edges,c,n)
            top=max(counts[c],key=lambda a:counts[c][a]) if counts[c] else ''
            summary=dict(cluster=c,n_proteins=n,proteins_with_domains=with_domain[c],fraction_with_domains=with_domain[c]/n,
                         distinct_domains=len(counts[c]),domain_copies=sum(copies[c].values()),cooccurrence_pairs=len(pairs[c]),
                         top_domain=top,top_domain_fraction=counts[c][top]/n if top else 0.)
            if standalone:
                summary.pop('cluster');summary.pop('n_proteins')
                summary.update(dataset=c,total_proteins=n,proteins_without_domains=n-with_domain[c])
            sw.writerow(summary);rows_for_html.append(summary)
    finally:sf.close();lf.close()
    if standalone:
        from domain_reports import draw_prevalence,render_collection_report
        if not args.no_plots:draw_prevalence(collection_rows,out/'domain_prevalence.svg')
        render_collection_report(out,rows_for_html[0],collection_rows)
        return dict(proteins=total,proteins_with_domains=sum(with_domain.values()),proteins_without_domains=total-sum(with_domain.values()),distinct_domains=len(global_counts))
    text='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Protein cluster domain profiles</title><style>body{font:15px system-ui;margin:30px;color:#244}td,th{padding:8px;border:1px solid #ccd}table{border-collapse:collapse}input{padding:10px;margin:10px}a{color:#08796c}</style><h1>Protein cluster domain profiles</h1><p>Prevalence uses all cluster proteins, including proteins with no detected domains. Edges mean domains in the same protein; repeats do not multiply edges. These are descriptive associations, not evidence of physical interaction, operons or EET activity.</p><input id="q" placeholder="Search clusters"><table><thead><tr>'''
    fields=list(rows_for_html[0]) if rows_for_html else []
    text+=''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'<th>Files</th></tr></thead><tbody>'
    from urllib.parse import quote
    for row in rows_for_html:
        c=quote(row['cluster'],safe='');text+='<tr>'+''.join('<td>'+html.escape(str(row[k]))+'</td>' for k in fields)
        text+=f'<td><a href="clusters/{c}/domain_prevalence.tsv">Domains</a> | <a href="clusters/{c}/cooccurrence_edges.tsv">Edges</a> | <a href="clusters/{c}/domain_network.graphml">GraphML</a>'
        if not args.no_plots:text+=f' | <a href="clusters/{c}/domain_network.svg">Network SVG</a>'
        text+='</td></tr>'
    text+='''</tbody></table><script>document.getElementById('q').addEventListener('input',e=>document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(e.target.value.toLowerCase())));</script></html>'''
    (out/'index.html').write_text(text)
    return dict(clusters=len(clusters),proteins=total,proteins_with_domains=sum(with_domain.values()),distinct_domains=len(global_counts))


def run(args):
    if not all(math.isfinite(x) for x in [args.min_model_coverage,args.overlap_fraction]+([args.domain_evalue] if args.domain_evalue is not None else [])):raise ValueError('Numeric parameters must be finite')
    if args.threads<1 or args.model_batch_size<1 or not 0<=args.min_model_coverage<=1 or not 0<args.overlap_fraction<=1 or args.min_pair_count<1 or args.max_plot_nodes<1 or args.max_plot_edges<0 or (args.domain_evalue is not None and args.domain_evalue<=0):
        raise ValueError('Invalid threads/batch-size/coverage/overlap/network limits or E-value')
    input_mode,files=resolve_input(args.input)
    if any(args.output.resolve() in path.resolve().parents for _,path in files):
        raise ValueError('Choose an output directory that does not contain the input protein files')
    log=prepare_output(args.output,args.resume,args.force);db=None
    try:
        import pyhmmer
        manifest=base_manifest(args);manifest['software']=versions(['biopython','pyhmmer','networkx','matplotlib'])
        manifest['input_mode']=input_mode
        manifest['protein_collections']=[name for name,_ in files]
        manifest['input_sha256']={str(p.resolve()):sha256(p) for _,p in files}
        manifest['hmm_sha256']=sha256(args.hmm);manifest['search_options']={}
        work=args.output/'work';work.mkdir(exist_ok=True)
        combined=work/'proteins.faa';index=work/'proteins.sqlite'
        # Keep the domain DB separate: normalized mapping can be regenerated safely.
        db,total=normalize(files,combined,index,args.output/'id_mapping.tsv',include_cluster=input_mode=='clusters')
        search_db=work/'search.sqlite';raw=work/'domains.domtblout';stamp=work/'search.stage.json'
        key,legacy_key=search_cache_keys(files,combined,manifest['hmm_sha256'],pyhmmer.__version__,args.domain_evalue)
        reuse=stage_reusable(stamp,key,[search_db,raw],args.resume and not args.force)
        legacy_reuse=not reuse and stage_reusable(stamp,legacy_key,[search_db,raw],args.resume and not args.force)
        reuse=reuse or legacy_reuse
        if reuse:
            LOG.info('Reusing completed HMM search')
            source=sqlite3.connect(search_db)
            source.backup(db)  # Restore mapping plus HMM results; verify membership below.
            source.close()
            if legacy_reuse:
                LOG.info('Migrating compatible v1.1.0 search cache; no HMM rescan needed')
                mark_stage(stamp,key,[search_db,raw])
            # Cluster labels/headers may have changed without changing normalized sequence bytes.
            # Such changes are part of the input signature and invalidate this cache.
        else:
            nmodels=model_inventory(args.hmm,db,args.domain_evalue)
            LOG.info('Searching %d HMMs against %d proteins',nmodels,total)
            dump(args.output/'manifest.inprogress.json',manifest)
            scan(db,args.hmm,combined,raw,args,manifest)
            target=sqlite3.connect(search_db);db.backup(target);target.close()
            mark_stage(stamp,key,[search_db,raw])
        # Persist the search definition even on reuse.
        manifest['search_definition']=dict(algorithm='PyHMMER hmmsearch on input proteins' if input_mode=='proteins' else 'PyHMMER hmmsearch on combined cluster proteins',
                                          ga_preferred=True,fallback_evalue=args.domain_evalue,total_proteins=total,
                                          models=db.execute('SELECT COUNT(*) FROM models').fetchone()[0],cache_reused=reuse,
                                          legacy_cache_migrated=legacy_reuse)
        mf,mw=writer(args.output/'model_inventory.tsv',['accession','name','model_length','threshold','mode'])
        for row in db.execute('SELECT * FROM models ORDER BY accession'):mw.writerow(dict(zip(['accession','name','model_length','threshold','mode'],row)))
        mf.close()
        manifest['counts']=summarize(db,args.output,args,input_mode=input_mode)
        manifest['interpretation']='Prevalence and within-protein domain co-occurrence; no statistical enrichment significance or physical interaction inferred'
        manifest['output_sha256']={p.name:sha256(p) for p in args.output.glob('*.tsv')}
        finish(args.output,manifest);LOG.info('Finished: %s',args.output/'index.html')
    finally:
        if db:db.close()
        LOG.removeHandler(log);log.close()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i','--input','--clusters',dest='input',type=Path,required=True,help='Protein FASTA (.faa/.faa.gz) or folder of cluster FASTAs; --clusters remains an alias')
    p.add_argument('--hmm',type=Path,required=True,help='Pfam-A.hmm or custom HMM database; decompress .gz first')
    p.add_argument('-o','--output',type=Path,required=True,help='Output profile directory')
    p.add_argument('-t','--threads',type=int,default=8,help='HMM search threads')
    p.add_argument('--domain-evalue',type=float,help='Explicit sequence/independent-domain fallback for models lacking GA; GA still preferred')
    p.add_argument('--min-model-coverage',type=float,default=0.,help='Minimum aligned model fraction; default lets GA govern short repeats')
    p.add_argument('--overlap-policy',choices=['best','all'],default='best',help='Best: resolve strongly overlapping different models before summaries')
    p.add_argument('--overlap-fraction',type=float,default=.5,help='Intersection/shorter-hit-length threshold for overlap conflicts')
    p.add_argument('--min-pair-count',type=int,default=1,help='Minimum same-protein pair count in GraphML/plots; edge TSV retains all pairs')
    p.add_argument('--max-plot-nodes',type=int,default=60,help='Plot most prevalent nodes; GraphML/TSVs retain all nodes')
    p.add_argument('--max-plot-edges',type=int,default=300,help='Plot strongest edges; GraphML/TSVs retain all qualifying edges')
    p.add_argument('--model-batch-size',type=int,default=128,help='Models per batch; bounded HMM memory')
    p.add_argument('--low-memory',action='store_true',help='Stream targets from disk instead of holding digital protein sequences in memory')
    p.add_argument('--no-plots',action='store_true',help='Skip SVG rendering (still writes GraphML, edge tables, HTML index)')
    p.add_argument('--resume',action='store_true',help='Reuse checksum-matching HMM searches, regenerate summaries/networks')
    p.add_argument('--force',action='store_true',help='Recompute into existing output')
    args=p.parse_args(argv);logging.basicConfig(level=logging.INFO,format='%(levelname)s %(message)s')
    try:run(args)
    except (ValueError,OSError,sqlite3.Error,ImportError) as e:LOG.error('%s. Install dependencies as described in README.md.',e);return 2
    return 0

if __name__=='__main__':sys.exit(main())
