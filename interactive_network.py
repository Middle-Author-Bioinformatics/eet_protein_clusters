"""Offline, bounded interactive network exports with explicit scientific labels."""
import html
import json
import math
import os
import re
from pathlib import Path
from collections import Counter, defaultdict
from common import writer

DEFAULT_NODES = 1000
DEFAULT_EDGES = 5000


def options(parser):
    parser.add_argument('--no-interactive', action='store_true', help='Skip interactive HTML (static plots and machine-readable analysis still available)')
    parser.add_argument('--network-max-nodes', type=int, default=DEFAULT_NODES,
                        help='Maximum nodes embedded per interactive HTML; complete analysis is retained')
    parser.add_argument('--network-max-edges', type=int, default=DEFAULT_EDGES,
                        help='Maximum strongest edges embedded per interactive HTML')
    parser.add_argument('--network-max-proteins', type=int, default=20000,
                        help='Maximum unique domain-carrying proteins embedded per domain HTML (default 20000; 0 disables lists)')


def validate(args):
    if args.network_max_nodes<1 or args.network_max_edges<0:
        raise ValueError('network-max-nodes must be >=1 and network-max-edges >=0')
    if getattr(args,'network_max_proteins',20000)<0:
        raise ValueError('network-max-proteins must be >=0')


def write_html(path, nodes, edges, title, kind, definition, total_nodes=None, total_edges=None, selection='All nodes', downloads=(), members=None):
    """All data, layout and JavaScript are embedded; file:// works without a CDN."""
    groups=defaultdict(list)
    for node in nodes: groups[str(node.get('group', 'ungrouped'))].append(node)
    for gi, (group, group_nodes) in enumerate(sorted(groups.items())):
        angle=gi*2.399963229728653;radius=math.sqrt(gi)*200
        cx,cy=radius*math.cos(angle),radius*math.sin(angle)
        for i,node in enumerate(group_nodes):
            a=i*2.399963229728653;r=math.sqrt(i)*28
            node.update(x=cx+r*math.cos(a),y=cy+r*math.sin(a))
    payload=dict(kind=kind,title=title,nodes=nodes,edges=edges,definition=definition,selection=selection,
                 total_nodes=len(nodes) if total_nodes is None else total_nodes,
                 total_edges=len(edges) if total_edges is None else total_edges)
    if kind=='domains':
        payload['members']=members or dict(available=False,reason='Protein IDs unavailable: aggregate graph data do not contain carriers. Regenerate with the original protein_architectures.tsv.')
    encoded=json.dumps(payload,ensure_ascii=False,allow_nan=False).replace('&','\\u0026').replace('<','\\u003c').replace('>','\\u003e').replace('\u2028','\\u2028').replace('\u2029','\\u2029')
    template=Path(__file__).with_name('network_viewer.html').read_text()
    script=Path(__file__).with_name('network_viewer.js').read_text()
    links=' | '.join(f'<a href="{html.escape(str(url),quote=True)}">{html.escape(str(label))}</a>' for label,url in downloads)
    parts=dict(TITLE=html.escape(title),DEFINITION=html.escape(definition),DOWNLOADS=links,DATA=encoded,SCRIPT=script)
    result=re.sub(r'@@(TITLE|DEFINITION|DOWNLOADS|DATA|SCRIPT)@@',lambda match:parts[match[1]],template)
    path.write_text(result)
    return dict(nodes=len(nodes),edges=len(edges),total_nodes=payload['total_nodes'],total_edges=payload['total_edges'],selection=selection)


def domain_communities(graph):
    """Weighted Louvain on the complete exported graph; isolates stay distinct."""
    import networkx as nx
    if not graph.number_of_edges(): communities=[{n} for n in graph]
    else: communities=nx.community.louvain_communities(graph,weight='jaccard',seed=17)
    ordered=sorted(communities,key=lambda c:(-len(c),tuple(sorted(c))))
    for i, members in enumerate(ordered,1):
        for node in members: graph.nodes[node]['community']=f'domain_community_{i:04d}'
    graph.graph['community_method']='Louvain; Jaccard edge weights; resolution=1; seed=17; complete exported graph'


def domain_network(graph, folder, args, stem='domain_network', title='Domain co-occurrence', members=None):
    domain_communities(graph)
    f,w=writer(folder/(stem+'_communities.tsv'),['accession','name','community','proteins_with_domain','domain_copies'])
    for accession, data in sorted(graph.nodes(data=True)):
        w.writerow(dict(data,accession=accession))
    f.close()
    if args.no_interactive:
        (folder/(stem+'.html')).unlink(missing_ok=True)
        return dict(skipped=True)
    chosen=sorted(graph,key=lambda n:(-graph.nodes[n]['proteins_with_domain'],n))[:args.network_max_nodes]
    nodes=[]
    for accession in chosen:
        d=dict(graph.nodes[accession]);name=d.get('name') or accession
        nodes.append(dict(d,id=accession,label=f'{name} ({accession})',accession=accession,group=d['community'],
                          size=6+min(18,3*math.log1p(d['proteins_with_domain']))))
    chosen=set(chosen)
    edges=[]
    for a,b,d in sorted(graph.edges(data=True),key=lambda e:(-e[2]['proteins_with_both'],e[0],e[1])):
        if a in chosen and b in chosen:
            edges.append(dict(d,source=a,target=b,value=d['jaccard'],
                              domain_a_name=graph.nodes[a]['name'],domain_b_name=graph.nodes[b]['name']))
            if len(edges)>=args.network_max_edges: break
    if args.network_max_edges==0:edges=[]
    links=[('Domain names and communities',stem+'_communities.tsv'),('GraphML',stem+'.graphml')]
    if members and members.available:
        links.append(('Complete protein architectures',os.path.relpath(members.source,folder)))
    clustering=getattr(args,'clustering',None)
    if clustering:
        links.append(('Protein identity / MCL network',os.path.relpath(Path(clustering)/'protein_network.html',folder)))
    carriers=members.payload(graph,chosen,getattr(args,'network_max_proteins',20000)) if members else None
    return write_html(folder/(stem+'.html'),nodes,edges,title,'domains',
                      'An edge means two different domains occur within the same protein. Repeats count once per protein. Node color: weighted Louvain community on the complete exported co-occurrence graph; not a physical interaction or an EET activity call.',
                      len(graph),graph.number_of_edges(),'Most prevalent domains; strongest shared-protein edges. Community assignments use the complete exported graph.',
                      links,members=carriers)


def protein_network(db, out, args, profiles=None):
    """Bounded HTML preview of the exact MCL input graph, retaining HSP identity."""
    wanted=getattr(args,'network_cluster',None)
    groups=defaultdict(list)
    for n, cluster in db.execute('SELECT p.n,m.cluster FROM proteins p JOIN membership m ON p.n=m.n ORDER BY p.n'):
        if not wanted or cluster==wanted:groups[cluster].append(n)
    if wanted and not groups: raise ValueError(f'Unknown --network-cluster {wanted}; see membership.tsv')
    ordered=sorted(groups,key=lambda g:(-len(groups[g]),g))[:args.network_max_nodes]
    selected=[];i=0
    while len(selected)<args.network_max_nodes:
        added=False
        for group in ordered:
            if i<len(groups[group]):selected.append(groups[group][i]);added=True
            if len(selected)>=args.network_max_nodes:break
        if not added:break
        i+=1
    db.execute('DROP TABLE IF EXISTS temp.network_selected')
    db.execute('CREATE TEMP TABLE network_selected(n INTEGER PRIMARY KEY)')
    db.executemany('INSERT INTO network_selected VALUES (?)',((n,) for n in selected))
    annotations={}
    chosen_ids={r[0] for r in db.execute('SELECT original FROM proteins WHERE n IN (SELECT n FROM network_selected)')}
    if profiles:
        import csv
        inventory=Path(profiles)/'model_inventory.tsv'
        name_map={}
        if inventory.exists():
            with open(inventory) as source:name_map={r['accession']:r['name'] for r in csv.DictReader(source,delimiter='\t')}
        with open(Path(profiles)/'protein_architectures.tsv') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['protein_id'] not in chosen_ids: continue
                annotations[r['protein_id']]=dict(domain_names=r.get('domain_names',';'.join(name_map.get(a,a) for a in r.get('domain_set','').split(';') if a)),
                                                 architecture_names=r.get('architecture_names',';'.join(f'{name_map.get(part.rsplit("@",1)[0],part.rsplit("@",1)[0])} ({part})' for part in r.get('architecture','').split(';') if part)))
    nodes=[]
    for n,internal,pid,header,length,cluster in db.execute('SELECT p.n,p.internal,p.original,p.header,p.length,m.cluster FROM proteins p JOIN membership m ON m.n=p.n WHERE p.n IN (SELECT n FROM network_selected) ORDER BY p.n'):
        nodes.append(dict(id=str(n),label=pid,protein_id=pid,original_header=header,length_aa=length,
                          group=cluster,cluster=cluster,size=6,**annotations.get(pid,{})))
    rows=db.execute('SELECT a,b,weight,identity_percent,coverage_a,coverage_b,evalue,bitscore FROM edges WHERE a IN (SELECT n FROM network_selected) AND b IN (SELECT n FROM network_selected) ORDER BY weight DESC,a,b LIMIT ?', (args.network_max_edges,))
    edges=[dict(source=str(a),target=str(b),weight=w,value=identity,identity_percent=identity,
                coverage_a=ca,coverage_b=cb,evalue=ev,bitscore=bits) for a,b,w,identity,ca,cb,ev,bits in rows]
    total_nodes=db.execute('SELECT COUNT(*) FROM proteins').fetchone()[0]
    total_edges=db.execute('SELECT COUNT(*) FROM edges').fetchone()[0]
    selection='Largest MCL clusters first; round-robin sampling of input-order members; strongest eligible edges. Filters change the view, not MCL membership.'
    if wanted: selection=f'Focused on {wanted}. '+selection
    links=[('Full membership','membership.tsv'),('Cluster summary','cluster_summary.tsv')]
    if profiles:
        domain_html=Path(profiles)/'all_domains_network.html'
        if not domain_html.exists():domain_html=Path(profiles)/'domain_network.html'
        links.append(('Domain co-occurrence network',os.path.relpath(domain_html,out)))
    summary=write_html(out/'protein_network.html',nodes,edges,'Protein amino-acid similarity and MCL clusters','proteins',
                       'Nodes are proteins; color is the existing MCL cluster. Edge identity is DIAMOND pairwise amino-acid identity for the HSP selected by the MCL graph weight. This is not genome-wide average amino-acid identity (AAI).',
                       total_nodes,total_edges,selection,links)
    f,w=writer(out/'protein_network_nodes.tsv',['id','protein_id','original_header','length_aa','cluster','domain_names','architecture_names']);w.writerows(nodes);f.close()
    f,w=writer(out/'protein_network_edges.tsv',['source','target','identity_percent','coverage_a','coverage_b','evalue','bitscore','weight']);w.writerows(edges);f.close()
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Protein clusters</title><h1>Protein clusters</h1><p><a href="protein_network.html">Interactive amino-acid similarity network</a></p><p><a href="membership.tsv">All protein memberships</a> | <a href="cluster_summary.tsv">All cluster sizes</a></p>')
    return summary
