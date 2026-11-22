#!/usr/bin/env python3
"""Regenerate interactive reports from completed clustering/domain outputs; no searches."""
import argparse
import csv
import json
import logging
from pathlib import Path
import sqlite3
import sys
from common import LOG, dump, sha256, base_manifest, versions
from interactive_network import options, validate, protein_network, domain_network


def refresh_proteins(folder, args, profiles=None):
    dbpath=folder/'work/proteins.sqlite'
    if not dbpath.is_file():raise ValueError(f'Missing saved protein graph: {dbpath}; rerun cluster_proteins.py --resume')
    with sqlite3.connect(dbpath) as db:
        fields={r[1] for r in db.execute('PRAGMA table_info(edges)')}
        if 'identity_percent' not in fields:
            # Older releases retain raw hits and memberships, so do not rerun DIAMOND.
            from argparse import Namespace
            from cluster_proteins import build_graph
            manifest=json.loads((folder/'manifest.json').read_text());params=manifest['parameters']
            hits=Path(params['hits']) if params.get('hits') else folder/'work/diamond_hits.tsv'
            if not hits.is_file():raise ValueError(f'Older graph lacks identity metadata and saved hits are missing: {hits}')
            expected=manifest.get('hits_sha256')
            stamp=folder/'work/diamond.stage.json'
            if not expected and stamp.exists():expected=json.loads(stamp.read_text()).get('sha256',{}).get(hits.name)
            if not expected or sha256(hits)!=expected:raise ValueError('Saved alignment checksum unavailable or mismatched; rerun cluster_proteins.py --resume to validate search inputs')
            build_graph(db,hits,folder/'work/edges.abc',Namespace(**params),internal=not params.get('hits'))
        db.execute('CREATE TABLE IF NOT EXISTS membership(n INTEGER PRIMARY KEY,cluster TEXT)')
        if not db.execute('SELECT COUNT(*) FROM membership').fetchone()[0]:
            with open(folder/'membership.tsv') as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    n=db.execute('SELECT n FROM proteins WHERE original=?',(row['protein_id'],)).fetchone()
                    if not n:raise ValueError('Saved membership does not match protein mapping')
                    db.execute('INSERT INTO membership VALUES (?,?)',(n[0],row['cluster']))
        if db.execute('SELECT COUNT(*) FROM membership').fetchone()[0]!=db.execute('SELECT COUNT(*) FROM proteins').fetchone()[0]:
            raise ValueError('Saved membership is incomplete')
        db.commit()
        return protein_network(db,folder,args,profiles)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--clustering',type=Path,help='Completed cluster_proteins output folder')
    p.add_argument('--domains',type=Path,help='Completed cluster_domains output folder; also annotates protein nodes')
    p.add_argument('--cluster',dest='network_cluster',help='Focus protein preview on one MCL cluster without changing clustering')
    options(p);args=p.parse_args(argv);logging.basicConfig(level=logging.INFO,format='%(levelname)s %(message)s')
    try:
        validate(args)
        if not args.clustering and not args.domains:raise ValueError('Supply --clustering and/or --domains')
        results={};manifest=base_manifest(args);manifest.update(mode='interactive_reports_only',software=versions(['networkx']))
        if args.clustering:
            results['proteins']=refresh_proteins(args.clustering,args,args.domains)
            LOG.info('Protein network: %s',args.clustering/'protein_network.html')
        if args.domains:
            import networkx as nx
            graphs=sorted(args.domains.rglob('domain_network.graphml'))
            pooled=args.domains/'all_domains_network.graphml'
            if pooled.exists():graphs.insert(0,pooled)
            if not graphs:raise ValueError('No saved domain GraphML files found; run cluster_domains.py --resume')
            for path in graphs:
                graph=nx.read_graphml(path);stem=path.stem
                domain_network(graph,path.parent,args,stem=stem,title='Domain co-occurrence: '+str(graph.graph.get('dataset',graph.graph.get('cluster','all proteins'))))
                nx.write_graphml(graph,path)
            results['domain_graphs']=len(graphs)
            LOG.info('Refreshed %d domain networks',len(graphs))
        manifest['results']=results
        for folder in [args.clustering,args.domains]:
            if folder:
                manifest['source_manifest_sha256']=sha256(folder/'manifest.json') if (folder/'manifest.json').exists() else 'unavailable'
                dump(folder/'network_report_manifest.json',manifest)
    except (ValueError,OSError,sqlite3.Error,KeyError) as exc:LOG.error('%s',exc);return 2
    return 0


if __name__=='__main__':sys.exit(main())
