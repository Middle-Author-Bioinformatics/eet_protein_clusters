#!/usr/bin/env python3
"""DIAMOND all-versus-all similarity -> MCL -> a FASTA for every cluster."""
import argparse
import csv
import json
import logging
import math
import os
import shutil
import tempfile
from interactive_network import options as network_options, validate as validate_network, protein_network
from pathlib import Path
import sqlite3
import sys
from common import (LOG,normalize,resolve,sha256,signature,writer,dump,versions,
                    tool_version,command,prepare_output,base_manifest,finish,
                    stage_reusable,mark_stage)

FIELDS=['qseqid','sseqid','pident','length','qstart','qend','sstart','send','evalue','bitscore','qlen','slen']


def build_graph(db,hits,output,args,internal=False):
    """Deduplicate a symmetric graph on disk, max weight per unordered pair."""
    db.execute('DROP TABLE IF EXISTS edges')
    db.execute('CREATE TABLE edges (a INTEGER,b INTEGER,weight REAL,identity_percent REAL,coverage_a REAL,coverage_b REAL,evalue REAL,bitscore REAL,PRIMARY KEY(a,b)) WITHOUT ROWID')
    stats=dict(hit_rows=0,self_hits=0,filtered_rows=0,qualifying_rows=0)
    # A small ID cache avoids a SQL lookup per repeated hit without O(N^2) memory.
    from functools import lru_cache
    lookup=lru_cache(maxsize=100000)(lambda s:resolve(db,s,internal))
    with open(hits) as f:
        for line_num,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'):continue
            c=line.split()
            if len(c)!=len(FIELDS):raise ValueError(f'{hits}:{line_num}: expected 12 fields: {" ".join(FIELDS)}')
            try:
                qi,qlen=lookup(c[0]);si,slen=lookup(c[1])
                identity=float(c[2]);length=int(c[3]);qs,qe,ss,se=map(int,c[4:8])
                ev,bit=float(c[8]),float(c[9]);reported_qlen,reported_slen=int(c[10]),int(c[11])
            except (ValueError,OverflowError) as e:raise ValueError(f'{hits}:{line_num}: {e}') from e
            if not all(math.isfinite(x) for x in [identity,ev,bit]) or not 0<=identity<=100 or ev<0 or bit<=0 or length<1:
                raise ValueError(f'{hits}:{line_num}: invalid numeric alignment fields')
            if reported_qlen!=qlen or reported_slen!=slen or not (1<=qs<=qe<=qlen and 1<=ss<=se<=slen):
                raise ValueError(f'{hits}:{line_num}: alignment coordinates/lengths disagree with FASTA')
            stats['hit_rows']+=1
            if qi==si:stats['self_hits']+=1;continue
            qcov=(qe-qs+1)/qlen;scov=(se-ss+1)/slen
            if ev>args.evalue or identity<args.min_identity or min(qcov,scov)<args.min_coverage:
                stats['filtered_rows']+=1;continue
            weight=min(args.evalue_cap,-math.log10(max(ev,10**(-args.evalue_cap)))) if args.weight=='neglog10-evalue' else bit
            if weight<=0:stats['filtered_rows']+=1;continue
            a,b=sorted((qi,si))
            ca,cb=(qcov,scov) if qi==a else (scov,qcov)
            db.execute('INSERT INTO edges VALUES(?,?,?,?,?,?,?,?) ON CONFLICT(a,b) DO UPDATE SET weight=excluded.weight,identity_percent=excluded.identity_percent,coverage_a=excluded.coverage_a,coverage_b=excluded.coverage_b,evalue=excluded.evalue,bitscore=excluded.bitscore WHERE excluded.weight>edges.weight OR (excluded.weight=edges.weight AND excluded.bitscore>edges.bitscore)',(a,b,weight,identity,ca,cb,ev,bit))
            stats['qualifying_rows']+=1
            if stats['hit_rows']%1000000==0:
                db.commit();LOG.info('Processed %d alignment rows',stats['hit_rows'])
    db.commit()
    stats['undirected_edges']=db.execute('SELECT COUNT(*) FROM edges').fetchone()[0]
    with open(output,'w') as out:
        for a,b,w in db.execute('SELECT a,b,weight FROM edges ORDER BY a,b'):
            out.write(f'p{a:010d}\tp{b:010d}\t{w:.8g}\n')
    return stats


def memberships(db,path):
    """Read MCL clusters, restore isolates and verify a complete partition."""
    groups=[];seen=set()
    if path.exists():
        for line in path.read_text().splitlines():
            if not line.strip():continue
            members=[]
            for label in line.split():
                row=resolve(db,label,True);n=row[0]
                if n in seen:raise ValueError(f'MCL output repeats protein {label}')
                seen.add(n);members.append(n)
            groups.append(sorted(members))
    active={n for (n,) in db.execute('SELECT a FROM edges UNION SELECT b FROM edges')}
    if active-seen:raise ValueError('MCL output omitted proteins with qualifying graph edges; inspect work/mcl.log')
    for (n,) in db.execute('SELECT n FROM proteins ORDER BY n'):
        if n not in seen:groups.append([n])
    return sorted(groups,key=lambda xs:(-len(xs),xs[0]))


def write_clusters(db,normalized_fasta,groups,out):
    from Bio import SeqIO
    db.execute('DROP TABLE IF EXISTS membership')
    db.execute('CREATE TABLE membership(n INTEGER PRIMARY KEY,cluster TEXT)')
    # SeqIO.index uses file offsets, not the protein sequences in RAM.
    index=SeqIO.index(str(normalized_fasta),'fasta')
    cluster_dir=out/'clusters';cluster_dir.mkdir(exist_ok=True)
    for old in cluster_dir.glob('cluster_*.faa'):old.unlink()
    sf,sw=writer(out/'cluster_summary.tsv',['cluster','n_proteins','min_length','mean_length','max_length','singleton','fasta'])
    mf,mw=writer(out/'membership.tsv',['cluster','protein_id','internal_id','length_aa'])
    try:
        for i,members in enumerate(groups,1):
            name=f'cluster_{i:06d}';lengths=[]
            with open(cluster_dir/(name+'.faa'),'w') as fasta:
                for n in members:
                    db.execute('INSERT INTO membership VALUES (?,?)',(n,name))
                    internal,original,header,length=db.execute('SELECT internal,original,header,length FROM proteins WHERE n=?',(n,)).fetchone()
                    fasta.write(f'>{header}\n{str(index[internal].seq)}\n');lengths.append(length)
                    mw.writerow(dict(cluster=name,protein_id=original,internal_id=internal,length_aa=length))
            sw.writerow(dict(cluster=name,n_proteins=len(members),min_length=min(lengths),mean_length=round(sum(lengths)/len(lengths),2),max_length=max(lengths),singleton=len(members)==1,fasta=f'clusters/{name}.faa'))
    finally:sf.close();mf.close();index.close()
    db.commit()



def storage_info(path):
    path=Path(path).resolve()
    usage=shutil.disk_usage(path)
    info=dict(path=str(path),total_bytes=usage.total,free_bytes=usage.free)
    if hasattr(os,'statvfs'):
        fs=os.statvfs(path)
        info.update(free_inodes=fs.f_favail,total_inodes=fs.f_files)
    return info


def diamond_search(args,normalized,work,hits,manifest):
    """Use an isolated scratch directory, cleaning only files owned by this run."""
    requested=getattr(args,'tmpdir',None)
    parent=Path(requested or os.environ.get('TMPDIR') or work).expanduser().resolve()
    parent.mkdir(parents=True,exist_ok=True)
    manifest['storage']=dict(output=storage_info(work),scratch_parent=storage_info(parent),
                             scratch_source='--tmpdir' if requested else ('TMPDIR' if os.environ.get('TMPDIR') else 'output/work'))
    for label in ['output','scratch_parent']:
        item=manifest['storage'][label]
        LOG.info('%s: %s; %.1f GiB free',label,item['path'],item['free_bytes']/1024**3)
    with tempfile.TemporaryDirectory(prefix='eet-diamond-',dir=parent) as scratch:
        manifest['storage'].update(scratch_run_directory=scratch,scratch_cleanup='pending')
        dump(args.output/'manifest.inprogress.json',manifest)
        try:
            command(['diamond','makedb','--in',normalized,'--db',work/'proteins.dmnd','--threads',args.threads],work/'diamond_makedb.log',manifest)
            command(['diamond','blastp','--db',work/'proteins.dmnd','--query',normalized,'--out',hits,
                     '--threads',args.threads,'--'+args.sensitivity,'--evalue',args.evalue,
                     '--max-target-seqs','0','--max-hsps','1','--block-size',args.block_size,
                     '--tmpdir',scratch,'--outfmt','6',*FIELDS],work/'diamond_blastp.log',manifest)
        except (ValueError,OSError) as exc:
            detail=str(exc).lower()
            if any(term in detail for term in ['no space left','disk quota exceeded']):
                raise ValueError(f'DIAMOND exhausted disk space, inodes or quota. Scratch parent: {parent}; output/work: {work.resolve()}. '
                                 'Check df -h and df -i on BOTH paths (and your quota). Use --tmpdir on a larger filesystem and put --output there too if necessary. '
                                 'A smaller --block-size (e.g. 0.4) reduces temporary space but can change alignments slightly. '
                                 'Rerun with --resume; this failed DIAMOND search restarts, while completed valid searches can be reused. '
                                 f'Original error: {exc}') from exc
            raise
        finally:
            manifest['storage']['scratch_cleanup']='managed_by_TemporaryDirectory'
    manifest['storage']['scratch_cleanup']='removed'


def run(args):
    validate_network(args)
    if not all(math.isfinite(x) for x in [args.inflation,args.evalue,args.min_identity,args.min_coverage,args.evalue_cap,args.block_size]):raise ValueError('Numeric parameters must be finite')
    if args.threads<1 or args.inflation<=1 or not 0<args.evalue<1 or not 0<=args.min_identity<=100 or not 0<=args.min_coverage<=1 or not 1<=args.evalue_cap<=300 or args.block_size<=0:
        raise ValueError('threads >=1, inflation >1, 0<evalue<1, identity 0..100, coverage 0..1, E-value cap 1..300, block-size>0 required')
    log=prepare_output(args.output,args.resume,args.force)
    db=None
    manifest=base_manifest(args)
    try:
        manifest['input_sha256']=sha256(args.input)
        manifest['software']=versions(['biopython'])
        manifest['software']['mcl']=tool_version('mcl')
        if not args.hits:manifest['software']['mcxload']=tool_version('mcxload')
        dump(args.output/'manifest.inprogress.json',manifest)
        work=args.output/'work';work.mkdir(exist_ok=True)
        normalized=work/'proteins.faa'
        db,total=normalize([('',args.input)],normalized,work/'proteins.sqlite',args.output/'id_mapping.tsv')
        if args.hits:
            hits=args.hits;manifest['hits_sha256']=sha256(hits)
        else:
            manifest['software']['diamond']=tool_version('diamond','version')
            hits=work/'diamond_hits.tsv'
            key=signature(dict(input=sha256(normalized),diamond=manifest['software']['diamond'],evalue=args.evalue,
                               sensitivity=args.sensitivity,threads=args.threads,block_size=args.block_size,max_targets=0,max_hsps=1))
            stamp=work/'diamond.stage.json'
            if stage_reusable(stamp,key,[hits],args.resume and not args.force):
                LOG.info('Reusing completed DIAMOND search')
                manifest['diamond_search_cache_reused']=True
                manifest['storage']=dict(output=storage_info(work),scratch_not_used='completed search reused')
            else:
                manifest['diamond_search_cache_reused']=False
                stamp.unlink(missing_ok=True)
                diamond_search(args,normalized,work,hits,manifest)
                mark_stage(stamp,key,[hits])
        edges=work/'edges.abc'
        stats=build_graph(db,hits,edges,args,internal=not args.hits)
        LOG.info('%d proteins; %d qualifying undirected edges',total,stats['undirected_edges'])
        mcl_output=work/'mcl_clusters.tsv'
        if stats['undirected_edges']:
            # Always use the native sparse matrix path for large networks.
            manifest['software']['mcxload']=tool_version('mcxload')
            matrix=work/'graph.mci';tab=work/'graph.tab'
            command(['mcxload','-abc',edges,'--stream-mirror','-re','max','-o',matrix,'-write-tab',tab],work/'mcxload.log',manifest)
            command(['mcl',matrix,'-I',args.inflation,'-te',args.threads,'-use-tab',tab,'-o',mcl_output],work/'mcl.log',manifest)
        else:
            mcl_output.write_text('');manifest['warnings'].append('No qualifying similarity edges; all proteins retained as singleton clusters.')
        groups=memberships(db,mcl_output)
        write_clusters(db,normalized,groups,args.output)
        if not args.no_interactive:
            manifest['interactive_network']=protein_network(db,args.output,args,args.domain_profiles)
        else:
            for name in ['protein_network.html','protein_network_nodes.tsv','protein_network_edges.tsv','index.html']:(args.output/name).unlink(missing_ok=True)
        manifest['counts']=dict(proteins=total,clusters=len(groups),singletons=sum(len(g)==1 for g in groups),**stats)
        manifest['graph_definition']='Undirected max over qualifying directional hits; self-hits excluded; isolates restored as singletons'
        manifest['output_sha256']={name:sha256(args.output/name) for name in ['membership.tsv','cluster_summary.tsv','id_mapping.tsv']}
        finish(args.output,manifest)
        LOG.info('Wrote %d cluster FASTAs to %s',len(groups),args.output/'clusters')
    except Exception as exc:
        scratch=manifest.get('storage',{}).get('scratch_run_directory')
        if scratch: manifest['storage']['scratch_cleanup']='removed' if not Path(scratch).exists() else 'incomplete'
        manifest['error']=str(exc)
        try: dump(args.output/'manifest.inprogress.json',manifest)
        except OSError as report_error: LOG.error('Could not save failure manifest: %s',report_error)
        raise
    finally:
        if db:db.close()
        LOG.removeHandler(log);log.close()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i','--input',type=Path,required=True,help='Protein .faa or .faa.gz; globally unique first-token IDs')
    p.add_argument('-o','--output',type=Path,required=True,help='Output directory; cluster FASTAs in clusters/')
    p.add_argument('-t','--threads',type=int,default=8,help='DIAMOND and MCL threads')
    p.add_argument('--inflation',type=float,default=2.,help='MCL granularity; higher usually produces finer clusters')
    p.add_argument('--evalue',type=float,default=1e-5,help='Maximum similarity-search E-value')
    p.add_argument('--min-identity',type=float,default=30.,help='Minimum percent identity of best HSP')
    p.add_argument('--min-coverage',type=float,default=.5,help='Minimum fraction covered on BOTH proteins, from aligned coordinate spans')
    p.add_argument('--sensitivity',choices=['sensitive','very-sensitive','ultra-sensitive'],default='very-sensitive',help='DIAMOND mode')
    p.add_argument('--weight',choices=['neglog10-evalue','bitscore'],default='neglog10-evalue',help='MCL edge weight')
    p.add_argument('--evalue-cap',type=float,default=200.,help='Cap -log10(E); handles zero E-values')
    p.add_argument('--block-size',type=float,default=2.,help='DIAMOND block size in billions of letters; lower reduces memory')
    p.add_argument('--tmpdir','--diamond-tmpdir',dest='tmpdir',type=Path,help='DIAMOND scratch parent; default $TMPDIR or OUTPUT/work. Final results still use --output')
    network_options(p)
    p.add_argument('--network-cluster',help='Focus the interactive protein view on one MCL cluster; clustering itself is unchanged')
    p.add_argument('--domain-profiles',type=Path,help='Optional cluster_domains output folder; annotate protein nodes with domain names')
    p.add_argument('--hits',type=Path,help='Reuse custom 12-column similarity TSV with ORIGINAL IDs; see README')
    p.add_argument('--resume',action='store_true',help='Reuse checksum-matching DIAMOND output; regenerate graph and clusters')
    p.add_argument('--force',action='store_true',help='Recompute into an existing output directory')
    args=p.parse_args(argv)
    logging.basicConfig(level=logging.INFO,format='%(levelname)s %(message)s')
    try:run(args)
    except (ValueError,OSError,sqlite3.Error,ImportError) as e:LOG.error('%s',e);return 2
    return 0

if __name__=='__main__':sys.exit(main())
