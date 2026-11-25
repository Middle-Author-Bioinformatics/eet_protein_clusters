#!/usr/bin/env python3
"""CDD domain profiles and networks for a protein FASTA or folder of clusters."""
import argparse
import csv
import json
import logging
import math
from pathlib import Path
import re
import shutil
import sqlite3
import sys
from common import (LOG,normalize,sha256,signature,writer,dump,versions,prepare_output,
                    base_manifest,finish,stage_reusable,mark_stage,command,tool_version)
from cluster_domains import resolve_input,summarize
from download_cdd import checked_snapshot
from interactive_network import options,validate

RAW_FORMAT='6 qseqid sseqid qstart qend sstart send evalue bitscore qlen slen'
HIT_FIELDS=['protein_id','internal_id','cluster','accession','model_name','pssm_id','reported_pssm_id','source_database',
            'hit_type','start','end','model_start','model_end','model_length','bitscore','evalue',
            'model_coverage','protein_coverage','incomplete','superfamily_pssm_id',
            'specific_bitscore_cutoff','threshold','coverage_status']
HIT_TYPES={'Specific','Non-Specific','Non-specific','Superfamily','Multidom'}


def source_database(accession):
    for prefix,name in [('pfam','Pfam_in_CDD'),('cd','NCBI_CDD'),('cl','CDD_superfamily'),
                        ('COG','COG'),('KOG','KOG'),('smart','SMART'),('PRK','PRK'),
                        ('TIGR','TIGRFAM'),('NF','NCBIFAM')]:
        if accession.startswith(prefix):return name
    return 'other_CDD_source'


def metadata(db,data_dir,evalue):
    db.execute('CREATE TABLE cdd_models (pssm_id TEXT PRIMARY KEY,accession TEXT UNIQUE,name TEXT,description TEXT,length INTEGER,cutoff REAL)')
    cutoffs={}
    with open(data_dir/'bitscore_specific.txt') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            c=line.rstrip('\n').split('\t')
            if len(c)!=3:raise ValueError('Malformed bitscore_specific.txt')
            value=float(c[2])
            if not c[0].isdigit() or not math.isfinite(value):raise ValueError('Invalid specific bit-score threshold')
            cutoffs[c[0]]=(c[1],value)
    with open(data_dir/'cddid.tbl') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            c=line.rstrip('\n').split('\t')
            if len(c)!=5 or not c[0].isdigit() or not c[4].isdigit() or (int(c[4])==0 and not c[1].startswith('cl')):
                raise ValueError('Malformed cddid.tbl; expected PSSM ID, accession, name, description, length')
            cutoff=cutoffs.get(c[0])
            if cutoff and cutoff[0]!=c[1]:raise ValueError('Specific thresholds disagree with inventory; use a consistent CDD snapshot')
            db.execute('INSERT INTO cdd_models VALUES (?,?,?,?,?,?)',(*c[:4],int(c[4]),cutoff[1] if cutoff else None))
    db.execute('CREATE TABLE models (accession TEXT PRIMARY KEY,name TEXT,length INTEGER,threshold TEXT,mode TEXT)')
    for acc,name,length,cutoff in db.execute('SELECT accession,name,length,cutoff FROM cdd_models'):
        db.execute('INSERT INTO models VALUES (?,?,?,?,?)',(acc,name,length,f'RPS-BLAST E-value <= {evalue}; specific cutoff={cutoff}','CDD_PSSM'))
    db.commit()


def read_alignments(db,path):
    db.execute('CREATE TABLE rps_alignments (protein_n INTEGER,pssm_id TEXT,start INTEGER,end INTEGER,model_start INTEGER,model_end INTEGER,evalue REAL,bitscore REAL,model_length INTEGER)')
    db.execute('CREATE INDEX rps_lookup ON rps_alignments(protein_n,pssm_id,start,end)')
    with open(path) as stream:
        for line_no,line in enumerate(stream,1):
            if not line.strip() or line.startswith('#'):continue
            c=line.rstrip('\r\n').split('\t')
            if len(c)!=10:raise ValueError(f'{path}:{line_no}: expected ten columns ({RAW_FORMAT})')
            protein=db.execute('SELECT n,length FROM proteins WHERE internal=?',(c[0],)).fetchone()
            match=re.fullmatch(r'(?:gnl\|CDD\||CDD:)?(\d+)',c[1])
            if not protein or not match:raise ValueError(f'{path}:{line_no}: unknown query or invalid CDD subject ID')
            qs,qe,ss,se=map(int,c[2:6]);ev,score=map(float,c[6:8]);ql,sl=map(int,c[8:10])
            if not (ql==protein[1] and 1<=qs<=qe<=ql and 1<=ss<=se<=sl and math.isfinite(ev) and ev>=0 and math.isfinite(score) and score>0):
                raise ValueError(f'{path}:{line_no}: invalid coordinates, lengths, or scores')
            model=db.execute('SELECT length FROM cdd_models WHERE pssm_id=?',(match[1],)).fetchone()
            if not model or model[0]!=sl:raise ValueError('BLAST model ID/length differs from metadata; use the matching CDD snapshot')
            db.execute('INSERT INTO rps_alignments VALUES (?,?,?,?,?,?,?,?,?)',(protein[0],match[1],qs,qe,ss,se,ev,score,sl))
    db.commit()


def parse_native(db,path,evalue):
    """Use native QUERY defline IDs, never guesses from ordinal Query_N labels."""
    db.execute('CREATE TABLE cdd_hits (protein_n INTEGER,pssm_id TEXT,start INTEGER,end INTEGER,evalue REAL,bitscore REAL,hit_type TEXT,incomplete TEXT,superfamily_id TEXT,model_start INTEGER,model_end INTEGER,reported_pssm_id TEXT)')
    data=False;ended=False;session=None;query=None;domains=False;seen=set();count=0
    with open(path) as stream:
        for line_no,line in enumerate(stream,1):
            if not line.strip() or line.startswith('#'):continue
            c=line.rstrip('\r\n').split('\t');token=c[0]
            def error(message):raise ValueError(f'{path}:{line_no}: {message}')
            if token=='DATA':
                if data or ended:error('duplicate DATA')
                data=True
            elif token=='ENDDATA':
                if not data or query or session or domains:error('unclosed block at ENDDATA')
                data=False;ended=True
            elif not data:error('expected DATA/ENDDATA around native output')
            elif token=='SESSION':
                if session or len(c)<3 or c[2].upper().split()[0] not in ['RPSBLAST','BLASTP']:error('expected protein RPSBLAST session')
                session=c[1]
            elif token=='ENDSESSION':
                if not session or query or len(c)<2 or c[1]!=session:error('invalid ENDSESSION')
                session=None
            elif token=='QUERY':
                if not session or query or len(c)<5 or c[2]!='Peptide':error('invalid peptide QUERY header')
                internal=c[4].lstrip('>').split()[0] if c[4].strip() else ''
                row=db.execute('SELECT n,length FROM proteins WHERE internal=?',(internal,)).fetchone()
                if not row or row[1]!=int(c[3]) or internal in seen:error('unknown, duplicate, or length-mismatched query definition')
                seen.add(internal);query=(c[1],row[0],row[1])
            elif token=='ENDQUERY':
                if not query or domains:error('invalid ENDQUERY')
                query=None
            elif token=='DOMAINS':
                if not query or domains:error('invalid DOMAINS')
                domains=True
            elif token=='ENDDOMAINS':
                if not domains:error('invalid ENDDOMAINS')
                domains=False
            elif domains:
                if len(c)!=12 or c[0]!=session or c[1]!=query[0] or c[2] not in HIT_TYPES:error('invalid domain row')
                uid=c[3];start,end=int(c[4]),int(c[5]);ev,score=float(c[6]),float(c[7])
                if not uid.isdigit() or not 1<=start<=end<=query[2] or not math.isfinite(ev) or ev<0 or not math.isfinite(score) or score<=0:error('invalid coordinates or scores')
                if c[10] not in ['-','N','C','NC'] or not (c[11]=='-' or c[11].isdigit()):error('invalid incompleteness/superfamily')
                reported_uid=uid
                if uid=='0' and c[2]=='Superfamily':
                    parent=db.execute('SELECT pssm_id FROM cdd_models WHERE accession=?',(c[8],)).fetchone()
                    if not parent:error('derived superfamily accession missing from snapshot')
                    uid=parent[0]
                model=db.execute('SELECT accession FROM cdd_models WHERE pssm_id=?',(uid,)).fetchone()
                if not model or model[0]!=c[8]:error('PSSM/accession missing from this snapshot')
                if ev>evalue:continue
                raw=db.execute('SELECT model_start,model_end FROM rps_alignments WHERE protein_n=? AND pssm_id=? AND start=? AND end=? ORDER BY bitscore DESC,evalue LIMIT 1',(query[1],uid,start,end)).fetchone()
                db.execute('INSERT INTO cdd_hits VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',(query[1],uid,start,end,ev,score,c[2],c[10],c[11],*(raw or (None,None)),reported_uid));count+=1
            else:error(f'unexpected token {token!r}; request rpsbproc -t doms')
    if not ended or data or session or query or domains:raise ValueError('Incomplete native rpsbproc output')
    if len(seen)!=db.execute('SELECT COUNT(*) FROM proteins').fetchone()[0]:
        raise ValueError('Native output does not assess every input protein; do not treat missing queries as domain-negative')
    db.commit();return count


def hit_rows(db):
    sql='''SELECT p.n,p.original,p.internal,p.cluster,p.length,m.accession,m.name,m.pssm_id,m.length,m.cutoff,
           h.start,h.end,h.model_start,h.model_end,h.evalue,h.bitscore,h.hit_type,h.incomplete,h.superfamily_id,h.reported_pssm_id
           FROM cdd_hits h JOIN proteins p ON p.n=h.protein_n JOIN cdd_models m ON m.pssm_id=h.pssm_id
           ORDER BY p.n,h.start,h.end,m.accession'''
    for n,pid,internal,cluster,plen,acc,name,uid,mlen,cutoff,start,end,ms,me,ev,score,kind,incomplete,parent,reported_uid in db.execute(sql):
        yield n,dict(protein_id=pid,internal_id=internal,cluster=cluster,accession=acc,model_name=name,pssm_id=uid,
                     reported_pssm_id=reported_uid,source_database=source_database(acc),hit_type=kind,start=start,end=end,model_start=ms if ms is not None else '',
                     model_end=me if me is not None else '',model_length=mlen,bitscore=score,evalue=ev,
                     model_coverage=(me-ms+1)/mlen if ms is not None and mlen else '',protein_coverage=(end-start+1)/plen,
                     incomplete=incomplete,superfamily_pssm_id=parent,specific_bitscore_cutoff=cutoff if cutoff is not None else '',
                     threshold='E-value filter; Specific assignments use rpsbproc model-specific criteria',
                     coverage_status='matched_raw_alignment' if ms is not None else 'unavailable_for_derived_or_unmatched_hit')


def select_hits(hits,policy='best',overlap_fraction=.5,min_coverage=0.):
    """Avoid hierarchy and overlapping-repeat inflation; keep disjoint repeats."""
    kept=[];decisions=[]
    def overlap(a,b):return max(0,min(a['end'],b['end'])-max(a['start'],b['start'])+1)/min(a['end']-a['start']+1,b['end']-b['start']+1)
    for hit in sorted(hits,key=lambda h:({'Specific':0,'Non-Specific':1,'Non-specific':1,'Superfamily':2,'Multidom':3}[h['hit_type']],h['evalue'],-h['bitscore'],h['accession'],h['start'],h['end'])):
        reason='accepted'
        if hit['hit_type']=='Multidom':reason='whole_protein_multidomain_model_not_a_domain'
        elif min_coverage and hit['model_coverage']=='':reason='model_coverage_unavailable'
        elif hit['model_coverage']!='' and hit['model_coverage']<min_coverage:reason='below_min_model_coverage'
        else:
            for old in kept:
                if overlap(hit,old)<overlap_fraction:continue
                if hit['pssm_id']==old['pssm_id']:reason='overlapping_duplicate_model';break
                if hit['hit_type']=='Superfamily' and old['superfamily_pssm_id']==hit['pssm_id']:
                    reason='superfamily_redundant_with_selected_domain';break
        if reason=='accepted':kept.append(hit)
        decisions.append(dict(hit,selected=reason=='accepted',selection_reason=reason))
    return sorted(kept,key=lambda h:(h['start'],h['end'],h['accession'])),decisions


def run(args):
    validate(args)
    if args.threads<1 or not math.isfinite(args.evalue) or args.evalue<=0 or not math.isfinite(args.min_model_coverage) or not 0<=args.min_model_coverage<=1 or not 0<args.overlap_fraction<=1 or args.min_pair_count<1 or args.max_plot_nodes<1 or args.max_plot_edges<0:
        raise ValueError('Invalid threads, E-value, coverage, overlap, or network parameters')
    mode,files=resolve_input(args.input)
    if any(args.output.resolve() in p.resolve().parents for _,p in files):raise ValueError('Output must not contain input FASTAs')
    if args.cdd.resolve()==args.output.resolve() or args.cdd.resolve() in args.output.resolve().parents or args.output.resolve() in args.cdd.resolve().parents:
        raise ValueError('Use separate, non-nested CDD database and analysis output folders')
    if bool(args.rpsbproc_output)!=bool(args.raw_hits):raise ValueError('Imports require both --rpsbproc-output and --raw-hits from the same normalized queries/database')
    info=checked_snapshot(args.cdd);prefix=(args.cdd/info['db_prefix']).resolve();data_dir=(args.cdd/info['data_dir']).resolve()
    for name in ['manifest.json','manifest.inprogress.json']:
        previous=args.output/name
        if previous.is_file() and json.loads(previous.read_text()).get('domain_source')!='CDD':
            raise ValueError('Existing output belongs to another analysis; choose a separate CDD output folder')
    log=prepare_output(args.output,args.resume,args.force);db=None
    manifest=base_manifest(args);manifest.update(input_mode=mode,domain_source='CDD',software=versions(['biopython','networkx','matplotlib']),database=info)
    try:
        work=args.output/'work';work.mkdir(exist_ok=True)
        combined=work/'proteins.faa';index=work/'proteins.sqlite'
        db,total=normalize(files,combined,index,args.output/'id_mapping.tsv',include_cluster=mode=='clusters')
        manifest['input_sha256']={str(p.resolve()):sha256(p) for _,p in files}
        imported=bool(args.rpsbproc_output)
        executables={name:shutil.which(name) for name in ['rpsblast','blast_formatter','rpsbproc']}
        executables['rpsblast']=executables['rpsblast'] or shutil.which('rpsblast+')
        if not imported and not all(executables.values()):
            raise ValueError('CDD needs BLAST+ and rpsbproc. Install: conda install -n eet-clusters -c conda-forge -c bioconda blast rpsbproc')
        tool_versions={} if imported else {name:tool_version(executable,'-version') for name,executable in executables.items()}
        manifest['software'].update(tool_versions)
        native=work/'rpsbproc.out';raw=work/'rpsblast.tsv';archive=work/'rpsblast.asn'
        if imported:
            if any(p.resolve() in [native.resolve(),raw.resolve()] for p in [args.rpsbproc_output,args.raw_hits]):raise ValueError('Keep imported files outside generated work paths')
            manifest['import_sha256']={str(p.resolve()):sha256(p) for p in [args.rpsbproc_output,args.raw_hits]}
        key=signature(dict(schema=1,input=manifest['input_sha256'],sequence=sha256(combined),database=sha256(args.cdd/'cdd_manifest.json'),tools=tool_versions,
                           evalue=args.evalue,threads=args.threads,imported=manifest.get('import_sha256'),search=dict(comp_based_stats=1,seg='no',max_target_seqs=1000000,rpsbproc_mode='rep',show_families=False)))
        stamp=work/'cdd.stage.json';reuse=stage_reusable(stamp,key,[native,raw],args.resume and not args.force)
        if reuse:LOG.info('Reusing completed RPS-BLAST/rpsbproc search')
        else:
            stamp.unlink(missing_ok=True);dump(args.output/'manifest.inprogress.json',manifest)
            if imported:shutil.copyfile(args.rpsbproc_output,native);shutil.copyfile(args.raw_hits,raw)
            else:
                command([executables['rpsblast'],'-query',combined.resolve(),'-db',prefix,'-evalue',args.evalue,'-num_threads',args.threads,'-seg','no','-comp_based_stats','1','-max_target_seqs','1000000','-outfmt','11','-out',archive.resolve()],work/'rpsblast.log',manifest)
                command([executables['blast_formatter'],'-archive',archive.resolve(),'-outfmt',RAW_FORMAT,'-out',raw.resolve()],work/'blast_formatter.log',manifest)
                command([executables['rpsbproc'],'-i',archive.resolve(),'-o',native.resolve(),'-d',data_dir,'-e',args.evalue,'-m','rep','-t','doms'],work/'rpsbproc.log',manifest)
        metadata(db,data_dir,args.evalue);read_alignments(db,raw);count=parse_native(db,native,args.evalue)
        mark_stage(stamp,key,[native,raw])
        f,w=writer(args.output/'model_inventory.tsv',['pssm_id','accession','name','description','model_length','specific_bitscore_cutoff','source_database','evalue_threshold'])
        nodes={}
        for uid,acc,name,desc,length,cutoff in db.execute('SELECT * FROM cdd_models ORDER BY accession'):
            w.writerow(dict(pssm_id=uid,accession=acc,name=name,description=desc,model_length=length,specific_bitscore_cutoff=cutoff if cutoff is not None else '',source_database=source_database(acc),evalue_threshold=args.evalue))
            nodes[acc]=dict(domain_source='CDD',cdd_pssm_id=uid,source_database=source_database(acc),description=desc)
        f.close()
        manifest['counts']=summarize(db,args.output,args,input_mode=mode,hit_source=hit_rows,domain_selector=select_hits,hit_fields=HIT_FIELDS,node_metadata=nodes)
        manifest.update(search_definition=dict(algorithm='RPS-BLAST + rpsbproc representative domain hits',evalue=args.evalue,composition_correction=True,seg=False,cache_reused=reuse,native_hits=count,total_proteins=total),
                        interpretation='Within-protein representative CDD domains; whole-protein Multidom models excluded and overlapping hierarchy duplicates suppressed. Imported Pfam models are not independent evidence of a Pfam hit.',
                        skipped=['CDD conserved-site annotation (domain-only run)','localization predictions','FeGenie/FEET searches'])
        manifest['output_sha256']={p.name:sha256(p) for p in args.output.glob('*.tsv')}
        finish(args.output,manifest);LOG.info('Finished: %s',args.output/'index.html')
    except Exception as exc:
        manifest['error']=str(exc);dump(args.output/'manifest.inprogress.json',manifest);raise
    finally:
        if db:db.close()
        LOG.removeHandler(log);log.close()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i','--input','--clusters',type=Path,required=True)
    p.add_argument('-o','--output',type=Path,required=True)
    p.add_argument('--cdd',type=Path,required=True,help='Database directory prepared by download_cdd.py')
    p.add_argument('-t','--threads',type=int,default=8)
    p.add_argument('--evalue',type=float,default=.01)
    p.add_argument('--min-model-coverage',type=float,default=0.)
    p.add_argument('--overlap-fraction',type=float,default=.5,help='Overlap/shorter hit for duplicate model/hierarchy suppression')
    p.add_argument('--min-pair-count',type=int,default=1)
    p.add_argument('--max-plot-nodes',type=int,default=60)
    p.add_argument('--max-plot-edges',type=int,default=300)
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--clustering',type=Path,help='Link to an existing MCL protein network')
    p.add_argument('--rpsbproc-output',type=Path,help='Import native representative domain output using generated pNNNNNNNNNN query IDs')
    p.add_argument('--raw-hits',type=Path,help='Paired imported ten-column blast_formatter table; see CDD.md')
    p.add_argument('--resume',action='store_true');p.add_argument('--force',action='store_true')
    options(p);args=p.parse_args(argv);args.overlap_policy='representative';args.domain_source='CDD'
    logging.basicConfig(level=logging.INFO,format='%(levelname)s %(message)s')
    try:run(args)
    except (ValueError,OSError,sqlite3.Error,KeyError,ImportError) as e:LOG.error('%s',e);return 2
    return 0


if __name__=='__main__':sys.exit(main())
