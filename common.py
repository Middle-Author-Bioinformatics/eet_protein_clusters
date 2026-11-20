"""Shared I/O and provenance for the two standalone protein-analysis scripts."""
import csv
import gzip
import hashlib
import importlib.metadata
import json
import logging
import re
import shutil
import sqlite3
import subprocess
from pathlib import Path
from datetime import datetime, timezone

VERSION = '1.1.0'
LOG = logging.getLogger('eet-clusters')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def signature(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()


def dump(path, obj):
    tmp=Path(str(path)+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    tmp.replace(path)


def writer(path, fields):
    stream=open(path,'w',newline='')
    w=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',extrasaction='ignore')
    w.writeheader()
    return stream,w


def versions(packages=()):
    import platform
    result={'python':platform.python_version(),'scripts':VERSION}
    for package in packages:
        try:result[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:result[package]='not_installed'
    return result


def executable(name):
    path=shutil.which(name)
    if not path:raise ValueError(f'{name} not found in PATH. Activate the environment described in README.md.')
    return path


def tool_version(name, flag='--version'):
    p=subprocess.run([executable(name),flag],text=True,capture_output=True)
    if p.returncode:raise ValueError(f'Could not run {name} {flag}: {p.stderr.strip()}')
    return (p.stdout+'\n'+p.stderr).strip()


def command(args, logfile, manifest):
    LOG.info('Running %s',args[0])
    manifest['commands'].append([str(a) for a in args])
    with open(logfile,'w') as out:
        result=subprocess.run([str(a) for a in args],stdout=out,stderr=subprocess.STDOUT)
    if result.returncode:
        tail='\n'.join(Path(logfile).read_text(errors='replace').splitlines()[-12:])
        raise ValueError(f'{args[0]} exited {result.returncode}. See {logfile}\n{tail}')


def stage_reusable(stamp, key, paths, resume):
    if not resume or not stamp.exists() or not all(p.is_file() for p in paths):return False
    try:
        old=json.loads(stamp.read_text())
        return old['key']==key and old['sha256']=={str(p.name):sha256(p) for p in paths}
    except (KeyError,ValueError):return False


def mark_stage(stamp,key,paths):
    dump(stamp,{'key':key,'sha256':{str(p.name):sha256(p) for p in paths}})


def prepare_output(output,resume=False,force=False):
    output.mkdir(parents=True,exist_ok=True)
    if any(output.iterdir()) and not (resume or force):
        raise ValueError(f'{output} is not empty; choose a new output directory or use --resume/--force.')
    log=logging.FileHandler(output/'run.log')
    log.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    LOG.addHandler(log)
    (output/'manifest.json').unlink(missing_ok=True)
    return log


def base_manifest(args):
    return {'version':VERSION,'started_utc':datetime.now(timezone.utc).isoformat(),
            'parameters':{k:str(v.resolve()) if isinstance(v,Path) else v for k,v in vars(args).items()},
            'commands':[],'warnings':[], 'script_sha256':{p.name:sha256(p) for p in Path(__file__).parent.glob('*.py')}}


def finish(output,manifest):
    manifest['completed_utc']=datetime.now(timezone.utc).isoformat()
    dump(output/'manifest.json',manifest)
    (output/'manifest.inprogress.json').unlink(missing_ok=True)


def open_fasta(path):
    return gzip.open(path,'rt') if str(path).lower().endswith('.gz') else open(path)


def normalize(files, fasta_path, database_path, mapping_path):
    """Stream sequences to FASTA and disk-backed index; IDs unique across all files.

    files is a list of (cluster_label, path). Original headers are retained in SQL.
    Synthetic IDs keep DIAMOND/MCL/HMMER parsing independent of user ID punctuation.
    """
    from Bio import SeqIO
    database_path.unlink(missing_ok=True)
    db=sqlite3.connect(database_path)
    db.execute('CREATE TABLE proteins (n INTEGER PRIMARY KEY, internal TEXT UNIQUE, original TEXT UNIQUE, header TEXT, length INTEGER, cluster TEXT)')
    n=0
    mf,mw=writer(mapping_path,['internal_id','protein_id','cluster','length_aa','original_header'])
    try:
        with open(fasta_path,'w') as dest:
            for cluster,path in files:
                local_count=0
                with open_fasta(path) as f:
                    for rec in SeqIO.parse(f,'fasta'):
                        seq=str(rec.seq).upper().rstrip('*')
                        if not seq or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+',seq):
                            raise ValueError(f'{path}: invalid sequence/internal stop in {rec.id}. Provide protein FASTA.')
                        n+=1;local_count+=1;internal=f'p{n:010d}'
                        try:db.execute('INSERT INTO proteins VALUES (?,?,?,?,?,?)',(n,internal,rec.id,rec.description,len(seq),cluster))
                        except sqlite3.IntegrityError as e:
                            raise ValueError(f'Duplicate protein ID {rec.id!r}; make IDs globally unique (e.g. genome|locus) before analysis.') from e
                        dest.write(f'>{internal}\n{seq}\n')
                        mw.writerow(dict(internal_id=internal,protein_id=rec.id,cluster=cluster,length_aa=len(seq),original_header=rec.description))
                if not local_count:raise ValueError(f'{path}: empty protein FASTA')
        db.commit()
    except Exception:
        db.close();raise
    finally:mf.close()
    if not n:
        db.close();raise ValueError('No protein sequences found')
    return db,n


def resolve(db,label,internal=False):
    column='internal' if internal else 'original'
    row=db.execute(f'SELECT n,length FROM proteins WHERE {column}=?',(label,)).fetchone()
    if row is None:raise ValueError(f'Unknown protein ID in similarity hits: {label!r}')
    return row
