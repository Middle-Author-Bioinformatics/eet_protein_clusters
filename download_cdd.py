#!/usr/bin/env python3
"""Download an isolated, checksummed NCBI CDD RPS-BLAST database snapshot."""
import argparse
from datetime import datetime, timezone
import gzip
import json
import logging
from pathlib import Path, PurePosixPath
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from common import LOG, dump, sha256

BASE='https://ftp.ncbi.nih.gov/pub/mmdb/cdd/'
ANNOTATIONS=['cddid.tbl.gz','cdtrack.txt','family_superfamily_links',
             'cddannot.dat.gz','cddannot_generic.dat.gz','bitscore_specific.txt']


def download(url,path):
    """Do not promote partial transfers; record HTTP metadata and local digests."""
    partial=path.with_name(path.name+'.part')
    LOG.info('Downloading %s',url)
    try:
        with urllib.request.urlopen(url,timeout=120) as src,partial.open('wb') as dst:
            expected=src.headers.get('Content-Length');size=0
            meta={k:src.headers.get(k) for k in ['ETag','Last-Modified','Content-Length']}
            while chunk:=src.read(1024*1024):dst.write(chunk);size+=len(chunk)
        if expected is not None and size!=int(expected):raise ValueError(f'Incomplete download: {url}')
        partial.replace(path)
    finally:partial.unlink(missing_ok=True)
    return dict(url=url,sha256=sha256(path),bytes=size,http=meta)


def extract(archive,destination):
    """Only regular files/directories; no links, devices, or path traversal."""
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            name=PurePosixPath(member.name)
            if name.is_absolute() or '..' in name.parts or '\\' in member.name:
                raise ValueError(f'Unsafe archive path: {member.name}')
            target=destination.joinpath(*name.parts)
            if member.isdir():target.mkdir(parents=True,exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True,exist_ok=True)
                with tar.extractfile(member) as source,target.open('wb') as dest:shutil.copyfileobj(source,dest)
            else:raise ValueError(f'Unsupported archive member: {member.name}')


def database_prefix(folder):
    aliases=sorted(folder.rglob('*.pal'))
    if len(aliases)==1:return aliases[0].with_suffix('')
    rps=sorted(folder.rglob('*.rps'))
    if len(rps)==1:return rps[0].with_suffix('')
    raise ValueError('Expected one RPS database or alias in the archive; none or ambiguous')


def checked_snapshot(folder):
    folder=Path(folder)
    path=folder/'cdd_manifest.json'
    if not path.is_file():raise ValueError(f'Missing {path}; run download_cdd.py first')
    data=json.loads(path.read_text())
    if data.get('schema')!=1:raise ValueError('Unsupported CDD snapshot manifest')
    for name in [*data['files'],data['db_prefix'],data['data_dir']]:
        if folder.resolve() not in (folder/name).resolve().parents:raise ValueError('CDD snapshot paths must stay inside its directory')
    for name,digest in data['files'].items():
        p=folder/name
        if not p.is_file() or sha256(p)!=digest:
            raise ValueError(f'CDD snapshot missing/changed file: {name}; use download --force to create a new snapshot')
    return data


def setup(args):
    args.output.mkdir(parents=True,exist_ok=True)
    manifest_path=args.output/'cdd_manifest.json'
    if manifest_path.exists() and not args.force:
        if not args.resume:raise ValueError('CDD already installed; use --resume to verify/reuse or --force for a new snapshot')
        old=checked_snapshot(args.output)
        if old['subset']!=args.subset:raise ValueError('Requested subset differs; choose another output folder or --force')
        LOG.info('Verified existing snapshot: %s',old['release_info'].splitlines()[0]);return old
    lock=args.output/'.download.lock'
    try:handle=lock.open('x')
    except FileExistsError:raise ValueError(f'Another setup may be active: {lock}; remove it only after checking no setup is running')
    snapshot=None
    try:
        handle.close();snapshots=args.output/'snapshots';snapshots.mkdir(exist_ok=True)
        snapshot=Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ-'),dir=snapshots))
        downloads=snapshot/'downloads';downloads.mkdir()
        db=snapshot/'db';db.mkdir();ann=snapshot/'data';ann.mkdir();sources={}
        sources['release_info']=download(BASE+'cdd.info',snapshot/'cdd.info')
        info=(snapshot/'cdd.info').read_text()
        archive='Cdd_LE.tar.gz' if args.subset=='full' else 'Cdd_NCBI_LE.tar.gz'
        sources['database']=download(BASE+'little_endian/'+archive,downloads/archive)
        extract(downloads/archive,db)
        for name in ANNOTATIONS:
            sources[name]=download(BASE+name,downloads/name)
            if name.endswith('.gz'):
                with gzip.open(downloads/name,'rb') as src,(ann/name[:-3]).open('wb') as dst:shutil.copyfileobj(src,dst)
            else:shutil.copyfile(downloads/name,ann/name)
        for name in ['README','rpsbproc/README']:
            target=snapshot/('CDD_README.txt' if name=='README' else 'rpsbproc_README.txt')
            sources[name]=download(BASE+name,target)
        sources['release_info_end']=download(BASE+'cdd.info',snapshot/'cdd.info.end')
        if (snapshot/'cdd.info.end').read_bytes()!=(snapshot/'cdd.info').read_bytes():
            raise ValueError('CDD release changed during download; rerun setup')
        prefix=database_prefix(db)
        if not args.keep_archive:shutil.rmtree(downloads)
        files={str(p.relative_to(args.output)):sha256(p) for p in sorted(snapshot.rglob('*')) if p.is_file()}
        result=dict(schema=1,subset=args.subset,created_utc=datetime.now(timezone.utc).isoformat(),release_info=info,
                    db_prefix=str(prefix.relative_to(args.output)),data_dir=str(ann.relative_to(args.output)),files=files,sources=sources,
                    note='Local SHA256 provenance, not signed upstream verification. Some annotation files update independently of cdd.info; preserve this entire snapshot.')
        dump(manifest_path,result);LOG.info('CDD ready: %s',prefix);return result
    except Exception:
        if snapshot and snapshot.exists():shutil.rmtree(snapshot)
        raise
    finally:lock.unlink(missing_ok=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('-o','--output',type=Path,required=True)
    p.add_argument('--subset',choices=['full','ncbi'],default='full',help='Full CDD or NCBI-curated subset; full includes imported Pfam and other sources')
    p.add_argument('--resume',action='store_true',help='Verify/reuse a complete snapshot; interrupted file downloads restart')
    p.add_argument('--force',action='store_true',help='Create a new snapshot, retaining old complete snapshots')
    p.add_argument('--keep-archive',action='store_true')
    args=p.parse_args(argv);logging.basicConfig(level=logging.INFO,format='%(levelname)s %(message)s')
    try:setup(args)
    except (ValueError,OSError,tarfile.TarError) as e:LOG.error('%s',e);return 2
    return 0


if __name__=='__main__':sys.exit(main())
