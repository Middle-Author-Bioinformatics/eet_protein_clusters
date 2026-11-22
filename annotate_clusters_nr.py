#!/usr/bin/env python3
"""Append top DIAMOND/NR hit annotations to a folder of protein FASTAs.

Standalone: Python >=3.9 standard library plus an installed DIAMOND executable.
Input files are never edited. See NR_ANNOTATION.md for ranking and resume rules.
"""
import argparse
from collections import Counter, deque
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import platform
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile

TOOL = 'annotate_clusters_nr'
VERSION = '1.0.0'
FIELDS = 'qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen stitle'.split()
LOG = logging.getLogger(TOOL)
SUFFIXES = ('.faa', '.fa', '.fasta', '.faa.gz', '.fa.gz', '.fasta.gz')


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def fasta_open(path, mode='rt'):
    return (gzip.open if path.name.lower().endswith('.gz') else open)(path, mode, newline='')


def records(path):
    """Normalize search sequences only; output later copies original sequence lines."""
    header = None
    parts = []
    with fasta_open(path) as source:
        for line_no, line in enumerate(source, 1):
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header = line[1:].rstrip('\r\n')
                if not header.strip():
                    raise ValueError(f'{path}:{line_no}: empty FASTA header')
                parts = []
            elif line.strip():
                if header is None:
                    raise ValueError(f'{path}:{line_no}: sequence before the first FASTA header')
                parts.append(''.join(line.split()))
    if header is not None:
        yield header, ''.join(parts)


def prepare(files, work):
    index = work / 'queries.sqlite'
    index.unlink(missing_ok=True)
    db = sqlite3.connect(index)
    db.execute('CREATE TABLE queries(qid TEXT PRIMARY KEY, filename TEXT, ordinal INTEGER, protein_id TEXT, header TEXT, length INTEGER)')
    db.execute('CREATE UNIQUE INDEX file_record ON queries(filename,ordinal)')
    db.execute('CREATE TABLE hits(qid TEXT PRIMARY KEY,payload TEXT)')
    counts = Counter()
    total = 0
    try:
        with open(work / 'queries.faa', 'w') as query:
            for path in files:
                counts[path.name] = 0
                for ordinal, (header, sequence) in enumerate(records(path), 1):
                    normalized = sequence.upper().rstrip('*')
                    if not normalized or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', normalized):
                        raise ValueError(f'{path}: {header.split()[0]} has an empty/invalid protein sequence or internal stop')
                    total += 1
                    qid = f'q{total:012d}'
                    counts[path.name] += 1
                    db.execute('INSERT INTO queries VALUES (?,?,?,?,?,?)', (qid, path.name, ordinal, header.split()[0], header, len(normalized)))
                    query.write(f'>{qid}\n{normalized}\n')
        if not sum(counts.values()):
            raise ValueError('No protein records found in the input FASTA files')
        db.commit()
    except BaseException:
        db.close()
        raise
    return db, counts


def clean_title(title, subject):
    title = ' '.join(''.join(c if ord(c) >= 32 else ' ' for c in title).split())
    if title == subject:
        return ''
    if title.startswith(subject + ' '):
        title = title[len(subject):].lstrip()
    return title


def parse_hits(db, path, args):
    """Select highest bit score; resolve ties deterministically, without name bias."""
    reported = 0
    with open(path) as source:
        for line_no, line in enumerate(source, 1):
            if not line.strip():
                continue
            values = line.rstrip('\r\n').split('\t', len(FIELDS) - 1)
            if len(values) != len(FIELDS):
                raise ValueError(f'{path}:{line_no}: expected {len(FIELDS)} tab-separated columns, including stitle')
            qid, subject, identity, length, qs, qe, ss, se, evalue, bits, qlen, slen, title = values
            query = db.execute('SELECT length FROM queries WHERE qid=?', (qid,)).fetchone()
            if query is None:
                raise ValueError(f'{path}:{line_no}: unknown query {qid}')
            try:
                identity, evalue, bits = float(identity), float(evalue), float(bits)
                length, qs, qe, ss, se, qlen, slen = map(int, (length, qs, qe, ss, se, qlen, slen))
            except ValueError as exc:
                raise ValueError(f'{path}:{line_no}: invalid numeric alignment field') from exc
            if not all(map(math.isfinite, (identity, evalue, bits))) or not 0 <= identity <= 100 or evalue < 0 or bits < 0 or length < 1 or qlen != query[0] or not (1 <= qs <= qe <= qlen and 1 <= ss <= se <= slen) or not subject:
                raise ValueError(f'{path}:{line_no}: inconsistent alignment coordinates, length or score')
            qcov, scov = 100 * (qe - qs + 1) / qlen, 100 * (se - ss + 1) / slen
            if evalue > args.evalue or identity < args.min_identity or qcov + 1e-8 < args.query_cover or scov + 1e-8 < args.subject_cover:
                continue
            hit = dict(nr_hit=subject, nr_title=title, annotation=clean_title(title, subject), identity_percent=identity,
                       alignment_length=length, query_start=qs, query_end=qe, subject_start=ss, subject_end=se,
                       evalue=evalue, bitscore=bits, query_length=qlen, subject_length=slen,
                       query_coverage_percent=qcov, subject_coverage_percent=scov, best_score_ties=1)
            previous = db.execute('SELECT payload FROM hits WHERE qid=?', (qid,)).fetchone()
            if previous:
                old = json.loads(previous[0])
                if bits < old['bitscore']:
                    continue
                if bits == old['bitscore']:
                    hit['best_score_ties'] = old['best_score_ties'] + 1
                    rank = lambda h: (-h['evalue'], h['query_coverage_percent'], h['subject_coverage_percent'], h['identity_percent'])
                    if rank(hit) < rank(old) or (rank(hit) == rank(old) and (subject, title) >= (old['nr_hit'], old['nr_title'])):
                        old['best_score_ties'] = hit['best_score_ties']
                        hit = old
            db.execute('INSERT OR REPLACE INTO hits VALUES (?,?)', (qid, json.dumps(hit)))
            reported += 1
    db.commit()
    return reported


def emit(db, files, output):
    fields = ['file', 'record', 'protein_id', 'status', 'nr_hit', 'annotation', 'identity_percent',
              'query_coverage_percent', 'subject_coverage_percent', 'evalue', 'bitscore', 'best_score_ties',
              'alignment_length', 'query_start', 'query_end', 'subject_start', 'subject_end',
              'query_length', 'subject_length', 'nr_title', 'original_header']
    summaries = []
    with open(output / 'annotations.tsv', 'w', newline='') as table:
        writer = csv.DictWriter(table, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for path in files:
            rows = iter(db.execute('SELECT q.ordinal,q.protein_id,q.header,h.payload FROM queries q LEFT JOIN hits h USING(qid) WHERE q.filename=? ORDER BY q.ordinal', (path.name,)))
            total, matches, missing_titles = 0, 0, 0
            fd, temporary = tempfile.mkstemp(prefix='.nr-', suffix='.gz' if path.name.lower().endswith('.gz') else '.faa', dir=output)
            os.close(fd)
            try:
                with fasta_open(path) as source, fasta_open(Path(temporary), 'wt') as dest:
                    for line in source:
                        if not line.startswith('>'):
                            dest.write(line)
                            continue
                        row = next(rows, None)
                        if row is None:
                            raise ValueError(f'{path}: input record count changed during annotation')
                        ordinal, protein_id, header, payload = row
                        if header != line[1:].rstrip('\r\n'):
                            raise ValueError(f'{path}: input changed during annotation')
                        hit = json.loads(payload) if payload else {}
                        status = 'hit' if hit.get('annotation') else ('hit_title_unavailable' if payload else 'no_qualifying_hit')
                        total += 1
                        matches += bool(payload)
                        missing_titles += status == 'hit_title_unavailable'
                        if payload:
                            annotation = hit['annotation'] or '[title unavailable in DIAMOND database]'
                            suffix = f" | nr_hit={hit['nr_hit']} | nr_annotation={annotation}"
                        else:
                            suffix = ' | nr_status=no_qualifying_hit'
                        newline = '\r\n' if line.endswith('\r\n') else '\n'
                        dest.write('>' + header + suffix + newline)
                        writer.writerow(dict(hit, file=path.name, record=ordinal, protein_id=protein_id,
                                             original_header=header, status=status))
                if next(rows, None) is not None:
                    raise ValueError(f'{path}: input record count changed during annotation')
                os.replace(temporary, output / path.name)
            finally:
                Path(temporary).unlink(missing_ok=True)
            summaries.append(dict(file=path.name, proteins=total, with_hit=matches, without_hit=total - matches,
                                  hit_title_unavailable=missing_titles))
    with open(output / 'cluster_summary.tsv', 'w', newline='') as table:
        writer = csv.DictWriter(table, list(summaries[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(summaries)
    return summaries


def execute(command, log, manifest):
    manifest['commands'].append([str(a) for a in command])
    LOG.info('Running DIAMOND; log: %s', log)
    with open(log, 'w') as stream:
        process = subprocess.run([str(a) for a in command], stdout=stream, stderr=subprocess.STDOUT)
    if process.returncode:
        with open(log, errors='replace') as stream:
            tail = ''.join(deque(stream, maxlen=15))
        raise RuntimeError(f'DIAMOND exited {process.returncode}; see {log}\n{tail}')
    return log.read_text().strip() if log.stat().st_size < 100000 else ''


def run(args):
    args.input = args.input.resolve()
    args.output = args.output.resolve()
    args.db = args.db.resolve()
    if not args.db.is_file() and Path(str(args.db) + '.dmnd').is_file():
        args.db = Path(str(args.db) + '.dmnd')
    if not args.input.is_dir():
        raise ValueError(f'Input folder not found: {args.input}')
    if args.input == args.output or args.input in args.output.parents or args.output in args.input.parents:
        raise ValueError('Use a separate output folder, outside the input folder; source files are never overwritten')
    if not args.db.is_file():
        raise ValueError(f'DIAMOND database not found: {args.db}; use --db /path/to/nr.dmnd')
    if args.threads < 1 or args.max_target_seqs < 0 or not math.isfinite(args.block_size) or args.block_size <= 0 or not math.isfinite(args.evalue) or args.evalue <= 0 or any(not math.isfinite(v) or not 0 <= v <= 100 for v in (args.query_cover, args.subject_cover, args.min_identity)):
        raise ValueError('Threads/block-size/E-value must be positive, max-target-seqs >=0; coverage and identity must be percentages from 0 to 100')
    files = sorted(p for p in args.input.iterdir() if p.is_file() and p.name.lower().endswith(SUFFIXES))
    if not files:
        raise ValueError(f'No .faa/.fa/.fasta (optionally .gz) files in {args.input}')
    if any((args.output / path.name).resolve() == path.resolve() for path in files):
        raise ValueError('An output filename resolves to an input file; use a separate output without input symlinks')
    binary = shutil.which(args.diamond)
    if not binary:
        raise ValueError('DIAMOND not found. Activate eet-clusters or run: conda install -c conda-forge -c bioconda diamond')
    previous = {}
    if args.output.exists() and any(args.output.iterdir()):
        try:
            previous = json.loads((args.output / 'manifest.json').read_text())
        except (OSError, ValueError):
            pass
        if not (args.resume or args.force) or previous.get('tool') != TOOL:
            raise ValueError('Output folder is not empty. Choose a new folder, or --resume/--force for an existing output from this script')
    args.output.mkdir(parents=True, exist_ok=True)
    work = args.output / 'work'
    work.mkdir(exist_ok=True)
    handler = logging.FileHandler(args.output / 'run.log')
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    LOG.addHandler(handler)
    manifest = dict(tool=TOOL, version=VERSION, status='running', started_utc=datetime.now(timezone.utc).isoformat(),
                    parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}, commands=[],
                    software=dict(python=platform.python_version()), script_sha256=sha256(Path(__file__)),
                    generated_fastas=sorted(set(previous.get('generated_fastas', [])) | {p.name for p in files}))
    manifest_path = args.output / 'manifest.json'
    db = None
    try:
        write_json(manifest_path, manifest)
        manifest['input_sha256'] = {p.name: sha256(p) for p in files}
        before = args.db.stat()
        database = dict(path=str(args.db), size_bytes=before.st_size, mtime_ns=before.st_mtime_ns, version_label=args.db_version)
        if args.hash_db:
            LOG.info('Hashing full DIAMOND database (optional; may take time for NR)')
            database['sha256'] = sha256(args.db)
        manifest['database'] = database
        version = execute([binary, 'version'], work / 'diamond_version.log', manifest)
        manifest['software']['diamond'] = version
        execute([binary, 'dbinfo', '--db', args.db], work / 'diamond_dbinfo.log', manifest)
        db, counts = prepare(files, work)
        LOG.info('Prepared %d proteins from %d files', sum(counts.values()), len(files))
        key = dict(schema=1, input_sha256=manifest['input_sha256'], query_sha256=sha256(work / 'queries.faa'),
                   database=database, diamond_version=version, diamond_executable=str(Path(binary).resolve()),
                   sensitivity=args.sensitivity, threads=args.threads, block_size=args.block_size, evalue=args.evalue,
                   min_identity=args.min_identity, query_cover=args.query_cover, subject_cover=args.subject_cover,
                   max_target_seqs=args.max_target_seqs, max_hsps=1, outfmt=FIELDS)
        signature = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
        stamp, hits = work / 'diamond.stage.json', work / 'nr_hits.tsv'
        reusable = False
        if args.resume and stamp.is_file() and hits.is_file():
            try:
                cached = json.loads(stamp.read_text())
                reusable = cached['signature'] == signature and cached['hits_sha256'] == sha256(hits)
            except (ValueError, KeyError):
                pass
        manifest['search_cache_reused'] = reusable
        manifest['search_definition'] = key
        if reusable:
            LOG.info('Reusing completed, checksum-verified DIAMOND search')
        else:
            stamp.unlink(missing_ok=True)
            parent = Path(args.tmpdir or os.environ.get('TMPDIR') or work).resolve()
            parent.mkdir(parents=True, exist_ok=True)
            partial = work / 'nr_hits.tsv.partial'
            try:
                with tempfile.TemporaryDirectory(prefix='nr-annotation-', dir=parent) as scratch:
                    manifest['scratch_directory'] = scratch
                    manifest['free_bytes'] = dict(output=shutil.disk_usage(work).free, scratch=shutil.disk_usage(scratch).free)
                    cmd = [binary, 'blastp', '--query', work / 'queries.faa', '--db', args.db, '--out', partial,
                           '--outfmt', '6', *FIELDS, '--threads', args.threads, '--evalue', args.evalue,
                           '--max-target-seqs', args.max_target_seqs, '--max-hsps', '1', '--id', args.min_identity,
                           '--query-cover', args.query_cover, '--subject-cover', args.subject_cover,
                           '--block-size', args.block_size, '--tmpdir', scratch, '--unal', '0']
                    if args.sensitivity != 'default':
                        cmd += ['--' + args.sensitivity]
                    write_json(manifest_path, manifest)
                    execute(cmd, work / 'diamond_blastp.log', manifest)
                    partial.replace(hits)
                manifest['scratch_cleanup'] = 'removed'
            finally:
                partial.unlink(missing_ok=True)
            write_json(stamp, dict(signature=signature, hits_sha256=sha256(hits)))
        parse_hits(db, hits, args)
        if {p.name: sha256(p) for p in files} != manifest['input_sha256']:
            raise ValueError('Input files changed during the search; rerun with stable inputs')
        after = args.db.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            stamp.unlink(missing_ok=True)
            raise ValueError('The DIAMOND database changed during this run; rerun with a stable database')
        manifest['summary'] = emit(db, files, args.output)
        current = {p.name for p in files}
        for name in set(previous.get('generated_fastas', [])) - current:
            if Path(name).name == name and name.lower().endswith(SUFFIXES):
                (args.output / name).unlink(missing_ok=True)
        manifest['generated_fastas'] = sorted(current)
        manifest['output_sha256'] = {name: sha256(args.output / name) for name in [*sorted(current), 'annotations.tsv', 'cluster_summary.tsv']}
        manifest['status'] = 'complete'
        manifest['completed_utc'] = datetime.now(timezone.utc).isoformat()
        write_json(manifest_path, manifest)
        LOG.info('Annotated FASTAs: %s; complete hit table: annotations.tsv', args.output)
    except (Exception, KeyboardInterrupt) as exc:
        manifest.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed', error=str(exc))
        if manifest.get('scratch_directory'):
            manifest['scratch_cleanup'] = 'removed' if not Path(manifest['scratch_directory']).exists() else 'incomplete'
        try:
            write_json(manifest_path, manifest)
        except OSError:
            LOG.error('Could not update manifest after failure')
        raise
    finally:
        if db is not None:
            db.close()
        LOG.removeHandler(handler)
        handler.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=Path, required=True, help='Folder of cluster protein FASTAs; top-level files only')
    parser.add_argument('-o', '--output', type=Path, required=True, help='Separate folder for annotated FASTAs and tables')
    parser.add_argument('--db', type=Path, default=Path('/home/ark/databases/nr.dmnd'), help='Existing DIAMOND protein database')
    parser.add_argument('-t', '--threads', type=int, default=32)
    parser.add_argument('--evalue', type=float, default=1e-5)
    parser.add_argument('--min-identity', type=float, default=0, help='Minimum identity (percent)')
    parser.add_argument('--query-cover', type=float, default=0, help='Minimum single-HSP query coverage (percent); 0 disables')
    parser.add_argument('--subject-cover', type=float, default=0, help='Minimum single-HSP subject coverage (percent); 0 disables')
    parser.add_argument('--max-target-seqs', type=int, default=25, help='DIAMOND target candidate limit; best qualifying reported hit is retained (0 reports all targets)')
    parser.add_argument('--sensitivity', choices=['default', 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'], default='sensitive')
    parser.add_argument('--block-size', type=float, default=0.4, help='DIAMOND sequence block size, in billions of letters')
    parser.add_argument('--tmpdir', '--diamond-tmpdir', type=Path, help='Scratch parent; defaults to TMPDIR, otherwise OUTPUT/work')
    parser.add_argument('--diamond', default='diamond', help='DIAMOND executable or its full path')
    parser.add_argument('--db-version', default='', help='Optional NR release/build label, recorded in the manifest/cache key')
    parser.add_argument('--hash-db', action='store_true', help='Hash the entire database for stronger provenance (can be slow for NR)')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--resume', action='store_true', help='Reuse a completed search when inputs, database and search settings match')
    group.add_argument('--force', action='store_true', help='Rerun search in an existing output owned by this script')
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    try:
        run(args)
    except (OSError, ValueError, RuntimeError, sqlite3.Error) as exc:
        LOG.error('%s', exc)
        return 2
    except KeyboardInterrupt:
        LOG.error('Interrupted; rerun with --resume to reuse completed search results')
        return 130
    return 0


if __name__ == '__main__':
    sys.exit(main())
