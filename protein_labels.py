"""Import display annotations from FASTA headers without changing protein IDs."""
import gzip
import json
import os
from pathlib import Path
import sqlite3
import tempfile
from common import LOG, sha256

SUFFIXES = ('.faa', '.fa', '.fasta', '.faa.gz', '.fa.gz', '.fasta.gz')
LABEL_FIELDS = ['annotation', 'nr_hit', 'annotation_status', 'annotated_header']


def header_label(header):
    if not header.strip():
        raise ValueError('Empty protein FASTA header')
    protein_id = header.split()[0]
    annotation, nr_hit, status = '', '', 'unannotated'
    if header.endswith(' | nr_status=no_qualifying_hit'):
        status = 'no_qualifying_hit'
    elif ' | nr_annotation=' in header:
        prefix, annotation = header.rsplit(' | nr_annotation=', 1)
        hit = prefix.rsplit(' | nr_hit=', 1)
        if len(hit) != 2 or not hit[1].strip():
            raise ValueError(f'{protein_id}: nr_annotation has no matching nr_hit')
        nr_hit, annotation, status = hit[1].strip(), annotation.strip(), 'nr_top_hit'
        if annotation == '[title unavailable in DIAMOND database]':
            annotation, status = '', 'title_unavailable'
    else:
        parts = header.split(None, 1)
        if len(parts) == 2:
            annotation, status = parts[1].strip(), 'header_description'
    return protein_id, dict(annotation=annotation, nr_hit=nr_hit, annotation_status=status, annotated_header=header)


class ProteinLabels:
    def __init__(self, source, work):
        source, work = Path(source).resolve(), Path(work)
        if source.is_dir():
            files = sorted(p for p in source.iterdir() if p.is_file() and p.name.lower().endswith(SUFFIXES))
        elif source.is_file() and source.name.lower().endswith(SUFFIXES):
            files = [source]
        else:
            files = []
        if not files:
            raise ValueError(f'No protein FASTAs found at {source}')
        self.provenance = dict(source=str(source), files={str(p): sha256(p) for p in files}, matches={})
        key = json.dumps(['protein-labels-v1', self.provenance['files']], sort_keys=True)
        work.mkdir(parents=True, exist_ok=True)
        self.path = work / 'protein_labels.sqlite'
        if self.path.exists():
            try:
                with sqlite3.connect(self.path) as db:
                    if db.execute('SELECT signature FROM metadata').fetchone() == (key,):
                        self.provenance['proteins'] = db.execute('SELECT COUNT(*) FROM labels').fetchone()[0]
                        return
            except sqlite3.DatabaseError:
                pass
        fd, name = tempfile.mkstemp(prefix='protein-labels-', suffix='.sqlite', dir=work)
        os.close(fd)
        try:
            with sqlite3.connect(name) as db:
                db.executescript('CREATE TABLE labels(protein_id TEXT PRIMARY KEY,payload TEXT); CREATE TABLE metadata(signature TEXT);')
                for file in files:
                    opener = gzip.open if file.name.lower().endswith('.gz') else open
                    seen_header = False
                    with opener(file, 'rt') as stream:
                        for line in stream:
                            if line.startswith('>'):
                                seen_header = True
                                pid, label = header_label(line[1:].rstrip('\r\n'))
                                try:
                                    db.execute('INSERT INTO labels VALUES (?,?)', (pid, json.dumps(label)))
                                except sqlite3.IntegrityError as exc:
                                    raise ValueError(f'Duplicate protein ID {pid!r} in supplied FASTAs; IDs must be globally unique for the saved networks') from exc
                            elif line.strip() and not seen_header:
                                raise ValueError(f'{file}: sequence/text before its first FASTA header')
                self.provenance['proteins'] = db.execute('SELECT COUNT(*) FROM labels').fetchone()[0]
                if not self.provenance['proteins']:
                    raise ValueError('No protein headers found in supplied FASTAs')
                db.execute('INSERT INTO metadata VALUES (?)', (key,))
            os.replace(name, self.path)
        finally:
            Path(name).unlink(missing_ok=True)

    def validate_ids(self, database, column, dataset):
        if column not in {'protein_id', 'original'}:
            raise ValueError('Invalid protein index column')
        with sqlite3.connect(database) as db:
            db.execute('ATTACH DATABASE ? AS display_labels', (str(self.path.resolve()),))
            total = db.execute('SELECT COUNT(*) FROM proteins').fetchone()[0]
            matched = db.execute(f'SELECT COUNT(*) FROM proteins p JOIN display_labels.labels l ON p.{column}=l.protein_id').fetchone()[0]
        self.provenance['matches'][dataset] = dict(matched=matched, total=total, without_imported_header=total-matched)
        if total and not matched:
            raise ValueError('No FASTA protein IDs match the saved analysis. Keep the original first whitespace-delimited IDs; append annotations after them. Use FASTAs from this analysis.')
        LOG.info('Imported display headers for %d/%d proteins (%s)', matched, total, dataset)

    def decorate(self, rows):
        if not rows:
            return rows
        with sqlite3.connect(self.path) as db:
            db.execute('CREATE TEMP TABLE requested(protein_id TEXT PRIMARY KEY)')
            db.executemany('INSERT OR IGNORE INTO requested VALUES (?)', ((r['protein_id'],) for r in rows))
            found = {pid: json.loads(payload) for pid, payload in db.execute('SELECT l.protein_id,l.payload FROM labels l JOIN requested r USING(protein_id)')}
        for row in rows:
            row.update(found.get(row['protein_id'], {}))
        return rows
