"""Disk-backed carrier lookup from selected architectures, never aggregate edges."""
import csv
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
from common import sha256, writer


class DomainMembers:
    def __init__(self, folder):
        self.folder = Path(folder)
        self.source = self.folder / 'protein_architectures.tsv'
        self.available = self.source.is_file()
        self.reason = 'Protein IDs unavailable: aggregate GraphML does not contain carriers. Supply the original protein_architectures.tsv in the domain results folder.'
        if not self.available:
            return
        work = self.folder / 'work'
        work.mkdir(exist_ok=True)
        self.path = work / 'network_members.sqlite'
        self.checksum = sha256(self.source)
        key = json.dumps(['domain-members-v1', self.checksum])
        if self.path.exists():
            try:
                with sqlite3.connect(self.path) as db:
                    if db.execute('SELECT value FROM metadata WHERE key="signature"').fetchone() == (key,):
                        self.has_clusters = bool(db.execute('SELECT value FROM metadata WHERE key="has_clusters"').fetchone()[0])
                        return
            except sqlite3.DatabaseError:
                pass
        fd, name = tempfile.mkstemp(prefix='network_members-', suffix='.sqlite', dir=work)
        os.close(fd)
        try:
            with sqlite3.connect(name) as db, open(self.source) as source:
                rows = csv.DictReader(source, delimiter='\t')
                if not {'protein_id', 'domain_set'}.issubset(rows.fieldnames or []):
                    raise ValueError(f'{self.source}: requires protein_id and domain_set columns')
                self.has_clusters = 'cluster' in rows.fieldnames
                db.executescript('''CREATE TABLE metadata(key TEXT PRIMARY KEY,value);
                    CREATE TABLE proteins(n INTEGER PRIMARY KEY,protein_id TEXT UNIQUE,cluster TEXT,domain_names TEXT,architecture_names TEXT);
                    CREATE TABLE carriers(n INTEGER,accession TEXT,PRIMARY KEY(n,accession)) WITHOUT ROWID;''')
                for n, row in enumerate(rows):
                    if not row['protein_id']:
                        raise ValueError(f'{self.source}: empty protein_id')
                    db.execute('INSERT INTO proteins VALUES (?,?,?,?,?)', (n, row['protein_id'], row.get('cluster', ''), row.get('domain_names', row['domain_set']), row.get('architecture_names', row.get('architecture', ''))))
                    db.executemany('INSERT INTO carriers VALUES (?,?)', ((n, a) for a in sorted(set(row['domain_set'].split(';'))) if a))
                db.execute('CREATE INDEX carrier_domain ON carriers(accession,n)')
                db.executemany('INSERT INTO metadata VALUES (?,?)', [('signature', key), ('has_clusters', int(self.has_clusters))])
            os.replace(name, self.path)
        finally:
            Path(name).unlink(missing_ok=True)

    def payload(self, graph, chosen, limit):
        if not self.available:
            return dict(available=False, reason=self.reason)
        scope = graph.graph.get('cluster') if self.has_clusters and graph.graph.get('input_mode') != 'proteins' else None
        with sqlite3.connect(self.path) as db:
            db.execute('CREATE TEMP TABLE selected(accession TEXT PRIMARY KEY)')
            db.executemany('INSERT INTO selected VALUES (?)', ((a,) for a in chosen))
            scope_sql = ' AND p.cluster=?' if scope is not None else ''
            params = (scope,) if scope is not None else ()
            # Refuse stale/mismatched architecture tables rather than inventing lists.
            counts = dict(db.execute('SELECT c.accession,COUNT(*) FROM carriers c JOIN selected s USING(accession) JOIN proteins p USING(n) WHERE 1=1' + scope_sql + ' GROUP BY c.accession', params))
            if any(counts.get(a, 0) != graph.nodes[a]['proteins_with_domain'] for a in chosen):
                return dict(available=False, reason='Protein IDs unavailable: architecture carrier counts do not match this graph. Regenerate both from the same domain analysis.')
            query = 'SELECT p.n,p.protein_id,p.cluster,p.domain_names,p.architecture_names FROM proteins p WHERE EXISTS (SELECT 1 FROM carriers c JOIN selected s USING(accession) WHERE c.n=p.n)' + scope_sql
            total = db.execute('SELECT COUNT(*) FROM (' + query + ')', params).fetchone()[0]
            rows = list(db.execute(query + ' ORDER BY p.protein_id LIMIT ?', params + (limit,)))
            db.execute('CREATE TEMP TABLE included(n INTEGER PRIMARY KEY,i INTEGER)')
            db.executemany('INSERT INTO included VALUES (?,?)', ((r[0], i) for i, r in enumerate(rows)))
            by_domain = {a: [] for a in chosen}
            for a, i in db.execute('SELECT c.accession,i.i FROM carriers c JOIN selected s USING(accession) JOIN included i USING(n) ORDER BY i.i'):
                by_domain[a].append(i)
            # Also check full pair counts; marginal counts alone can conceal mismatched inputs.
            for a, b, d in graph.edges(data=True):
                if a not in by_domain or b not in by_domain:
                    continue
                count = db.execute('SELECT COUNT(*) FROM carriers a JOIN carriers b USING(n) JOIN proteins p USING(n) WHERE a.accession=? AND b.accession=?' + scope_sql, (a, b) + params).fetchone()[0]
                if count != d['proteins_with_both']:
                    return dict(available=False, reason='Protein IDs unavailable: architecture pair counts do not match this graph. Regenerate both from the same domain analysis.')
        return dict(available=True, source_sha256=self.checksum, matching_proteins=total, embedded_count=len(rows), truncated=len(rows) < total,
                    proteins=[dict(protein_id=r[1], cluster=r[2], domain_names=r[3], architecture_names=r[4]) for r in rows], by_domain=by_domain)

    def export_pair(self, accessions, path, cluster=None):
        if not self.available:
            raise ValueError(self.reason)
        with sqlite3.connect(self.path) as db:
            resolved = []
            for accession in accessions:
                if db.execute('SELECT 1 FROM carriers WHERE accession=? LIMIT 1', (accession,)).fetchone():
                    resolved.append(accession)
                    continue
                matches = [r[0] for r in db.execute('SELECT DISTINCT accession FROM carriers') if re.fullmatch(re.escape(accession) + r'\.\d+', r[0])]
                if len(matches) != 1:
                    raise ValueError(f'Unknown or ambiguous domain {accession}; use an exact accession from protein_architectures.tsv')
                resolved.extend(matches)
            if resolved[0] == resolved[1]:
                raise ValueError('Supply two different domains for a co-occurrence query')
            if cluster and not self.has_clusters:
                raise ValueError('--cluster requires cluster annotations; this is a single-FASTA domain analysis')
            if cluster and not db.execute('SELECT 1 FROM proteins WHERE cluster=? LIMIT 1', (cluster,)).fetchone():
                raise ValueError(f'Unknown cluster {cluster}')
            query = 'SELECT p.protein_id,p.cluster,p.domain_names,p.architecture_names FROM proteins p JOIN carriers a USING(n) JOIN carriers b USING(n) WHERE a.accession=? AND b.accession=?'
            params = tuple(resolved)
            if cluster:
                query += ' AND p.cluster=?'
                params += (cluster,)
            fields = ['protein_id'] + (['cluster'] if self.has_clusters else []) + ['domain_a', 'domain_b', 'domain_names', 'architecture_names']
            f, w = writer(path, fields)
            count = 0
            try:
                for pid, group, names, arch in db.execute(query + ' ORDER BY p.protein_id', params):
                    w.writerow(dict(protein_id=pid, cluster=group, domain_a=resolved[0], domain_b=resolved[1], domain_names=names, architecture_names=arch))
                    count += 1
            finally:
                f.close()
            return count
