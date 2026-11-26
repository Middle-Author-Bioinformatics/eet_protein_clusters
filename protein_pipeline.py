#!/usr/bin/env python3
"""Protein FASTA -> MCL families -> NR/domains -> tables, curated HMMs, and trees."""
import argparse
from collections import defaultdict
import html
import json
import logging
import math
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys

from Bio import SeqIO
from common import LOG, base_manifest, command, dump, finish, open_fasta, prepare_output, sha256, versions
from family_analysis import (DEFAULTS, SEARCH, STYLE, analyze_family, cluster_table,
                             family_report, read_tsv, write_tsv)

ROOT = Path(__file__).resolve().parent


def program(name, alternatives=()):
    for candidate in [name, *alternatives]:
        found = shutil.which(candidate)
        if found:
            return found
    raise ValueError(f'{name} not found. Install the workflow dependencies: '
                     'conda install -c conda-forge -c bioconda diamond mcl mafft fasttree; '
                     'CDD additionally needs blast and rpsbproc. See PROTEIN_WORKFLOW.md.')


def version(path, flag):
    result = subprocess.run([path, flag], capture_output=True, text=True)
    text = (result.stdout + '\n' + result.stderr).strip()
    if result.returncode or not text:
        raise ValueError(f'Could not run {path} {flag}: {text[-1000:]}')
    return text


def preflight(args):
    if not args.input.is_file():
        raise ValueError(f'Protein FASTA not found: {args.input}')
    for protected in [args.input, args.nr_db, args.hmm if args.domain_db == 'pfam' else args.cdd]:
        if protected is not None and (protected.resolve() == args.output.resolve()
                                       or args.output.resolve() in protected.resolve().parents):
            raise ValueError('Output must not contain the input FASTA or a reference database')
    if args.output.is_dir() and any(args.output.iterdir()):
        previous = next((args.output/p for p in ['manifest.json', 'manifest.inprogress.json', 'manifest.failed.json']
                         if (args.output/p).is_file()), None)
        if previous is None or json.loads(previous.read_text()).get('workflow') != 'protein_families':
            raise ValueError('Choose a new output directory; this folder is not a protein_pipeline run')
        if json.loads(previous.read_text())['parameters']['domain_db'] != args.domain_db:
            raise ValueError('Use a separate output directory when switching between Pfam and CDD')
    if args.threads < 1 or args.min_hmm_sequences < 2 or args.min_hmm_columns < 1 or not 1 <= args.seed <= 2147483647:
        raise ValueError('Positive threads/columns, >=2 HMM sequences, and seed 1..2147483647 are required')
    positive = [args.inflation, args.evalue, args.nr_evalue, args.block_size, args.cdd_evalue]
    fractions = [args.min_coverage, args.min_model_coverage, args.min_column_occupancy, args.min_sequence_coverage]
    if (not all(math.isfinite(v) and v > 0 for v in positive)
            or args.inflation <= 1
            or not all(math.isfinite(v) and 0 <= v <= 1 for v in fractions)
            or args.min_column_occupancy == 0
            or not all(math.isfinite(v) and 0 <= v <= 100 for v in [args.min_identity, args.nr_min_identity, args.nr_query_cover, args.nr_subject_cover])
            or args.nr_max_target_seqs < 0
            or (args.domain_evalue is not None and (not math.isfinite(args.domain_evalue) or args.domain_evalue <= 0))):
        raise ValueError('Invalid thresholds: inflation must be >1, E-values/block size positive, '
                         'coverage/occupancy fractions 0..1, identity/NR coverage percentages 0..100')
    if not args.nr_db.is_file():
        raise ValueError(f'DIAMOND NR database not found: {args.nr_db}; supply --nr-db /path/to/nr.dmnd')
    # Fail before expensive searches if the FASTA or model thresholds are invalid.
    seen = set()
    with open_fasta(args.input) as stream:
        for record in SeqIO.parse(stream, 'fasta'):
            seq = str(record.seq).upper().rstrip('*')
            if record.id in seen:
                raise ValueError(f'Duplicate protein ID {record.id!r}; IDs must be globally unique')
            if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', seq):
                raise ValueError(f'{record.id}: invalid protein sequence or internal stop')
            seen.add(record.id)
    if not seen:
        raise ValueError('Empty protein FASTA')
    tools = {n: program(n) for n in ['diamond', 'mcl', 'mcxload']}
    tools['mafft'] = program(args.mafft)
    if not args.skip_tree:
        tools['fasttree'] = program(args.fasttree, ['fasttree'] if args.fasttree == 'FastTree' else [])
    if args.domain_db == 'pfam':
        if not args.hmm.is_file():
            raise ValueError(f'HMM database not found: {args.hmm}; supply --hmm /path/to/Pfam-A.hmm')
        from cluster_domains import model_inventory
        with sqlite3.connect(':memory:') as db:
            model_inventory(args.hmm, db, args.domain_evalue)
    else:
        from download_cdd import checked_snapshot
        checked_snapshot(args.cdd)
        for name, alternatives in [('rpsblast', ['rpsblast+']), ('blast_formatter', []), ('rpsbproc', [])]:
            tools[name] = program(name, alternatives)
    return tools, len(seen)


def child(script, argv, output, args, manifest):
    cmd = [sys.executable, str(ROOT/script), *map(str, argv)]
    cmd += ['--force'] if args.force else (['--resume'] if args.resume else [])
    log = args.output/'logs'/f'{Path(script).stem}.log'
    command(cmd, log, manifest)
    path = output/'manifest.json'
    manifest['stages'][script] = dict(manifest=str(path.resolve()), sha256=sha256(path),
                                    provenance=json.loads(path.read_text()))
    dump(args.output/'manifest.inprogress.json', manifest)


def grouped(rows, key, unique):
    result = defaultdict(dict)
    for row in rows:
        cluster, pid = row[key], row[unique]
        if pid in result[cluster]:
            raise ValueError(f'Duplicate {unique}: {pid} in {cluster}')
        result[cluster][pid] = row
    return result


def reports(args, tools, manifest):
    out = args.output
    members = grouped(read_tsv(out/'clustering/membership.tsv'), 'cluster', 'protein_id')
    nr = grouped(read_tsv(out/'nr_annotated/annotations.tsv'), 'file', 'protein_id')
    domains = grouped(read_tsv(out/'domains/protein_architectures.tsv'), 'cluster', 'protein_id')
    inventory = {r['accession']: r['name'] for r in read_tsv(out/'domains/model_inventory.tsv')}
    if set(nr) != {c+'.faa' for c in members} or set(domains) != set(members):
        raise ValueError('NR/domain files do not match the current MCL membership')
    families = out/'families'
    families.mkdir(exist_ok=True)
    for stale in families.glob('cluster_*'):
        if stale.is_dir() and stale.name not in members and re.fullmatch('cluster_[0-9]+', stale.name):
            shutil.rmtree(stale)
    settings = {k: getattr(args, k) for k in DEFAULTS}
    sequence_software = {k: v for k, v in manifest['software'].items()
                         if k in ['mafft', 'fasttree', 'pyhmmer', 'biopython']}
    summaries = []
    combined = out/'family_profiles.hmm'
    combined.unlink(missing_ok=True)
    for i, (cluster, membership) in enumerate(members.items(), 1):
        LOG.info('Family %d/%d: %s (%d proteins)', i, len(members), cluster, len(membership))
        folder = families/cluster
        folder.mkdir(exist_ok=True)
        records = {r.id: str(r.seq) for r in SeqIO.parse(out/'clustering/clusters'/f'{cluster}.faa', 'fasta')}
        if set(records) != set(membership):
            raise ValueError(f'{cluster}: FASTA and membership differ')
        rows, columns, counts = cluster_table(cluster, records, nr[cluster+'.faa'], domains[cluster], inventory, folder)
        result = analyze_family(cluster, records, folder, settings, sequence_software, tools['mafft'],
                                tools.get('fasttree'), args.threads, args.resume and not args.force,
                                not args.skip_hmm, not args.skip_tree)
        result.update(counts)
        if result['distinct_domain_architectures'] > 1:
            result['review_flags'] = 'heterogeneous_domain_architectures'
        else:
            result['review_flags'] = ''
        if result['hmm_status'] != 'built_uncalibrated':
            result['review_flags'] += (';' if result['review_flags'] else '') + result['hmm_status']
        if result['tree_status'] == 'unresolved_membership_only':
            result['review_flags'] += (';' if result['review_flags'] else '') + 'tree_not_inferred'
        family_report(folder, cluster, rows, columns, result, args.domain_db)
        summaries.append(result)
        if (folder/'profile.hmm').is_file():
            with combined.open('ab') as dest, (folder/'profile.hmm').open('rb') as source:
                shutil.copyfileobj(source, dest)
    write_tsv(out/'cluster_summary.tsv', list(summaries[0]), summaries)
    manifest['family_analysis'] = dict(summary='cluster_summary.tsv', families=len(summaries),
        proteins=sum(r['proteins'] for r in summaries),
        hmms_built=sum(r['hmm_status'] == 'built_uncalibrated' for r in summaries),
        trees_inferred=sum(r['tree_status'] == 'fasttree_approximate_ml_unrooted' for r in summaries),
        sequence_caches_reused=sum(r['cache_reused'] for r in summaries),
        manifests={r['cluster']: sha256(families/r['cluster']/'sequence_manifest.json') for r in summaries})
    manifest['skipped_analyses'] = [dict(cluster=r['cluster'], hmm=r['hmm_status'], tree=r['tree_status'])
                                  for r in summaries if r['hmm_status'] != 'built_uncalibrated'
                                  or r['tree_status'] != 'fasttree_approximate_ml_unrooted']
    header = ['Cluster', 'Proteins', 'Domains', 'Curated sequences', 'Curated columns', 'HMM', 'Tree', 'Review flags']
    body = ''.join('<tr><td><a href="families/'+r['cluster']+'/index.html">'+r['cluster']+'</a></td>'
                   + ''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in
                             ['proteins', 'domains', 'hmm_training_sequences', 'curated_columns',
                              'hmm_status', 'tree_status', 'review_flags'])+'</tr>' for r in summaries)
    (out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Protein families</title>'
        + STYLE + '<h1>Protein family analysis</h1><p><a href="cluster_summary.tsv">Cluster summary</a> | '
        '<a href="nr_annotated/annotations.tsv">Detailed NR hits</a> | <a href="domains/index.html">Domain profiles</a> | '
        '<a href="clustering/index.html">Protein similarity networks</a> | <a href="manifest.json">Run manifest</a>'
        + (' | <a href="family_profiles.hmm">Combined uncalibrated family HMMs</a>' if combined.is_file() else '')
        + '</p><p class="note">Click a cluster for its protein-by-domain table, NR descriptions, alignments, curation audit, '
        'HMM and tree. Automated curation needs review; small or degenerate families have explicitly unresolved trees '
        'and may lack an HMM. Protein similarity, domain co-occurrence and NR labels do not establish EET function.</p>'
        '<input aria-label="Search families" placeholder="Search clusters or analysis status"><div class="scroll"><table><thead><tr>'
        + ''.join('<th>'+h+'</th>' for h in header)+'</tr></thead><tbody>'+body+'</tbody></table></div>'+SEARCH+'</html>')


def run(args):
    LOG.info('Checking protein FASTA, reference databases, and executable dependencies')
    tools, total = preflight(args)
    log = prepare_output(args.output, args.resume, args.force)
    manifest = base_manifest(args)
    manifest.update(workflow='protein_families', stages={}, input_sha256=sha256(args.input),
                    input_proteins=total, software=versions(['biopython', 'pyhmmer']),
                    skipped_analyses=[], tool_paths=tools)
    try:
        for name, flag in [('diamond', 'version'), ('mcl', '--version'), ('mcxload', '--version'),
                           ('mafft', '--version'), ('fasttree', '-help')]:
            if name in tools:
                manifest['software'][name] = version(tools[name], flag).splitlines()[0]
        git = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True)
        manifest['repository_commit'] = git.stdout.strip() if git.returncode == 0 else 'unavailable'
        manifest['interpretation'] = 'NR best-reported-hit descriptions and detected domains are evidence, not confirmed function. HMM profiles are uncalibrated family models. Trees use full alignments; HMM training uses separately audited curation.'
        (args.output/'logs').mkdir(exist_ok=True)
        dump(args.output/'manifest.inprogress.json', manifest)
        (args.output/'manifest.failed.json').unlink(missing_ok=True)
        clustering, annotated, domains = [args.output/p for p in ['clustering', 'nr_annotated', 'domains']]
        scratch = ['--tmpdir', args.tmpdir] if args.tmpdir else []
        child('cluster_proteins.py', ['-i', args.input, '-o', clustering, '--threads', args.threads,
              '--inflation', args.inflation, '--min-identity', args.min_identity, '--min-coverage', args.min_coverage,
              '--evalue', args.evalue, '--block-size', args.block_size, *scratch], clustering, args, manifest)
        child('annotate_clusters_nr.py', ['--input', clustering/'clusters', '--output', annotated,
              '--db', args.nr_db, '--threads', args.threads, '--evalue', args.nr_evalue,
              '--min-identity', args.nr_min_identity, '--query-cover', args.nr_query_cover,
              '--subject-cover', args.nr_subject_cover, '--max-target-seqs', args.nr_max_target_seqs,
              '--block-size', args.block_size, *scratch, *(['--hash-db'] if args.hash_nr else [])],
              annotated, args, manifest)
        domain_args = ['--input', clustering/'clusters', '--output', domains, '--threads', args.threads,
                       '--min-model-coverage', args.min_model_coverage, '--clustering', clustering]
        if args.domain_db == 'pfam':
            child('cluster_domains.py', [*domain_args, '--hmm', args.hmm,
                  *(['--domain-evalue', args.domain_evalue] if args.domain_evalue is not None else [])],
                  domains, args, manifest)
        else:
            child('cluster_cdd.py', [*domain_args, '--cdd', args.cdd, '--evalue', args.cdd_evalue], domains, args, manifest)
        command([sys.executable, ROOT/'network_reports.py', '--clustering', clustering, '--domains', domains,
                 '--protein-fastas', annotated], args.output/'logs/network_reports.log', manifest)
        # Search manifests remain immutable; record the later presentation refresh separately.
        manifest['network_reports'] = {
            p.name: dict(manifest=str((p/'network_report_manifest.json').resolve()),
                         sha256=sha256(p/'network_report_manifest.json'),
                         current_html_sha256={str(h.relative_to(p)): sha256(h) for h in p.rglob('*.html')})
            for p in [clustering, domains]}
        reports(args, tools, manifest)
        manifest['output_sha256'] = {str(p.relative_to(args.output)): sha256(p) for p in
                                    [args.output/'cluster_summary.tsv', args.output/'index.html',
                                     *sorted((args.output/'families').glob('*/proteins.tsv')),
                                     *sorted(args.output.glob('family_profiles.hmm'))]}
        finish(args.output, manifest)
        LOG.info('Complete: %s', args.output/'index.html')
    except Exception as exc:
        manifest['error'] = str(exc)
        dump(args.output/'manifest.failed.json', manifest)
        raise
    finally:
        LOG.removeHandler(log)
        log.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', required=True, type=Path, help='Protein FASTA (.faa or .faa.gz); globally unique IDs')
    p.add_argument('-o', '--output', required=True, type=Path)
    p.add_argument('-t', '--threads', type=int, default=8)
    p.add_argument('--inflation', type=float, default=2.0)
    p.add_argument('--min-identity', type=float, default=35, help='Clustering alignment identity percent')
    p.add_argument('--min-coverage', type=float, default=0.5, help='Clustering coverage fraction on BOTH proteins')
    p.add_argument('--evalue', type=float, default=1e-5, help='Clustering search E-value')
    p.add_argument('--block-size', type=float, default=0.4, help='DIAMOND block size in billions of letters, both searches')
    p.add_argument('--tmpdir', '--diamond-tmpdir', type=Path, help='DIAMOND scratch on a filesystem with sufficient free space')
    p.add_argument('--nr-db', type=Path, default=Path('/home/ark/databases/nr.dmnd'))
    p.add_argument('--nr-evalue', type=float, default=1e-5)
    p.add_argument('--nr-min-identity', type=float, default=0)
    p.add_argument('--nr-query-cover', type=float, default=0, help='Percent, independent of clustering coverage')
    p.add_argument('--nr-subject-cover', type=float, default=0, help='Percent')
    p.add_argument('--nr-max-target-seqs', type=int, default=25, help='Reported targets; 0 for all. Best qualifying reported hit is retained')
    p.add_argument('--hash-nr', action='store_true', help='Hash the entire NR database instead of using the default stat fingerprint')
    p.add_argument('--domain-db', choices=['pfam', 'cdd'], default='pfam')
    p.add_argument('--hmm', type=Path, default=Path('/home/ark/databases/HMMs/Pfam-A.hmm'))
    p.add_argument('--cdd', type=Path, default=Path('/home/ark/databases/CDD'))
    p.add_argument('--domain-evalue', type=float, help='Explicit fallback only for Pfam/custom models without gathering cutoffs')
    p.add_argument('--cdd-evalue', type=float, default=0.01)
    p.add_argument('--min-model-coverage', type=float, default=0)
    p.add_argument('--mafft', default='mafft', help='Executable or path')
    p.add_argument('--fasttree', default='FastTree', help='Executable or path; fasttree fallback is recognized')
    p.add_argument('--min-column-occupancy', type=float, default=DEFAULTS['min_column_occupancy'], help='Standard-AA fraction for HMM columns, after exact-sequence deduplication')
    p.add_argument('--min-sequence-coverage', type=float, default=DEFAULTS['min_sequence_coverage'], help='Standard-AA coverage of the current HMM core for each training sequence')
    p.add_argument('--min-hmm-sequences', type=int, default=DEFAULTS['min_hmm_sequences'], help='Minimum unique retained sequences')
    p.add_argument('--min-hmm-columns', type=int, default=DEFAULTS['min_hmm_columns'])
    p.add_argument('--seed', type=int, default=DEFAULTS['seed'])
    p.add_argument('--skip-hmm', action='store_true', help='Still write and audit the curated alignment')
    p.add_argument('--skip-tree', action='store_true', help='FastTree is then optional')
    g = p.add_mutually_exclusive_group()
    g.add_argument('--resume', action='store_true', help='Reuse matching searches and family sequence products after checksum checks')
    g.add_argument('--force', action='store_true', help='Recompute generated products in an existing output directory')
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    try:
        run(args)
    except (ValueError, OSError, KeyError) as exc:
        LOG.error('%s', exc)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
