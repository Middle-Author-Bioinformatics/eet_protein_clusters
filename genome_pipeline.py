#!/usr/bin/env python3
"""Genomes -> whole-proteome MCL clusters -> domain profiles -> localization and motifs."""
import argparse
from collections import Counter
import csv
import html
import json
import logging
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
from Bio import SeqIO
from common import (LOG, base_manifest, command, dump, finish, mark_stage, prepare_output,
                    sha256, signature, stage_reusable, versions, writer)
from genome_inputs import MAP_FIELDS, discover, gram_assignments, prepare
from localization import LOCALIZATION_FIELDS, TOOLS, collect, infer, EXTENSIONS, EXECUTABLES
from interactive_network import options as network_options, validate as validate_network
from predictor_setup import configure
from protein_motifs import COUNT_FIELDS, PATTERNS, SOURCES, characterize

ROOT = Path(__file__).resolve().parent


def read_tsv(path):
    with open(path) as f: yield from csv.DictReader(f, delimiter='\t')


def load_config(path):
    if not path: return {}
    value = json.loads(path.read_text())
    if not isinstance(value, dict) or set(value)-set(TOOLS): raise ValueError('Localization config keys: psortb, signalp, phobius, lipop')
    for name, settings in value.items():
        if not isinstance(settings, dict) or set(settings)-{'argv', 'output', 'version', 'model_files', 'format', 'version_argv', 'stdin'}:
            raise ValueError(f'{name}: allowed config fields: argv, output, version, model_files, format, version_argv, stdin')
        if 'argv' in settings and (not isinstance(settings['argv'], list) or not settings['argv'] or not all(isinstance(a, str) for a in settings['argv'])):
            raise ValueError(f'{name}: argv must be a nonempty list of strings')
        if 'format' in settings and (name!='signalp' or settings['format'] not in ['signalp4', 'signalp6']):
            raise ValueError('format applies only to signalp and must be signalp4 or signalp6')
        if 'stdin' in settings and not isinstance(settings['stdin'], bool): raise ValueError(f'{name}: stdin must be boolean')
        if 'version_argv' in settings and (not isinstance(settings['version_argv'], list) or not settings['version_argv'] or not all(isinstance(a,str) for a in settings['version_argv'])):
            raise ValueError(f'{name}: version_argv must be a nonempty argument list')
        for field in ['output', 'version']:
            if field in settings and not isinstance(settings[field], str): raise ValueError(f'{name}: {field} must be a string')
        if not isinstance(settings.get('model_files', []), list) or not all(isinstance(p, str) for p in settings.get('model_files', [])):
            raise ValueError(f'{name}: model_files must be a list of file paths')
        settings['model_files'] = [str((path.parent/Path(p)).resolve()) for p in settings.get('model_files', [])]
    return value


def preflight(args, files, config):
    validate_network(args)
    import math
    if args.hmm_metadata and not args.hmm_metadata.is_file(): raise ValueError(f'Missing HMM metadata file: {args.hmm_metadata}')
    if (args.threads < 1 or not all(math.isfinite(x) for x in [args.inflation, args.min_identity, args.min_coverage, args.evalue, args.block_size])
        or args.block_size <= 0 or args.inflation <= 1 or not 0 <= args.min_coverage <= 1 or not 0 <= args.min_identity <= 100 or not 0 < args.evalue < 1):
        raise ValueError('Require threads>=1, inflation>1, coverage 0..1, identity 0..100 and 0<evalue<1')
    if args.domain_evalue is not None and (not math.isfinite(args.domain_evalue) or args.domain_evalue <= 0):
        raise ValueError('--domain-evalue must be positive and finite')
    if not math.isfinite(args.min_model_coverage) or not 0 <= args.min_model_coverage <= 1:
        raise ValueError('--min-model-coverage must be 0..1')
    for tool in ['diamond', 'mcl', 'mcxload']:
        if not shutil.which(tool): raise ValueError(f'{tool} not found in PATH. conda env create -f environment.yml; conda activate eet-clusters')
    if any(kind == 'assembly' for _, _, kind in files):
        try: import pyrodigal
        except ImportError as exc: raise ValueError('Install gene prediction: python -m pip install -r requirements-genomes.txt') from exc
    # Validate all model thresholds before spending time on all-versus-all alignment.
    from cluster_domains import model_inventory
    with sqlite3.connect(':memory:') as db: model_inventory(args.hmm, db, args.domain_evalue)
    for tool in TOOLS:
        if getattr(args, tool) != 'run': continue
        all_imported = args.localization_import and all((args.localization_import/f'{g}.{EXTENSIONS[tool]}').is_file() for g, _, _ in files)
        if all_imported: continue
        if tool == 'lipop' and 'argv' not in config.get(tool, {}):
            raise ValueError('LipoP run requires installation-specific argv in --localization-config; see GENOME_WORKFLOW.md')
        default_executable = 'signalp' if tool=='signalp' and config.get(tool, {}).get('format')=='signalp4' else EXECUTABLES[tool]
        executable = config.get(tool, {}).get('argv', [default_executable])[0]
        if not shutil.which(executable):
            raise ValueError(f'{tool}: {executable} not found. Install/configure it, import native results, or select --{tool} skip. See GENOME_WORKFLOW.md.')


def child(script, argv, output, args, manifest):
    cmd = [sys.executable, str(ROOT/script), *map(str, argv)]
    if args.force: cmd.append('--force')
    elif args.resume: cmd.append('--resume')
    command(cmd, output, manifest)
    path = Path(argv[argv.index('-o')+1])/'manifest.json'
    manifest['stages'][script] = dict(manifest=str(path.resolve()), sha256=sha256(path),
                                    provenance=json.loads(path.read_text()))


def assemble_annotations(output, genomes, args, config, manifest):
    work = output/'work'; work.mkdir(exist_ok=True)
    db_path = work/'annotations.sqlite'; db_path.unlink(missing_ok=True)
    db = sqlite3.connect(db_path)
    try:
        db.execute('CREATE TABLE annotation (protein_id TEXT PRIMARY KEY, genome TEXT, cluster TEXT, metadata TEXT, domains TEXT, localization TEXT)')
        db.executemany('INSERT INTO annotation VALUES (?,?,NULL,?,NULL,NULL)',
                       ((r['protein_id'], r['genome'], json.dumps(r)) for r in read_tsv(output/'inputs/protein_mapping.tsv')))
        db.execute('CREATE INDEX genome_idx ON annotation(genome)')
        for row in read_tsv(output/'clustering/membership.tsv'):
            if db.execute('UPDATE annotation SET cluster=? WHERE protein_id=?', (row['cluster'], row['protein_id'])).rowcount != 1:
                raise ValueError('Clustering contains protein absent from genome mapping')
        for row in read_tsv(output/'domains/protein_architectures.tsv'):
            if db.execute('UPDATE annotation SET domains=? WHERE protein_id=?', (json.dumps(row), row['protein_id'])).rowcount != 1:
                raise ValueError('Domain output contains protein absent from genome mapping')
        db.commit()
        for genome in genomes:
            mapping = [json.loads(r[0]) for r in db.execute('SELECT metadata FROM annotation WHERE genome=?', (genome['genome'],))]
            evidence = collect(genome, mapping, args, config, output/'localization'/Path(genome['fasta']).stem, manifest)
            native_path = output/'localization'/Path(genome['fasta']).stem/'evidence.json'
            native_path.parent.mkdir(parents=True, exist_ok=True)
            dump(native_path, evidence)
            for row in mapping:
                pid = row['protein_id']
                raw = {tool: evidence[tool].get(pid, {}) for tool in TOOLS}
                db.execute('UPDATE annotation SET localization=? WHERE protein_id=?', (json.dumps(infer(raw, genome['gram'])), pid))
            db.commit()
        if db.execute('SELECT COUNT(*) FROM annotation WHERE cluster IS NULL OR domains IS NULL OR localization IS NULL').fetchone()[0]:
            raise ValueError('Incomplete cluster/domain/localization mapping; inspect stage outputs')
        return write_annotations(db, output, genomes, manifest)
    finally: db.close()


def write_annotations(db, output, genomes, manifest):
    fields = MAP_FIELDS + ['cluster'] + COUNT_FIELDS + ['n_domain_hits', 'n_distinct_domains', 'architecture', 'domain_set', 'architecture_names', 'domain_names'] + LOCALIZATION_FIELDS
    f, w = writer(output/'protein_annotations.tsv', fields)
    m, mw = writer(output/'motif_hits.tsv', ['protein_id', 'genome', 'cluster', 'motif', 'regex', 'start', 'end', 'match'])
    counts = {g['genome']: Counter() for g in genomes}; cluster_counts = {}; preview = []; distributions = Counter()
    try:
        for rec in SeqIO.parse(output/'inputs/proteins.faa', 'fasta'):
            genome, cluster, meta, domains, loc = db.execute('SELECT genome,cluster,metadata,domains,localization FROM annotation WHERE protein_id=?', (rec.id,)).fetchone()
            row = json.loads(meta); row.update(json.loads(domains)); row.update(json.loads(loc)); row['cluster'] = cluster
            motifs, hits = characterize(rec.seq); row.update(motifs); w.writerow(row)
            for hit in hits: mw.writerow(dict(hit, protein_id=rec.id, genome=genome, cluster=cluster))
            for c in [counts[genome], cluster_counts.setdefault(cluster, Counter())]:
                c['proteins'] += 1
                c['proteins_with_domains'] += int(int(row['n_domain_hits']) > 0)
                c['canonical_mhc_ge2'] += int(row['cxxch_count'] >= 2)
                c['hematite_motif_proteins'] += int(row['hematite_count'] > 0)
                c['canonical_heme_motifs'] += row['cxxch_count']; c['hematite_motifs'] += row['hematite_count']
                c['psortb_assessed'] += row['psortb_status']=='assessed'
                c['psortb_unknown'] += row['psortb_localization'].lower()=='unknown'
                c['psortb_unavailable'] += row['psortb_status']=='unavailable'
                c['unknown_with_feature_refinement'] += row['psortb_localization'].lower()=='unknown' and row['localization_interpretation']!='Unknown'
                c['localization_conflicts'] += bool(row['localization_conflict'])
                c['signalp_assessed'] += row['signalp_status']=='assessed'
                c['phobius_assessed'] += row['phobius_status']=='assessed'
                c['lipop_assessed'] += row['lipop_status']=='assessed'
            for scope, group in [('genome', genome), ('cluster', cluster)]:
                for source, field in [('SignalP', 'signalp_prediction'), ('Phobius_SP', 'phobius_signal_peptide'), ('Phobius_TM_count', 'phobius_tm_count'), ('feature_interpretation', 'localization_interpretation')]+([('PSORTb', 'psortb_localization')] if row['psortb_status']=='assessed' else []):
                    distributions[scope, group, source, row[field]] += 1
            if len(preview)<1000: preview.append(row)
    finally: f.close(); m.close()
    cf = ['proteins', 'proteins_with_domains', 'canonical_mhc_ge2', 'canonical_heme_motifs', 'hematite_motif_proteins', 'hematite_motifs',
          'psortb_assessed', 'psortb_unknown', 'psortb_unavailable', 'unknown_with_feature_refinement',
          'localization_conflicts', 'signalp_assessed', 'phobius_assessed', 'lipop_assessed']
    f, w = writer(output/'genome_summary.tsv', ['genome', 'input_type', 'gram']+cf)
    for genome in genomes: w.writerow(dict(genome, **counts[genome['genome']]))
    f.close()
    f, w = writer(output/'cluster_evidence_summary.tsv', ['cluster']+cf)
    for cluster, count in cluster_counts.items(): w.writerow(dict(count, cluster=cluster))
    f.close()
    f, w = writer(output/'localization_distribution.tsv', ['scope', 'group_id', 'source', 'category', 'proteins', 'total_proteins', 'fraction_proteins'])
    for (scope, group, source, category), count in sorted(distributions.items()):
        total = (counts if scope=='genome' else cluster_counts)[group]['proteins']
        w.writerow(dict(scope=scope, group_id=group, source=source, category=category, proteins=count, total_proteins=total, fraction_proteins=count/total))
    f.close()
    manifest['counts'] = dict(genomes=len(genomes), proteins=sum(c['proteins'] for c in counts.values()),
                              clusters=len(cluster_counts), canonical_mhc_ge2=sum(c['canonical_mhc_ge2'] for c in counts.values()),
                              hematite_motif_proteins=sum(c['hematite_motif_proteins'] for c in counts.values()))
    return preview


def report(output, preview, manifest):
    esc = lambda v: html.escape(str(v), quote=True)
    text = '<!doctype html><html lang="en"><meta charset="utf-8"><title>Genome protein analysis</title>'
    text += '<style>body{font:15px system-ui;margin:2em;color:#173344}table{border-collapse:collapse;font-size:13px}td,th{padding:6px;border:1px solid #ccd}th{background:#e8f0f2}input{padding:8px;width:50%}.scroll{overflow:auto;max-height:65vh}a{color:#006b7a}</style>'
    text += '<h1>Genome protein analysis</h1><p>'+esc(manifest['counts'])+'</p>'
    text += '<p>Whole proteomes were clustered. Sequence patterns are candidate features; motif counts are not measured heme stoichiometry, binding or EET activity. Signal peptides support export, not a final destination or extracellular exposure. SignalP and Phobius results are preserved separately; interpretations describe export and membrane attachment. Negative and unassessed predictions remain distinct.</p>'
    links = [('protein_annotations.tsv', 'All protein annotations'), ('genome_summary.tsv', 'Genome summary'),
             ('cluster_evidence_summary.tsv', 'Cluster evidence summary'), ('localization_distribution.tsv', 'Localization distributions'), ('motif_hits.tsv', 'Motif coordinates'),
             ('inputs/protein_mapping.tsv', 'Original CDS mappings'), ('inputs/proteins.faa', 'All proteins'),
             ('clustering/cluster_summary.tsv', 'Cluster FASTA index'), ('clustering/membership.tsv', 'Cluster membership'),
             ('domains/index.html', 'Per-cluster domains and networks'), ('domains/domain_hits.tsv', 'Domain hits'),
             ('domains/domain_prevalence.tsv', 'Domain prevalence'), ('manifest.json', 'Run manifest')]
    if (output/'clustering/protein_network.html').exists():links.insert(0,('clustering/protein_network.html','Interactive protein similarity network'))
    if (output/'domains/all_domains_network.html').exists():links.insert(1,('domains/all_domains_network.html','Interactive domain co-occurrence network'))
    text += '<p>'+' | '.join(f'<a href="{p}">{label}</a>' for p, label in links)+'</p>'
    skipped = manifest['skipped_analyses']
    if skipped: text += '<h2>Unavailable analyses</h2><ul>'+''.join('<li>'+esc(f'{s["genome"]}: {s["tool"]}: {s["reason"]}')+'</li>' for s in skipped)+'</ul>'
    def table(title, rows, columns, label):
        rows = list(rows)
        return '<h2>'+esc(title)+'</h2><p>'+esc(label)+'</p><input placeholder="Search this table" oninput="filter(this)"><div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(c)+'</th>' for c in columns)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(row.get(c, ''))+'</td>' for c in columns)+'</tr>' for row in rows)+'</tbody></table></div>'
    text += table('Genomes', read_tsv(output/'genome_summary.tsv'), ['genome', 'gram', 'proteins', 'canonical_mhc_ge2', 'hematite_motif_proteins', 'signalp_assessed', 'phobius_assessed', 'localization_conflicts'], 'All input genomes')
    text += table('Proteins', preview, ['protein_id', 'genome', 'original_id', 'cluster', 'cxxch_count', 'mhcscan_hemes', 'fegenie_heme_c_total', 'heme_b_count', 'hematite_count', 'domain_names', 'architecture_names', 'signalp_prediction', 'signalp_generation', 'signalp_cleavage', 'signalp_d_score', 'phobius_signal_peptide', 'phobius_tm_count', 'localization_interpretation', 'localization_conflict'], 'Preview: first 1,000 proteins. Download protein_annotations.tsv for every protein; links above contain complete cluster/domain results.')
    text += '<script>function filter(i){let q=i.value.toLowerCase();i.nextElementSibling.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q))}</script></html>'
    (output/'index.html').write_text(text)


def run(args):
    args.input = args.input.resolve(); args.output = args.output.resolve(); args.hmm = args.hmm.resolve()
    if not args.input.is_dir(): raise ValueError('--input must be a folder of genomes')
    if args.output == args.input or args.output in args.input.parents or args.input in args.output.parents:
        raise ValueError('Keep input and output directories separate (neither may contain the other)')
    files = discover(args.input, args.input_type); grams = gram_assignments(files, args.gram, args.genome_metadata)
    config = configure(args, load_config(args.localization_config)); preflight(args, files, config)
    log = prepare_output(args.output, args.resume, args.force)
    manifest = base_manifest(args)
    manifest.update(software=versions(['biopython', 'pyhmmer', 'pyrodigal', 'networkx', 'matplotlib']),
                    motif_sources=SOURCES, motif_patterns=PATTERNS, motif_counting='nonoverlapping_per_pattern; positions 1-based inclusive',
                    predictor_configuration=config, localization=[], skipped_analyses=[], stages={}, input_sha256={str(p):sha256(p) for _, p, _ in files})
    commit = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    manifest['hmm_database'] = dict(path=str(args.hmm), sha256=sha256(args.hmm), model_inventory='domains/model_inventory.tsv')
    if args.hmm_metadata:
        manifest['hmm_database']['release_metadata'] = dict(path=str(args.hmm_metadata.resolve()), sha256=sha256(args.hmm_metadata), text=args.hmm_metadata.read_text())
    dirty = subprocess.run(['git', '-C', str(ROOT), 'status', '--porcelain'], capture_output=True, text=True)
    manifest['repository_dirty'] = bool(dirty.stdout.strip()) if dirty.returncode==0 else 'unavailable'
    manifest['repository_commit'] = commit.stdout.strip() if commit.returncode==0 else 'unavailable'
    try:
        dump(args.output/'manifest.inprogress.json', manifest)
        inputs = args.output/'inputs'; inputs.mkdir(exist_ok=True)
        key = signature(dict(files=[(g, str(p), sha256(p), k) for g,p,k in files], grams=grams,
                             table=args.translation_table, code=sha256(ROOT/'genome_inputs.py'), software=manifest['software']['pyrodigal'],
                             biopython=manifest['software']['biopython']))
        stamp = inputs/'prepare.stage.json'; index = inputs/'genomes.json'
        paths = [inputs/'proteins.faa', inputs/'protein_mapping.tsv', index]
        old_genomes = json.loads(index.read_text()) if index.exists() else []
        paths += [Path(g['fasta']) for g in old_genomes]
        reused = stage_reusable(stamp, key, paths, args.resume and not args.force)
        if reused: genomes = old_genomes; LOG.info('Reusing extracted proteomes')
        else:
            LOG.info('Extracting whole proteomes from %d genomes', len(files))
            genomes = prepare(files, grams, inputs, args.translation_table); dump(index, genomes)
            paths = [inputs/'proteins.faa', inputs/'protein_mapping.tsv', index]+[Path(g['fasta']) for g in genomes]
            mark_stage(stamp, key, paths)
        manifest['stages']['gene_extraction'] = dict(cache_reused=reused, signature=key, genomes=genomes)
        network_extra=['--network-max-nodes',args.network_max_nodes,'--network-max-edges',args.network_max_edges,'--network-max-proteins',args.network_max_proteins]+(['--no-interactive'] if args.no_interactive else [])
        cluster_extra=(['--tmpdir',args.tmpdir] if args.tmpdir else [])+network_extra
        child('cluster_proteins.py', ['-i', inputs/'proteins.faa', '-o', args.output/'clustering', '-t', args.threads,
              '--inflation', args.inflation, '--min-identity', args.min_identity, '--min-coverage', args.min_coverage,
              '--evalue', args.evalue, '--sensitivity', args.sensitivity, '--block-size', args.block_size]+cluster_extra, args.output/'clustering.log', args, manifest)
        domain_args = ['-i', args.output/'clustering/clusters', '-o', args.output/'domains', '--hmm', args.hmm,
                       '-t', args.threads, '--min-model-coverage', args.min_model_coverage, '--overlap-policy', args.overlap_policy]
        domain_args+=network_extra+['--clustering',args.output/'clustering']
        if args.domain_evalue is not None: domain_args += ['--domain-evalue', args.domain_evalue]
        if args.no_plots: domain_args.append('--no-plots')
        if args.low_memory: domain_args.append('--low-memory')
        child('cluster_domains.py', domain_args, args.output/'domains.log', args, manifest)
        if not args.no_interactive:
            from network_reports import refresh_proteins
            manifest['interactive_protein_network']=refresh_proteins(args.output/'clustering',args,args.output/'domains')
        dump(args.output/'manifest.inprogress.json', manifest)
        preview = assemble_annotations(args.output, genomes, args, config, manifest)
        report(args.output, preview, manifest)
        manifest['output_sha256'] = {p:sha256(args.output/p) for p in ['protein_annotations.tsv', 'motif_hits.tsv', 'genome_summary.tsv', 'cluster_evidence_summary.tsv', 'localization_distribution.tsv', 'index.html']}
        finish(args.output, manifest); LOG.info('Finished: %s', args.output/'index.html')
    except Exception as exc:
        manifest['error'] = str(exc); dump(args.output/'manifest.inprogress.json', manifest)
        raise
    finally: LOG.removeHandler(log); log.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--input', type=Path, required=True, help='Genome folder (one genome per file)')
    p.add_argument('-o', '--output', type=Path, required=True)
    p.add_argument('--hmm', type=Path, required=True, help='Uncompressed Pfam-A/custom HMM database')
    p.add_argument('--hmm-metadata', type=Path, help='Optional Pfam.version.txt or custom database release metadata, recorded verbatim with checksum')
    p.add_argument('--input-type', choices=['auto', 'assembly', 'genbank', 'proteins'], default='auto')
    p.add_argument('--gram', choices=['positive', 'negative'], help='Gram group for SignalP 4 (and optional PSORTb); metadata can override per genome')
    p.add_argument('--genome-metadata', type=Path, help='TSV genome (filename stem), gram (positive/negative)')
    p.add_argument('--translation-table', type=int, choices=[4,11], default=11, help='Assembly prediction table; default for GenBank CDS without transl_table')
    p.add_argument('-t', '--threads', type=int, default=8)
    p.add_argument('--inflation', type=float, default=2.)
    p.add_argument('--min-identity', type=float, default=30.)
    p.add_argument('--min-coverage', type=float, default=.5)
    p.add_argument('--evalue', type=float, default=1e-5)
    p.add_argument('--sensitivity', choices=['sensitive', 'very-sensitive', 'ultra-sensitive'], default='very-sensitive')
    p.add_argument('--block-size', type=float, default=2.)
    p.add_argument('--tmpdir','--diamond-tmpdir',dest='tmpdir',type=Path,help='DIAMOND scratch parent, forwarded to clustering; default $TMPDIR or clustering/work')
    p.add_argument('--domain-evalue', type=float, help='Explicit fallback for HMMs without GA; Pfam GA preferred')
    p.add_argument('--min-model-coverage', type=float, default=0.)
    p.add_argument('--overlap-policy', choices=['best', 'all'], default='best')
    network_options(p)
    p.add_argument('--low-memory', action='store_true')
    p.add_argument('--no-plots', action='store_true')
    for tool in TOOLS:
        p.add_argument('--'+tool, choices=['run', 'auto', 'skip'], default='skip' if tool in ['psortb', 'lipop'] else 'auto',
                       help='run: require executable or import; auto: run if configured/installed, otherwise mark unavailable; skip: omit explicitly')
    p.add_argument('--localization-import', type=Path, help='Native result folder, <genome>.psortb.tsv/.signalp.tsv/.phobius.txt/.lipop.txt')
    p.add_argument('--localization-config', type=Path, help='JSON argv/version/model_files overrides; no shell execution')
    p.add_argument('--predictor-cache-tag', default='', help='Change after updating models not covered by config model_files')
    p.add_argument('--signalp-version', choices=['4', '6'], default='4', help='SignalP generation for default execution; imports are detected from native format')
    p.add_argument('--signalp4-env', help='Existing SignalP 4 Conda environment name or prefix; executed with conda run, no activation needed')
    p.add_argument('--conda-executable', default='conda', help='Conda binary for --signalp4-env')
    p.add_argument('--phobius-dir', type=Path, help='Existing MHCscan bin folder containing phobius.pl, decoder, model and options')
    p.add_argument('--signalp-mode', choices=['fast', 'slow-sequential'], default='fast')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--force', action='store_true')
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    try: run(args)
    except (ValueError, OSError, sqlite3.Error, ImportError, subprocess.SubprocessError) as exc:
        LOG.error('%s', exc); return 2
    return 0


if __name__ == '__main__': sys.exit(main())
