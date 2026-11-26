"""Auditable protein alignments, conservative profile-HMM curation, and trees.

HMM training and phylogeny deliberately use different sequence sets: the full
alignment retains every member, whereas HMM training excludes exact duplicates
and low-core-coverage fragments. No EET function or detection cutoff is inferred.
"""
import csv
from collections import Counter
import html
import json
import os
from pathlib import Path
import re
import subprocess

from Bio import Phylo, SeqIO
from common import LOG, dump, sha256, signature, writer

AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
DEFAULTS = dict(min_column_occupancy=0.5, min_sequence_coverage=0.5,
                min_hmm_sequences=3, min_hmm_columns=20, seed=17)


def read_tsv(path):
    with open(path, newline='') as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def write_tsv(path, fields, rows):
    stream, table = writer(path, fields)
    try:
        table.writerows(rows)
    finally:
        stream.close()


def fasta(path, rows):
    with open(path, 'w') as stream:
        for name, seq in rows.items():
            stream.write(f'>{name}\n{seq}\n')


def validate_alignment(path, expected):
    rows = {}
    for record in SeqIO.parse(path, 'fasta'):
        if record.id in rows:
            raise ValueError(f'{path}: repeated aligned ID {record.id}')
        rows[record.id] = str(record.seq).upper()
    if set(rows) != set(expected) or len({len(s) for s in rows.values()}) != 1:
        raise ValueError(f'{path}: alignment IDs or lengths do not match the input')
    for name, seq in rows.items():
        if seq.replace('-', '') != expected[name].upper():
            raise ValueError(f'{path}: alignment changed residues for {name}')
    return {name: rows[name] for name in expected}


def curate(rows, settings):
    """Monotonically filter low-occupancy columns and low-core-coverage rows.

    Occupancy uses the 20 standard amino acids; ambiguous residues are missing
    evidence. Exact duplicates have one representative during HMM curation.
    All coordinates in the audit are 1-based columns of the full alignment.
    """
    length = len(next(iter(rows.values())))
    representatives, seen, duplicate_of = [], {}, {}
    for name, seq in rows.items():
        ungapped = seq.replace('-', '')
        if ungapped in seen:
            duplicate_of[name] = seen[ungapped]
        else:
            seen[ungapped] = name
            representatives.append(name)
    initial = {i: sum(rows[n][i] in AA for n in representatives) / len(representatives)
               for i in range(length)}
    members, columns = representatives[:], list(range(length))
    excluded = {}
    rounds = 0
    while members and columns:
        rounds += 1
        kept = [i for i in columns if sum(rows[n][i] in AA for n in members) / len(members)
                >= settings['min_column_occupancy']]
        if not kept:
            columns = []
            break
        survivors = []
        for name in members:
            coverage = sum(rows[name][i] in AA for i in kept) / len(kept)
            if coverage >= settings['min_sequence_coverage']:
                survivors.append(name)
            else:
                excluded[name] = dict(reason='low_core_coverage', coverage=coverage, round=rounds)
        changed = kept != columns or survivors != members
        columns, members = kept, survivors
        if not changed:
            break
    if not members:
        columns = []
    curated = {n: ''.join(rows[n][i] for i in columns) for n in members} if columns else {}
    seq_audit = []
    for name, seq in rows.items():
        representative = duplicate_of.get(name, name)
        state = ('exact_duplicate' if name in duplicate_of else
                 ('included' if name in curated else excluded.get(name, {}).get('reason', 'no_core_columns')))
        seq_audit.append(dict(protein_id=name, hmm_training_status=state,
                              representative=representative,
                              standard_aa=sum(c in AA for c in seq),
                              ambiguous_aa=sum(c not in AA and c != '-' for c in seq),
                              final_core_coverage=(sum(seq[i] in AA for i in columns) / len(columns)
                                                   if columns else ''),
                              coverage_at_exclusion=excluded.get(name, {}).get('coverage', ''),
                              exclusion_round=excluded.get(name, {}).get('round', '')))
    kept_columns = set(columns)
    col_audit = [dict(full_alignment_column=i+1, retained=int(i in kept_columns),
                      initial_unique_sequence_occupancy=initial[i],
                      final_training_sequence_occupancy=(sum(rows[n][i] in AA for n in members) / len(members)
                                                         if members else 0)) for i in range(length)]
    return curated, seq_audit, col_audit, rounds


def run_tool(argv, output, log, manifest, threads):
    manifest['commands'].append(list(map(str, argv)))
    LOG.info('Running %s for %s', Path(argv[0]).name, output.parent.name)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(output, 'w') as dest, open(log, 'w') as diagnostics:
        result = subprocess.run(list(map(str, argv)), stdout=dest, stderr=diagnostics, env=env)
    if result.returncode:
        tail = '\n'.join(log.read_text(errors='replace').splitlines()[-12:])
        raise ValueError(f'{argv[0]} exited {result.returncode}. See {log}\n{tail}')


def named_newick(tree, mapping):
    """Write only observed branch lengths; quote original IDs safely."""
    def clade(node):
        if node.clades:
            text = '(' + ','.join(clade(c) for c in node.clades) + ')'
            if node.confidence is not None:
                text += format(node.confidence, '.12g')
        else:
            text = "'" + mapping[node.name].replace("'", "''") + "'"
        if node.branch_length is not None:
            text += ':' + format(node.branch_length, '.12g')
        return text
    return clade(tree.root) + ';\n'


def check_tree(path, expected):
    tree = Phylo.read(path, 'newick')
    leaves = [n.name for n in tree.get_terminals()]
    if Counter(leaves) != Counter(list(expected)):
        raise ValueError(f'{path}: tree tips do not match all cluster protein IDs')
    return tree


def analyze_family(cluster, records, folder, settings, software, mafft, fasttree,
                   threads=1, resume=False, build_hmm=True, build_tree=True):
    """Build sequence-only products. Table/annotation changes do not invalidate them."""
    import pyhmmer as ph
    settings = dict(DEFAULTS, **settings)
    folder.mkdir(parents=True, exist_ok=True)
    stamp = folder/'sequence_manifest.json'
    key = signature(dict(cluster=cluster, sequences=list(records.items()), settings=settings,
                         software=software, module=sha256(Path(__file__)),
                         threads=threads, build_hmm=build_hmm, build_tree=build_tree))
    if resume and stamp.is_file():
        try:
            previous = json.loads(stamp.read_text())
            if previous.get('key') == key and all((folder/p).is_file() and sha256(folder/p) == digest
                                                 for p, digest in previous['output_sha256'].items()):
                return dict(previous['result'], cache_reused=True)
        except (ValueError, KeyError):
            LOG.info('Rebuilding %s after an invalid sequence cache manifest', cluster)
    stamp.unlink(missing_ok=True)
    # Only these generated products are replaced. Matrices and user-added files survive.
    products = ['sequences.safe.faa', 'id_mapping.tsv', 'alignment.safe.faa', 'alignment.full.faa',
                'alignment.curated.faa', 'curation_sequences.tsv', 'curation_columns.tsv',
                'profile.hmm', 'tree.safe.nwk', 'tree.nwk', 'mafft.log', 'fasttree.log']
    for name in products:
        (folder/name).unlink(missing_ok=True)
    manifest = dict(key=key, cluster=cluster, software=software, settings=settings, commands=[],
                    input_sequence_sha256=signature(list(records.items())), skipped_analyses=[])
    dump(folder/'sequence_manifest.inprogress.json', manifest)
    safe = {f't{i:07d}': seq for i, seq in enumerate(records.values(), 1)}
    mapping = dict(zip(safe, records))
    fasta(folder/'sequences.safe.faa', safe)
    write_tsv(folder/'id_mapping.tsv', ['tree_alignment_id', 'protein_id'],
              [dict(tree_alignment_id=k, protein_id=v) for k, v in mapping.items()])
    if len(records) == 1:
        fasta(folder/'alignment.safe.faa', safe)
        alignment_status = 'singleton_sequence_copy'
    else:
        run_tool([mafft, '--auto', '--amino', '--anysymbol', '--thread', str(threads),
                  '--threadit', '0', folder/'sequences.safe.faa'], folder/'alignment.safe.faa',
                 folder/'mafft.log', manifest, threads)
        alignment_status = 'mafft_auto'
    aligned_safe = validate_alignment(folder/'alignment.safe.faa', safe)
    aligned = {mapping[k]: seq for k, seq in aligned_safe.items()}
    fasta(folder/'alignment.full.faa', aligned)
    curated, seq_audit, col_audit, rounds = curate(aligned, settings)
    fasta(folder/'alignment.curated.faa', curated)
    write_tsv(folder/'curation_sequences.tsv', list(seq_audit[0]), seq_audit)
    write_tsv(folder/'curation_columns.tsv', list(col_audit[0]), col_audit)
    result = dict(cluster=cluster, proteins=len(records), unique_sequences=len(set(records.values())),
                  alignment_status=alignment_status, alignment_columns=len(next(iter(aligned.values()))),
                  hmm_training_sequences=len(curated), curated_columns=(len(next(iter(curated.values()))) if curated else 0),
                  curation_rounds=rounds, hmm_match_states=0, cache_reused=False)
    if not build_hmm:
        result['hmm_status'] = 'skipped_by_request'
    elif len(curated) < settings['min_hmm_sequences']:
        result['hmm_status'] = 'insufficient_unique_training_sequences'
    elif result['curated_columns'] < settings['min_hmm_columns']:
        result['hmm_status'] = 'insufficient_curated_columns'
    else:
        alphabet = ph.easel.Alphabet.amino()
        # Internal identifiers also isolate HMMER from original-ID punctuation.
        msa = ph.easel.TextMSA(name=cluster, sequences=[ph.easel.TextSequence(name=f'h{i}', sequence=s)
                                                       for i, s in enumerate(curated.values(), 1)])
        builder = ph.plan7.Builder(alphabet, weighting='pb', effective_number='entropy',
                                  symfrac=0.5, seed=settings['seed'])
        hmm, _, _ = builder.build_msa(msa.digitize(alphabet), ph.plan7.Background(alphabet))
        hmm.description = 'Uncalibrated protein-family profile; no EET specificity established'
        with open(folder/'profile.hmm', 'wb') as stream:
            hmm.write(stream)
        manifest['hmm_builder'] = dict(api='pyhmmer.plan7.Builder.build_msa', architecture='fast',
                                       weighting='pb', effective_number='entropy', symfrac=0.5,
                                       seed=settings['seed'], training_sequences=hmm.nseq,
                                       effective_sequences=hmm.nseq_effective, match_states=hmm.M,
                                       gathering_threshold=None)
        result.update(hmm_status='built_uncalibrated', hmm_match_states=hmm.M)
    if result['hmm_status'] != 'built_uncalibrated':
        manifest['skipped_analyses'].append(dict(analysis='HMM', reason=result['hmm_status']))
    if not build_tree:
        result['tree_status'] = 'skipped_by_request'
    else:
        variable = any(len({s[i] for s in aligned.values() if s[i] in AA}) > 1
                       for i in range(result['alignment_columns']))
        infer = len(records) >= 3 and result['unique_sequences'] >= 3 and variable
        if infer:
            run_tool([fasttree, '-lg', '-gamma', '-seed', str(settings['seed']), folder/'alignment.safe.faa'],
                     folder/'tree.safe.nwk', folder/'fasttree.log', manifest, threads)
            result['tree_status'] = 'fasttree_approximate_ml_unrooted'
        else:
            # This is a membership representation, without invented branch lengths or support.
            text = next(iter(safe)) if len(safe) == 1 else '(' + ','.join(safe) + ')'
            (folder/'tree.safe.nwk').write_text(text + ';\n')
            result['tree_status'] = 'unresolved_membership_only'
            manifest['skipped_analyses'].append(dict(analysis='tree_inference',
                reason='fewer_than_three_distinct_sequences_or_no_variable_columns'))
        tree = check_tree(folder/'tree.safe.nwk', safe)
        (folder/'tree.nwk').write_text(named_newick(tree, mapping))
        check_tree(folder/'tree.nwk', records)
    if not build_tree:
        manifest['skipped_analyses'].append(dict(analysis='tree', reason='skipped_by_request'))
    manifest.update(result=result, output_sha256={p: sha256(folder/p) for p in products if (folder/p).is_file()},
                    interpretation='Automated occupancy curation requires review. Tree uses all members in the full alignment; local SH-like support is not bootstrap support. FastTree optimizes under CAT, then rescales under Gamma20. Profiles have no calibrated detection cutoff or established EET specificity.')
    dump(stamp, manifest)
    (folder/'sequence_manifest.inprogress.json').unlink(missing_ok=True)
    return result


def cluster_table(cluster, proteins, annotations, architectures, domain_names, folder):
    """Exact-ID joins, binary domains, and column 2 reserved for the NR product."""
    if set(proteins) != set(annotations) or set(proteins) != set(architectures):
        raise ValueError(f'{cluster}: incomplete/extra NR or domain protein mapping')
    for row in architectures.values():
        if row.get('cluster') != cluster:
            raise ValueError(f'{cluster}: domain cluster does not match membership')
    domain_sets = {n: set(filter(None, architectures[n]['domain_set'].split(';'))) for n in proteins}
    accessions = sorted(set().union(*domain_sets.values()), key=lambda a: (domain_names.get(a, a), a))
    columns = [f'{domain_names.get(a, a)} [{a}]' for a in accessions]
    def annotation(name):
        row = annotations[name]
        return row['annotation'] or ('No qualifying NR hit' if row['status'] == 'no_qualifying_hit'
                                      else 'NR hit; description unavailable')
    rows = [dict(protein_id=n, annotation=annotation(n),
                 **{col: int(a in domain_sets[n]) for a, col in zip(accessions, columns)}) for n in proteins]
    write_tsv(folder/'proteins.tsv', ['protein_id', 'annotation', *columns], rows)
    write_tsv(folder/'domain_columns.tsv', ['column', 'accession', 'name'],
              [dict(column=c, accession=a, name=domain_names.get(a, a)) for a, c in zip(accessions, columns)])
    return rows, columns, dict(domains=len(accessions), proteins_with_domains=sum(bool(v) for v in domain_sets.values()),
                               distinct_domain_architectures=len({tuple(p.rsplit('@', 1)[0] for p in architectures[n]['architecture'].split(';') if p)
                                                                  for n in proteins}))


STYLE = '<style>body{font:15px system-ui;margin:2rem;color:#172638}a{color:#1563a5}table{border-collapse:collapse;width:100%}td,th{border:1px solid #d8dfe7;padding:.4rem;text-align:left}th{position:sticky;top:0;background:#edf2f7}.scroll{overflow:auto;max-height:70vh}input{padding:.6rem;width:25rem;max-width:80%}.note{background:#f3f5f7;padding:1rem}code{white-space:pre-wrap}</style>'
SEARCH = '<script>document.querySelector("input").oninput=function(){const q=this.value.toLowerCase();document.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));};</script>'


def family_report(folder, cluster, rows, columns, result, domain_db):
    links = [('proteins.tsv', 'Protein/domain table (TSV)'), ('domain_columns.tsv', 'Domain column key'),
             ('../../nr_annotated/'+cluster+'.faa', 'NR-annotated proteins'),
             ('alignment.full.faa', 'Full alignment'), ('alignment.curated.faa', 'Curated HMM alignment'),
             ('curation_sequences.tsv', 'Sequence curation audit'), ('curation_columns.tsv', 'Column curation audit'),
             ('profile.hmm', 'Uncalibrated family HMM'), ('tree.nwk', 'Protein tree (Newick)'),
             ('id_mapping.tsv', 'Safe tree/alignment ID mapping'), ('sequence_manifest.json', 'Sequence-analysis provenance'),
             ('../../domains/clusters/'+cluster+'/domain_network.html', 'Domain network')]
    downloads = ' | '.join(f'<a href="{html.escape(p, quote=True)}">{html.escape(label)}</a>'
                           for p, label in links if (folder/p).is_file())
    fields = ['protein_id', 'annotation', *columns]
    header = ''.join(f'<th>{html.escape(f)}</th>' for f in fields)
    body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(row[f]))}</td>' for f in fields) + '</tr>' for row in rows)
    status = '; '.join(f'{k}: {result[k]}' for k in ['alignment_status', 'hmm_status', 'tree_status'])
    text = f'<!doctype html><html lang="en"><meta charset="utf-8"><title>{html.escape(cluster)}</title>{STYLE}<h1>{html.escape(cluster)}</h1><p><a href="../../index.html">All clusters</a></p><p>{downloads}</p><p>{html.escape(status)}</p><p class="note">{domain_db.upper()} selected domains: 1 = detected, 0 = not detected under the configured search thresholds; repeats count once here. NR descriptions are similarity-based annotations. HMM curation is automated and needs review; trees use the full alignment. An unresolved tree is a membership representation only. None of these results establishes EET activity.</p><input aria-label="Search proteins" placeholder="Search IDs or annotations"><div class="scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>{SEARCH}</html>'
    (folder/'index.html').write_text(text)
