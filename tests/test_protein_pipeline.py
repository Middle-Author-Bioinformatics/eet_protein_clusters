"""Protein-family joins, auditable curation, real sequence tools, and orchestration."""
import json
from pathlib import Path
import random
import shutil
import subprocess

from Bio import Phylo, SeqIO
import pyhmmer as ph
import pytest

import family_analysis as fa
import protein_pipeline as pp
from test_workflow import make_hmm, read


def family_tools():
    mafft = shutil.which('mafft')
    fasttree = shutil.which('FastTree') or shutil.which('fasttree')
    if not mafft or not fasttree:
        pytest.skip('Real MAFFT/FastTree executables are not installed')
    return mafft, fasttree


def test_curation_fragments_duplicates_low_occupancy_and_ambiguity():
    rows = {'full1': 'ACDEFGHIKL---', 'full2': 'ACDEFGHIKM---',
            'full3': 'ACDEFGHIKN---', 'copy': 'ACDEFGHIKL---',
            'fragment': 'AC----------X'}
    curated, seqs, columns, rounds = fa.curate(rows, fa.DEFAULTS)
    assert list(curated) == ['full1', 'full2', 'full3']
    assert all(len(s) == 10 for s in curated.values())
    audit = {r['protein_id']: r for r in seqs}
    assert audit['copy']['hmm_training_status'] == 'exact_duplicate'
    assert audit['copy']['representative'] == 'full1'
    assert audit['fragment']['hmm_training_status'] == 'low_core_coverage'
    assert audit['fragment']['ambiguous_aa'] == 1
    assert [r['full_alignment_column'] for r in columns if r['retained']] == list(range(1, 11))
    assert rounds >= 2


def test_empty_core_is_audited():
    curated, seqs, cols, _ = fa.curate({'a': 'XX--', 'b': '--BB'}, fa.DEFAULTS)
    assert not curated and not any(r['retained'] for r in cols)
    assert {r['hmm_training_status'] for r in seqs} == {'no_core_columns'}


def test_matrix_keeps_all_proteins_and_counts_repeated_domains_once(tmp_path):
    proteins = {'A': 'MAC', 'negative': 'MFG', 'B': 'MLK'}
    nr = {p: dict(annotation='ATPase [Bacterium]' if p == 'A' else '', status='hit' if p == 'A' else 'no_qualifying_hit')
          for p in proteins}
    arch = {'A': dict(cluster='c1', domain_set='PF1.2;PF2.1', architecture='PF1.2@1-10;PF1.2@20-30;PF2.1@40-50'),
            'negative': dict(cluster='c1', domain_set='', architecture=''),
            'B': dict(cluster='c1', domain_set='PF1.2', architecture='PF1.2@1-10')}
    rows, columns, counts = fa.cluster_table('c1', proteins, nr, arch,
                                             {'PF1.2': 'Anchor', 'PF2.1': 'Cytochrome'}, tmp_path)
    table = read(tmp_path/'proteins.tsv')
    assert list(table[0]) == ['protein_id', 'annotation', 'Anchor [PF1.2]', 'Cytochrome [PF2.1]']
    assert table[0]['Anchor [PF1.2]'] == '1'
    assert list(r['protein_id'] for r in table) == list(proteins)
    assert table[1]['annotation'] == 'No qualifying NR hit'
    assert table[1]['Anchor [PF1.2]'] == table[1]['Cytochrome [PF2.1]'] == '0'
    assert counts['proteins_with_domains'] == 2 and counts['domains'] == 2
    with pytest.raises(ValueError, match='mapping'):
        fa.cluster_table('c1', proteins, {'A': nr['A']}, arch, {}, tmp_path)
    arch['B']['cluster'] = 'wrong'
    with pytest.raises(ValueError, match='cluster'):
        fa.cluster_table('c1', proteins, nr, arch, {}, tmp_path)


def test_domain_negative_matrix_has_only_requested_metadata_columns(tmp_path):
    fa.cluster_table('c1', ['p'], {'p': dict(annotation='', status='hit_title_unavailable')},
                     {'p': dict(cluster='c1', domain_set='', architecture='')}, {}, tmp_path)
    table = read(tmp_path/'proteins.tsv')
    assert list(table[0]) == ['protein_id', 'annotation']
    assert table[0]['annotation'] == 'NR hit; description unavailable'


def test_cdd_native_architectures_feed_the_same_named_matrix(tmp_path):
    import cluster_cdd
    from test_cdd import fixture
    database, proteins, native, raw = fixture(tmp_path)
    out = tmp_path/'cdd_out'
    assert cluster_cdd.main(['--input', str(proteins), '--cdd', str(database), '--rpsbproc-output', str(native),
                             '--raw-hits', str(raw), '--output', str(out), '--no-plots']) == 0
    inventory = {r['accession']: r['name'] for r in read(out/'model_inventory.tsv')}
    arch = {r['protein_id']: r for r in read(out/'protein_architectures.tsv') if r['cluster'] == 'c1'}
    matrix = tmp_path/'matrix'; matrix.mkdir()
    fa.cluster_table('c1', ['pair'], {'pair': dict(annotation='Product', status='hit')}, arch, inventory, matrix)
    assert read(matrix/'proteins.tsv') == [dict(protein_id='pair', annotation='Product',
                                               **{'Domain_A [cd00001]': '1', 'Domain_B [cd00002]': '1'})]


def test_full_alignment_validation_rejects_loss_and_changes(tmp_path):
    path = tmp_path/'alignment.faa'
    path.write_text('>a\nAC-D\n>b\nACED\n')
    assert fa.validate_alignment(path, {'a': 'ACD', 'b': 'ACED'})['a'] == 'AC-D'
    with pytest.raises(ValueError, match='changed residues'):
        fa.validate_alignment(path, {'a': 'ACF', 'b': 'ACED'})
    path.write_text('>a\nACD\n>a\nACD\n')
    with pytest.raises(ValueError, match='repeated'):
        fa.validate_alignment(path, {'a': 'ACD'})


def test_incomplete_tree_is_rejected(tmp_path):
    path = tmp_path/'tree.nwk'; path.write_text('(a,b);')
    with pytest.raises(ValueError, match='tips'):
        fa.check_tree(path, ['a', 'b', 'missing'])


@pytest.mark.parametrize('sequences', [{'p': 'ACDEFG'}, {'p': 'ACDEFG', 'q': 'ACDEFH'},
                                      {'p': 'ACDEFG', 'q': 'ACDEFG', 'r': 'ACDEFG'}])
def test_tiny_clusters_have_no_invented_tree_lengths_or_hmm(tmp_path, sequences):
    mafft, fasttree = family_tools()
    result = fa.analyze_family('c1', sequences, tmp_path, {}, {}, mafft, fasttree)
    assert result['hmm_status'] == 'insufficient_unique_training_sequences'
    assert result['tree_status'] == 'unresolved_membership_only'
    tree = fa.check_tree(tmp_path/'tree.nwk', sequences)
    assert all(n.branch_length is None and n.confidence is None for n in tree.find_clades())
    assert not (tmp_path/'profile.hmm').exists()
    assert len(list(SeqIO.parse(tmp_path/'alignment.full.faa', 'fasta'))) == len(sequences)


def test_real_alignment_hmm_tree_id_roundtrip_and_cache(tmp_path, monkeypatch):
    mafft, fasttree = family_tools()
    rng = random.Random(21)
    reference = ''.join(rng.choices(sorted(fa.AA), k=95))
    names = ['a|b', 'x:y', "quote'name", 'a(b)', 'duplicate,one']
    sequences = {n: (reference[:i*8] + 'W' + reference[i*8+1:]) for i, n in enumerate(names[:4])}
    sequences[names[-1]] = sequences[names[0]]
    result = fa.analyze_family('c1', sequences, tmp_path, {}, {}, mafft, fasttree, threads=2)
    assert result['hmm_status'] == 'built_uncalibrated'
    assert result['tree_status'] == 'fasttree_approximate_ml_unrooted'
    assert result['hmm_training_sequences'] == 4
    fa.check_tree(tmp_path/'tree.nwk', sequences)
    with ph.plan7.HMMFile(str(tmp_path/'profile.hmm')) as handle:
        model = next(handle)
        assert model.nseq == 4 and model.M >= 90 and model.cutoffs.gathering is None
    # A usable profile must be searchable, not merely a syntactically valid file.
    seqs = ph.easel.DigitalSequenceBlock(ph.easel.Alphabet.amino(),
        [ph.easel.TextSequence(name='test', sequence=reference).digitize(ph.easel.Alphabet.amino())])
    assert len(next(ph.hmmer.hmmsearch([model], seqs, cpus=1, E=10))) == 1
    original_run = fa.run_tool
    def forbidden(*a, **k):
        raise AssertionError('Matching cache should not execute external tools')
    monkeypatch.setattr(fa, 'run_tool', forbidden)
    cached = fa.analyze_family('c1', sequences, tmp_path, {}, {}, mafft, fasttree, threads=2, resume=True)
    assert cached['cache_reused']
    (tmp_path/'tree.nwk').write_text('corrupt')
    monkeypatch.setattr(fa, 'run_tool', original_run)
    assert not fa.analyze_family('c1', sequences, tmp_path, {}, {}, mafft, fasttree, threads=2, resume=True)['cache_reused']
    changed = dict(fa.DEFAULTS, min_column_occupancy=0.9)
    assert not fa.analyze_family('c1', sequences, tmp_path, changed, {}, mafft, fasttree, threads=2, resume=True)['cache_reused']


def pipeline_fixture(tmp_path):
    seqs = make_hmm(tmp_path/'pfam.hmm')
    a, b = seqs['A'], seqs['B']
    ref = a + 'GGGGGGGGGG' + a + 'GGGGGGGGGG' + b
    proteins = tmp_path/'input.faa'
    proteins.write_text('>protein|one original description\n'+ref+'\n>protein:two\n'+ref[:15]+'W'+ref[16:]
                       +'\n>protein(three)\n'+ref[:24]+'Y'+ref[25:]+'\n>duplicate\n'+ref+'\n>no_hit\nM'+'P'*90+'\n')
    nr = tmp_path/'nr.faa'
    nr.write_text('>WP_00001.1 Test cytochrome product [Synthetic bacterium]\n'+ref+'\n')
    result = subprocess.run(['diamond', 'makedb', '--in', str(nr), '--db', str(tmp_path/'nr.dmnd')], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return ['--input', str(proteins), '--output', str(tmp_path/'out'), '--hmm', str(tmp_path/'pfam.hmm'),
            '--nr-db', str(tmp_path/'nr.dmnd'), '--threads', '2']


def test_real_end_to_end_pfam_nr_mcl_matrices_hmms_trees_and_resume(tmp_path):
    family_tools()
    if not all(shutil.which(t) for t in ['diamond', 'mcl', 'mcxload']):
        pytest.skip('DIAMOND/MCL integration executables not installed')
    args = pipeline_fixture(tmp_path)
    assert pp.main(args) == 0
    out = tmp_path/'out'
    summary = read(out/'cluster_summary.tsv')
    assert len(summary) == 2 and sum(int(r['proteins']) for r in summary) == 5
    large = next(r for r in summary if r['proteins'] == '4')
    folder = out/'families'/large['cluster']
    table = read(folder/'proteins.tsv')
    assert list(table[0]) == ['protein_id', 'annotation', 'A [TEST_A.1]', 'B [TEST_B.1]']
    assert all(r['annotation'] == 'Test cytochrome product [Synthetic bacterium]' for r in table)
    assert all(r['A [TEST_A.1]'] == r['B [TEST_B.1]'] == '1' for r in table)
    assert large['hmm_status'] == 'built_uncalibrated'
    assert large['tree_status'] == 'fasttree_approximate_ml_unrooted'
    assert large['hmm_training_sequences'] == '3'
    fa.check_tree(folder/'tree.nwk', [r['protein_id'] for r in table])
    singleton = next(r for r in summary if r['proteins'] == '1')
    singleton_table = read(out/'families'/singleton['cluster']/'proteins.tsv')
    assert list(singleton_table[0]) == ['protein_id', 'annotation']
    assert singleton_table[0]['annotation'] == 'No qualifying NR hit'
    assert (folder/'index.html').exists() and (out/'index.html').exists()
    manifest = json.loads((out/'manifest.json').read_text())
    assert manifest['parameters']['inflation'] == 2 and manifest['parameters']['min_identity'] == 35
    assert manifest['family_analysis']['hmms_built'] == manifest['family_analysis']['trees_inferred'] == 1
    model_time = (folder/'profile.hmm').stat().st_mtime_ns
    assert pp.main(args+['--resume']) == 0
    assert (folder/'profile.hmm').stat().st_mtime_ns == model_time
    manifest = json.loads((out/'manifest.json').read_text())
    assert manifest['family_analysis']['sequence_caches_reused'] == 2
    assert manifest['stages']['cluster_domains.py']['provenance']['search_definition']['cache_reused']
    # Switching off optional analyses must remove their stale files and combined model DB.
    assert pp.main(args+['--resume', '--skip-tree', '--skip-hmm']) == 0
    assert not (folder/'profile.hmm').exists() and not (folder/'tree.nwk').exists()
    assert not (out/'family_profiles.hmm').exists()
    # Keep different domain databases out of the same output tree, even with force.
    assert pp.main(args+['--domain-db', 'cdd', '--force']) == 2
    assert (out/'manifest.json').exists()


def test_cdd_route_dispatches_existing_profiler_and_passes_its_parameters(tmp_path, monkeypatch):
    # Exercise root orchestration without a large external CDD installation.
    monkeypatch.setattr(pp, 'preflight', lambda a: ({}, 1))
    called = []
    def child(script, argv, output, options, manifest):
        called.append((script, list(map(str, argv))))
        output.mkdir(exist_ok=True)
        (output/'network_report_manifest.json').write_text('{}')
    monkeypatch.setattr(pp, 'child', child)
    monkeypatch.setattr(pp, 'command', lambda *a: None)
    def reports(options, *a):
        (options.output/'cluster_summary.tsv').write_text('cluster\n')
        (options.output/'index.html').write_text('<html></html>')
        (options.output/'families').mkdir()
    monkeypatch.setattr(pp, 'reports', reports)
    monkeypatch.setattr(pp, 'finish', lambda *a: None)
    fasta = tmp_path/'x.faa'; fasta.write_text('>p\nMACDE\n')
    out = tmp_path/'out'
    assert pp.main(['-i', str(fasta), '-o', str(out), '--domain-db', 'cdd', '--cdd', '/db/CDD',
                    '--cdd-evalue', '0.001', '--min-model-coverage', '0.4']) == 0
    cdd = next(a for script, a in called if script == 'cluster_cdd.py')
    assert cdd[cdd.index('--cdd')+1] == '/db/CDD'
    assert cdd[cdd.index('--evalue')+1] == '0.001'
    assert cdd[cdd.index('--min-model-coverage')+1] == '0.4'
    assert not any(script == 'cluster_domains.py' for script, _ in called)


def test_failed_child_has_actionable_error_and_no_success_manifest(tmp_path, monkeypatch):
    monkeypatch.setattr(pp, 'preflight', lambda a: ({}, 1))
    def fail(*a):
        raise ValueError('DIAMOND failed; see search.log')
    monkeypatch.setattr(pp, 'child', fail)
    fasta = tmp_path/'x.faa'; fasta.write_text('>p\nMACDE\n')
    out = tmp_path/'out'
    assert pp.main(['-i', str(fasta), '-o', str(out)]) == 2
    assert not (out/'manifest.json').exists()
    assert 'DIAMOND failed' in json.loads((out/'manifest.failed.json').read_text())['error']


def test_invalid_input_thresholds_and_protected_paths(tmp_path):
    fasta = tmp_path/'x.faa'; fasta.write_text('>p\nMACDE\n')
    database = tmp_path/'nr.dmnd'; database.write_text('dummy')
    base = ['-i', str(fasta), '--nr-db', str(database)]
    assert pp.main(base+['-o', str(tmp_path), '--force']) == 2
    assert fasta.read_text() == '>p\nMACDE\n'
    out = tmp_path/'out'
    assert pp.main(base+['-o', str(out), '--inflation', '1']) == 2
    assert not out.exists()
    fasta.write_text('>p\nMACDE\n>p\nMFGHI\n')
    assert pp.main(base+['-o', str(out)]) == 2
    assert not out.exists()
