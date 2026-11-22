"""Top-hit provenance, unchanged cluster membership/sequences, and safe resume."""
from argparse import Namespace
import csv
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import annotate_clusters_nr as nr


def table(path):
    with open(path) as source:
        return list(csv.DictReader(source, delimiter='\t'))


def inputs(tmp_path):
    folder = tmp_path/'clusters'; folder.mkdir()
    (folder/'cluster_000001.faa').write_bytes(b'>same|id original description\r\nMacDefGhIkLMNPQRSTVWY*\r\n>nohit\r\nMPPPPPPPPPPPPPPPPPPP\r\n')
    with gzip.open(folder/'cluster_000002.faa.gz','wt') as out:
        out.write('>same|id another cluster\nMKKKKKKKKKKKKKKKKKKK\n')
    (folder/'empty.faa').write_text('')
    (folder/'ignore.txt').write_text('not a FASTA')
    database = tmp_path/'nr.dmnd'; database.write_text('fake database')
    return folder, database


def fake_execute(monkeypatch, fail=False):
    calls = []
    def execute(command, log, manifest):
        command = [str(a) for a in command]; calls.append(command)
        manifest['commands'].append(command); log.write_text('test DIAMOND\n')
        if command[1] == 'version':
            return 'diamond version test'
        if command[1] == 'blastp':
            scratch = Path(command[command.index('--tmpdir')+1]); assert scratch.is_dir()
            (scratch/'owned').write_text('temporary')
            output = Path(command[command.index('--out')+1])
            output.write_text('q000000000001\tWP_123\t90\t21\t1\t21\t1\t21\t1e-20\t100\t21\t21\tWP_123 Cytochrome-like protein [Test bacterium]\n'
                              'q000000000003\tWP_BARE\t100\t20\t1\t20\t1\t20\t0\t120\t20\t20\tWP_BARE\n')
            if fail:
                raise RuntimeError('DIAMOND: No space left on device')
        return 'test database info'
    monkeypatch.setattr(nr, 'execute', execute)
    return calls


def command(folder, database, out):
    return ['-i',str(folder),'-o',str(out),'--db',str(database),'--diamond',sys.executable,'-t','1']


def test_preserves_ids_sequences_headers_files_and_no_hits(tmp_path, monkeypatch):
    folder, database = inputs(tmp_path); out = tmp_path/'annotated'
    before = {p.name: nr.sha256(p) for p in folder.iterdir()}
    calls = fake_execute(monkeypatch)
    scratch = tmp_path/'scratch with spaces'; scratch.mkdir(); (scratch/'another_job').write_text('keep')
    assert nr.main(command(folder,database,out)+['--tmpdir',str(scratch)]) == 0
    assert len([a for a in calls if a[1]=='blastp']) == 1
    for path in folder.iterdir():
        assert nr.sha256(path) == before[path.name]
        if path.name.endswith(nr.SUFFIXES):
            with nr.fasta_open(path) as a, nr.fasta_open(out/path.name) as b:
                aa, bb = list(a), list(b)
            assert [s for s in aa if not s.startswith('>')] == [s for s in bb if not s.startswith('>')]
            assert [s.split()[0] for s in aa if s.startswith('>')] == [s.split()[0] for s in bb if s.startswith('>')]
    text = (out/'cluster_000001.faa').read_text()
    assert '>same|id original description | nr_hit=WP_123 | nr_annotation=Cytochrome-like protein [Test bacterium]' in text
    assert '>nohit | nr_status=no_qualifying_hit' in text
    rows = table(out/'annotations.tsv'); assert len(rows) == 3
    assert [r['status'] for r in rows] == ['hit','no_qualifying_hit','hit_title_unavailable']
    assert rows[0]['query_coverage_percent'] == '100.0'
    assert rows[2]['protein_id'] == rows[0]['protein_id'] and rows[2]['file'] != rows[0]['file']
    assert table(out/'cluster_summary.tsv')[-1]['proteins'] == '0'
    manifest = json.loads((out/'manifest.json').read_text())
    assert manifest['status'] == 'complete' and manifest['scratch_cleanup'] == 'removed'
    assert list(scratch.iterdir()) == [scratch/'another_job']


def test_resume_checks_inputs_database_settings_and_raw_results(tmp_path, monkeypatch):
    folder, database = inputs(tmp_path); out = tmp_path/'out'; calls = fake_execute(monkeypatch)
    args = command(folder,database,out)
    assert nr.main(args) == 0
    hits = out/'work/nr_hits.tsv'; before = hits.stat().st_mtime_ns
    calls.clear()
    assert nr.main(args+['--resume','--tmpdir',str(tmp_path/'new_scratch')]) == 0
    assert not any(a[1]=='blastp' for a in calls) and hits.stat().st_mtime_ns == before
    assert not (tmp_path/'new_scratch').exists()
    assert json.loads((out/'manifest.json').read_text())['search_cache_reused']
    database.write_text('changed database')
    calls.clear(); assert nr.main(args+['--resume']) == 0
    assert any(a[1]=='blastp' for a in calls)
    calls.clear(); hits.write_text('corrupted raw file')
    assert nr.main(args+['--resume']) == 0 and any(a[1]=='blastp' for a in calls)
    calls.clear(); (folder/'cluster_000001.faa').write_text((folder/'cluster_000001.faa').read_text().replace('original description','new description'))
    assert nr.main(args+['--resume']) == 0 and any(a[1]=='blastp' for a in calls)
    assert 'new description | nr_hit=' in (out/'cluster_000001.faa').read_text()
    calls.clear(); assert nr.main(args+['--resume','--query-cover','50']) == 0
    assert any(a[1]=='blastp' for a in calls)


def test_failed_search_not_reused_and_cleanup_is_owned(tmp_path, monkeypatch):
    folder, database = inputs(tmp_path); out = tmp_path/'out'; args = command(folder,database,out)
    fake_execute(monkeypatch); assert nr.main(args) == 0
    fake_execute(monkeypatch,fail=True)
    assert nr.main(args+['--force']) == 2
    assert not (out/'work/diamond.stage.json').exists() and not (out/'work/nr_hits.tsv.partial').exists()
    manifest = json.loads((out/'manifest.json').read_text())
    assert manifest['status']=='failed' and manifest['scratch_cleanup']=='removed'
    calls = fake_execute(monkeypatch)
    assert nr.main(args+['--resume']) == 0 and any(a[1]=='blastp' for a in calls)


def test_stale_generated_fastas_removed_but_unrelated_files_kept(tmp_path, monkeypatch):
    folder, database = inputs(tmp_path); out = tmp_path/'out'; args = command(folder,database,out)
    fake_execute(monkeypatch); assert nr.main(args)==0
    (out/'user_notes.txt').write_text('keep'); (folder/'empty.faa').unlink()
    assert nr.main(args+['--resume'])==0
    assert not (out/'empty.faa').exists() and (out/'user_notes.txt').read_text()=='keep'


def test_ranking_is_by_score_and_ties_not_by_named_function(tmp_path):
    folder = tmp_path/'clusters'; folder.mkdir(); p = folder/'x.faa'; p.write_text('>protein\n'+'M'*100+'\n')
    work = tmp_path/'work'; work.mkdir(); db,_ = nr.prepare([p],work)
    raw = tmp_path/'hits'
    def line(subject, score, evalue='1e-30', title='hypothetical protein', start=1, end=100):
        return f'q000000000001\t{subject}\t80\t100\t{start}\t{end}\t1\t100\t{evalue}\t{score}\t100\t100\t{subject} {title}\n'
    raw.write_text(line('NAMED',90,'1e-80','known enzyme')+line('Z_TOP',100)+line('A_TOP',100)+line('LOW',80))
    args = Namespace(evalue=1e-5,min_identity=0,query_cover=0,subject_cover=0)
    nr.parse_hits(db,raw,args); hit = json.loads(db.execute('SELECT payload FROM hits').fetchone()[0])
    assert hit['nr_hit']=='A_TOP' and hit['annotation']=='hypothetical protein' and hit['best_score_ties']==2
    assert nr.clean_title('WP_1','WP_1')==''
    assert nr.clean_title('protein\x01description [species]', 'WP_1')=='protein description [species]'
    db.close()


@pytest.mark.parametrize('row',[
    'unknown\tS\t80\t100\t1\t100\t1\t100\t1e-20\t100\t100\t100\ttitle',
    'q000000000001\tS\tnan\t100\t1\t100\t1\t100\t1e-20\t100\t100\t100\ttitle',
    'q000000000001\tS\t80\t100\t1\t101\t1\t100\t1e-20\t100\t100\t100\ttitle',
    'only\ttwo',
])
def test_bad_alignment_records_rejected(tmp_path,row):
    folder = tmp_path/'in'; folder.mkdir(); p = folder/'x.faa'; p.write_text('>p\n'+'M'*100+'\n')
    work = tmp_path/'work'; work.mkdir(); db,_ = nr.prepare([p],work)
    raw = tmp_path/'hits'; raw.write_text(row+'\n')
    with pytest.raises(ValueError):nr.parse_hits(db,raw,Namespace(evalue=1e-5,min_identity=0,query_cover=0,subject_cover=0))
    db.close()


def test_safe_paths_and_invalid_sequences(tmp_path,monkeypatch):
    folder, database = inputs(tmp_path); args = command(folder,database,folder)
    assert nr.main(args)==2
    assert nr.main(command(folder,database,folder/'nested'))==2
    out=tmp_path/'unrelated'; out.mkdir(); (out/'userfile').write_text('keep')
    assert nr.main(command(folder,database,out)+['--resume'])==2
    assert (out/'userfile').read_text()=='keep'
    fake_execute(monkeypatch)
    (folder/'cluster_000001.faa').write_text('>bad\nMA*ACD\n')
    assert nr.main(command(folder,database,tmp_path/'bad'))==2
    assert json.loads((tmp_path/'bad/manifest.json').read_text())['status']=='failed'


@pytest.mark.skipif(not shutil.which('diamond'),reason='DIAMOND not installed')
def test_real_diamond_titles_ties_no_hits_and_resume(tmp_path):
    rng=random.Random(91); alphabet='ACDEFGHIKLMNPQRSTVWY'
    seq=''.join(rng.choices(alphabet,k=180)); other=''.join(rng.choices(alphabet,k=180))
    reference=tmp_path/'reference.faa'
    reference.write_text(f'>Z_REF Named enzyme [Test species]\n{seq}\n>A_REF hypothetical protein [Test species]\n{seq}\n>BARE\n{other}\n')
    database=tmp_path/'tiny.dmnd'
    subprocess.run(['diamond','makedb','--in',str(reference),'--db',str(database)],check=True,capture_output=True)
    folder=tmp_path/'clusters'; folder.mkdir()
    (folder/'cluster_000001.faa').write_text(f'>original|id original header\n{seq.lower()}*\n>unmatched\nM'+ 'P'*180+'\n')
    (folder/'cluster_000002.faa').write_text(f'>original|id second file\n{other}\n')
    out=tmp_path/'annotated'; scratch=tmp_path/'scratch'
    args=['-i',str(folder),'-o',str(out),'--db',str(database),'-t','1','--tmpdir',str(scratch),'--query-cover','90','--subject-cover','90']
    assert nr.main(args)==0
    rows=table(out/'annotations.tsv')
    assert rows[0]['nr_hit']=='A_REF' and rows[0]['annotation']=='hypothetical protein [Test species]'
    assert rows[0]['best_score_ties']=='2' and rows[0]['identity_percent']=='100.0'
    assert rows[1]['status']=='no_qualifying_hit' and rows[2]['status']=='hit_title_unavailable'
    assert seq.lower()+'*' in (out/'cluster_000001.faa').read_text()
    before=(out/'work/nr_hits.tsv').stat().st_mtime_ns
    assert nr.main(args+['--resume'])==0
    assert (out/'work/nr_hits.tsv').stat().st_mtime_ns==before
    manifest=json.loads((out/'manifest.json').read_text()); assert manifest['search_cache_reused']
    assert not list(scratch.iterdir())


@pytest.mark.skipif(not shutil.which('diamond'),reason='DIAMOND not installed')
def test_real_coverage_filter_selects_best_qualifying_hit(tmp_path):
    rng=random.Random(771); alphabet='ACDEFGHIKLMNPQRSTVWY'
    seq=''.join(rng.choices(alphabet,k=240))
    broad=''.join(rng.choice(alphabet.replace(aa,'')) if i%5 in (0,1) else aa for i,aa in enumerate(seq))
    reference=tmp_path/'reference.faa'
    reference.write_text(f'>SHORT Short exact domain [Test species]\n{seq[40:190]}\n>BROAD Full-length remote match [Test species]\n{broad}\n')
    database=tmp_path/'tiny.dmnd'
    subprocess.run(['diamond','makedb','--in',str(reference),'--db',str(database)],check=True,capture_output=True)
    folder=tmp_path/'clusters';folder.mkdir();(folder/'c.faa').write_text(f'>p\n{seq}\n')
    base=['-i',str(folder),'--db',str(database),'-t','1']
    assert nr.main(base+['-o',str(tmp_path/'any')])==0
    assert nr.main(base+['-o',str(tmp_path/'covered'),'--query-cover','90'])==0
    any_hit=table(tmp_path/'any/annotations.tsv')[0];covered_hit=table(tmp_path/'covered/annotations.tsv')[0]
    assert any_hit['nr_hit']=='SHORT' and covered_hit['nr_hit']=='BROAD'
    assert float(any_hit['query_coverage_percent'])<90<=float(covered_hit['query_coverage_percent'])
    assert float(any_hit['bitscore'])>float(covered_hit['bitscore'])
