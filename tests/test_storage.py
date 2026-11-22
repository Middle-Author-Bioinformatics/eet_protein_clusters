"""Scratch placement, safe cleanup and completed-search cache semantics."""
import json
from pathlib import Path
import shutil
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import cluster_proteins as cp


def input_fasta(tmp_path):
    inp=tmp_path/'all.faa';inp.write_text('>a\nMACDEFGHIKLMNPQRSTVWY\n>b\nMACDEFGHIKLMNPQRSTVWY\n')
    return inp


def fake_commands(monkeypatch, fail=False):
    """Inject a disk-full failure without consuming any disk space."""
    calls=[]
    def command(argv, logfile, manifest):
        argv=[str(a) for a in argv];calls.append(argv);manifest['commands'].append(argv)
        Path(logfile).write_text('test double\n')
        if argv[:2]==['diamond','blastp']:
            scratch=Path(argv[argv.index('--tmpdir')+1]);assert scratch.is_dir()
            (scratch/'owned-scratch-file').write_text('temporary')
            hits=Path(argv[argv.index('--out')+1]);hits.write_text('')
            if fail: raise ValueError('diamond exited 1: No space left on device')
    monkeypatch.setattr(cp,'command',command)
    monkeypatch.setattr(cp,'tool_version',lambda *args:'test-double-version')
    return calls


@pytest.mark.parametrize('source',['explicit','environment','default'])
def test_scratch_selection_and_cleanup(tmp_path,monkeypatch,source):
    inp=input_fasta(tmp_path);out=tmp_path/'out';calls=fake_commands(monkeypatch)
    monkeypatch.delenv('TMPDIR',raising=False)
    args=['-i',str(inp),'-o',str(out)]
    parent=(out/'work') if source=='default' else tmp_path/'scratch with spaces'
    parent.mkdir(parents=True);sentinel=parent/'unrelated-job';sentinel.write_text('keep')
    if source=='explicit':
        args+=['--tmpdir',str(parent)];monkeypatch.setenv('TMPDIR',str(tmp_path/'ignored'))
    elif source=='environment':monkeypatch.setenv('TMPDIR',str(parent))
    # prepare_output rejects the sentinel-containing output without --resume.
    if source=='default':args.append('--resume')
    assert cp.main(args)==0
    data=json.loads((out/'manifest.json').read_text())
    blast=next(a for a in calls if a[:2]==['diamond','blastp'])
    scratch=Path(blast[blast.index('--tmpdir')+1])
    assert scratch.parent==parent.resolve() and not scratch.exists()
    assert sentinel.read_text()=='keep' and data['storage']['scratch_cleanup']=='removed'
    assert data['storage']['scratch_parent']['path']==str(parent.resolve())
    assert data['storage']['scratch_parent']['free_bytes']>=0
    assert data['storage']['scratch_source']=={'explicit':'--tmpdir','environment':'TMPDIR','default':'output/work'}[source]


def test_disk_full_invalidates_old_stamp_and_restarts(tmp_path,monkeypatch):
    inp=input_fasta(tmp_path);out=tmp_path/'out';parent=tmp_path/'scratch';parent.mkdir()
    sentinel=parent/'another-job';sentinel.write_text('keep')
    args=['-i',str(inp),'-o',str(out),'--tmpdir',str(parent)]
    fake_commands(monkeypatch)
    assert cp.main(args)==0 and (out/'work/diamond.stage.json').exists()
    fake_commands(monkeypatch,fail=True)
    assert cp.main(args+['--force'])==2
    assert not (out/'manifest.json').exists() and not (out/'work/diamond.stage.json').exists()
    progress=json.loads((out/'manifest.inprogress.json').read_text())
    assert '--tmpdir' in progress['error'] and 'df -i' in progress['error']
    assert progress['storage']['scratch_cleanup']=='removed'
    assert not list(parent.glob('eet-diamond-*')) and sentinel.read_text()=='keep'
    calls=fake_commands(monkeypatch)
    assert cp.main(args+['--resume'])==0
    assert any(a[:2]==['diamond','blastp'] for a in calls)


def test_changing_scratch_preserves_completed_search(tmp_path,monkeypatch):
    inp=input_fasta(tmp_path);out=tmp_path/'out';parent=tmp_path/'scratch'
    calls=fake_commands(monkeypatch)
    args=['-i',str(inp),'-o',str(out)]
    assert cp.main(args+['--tmpdir',str(parent)])==0
    hits=out/'work/diamond_hits.tsv';before=hits.stat().st_mtime_ns
    calls.clear();new_parent=tmp_path/'unused-scratch'
    assert cp.main(args+['--resume','--tmpdir',str(new_parent)])==0
    assert not calls and hits.stat().st_mtime_ns==before and not new_parent.exists()
    assert json.loads((out/'manifest.json').read_text())['diamond_search_cache_reused']


def test_bad_scratch_path_is_visible_and_not_cached(tmp_path,monkeypatch):
    inp=input_fasta(tmp_path);out=tmp_path/'out';parent=tmp_path/'not-a-directory';parent.write_text('keep')
    fake_commands(monkeypatch)
    assert cp.main(['-i',str(inp),'-o',str(out),'--tmpdir',str(parent)])==2
    assert not (out/'work/diamond.stage.json').exists() and parent.read_text()=='keep'
    assert json.loads((out/'manifest.inprogress.json').read_text())['error']


@pytest.mark.skipif(not all(shutil.which(t) for t in ['diamond','mcl','mcxload']),reason='clustering binaries not installed')
def test_real_diamond_uses_custom_scratch(tmp_path):
    inp=input_fasta(tmp_path);out=tmp_path/'out';scratch=tmp_path/'custom scratch'
    assert cp.main(['-i',str(inp),'-o',str(out),'--tmpdir',str(scratch),'--block-size','0.4','-t','1'])==0
    manifest=json.loads((out/'manifest.json').read_text())
    blast=next(a for a in manifest['commands'] if a[:2]==['diamond','blastp'])
    assert Path(blast[blast.index('--tmpdir')+1]).parent==scratch.resolve()
    assert not list(scratch.iterdir()) and (out/'clusters/cluster_000001.faa').exists()


def test_disk_full_failure_is_not_masked_by_manifest_write(tmp_path,monkeypatch,caplog):
    inp=input_fasta(tmp_path);out=tmp_path/'out';fake_commands(monkeypatch,fail=True)
    original_dump=cp.dump
    def full_disk_dump(path,data):
        if 'error' in data:raise OSError(28,'No space left on device')
        return original_dump(path,data)
    monkeypatch.setattr(cp,'dump',full_disk_dump)
    assert cp.main(['-i',str(inp),'-o',str(out),'--tmpdir',str(tmp_path/'scratch')])==2
    assert '--tmpdir on a larger filesystem' in caplog.text
    assert 'Could not save failure manifest' in caplog.text
    assert not (out/'work/diamond.stage.json').exists()
