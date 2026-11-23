"""HTTP is mocked: no paid requests, private data transmission or real key needed."""
from argparse import Namespace
from io import BytesIO
import csv
import gzip
import json
from pathlib import Path
import sys
import urllib.error
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import summarize_clusters_llm as llm


def inputs(root):
    folder=root/'annotated';folder.mkdir()
    (folder/'cluster_1.faa').write_text(
        '>p1 | nr_hit=WP_1 | nr_annotation=Putative ATPase [Species one]\nMACDEFGHIK\n'
        '>p2 | nr_hit=WP_1 | nr_annotation=Putative ATPase [Species two]\nMACDEFGHIK\n'
        '>p3 | nr_hit=WP_3 | nr_annotation=hypothetical protein [Species one]\nMACDEFGHIK\n'
        '>p4 | nr_status=no_qualifying_hit\nMACDEFGHIK\n')
    with gzip.open(folder/'cluster_2.faa.gz','wt') as out:
        out.write('>p5 | nr_hit=WP_5 | nr_annotation=cytochrome c family protein [Species one]\nMCTTCHE\n')
    (folder/'unknown.faa').write_text('>u | nr_hit=WP_U | nr_annotation=uncharacterized protein [Species]\nMKKKKKK\n')
    (folder/'empty.faa').write_text('')
    return folder


def cli(folder, output):
    return ['-i',str(folder),'-o',str(output),'--pause','0']


def evidence(path):
    return llm.summarize(path,None)


def result_for(data):
    group=next(g for g in data['annotations'] if not g['uninformative'])
    return dict(annotation_prediction='Putative protein family',status='consensus',rationale='The supplied descriptions agree on a broad family.',
                supporting_annotation_ids=[group['annotation_id']],conflicting_annotation_ids=[],
                limitations=['Top-hit annotations are not experimental validation.'],recommended_checks=['Inspect Pfam domains and alignment coverage.'])


def response_for(data):
    return dict(id='resp_test',status='completed',model='test-snapshot',
                usage=dict(input_tokens=100,output_tokens=50,total_tokens=150),
                output=[dict(type='reasoning',content=[]),dict(type='message',content=[dict(type='output_text',text=json.dumps(result_for(data)))])])


def fake_http(monkeypatch):
    calls=[]
    def urlopen(request,timeout):
        assert request.full_url=='https://api.openai.com/v1/responses'
        body=json.loads(request.data);data=json.loads(body['input'][-1]['content']);calls.append(body)
        assert body['text']['format']['strict'] is True and body['store'] is False
        assert 'MACDEFGHIK' not in request.data.decode() and 'TEST_SECRET_KEY' not in request.data.decode()
        assert request.headers['Authorization']=='Bearer TEST_SECRET_KEY'
        return BytesIO(json.dumps(response_for(data)).encode())
    monkeypatch.setattr(llm.urllib.request,'urlopen',urlopen)
    monkeypatch.setenv('OPENAI_API_KEY','TEST_SECRET_KEY')
    monkeypatch.setattr(llm,'API_KEY','')
    return calls


def rows(path):
    with open(path) as source:return list(csv.DictReader(source,delimiter='\t'))


def test_headers_counts_species_pooling_and_unknowns(tmp_path):
    folder=inputs(tmp_path);full=evidence(folder/'cluster_1.faa')
    atp=next(g for g in full['annotations'] if g['annotation']=='Putative ATPase')
    assert full['n_proteins']==4 and full['n_informative']==2 and full['n_uninformative']==2
    assert atp['proteins']==2 and atp['distinct_nr_hits']==1 and len(atp['trailing_bracket_text'])==2
    assert atp['fraction_all_proteins']==.5
    assert llm.generic('Conserved hypothetical protein','nr_hit')
    assert not llm.generic('DUF4016 domain-containing protein','nr_hit')
    assert llm.parse_header('id original note')['annotation']=='original note'
    assert llm.parse_header('id')['source']=='unannotated'
    assert llm.parse_header('id | nr_hit=WP | nr_annotation=[title unavailable in DIAMOND database]')['source']=='missing_title'
    assert llm.parse_header('id | nr_hit=OLD | nr_annotation=old | nr_status=no_qualifying_hit')['source']=='no_hit'
    with pytest.raises(ValueError):llm.parse_header('id | nr_annotation=ATPase')


def test_dry_run_never_calls_api_and_omissions_are_explicit(tmp_path,monkeypatch):
    folder=inputs(tmp_path);out=tmp_path/'out'
    monkeypatch.delenv('OPENAI_API_KEY',raising=False)
    monkeypatch.setattr(llm.urllib.request,'urlopen',lambda *a,**k:pytest.fail('dry-run sent a request'))
    assert llm.main(cli(folder,out)+['--dry-run','--max-labels','1'])==0
    reports=rows(out/'cluster_annotations.tsv');assert len(reports)==4
    assert all(r['result_source']=='dry_run' and not r['annotation_prediction'] for r in reports)
    request=next(json.loads(p.read_text()) for p in (out/'work').glob('*.request.json') if 'cluster_1.faa' in p.read_text())
    sent=json.loads(request['input'][-1]['content'])
    assert sent['n_proteins']==4 and sent['omitted_annotation_groups']==2 and sent['omitted_proteins']==2
    assert 'not independent' in llm.SYSTEM_PROMPT and 'NOT\na pathway' in llm.SYSTEM_PROMPT


def test_mocked_api_results_resume_and_no_secret_saved(tmp_path,monkeypatch):
    folder=inputs(tmp_path);out=tmp_path/'out';calls=fake_http(monkeypatch)
    source_hash={p.name:llm.sha256(p) for p in folder.iterdir()}
    assert llm.main(cli(folder,out))==0 and len(calls)==2
    reports=rows(out/'cluster_annotations.tsv')
    assert len(reports)==4 and sum(r['result_source']=='local_insufficient_evidence' for r in reports)==2
    assert source_hash=={p.name:llm.sha256(p) for p in folder.iterdir()}
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['token_usage']['total_tokens']==300 and manifest['api_clusters_requested']==2
    for p in out.rglob('*'):
        if p.is_file():assert b'TEST_SECRET_KEY' not in p.read_bytes()
    calls.clear();assert llm.main(cli(folder,out)+['--resume'])==0 and not calls
    assert sum(r['result_source']=='cached_api' for r in rows(out/'cluster_annotations.tsv'))==2
    assert json.loads((out/'manifest.json').read_text())['token_usage']=={}
    calls.clear();assert llm.main(cli(folder,out)+['--resume','--model','different-model'])==0 and len(calls)==2
    response=next((out/'work').glob('*.response.json'));response.write_text('{}')
    calls.clear();assert llm.main(cli(folder,out)+['--resume','--model','different-model'])==0 and len(calls)==1


def test_manual_key_and_dry_run_to_live(tmp_path,monkeypatch):
    folder=inputs(tmp_path);out=tmp_path/'out';calls=fake_http(monkeypatch)
    assert llm.main(cli(folder,out)+['--dry-run','--max-clusters','1'])==0 and not calls
    monkeypatch.delenv('OPENAI_API_KEY');monkeypatch.setattr(llm,'API_KEY','TEST_SECRET_KEY')
    assert llm.main(cli(folder,out)+['--resume','--max-clusters','1'])==0 and len(calls)==1
    assert json.loads((out/'manifest.json').read_text())['files_not_selected']==3


def test_missing_key_is_actionable_and_no_hit_only_needs_no_key(tmp_path,monkeypatch):
    folder=inputs(tmp_path);monkeypatch.delenv('OPENAI_API_KEY',raising=False);monkeypatch.setattr(llm,'API_KEY','')
    monkeypatch.setattr(llm.urllib.request,'urlopen',lambda *a,**k:pytest.fail('missing-key request'))
    out=tmp_path/'out';assert llm.main(cli(folder,out))==2
    assert 'API_KEY' in json.loads((out/'manifest.json').read_text())['error']
    (folder/'cluster_1.faa').unlink();(folder/'cluster_2.faa.gz').unlink()
    assert llm.main(cli(folder,tmp_path/'unknowns'))==0


@pytest.mark.parametrize('kind',['refusal','incomplete','bad_json','unknown_id','generic_support','multiline','no_support','wrong_fields'])
def test_rejects_invalid_or_unsupported_model_outputs(tmp_path,kind):
    folder=inputs(tmp_path);data=evidence(folder/'cluster_1.faa');response=response_for(data)
    result=result_for(data)
    if kind=='refusal':response['output']=[dict(type='message',content=[dict(type='refusal',refusal='No')])]
    elif kind=='incomplete':response['status']='incomplete'
    elif kind=='bad_json':response['output'][1]['content'][0]['text']='not json'
    else:
        if kind=='unknown_id':result['supporting_annotation_ids']=['invented']
        if kind=='generic_support':result['supporting_annotation_ids']=[next(g['annotation_id'] for g in data['annotations'] if g['uninformative'])]
        if kind=='multiline':result['annotation_prediction']='line1\nline2'
        if kind=='no_support':result['supporting_annotation_ids']=[]
        if kind=='wrong_fields':result['extra']='bad'
        response['output'][1]['content'][0]['text']=json.dumps(result)
    with pytest.raises(ValueError):llm.decode_response(response,data)


def test_partial_api_failure_is_not_cached_and_resume_retries_only_failure(tmp_path,monkeypatch):
    folder=inputs(tmp_path);out=tmp_path/'out';calls=fake_http(monkeypatch)
    actual=llm.urllib.request.urlopen
    def broken(request,timeout):
        data=json.loads(json.loads(request.data)['input'][-1]['content'])
        if data['file']=='cluster_1.faa':
            bad=response_for(data);bad['status']='incomplete';return BytesIO(json.dumps(bad).encode())
        return actual(request,timeout)
    monkeypatch.setattr(llm.urllib.request,'urlopen',broken)
    assert llm.main(cli(folder,out))==2
    failed=next(r for r in rows(out/'cluster_annotations.tsv') if r['file']=='cluster_1.faa')
    assert failed['result_source']=='api_error' and not failed['annotation_prediction']
    assert len(list((out/'work').glob('*.cache.json')))==1
    calls.clear();monkeypatch.setattr(llm.urllib.request,'urlopen',actual)
    assert llm.main(cli(folder,out)+['--resume'])==0 and len(calls)==1


def test_http_retries_and_authentication_fail_fast(monkeypatch):
    attempts=[];sleeps=[]
    args=Namespace(retries=2,timeout=1)
    def flaky(request,timeout):
        attempts.append(1)
        if len(attempts)==1:raise urllib.error.HTTPError(llm.API_URL,429,'limit',{'Retry-After':'0'},BytesIO(b'{"error":{"code":"rate_limit_exceeded"}}'))
        return BytesIO(b'{"status":"completed"}')
    monkeypatch.setattr(llm.urllib.request,'urlopen',flaky);monkeypatch.setattr(llm.time,'sleep',sleeps.append)
    assert llm.call_api({},'test',args)['status']=='completed' and len(attempts)==2 and sleeps==[0]
    def unauthorized(*a,**k):raise urllib.error.HTTPError(llm.API_URL,401,'bad key',{},BytesIO(b'{"error":{"message":"SECRET"}}'))
    monkeypatch.setattr(llm.urllib.request,'urlopen',unauthorized)
    with pytest.raises(llm.APIError,match='401') as caught:llm.call_api({},'SECRET',args)
    assert caught.value.fatal and 'SECRET' not in str(caught.value)


def test_nr_alignment_table_join_and_mismatches(tmp_path,monkeypatch):
    folder=inputs(tmp_path);target=folder/'cluster_1.faa'
    fields=['file','record','protein_id','status','nr_hit','annotation','identity_percent','query_coverage_percent','subject_coverage_percent']
    rows_data=[]
    for n,line in enumerate([l[1:] for l in target.read_text().splitlines() if l.startswith('>')],1):
        h=llm.parse_header(line)
        rows_data.append(dict(file=target.name,record=n,protein_id=h['protein_id'],status='hit' if h['source']=='nr_hit' else 'no_qualifying_hit',nr_hit=h['nr_hit'],annotation=h['annotation'] if h['source']=='nr_hit' else '',identity_percent=80 if h['nr_hit'] else '',query_coverage_percent=30 if h['nr_hit'] else '',subject_coverage_percent=100 if h['nr_hit'] else ''))
    table=tmp_path/'annotations.tsv'
    with open(table,'w') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(rows_data)
    work=tmp_path/'work';work.mkdir();db=llm.nr_index(table,work)
    data=llm.summarize(target,db);atp=next(g for g in data['annotations'] if g['annotation']=='Putative ATPase')
    assert atp['alignment_statistics']['query_coverage_percent']['below_50_percent']==2
    assert atp['alignment_statistics']['identity_percent']['median']==80
    target.write_text(target.read_text().replace('nr_hit=WP_1','nr_hit=WRONG'))
    with pytest.raises(ValueError,match='mismatch'):llm.summarize(target,db)
    db.close()


def test_input_error_manifest_and_request_caps(tmp_path):
    folder=inputs(tmp_path);out=tmp_path/'out'
    table=tmp_path/'bad.tsv';table.write_text('wrong\tcolumns\n')
    assert llm.main(cli(folder,out)+['--dry-run','--nr-table',str(table)])==2
    assert json.loads((out/'manifest.json').read_text())['status']=='failed'
    assert llm.main(cli(folder,folder)+['--dry-run'])==2
    data=evidence(folder/'cluster_1.faa')
    bounded=llm.bounded_evidence(data,Namespace(max_labels=60,max_evidence_chars=3000))
    assert len(llm.packed(bounded))<=3000 and 'omitted_proteins' in bounded


def test_consumes_real_nr_writer_output_with_alignment_statistics(tmp_path):
    import annotate_clusters_nr as nr
    source=tmp_path/'original.faa';source.write_text('>protein original note\n'+'M'*100+'\n>unmatched\n'+'P'*100+'\n')
    work=tmp_path/'nr_work';work.mkdir();db,_=nr.prepare([source],work)
    raw=tmp_path/'hits.tsv';raw.write_text('q000000000001\tWP_REF\t75\t90\t1\t90\t1\t90\t1e-20\t150\t100\t100\tWP_REF Putative ATPase [Species]\n')
    nr.parse_hits(db,raw,Namespace(evalue=1e-5,min_identity=0,query_cover=0,subject_cover=0))
    annotated=tmp_path/'annotated';annotated.mkdir();nr.emit(db,[source],annotated);db.close()
    out=tmp_path/'llm';assert llm.main(cli(annotated,out)+['--dry-run'])==0
    sent=json.loads(json.loads(next((out/'work').glob('*.request.json')).read_text())['input'][-1]['content'])
    assert sent['n_proteins']==2 and sent['n_informative']==1 and sent['alignment_statistics_available']
    atp=next(g for g in sent['annotations'] if g['annotation']=='Putative ATPase')
    assert atp['alignment_statistics']['query_coverage_percent']['median']==90


def test_mixed_result_keeps_a_single_label_and_conflicting_evidence(tmp_path):
    folder=inputs(tmp_path);path=folder/'cluster_1.faa'
    with open(path,'a') as out:out.write('>conflict | nr_hit=WP_C | nr_annotation=Cytochrome c family protein [Species]\nMCTTCH\n')
    data=evidence(path);result=result_for(data)
    atp=next(g['annotation_id'] for g in data['annotations'] if g['annotation']=='Putative ATPase')
    heme=next(g['annotation_id'] for g in data['annotations'] if g['annotation']=='Cytochrome c family protein')
    result.update(annotation_prediction='Mixed protein cluster; function unresolved',status='mixed',
                  supporting_annotation_ids=[atp],conflicting_annotation_ids=[heme])
    assert llm.validate_result(result,data)['status']=='mixed'
    assert result['conflicting_annotation_ids']==[heme]
