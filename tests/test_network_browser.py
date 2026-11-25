"""Actual offline browser interactions; Playwright is a development dependency only."""
import json
import os
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from interactive_network import write_html


def test_offline_interactive_controls_and_downloads(tmp_path):
    executable=os.environ.get('EET_BROWSER_EXECUTABLE')
    if not executable:
        if os.environ.get('REQUIRE_BROWSER_TESTS'):pytest.fail('EET_BROWSER_EXECUTABLE is required for CI browser tests')
        pytest.skip('Set EET_BROWSER_EXECUTABLE to enable actual browser tests')
    pw=pytest.importorskip('playwright.sync_api')
    target=tmp_path/'network.html'
    nodes=[dict(id='A',accession='PF00395.29',label='SLH (PF00395)',name='S-layer homology',group='community_1',size=10,proteins_with_domain=3),
           dict(id='B',label='Heme binding (PF_TEST)',name='Heme binding',group='community_1',size=10,proteins_with_domain=2),
           dict(id='C',label='Isolated domain (PF_ISO)',name='Isolated domain',group='community_2',size=8,proteins_with_domain=0)]
    edge=dict(source='A',target='B',value=.4,jaccard=.4,lift=2,proteins_with_both=2)
    members=dict(available=True,matching_proteins=3,embedded_count=3,truncated=False,
                 proteins=[dict(protein_id='pair1',annotation='Putative cytochrome c protein',nr_hit='WP_1',cluster='c1',domain_names='SLH;Heme',architecture_names='SLH + Heme'),
                           dict(protein_id='pair2',annotation='Cell wall protein',nr_hit='WP_2',cluster='c2',domain_names='SLH;Heme',architecture_names='SLH + Heme'),
                           dict(protein_id='SLH_only',cluster='c1',domain_names='SLH',architecture_names='SLH')],by_domain=dict(A=[0,1,2],B=[0,1],C=[]),
                 cluster_index=dict(available=True,unassigned_proteins=0,source='protein_architectures.tsv',clusters=[
                     dict(id='c1',proteins=2,domains=dict(A=2,B=1),total_domains=2),
                     dict(id='c2',proteins=1,domains=dict(A=1,B=1),total_domains=2),
                     dict(id='c3',proteins=1,domains={},total_domains=0)]))
    write_html(target,nodes,[edge],'Domain network','domains','Same-protein co-occurrence',members=members)
    with pw.sync_playwright() as p:
        browser=p.chromium.launch(executable_path=executable,headless=True,args=['--no-sandbox'])
        context=browser.new_context(viewport={'width':1440,'height':1000},offline=True,accept_downloads=True)
        page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(target.as_uri());page.wait_for_function('window.networkView && document.querySelector("canvas").width>1')
        assert page.evaluate('networkView.getVisible().nodes.length')==3
        assert page.locator('#protein-cluster option').count()==4
        page.locator('#protein-cluster').select_option('c1')
        assert set(page.evaluate('networkView.getHighlight().domains'))=={'A','B'}
        assert page.evaluate('networkView.getVisible().nodes.length')==3
        assert page.evaluate('networkView.getVisible().edges.length')==1
        assert 'edges stay unchanged' in page.locator('#cluster-status').inner_text()
        page.locator('#protein-cluster').select_option('c3')
        assert page.evaluate('networkView.getHighlight().domains')==[]
        assert '0 / 0' in page.locator('#cluster-status').inner_text()
        page.locator('#reset').click()
        assert page.evaluate('networkView.getHighlight().cluster') is None
        colors=page.locator('#legend .swatch').evaluate_all('(nodes)=>nodes.map(n=>getComputedStyle(n).backgroundColor)')
        assert len(set(colors))==2
        page.locator('#search').fill('PF00395');page.locator('#include-neighbors').uncheck()
        assert page.evaluate('networkView.getVisible().nodes.length')==1
        page.locator('#include-neighbors').check();assert page.evaluate('networkView.getVisible().nodes.length')==2
        page.locator('#reset').click();page.locator('#threshold').fill('0.8');page.locator('#hide-isolated').check()
        assert page.evaluate('networkView.getVisible().nodes.length')==0
        page.locator('#reset').click();page.locator('#group').select_option('community_2')
        assert page.evaluate('networkView.getVisible().nodes.length')==1
        page.locator('#reset').click();page.locator('#browse-nodes summary').click();page.locator('#node-list button').first.click()
        assert 'SLH (PF00395)' in page.locator('#details').inner_text()
        assert '2 proteins' in page.locator('#neighbors').inner_text()
        assert not page.locator('#technical-details').evaluate('(el)=>el.open')
        assert 'proteins_with_domain' not in page.locator('aside').inner_text()
        assert page.locator('#pfam-link').get_attribute('href')=='https://www.ebi.ac.uk/interpro/entry/pfam/PF00395/'
        assert page.locator('#pfam-link').get_attribute('rel')=='noopener noreferrer'
        page.evaluate('networkView.data.nodes[0].domain_source="CDD";networkView.data.nodes[0].cdd_pssm_id="197609";networkView.selectNode("A")')
        assert page.locator('#cdd-link').get_attribute('href')=='https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=197609'
        assert page.locator('#cdd-link').get_attribute('rel')=='noopener noreferrer'
        assert page.locator('#pfam-link').count()==0
        page.evaluate('delete networkView.data.nodes[0].domain_source;networkView.selectNode("A")')
        assert '3 / 3' in page.locator('#member-status').inner_text()
        assert set(page.locator('#member-list .cluster-badge').all_inner_texts())=={'c1','c2'}
        page.locator('#protein-cluster').select_option('c2')
        assert page.locator('#member-list li').count()==3
        assert page.locator('#member-list li').first.locator('.member-id').inner_text()=='pair2'
        assert '1 in c2' in page.locator('#details').inner_text()
        page.locator('.shared-proteins').click()
        assert '2 / 2' in page.locator('#member-status').inner_text()
        assert 'SLH_only' not in page.locator('#member-list').inner_text()
        assert 'Putative cytochrome c protein' in page.locator('#member-list').inner_text()
        assert 'SLH + Heme' not in page.locator('#member-list').inner_text()
        page.locator('#member-list li').first.locator('summary').click()
        assert 'SLH + Heme' in page.locator('#member-list li').first.inner_text()
        page.locator('#member-search').fill('Cell wall')
        assert page.locator('#member-list li').count()==1
        page.locator('#member-search').fill('pair2')
        assert page.locator('#member-list li').count()==1
        with page.expect_download() as carriers:page.locator('#download-members').click()
        exported_carriers=Path(carriers.value.path()).read_text()
        assert 'pair1' in exported_carriers and 'pair2' in exported_carriers and 'SLH_only' not in exported_carriers
        assert 'Putative cytochrome c protein' in exported_carriers and 'WP_2' in exported_carriers
        assert 'c1' in exported_carriers and 'c2' in exported_carriers
        page.evaluate('networkView.selectNode("B")');assert page.locator('#pfam-link').count()==0
        page.evaluate('networkView.data.members.available=false;networkView.selectNode("A")')
        assert 'unavailable' in page.locator('#member-status').inner_text() and page.locator('#download-members').is_hidden()
        page.evaluate('networkView.data.members.available=true;networkView.data.members.by_domain.A=[0];networkView.selectNode("A")')
        assert 'Partial list' in page.locator('#member-status').inner_text()
        bbox=page.locator('#network').bounding_box()
        before=page.evaluate('networkView.getCamera().scale');page.mouse.move(bbox['x']+bbox['width']/2,bbox['y']+bbox['height']/2);page.mouse.wheel(0,-100)
        page.wait_for_function('(previous)=>networkView.getCamera().scale>previous',arg=before)
        # Canvas node drag changes its coordinates without mutating scientific metadata.
        point=page.evaluate('()=>{const n=networkView.data.nodes[0],c=networkView.getCamera();return {x:n.x*c.scale+c.x,y:n.y*c.scale+c.y,nx:n.x}}')
        page.mouse.move(bbox['x']+point['x'],bbox['y']+point['y']);page.mouse.down();page.mouse.move(bbox['x']+point['x']+30,bbox['y']+point['y']+10,steps=3);page.mouse.up()
        assert page.evaluate('networkView.data.nodes[0].x')!=point['nx']
        page.locator('#reset').click()
        with page.expect_download() as dl:page.locator('#download-data').click()
        exported=json.loads(Path(dl.value.path()).read_text());assert len(exported['nodes'])==3 and len(exported['edges'])==1
        with page.expect_download() as png:page.locator('#download-image').click()
        assert Path(png.value.path()).read_bytes().startswith(b'\x89PNG')
        if os.environ.get('NETWORK_TEST_ARTIFACTS'):
            page.evaluate('networkView.data.members.by_domain.A=[0,1,2];networkView.selectNode("A")')
            page.locator('#protein-cluster').select_option('c2')
            page.locator('.shared-proteins').click()
            folder=Path(os.environ['NETWORK_TEST_ARTIFACTS']);folder.mkdir(parents=True,exist_ok=True);page.screenshot(path=folder/'domain-network.png',full_page=True)
        with page.expect_download() as tab:page.locator('#download-table').click()
        assert 'SLH (PF00395)' in Path(tab.value.path()).read_text()
        # A highly connected domain starts with a short list; the rest is expandable.
        extras=[dict(id=f'N{i}',label=f'Domain {i}',name=f'Domain {i}',group='community_1',size=8,proteins_with_domain=1) for i in range(10)]
        more_edges=[dict(source='A',target=n['id'],value=.3,jaccard=.3,lift=1,proteins_with_both=1) for n in extras]
        many_members=json.loads(json.dumps(members))
        for n in extras:
            many_members['by_domain'][n['id']]=[0]
            many_members['cluster_index']['clusters'][0]['domains'][n['id']]=1
            many_members['cluster_index']['clusters'][0]['total_domains']+=1
        many=tmp_path/'many.html';write_html(many,nodes+extras,[edge]+more_edges,'Many connections','domains','Definition',members=many_members)
        page.goto(many.as_uri());page.wait_for_function('window.networkView');page.evaluate('networkView.selectNode("A")')
        assert page.locator('#neighbors .neighbor-target:visible').count()==8
        page.locator('#neighbors details summary').click()
        assert page.locator('#neighbors .neighbor-target:visible').count()==11
        page.locator('#technical-details summary').click()
        assert 'Jaccard' in page.locator('#statistics').inner_text()
        # Domain presence remains available when no carrier identities are embedded.
        capped=json.loads(json.dumps(members));capped.update(proteins=[],by_domain=dict(A=[]),embedded_count=0,truncated=True)
        for c in capped['cluster_index']['clusters']:c['domains']={k:v for k,v in c['domains'].items() if k=='A'}
        capfile=tmp_path/'capped.html';write_html(capfile,[nodes[0]],[],'Capped','domains','Definition',total_nodes=3,total_edges=1,members=capped)
        page.goto(capfile.as_uri());page.wait_for_function('window.networkView')
        page.locator('#protein-cluster').select_option('c2')
        assert page.evaluate('networkView.getHighlight().visible_domains')==['A']
        assert '1 / 2' in page.locator('#cluster-status').inner_text()
        page.evaluate('networkView.selectNode("A")')
        assert '0 / 3' in page.locator('#member-status').inner_text()
        legacy=json.loads(json.dumps(members));legacy.pop('cluster_index')
        legacyfile=tmp_path/'legacy.html';write_html(legacyfile,nodes,[edge],'Legacy','domains','Definition',members=legacy)
        page.goto(legacyfile.as_uri());page.wait_for_function('window.networkView')
        assert page.locator('#protein-cluster').is_disabled()
        assert not errors;browser.close()
