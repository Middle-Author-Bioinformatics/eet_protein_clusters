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
    nodes=[dict(id='A',label='SLH (PF00395)',name='S-layer homology',group='community_1',size=10,proteins_with_domain=3),
           dict(id='B',label='Heme binding (PF_TEST)',name='Heme binding',group='community_1',size=10),
           dict(id='C',label='Isolated domain (PF_ISO)',name='Isolated domain',group='community_2',size=8)]
    edge=dict(source='A',target='B',value=.4,jaccard=.4,lift=2,proteins_with_both=2)
    write_html(target,nodes,[edge],'Domain network','domains','Same-protein co-occurrence')
    with pw.sync_playwright() as p:
        browser=p.chromium.launch(executable_path=executable,headless=True,args=['--no-sandbox'])
        context=browser.new_context(viewport={'width':1440,'height':1000},offline=True,accept_downloads=True)
        page=context.new_page();errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(target.as_uri());page.wait_for_function('window.networkView && document.querySelector("canvas").width>1')
        assert page.evaluate('networkView.getVisible().nodes.length')==3
        page.locator('#search').fill('PF00395');page.locator('#include-neighbors').uncheck()
        assert page.evaluate('networkView.getVisible().nodes.length')==1
        page.locator('#include-neighbors').check();assert page.evaluate('networkView.getVisible().nodes.length')==2
        page.locator('#reset').click();page.locator('#threshold').fill('0.8');page.locator('#hide-isolated').check()
        assert page.evaluate('networkView.getVisible().nodes.length')==0
        page.locator('#reset').click();page.locator('#group').select_option('community_2')
        assert page.evaluate('networkView.getVisible().nodes.length')==1
        page.locator('#reset').click();page.locator('#node-list button').first.click()
        assert 'SLH (PF00395)' in page.locator('#details').inner_text()
        assert 'Shared proteins: 2' in page.locator('#neighbors').inner_text()
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
            folder=Path(os.environ['NETWORK_TEST_ARTIFACTS']);folder.mkdir(parents=True,exist_ok=True);page.screenshot(path=folder/'domain-network.png',full_page=True)
        assert not errors;browser.close()
