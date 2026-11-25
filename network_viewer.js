'use strict';
(() => {
const data=JSON.parse(document.getElementById('network-data').textContent), el=id=>document.getElementById(id);
const canvas=el('network'),ctx=canvas.getContext('2d'),plot=el('plot'),isDomain=data.kind==='domains';
const byId=new Map(data.nodes.map(n=>[n.id,n])),adj=new Map(data.nodes.map(n=>[n.id,new Set()]));
for(const e of data.edges){adj.get(e.source)?.add(e.target);adj.get(e.target)?.add(e.source);}
let visibleNodes=[],visibleEdges=[],selected=null,selectedPair=null,memberRows=[],camera={x:0,y:0,scale:1},width=1,height=1,drag=null,relaxing=false;
const hash=s=>[...String(s)].reduce((h,c)=>(Math.imul(h,31)+c.charCodeAt(0))|0,17)>>>0;
const color=g=>{const suffix=String(g).match(/(\d+)$/),hue=suffix?(Number(suffix[1])%1000000*137.507764+195)%360:hash(g)%360;return `hsl(${hue} 48% 44%)`;};
const groupCounts=new Map();for(const n of data.nodes)groupCounts.set(n.group,(groupCounts.get(n.group)||0)+1);
for(const g of [...groupCounts.keys()].sort()){const o=document.createElement('option');o.value=g;o.textContent=`${g} (${groupCounts.get(g)})`;el('group').append(o);}
el('threshold-label').textContent=isDomain?'Minimum Jaccard':'Minimum identity (%)';el('threshold').max=isDomain?'1':'100';el('threshold').step=isDomain?'0.01':'1';el('pair-label').hidden=!isDomain;
const proteinNoun=n=>Number(n)===1?'protein':'proteins';
const clusterIndex=data.members?.cluster_index,clusterById=new Map((clusterIndex?.clusters||[]).map(c=>[c.id,c]));
el('protein-cluster-control').hidden=!isDomain;
el('group').options[0].textContent=isDomain?'All domain communities':'All MCL clusters';
el('group').setAttribute('aria-label',isDomain?'Domain community':'MCL cluster');
if(isDomain&&clusterIndex?.available&&clusterById.size){
 for(const c of clusterById.values()){const option=document.createElement('option');option.value=c.id;option.textContent=`${c.id} (${c.proteins} ${proteinNoun(c.proteins)})`;el('protein-cluster').append(option);}
 el('protein-cluster').title='Highlight domains present in this protein cluster. Colors and pooled co-occurrence edges stay unchanged.';
}else{
 el('protein-cluster').disabled=true;el('protein-cluster').options[0].textContent='Cluster IDs unavailable';
 el('protein-cluster').title=clusterIndex?.reason||(!data.members?.available&&data.members?.reason)||'Regenerate with cluster assignments in protein_architectures.tsv, --clustering, or --membership.';
}
function activeCluster(){return isDomain?clusterById.get(el('protein-cluster').value):null;}
function highlightedDomains(){return new Set(Object.keys(activeCluster()?.domains||{}));}
function highlightStatus(){
 const cluster=activeCluster(),status=el('cluster-status');status.hidden=!cluster;
 if(cluster){const visible=visibleNodes.filter(n=>Object.hasOwn(cluster.domains,n.id)).length;status.textContent=`${cluster.id} · ${cluster.proteins} analyzed ${proteinNoun(cluster.proteins)} · ${visible} / ${cluster.total_domains} domains highlighted. Domain presence only; edges stay unchanged.`;}
}
el('scope').textContent=`Exported ${data.nodes.length.toLocaleString()} / ${data.total_nodes.toLocaleString()} nodes and ${data.edges.length.toLocaleString()} / ${data.total_edges.toLocaleString()} edges. ${data.selection}`;
function screen(n){return{x:n.x*camera.scale+camera.x,y:n.y*camera.scale+camera.y};}
function fit(){
 if(!visibleNodes.length)return;
 const xs=visibleNodes.map(n=>n.x),ys=visibleNodes.map(n=>n.y),xmin=Math.min(...xs)-40,xmax=Math.max(...xs)+40,ymin=Math.min(...ys)-40,ymax=Math.max(...ys)+40;
 camera.scale=Math.min(4,Math.max(.005,Math.min(width/(xmax-xmin),height/(ymax-ymin))*.88));camera.x=width/2-(xmin+xmax)/2*camera.scale;camera.y=height/2-(ymin+ymax)/2*camera.scale;draw();
}
function draw(){
 ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);ctx.fillStyle='white';ctx.fillRect(0,0,width,height);
 const neighbors=selected?adj.get(selected)||new Set():null,labelBoxes=[],cluster=activeCluster(),highlighted=highlightedDomains();
 for(const e of visibleEdges){const a=screen(byId.get(e.source)),b=screen(byId.get(e.target));ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.lineWidth=.4+Math.min(2.5,(isDomain?e.value:e.value/100)*2.2);ctx.strokeStyle=selected?(e.source===selected||e.target===selected?'#227c85':'#dfe9e960'):'#8aa5ad70';ctx.stroke();}
 for(const n of visibleNodes){const p=screen(n),r=Math.max(3,Math.min(22,(n.size||6)*Math.sqrt(camera.scale)));ctx.beginPath();ctx.arc(p.x,p.y,r,0,2*Math.PI);ctx.fillStyle=color(n.group);ctx.globalAlpha=(selected&&n.id!==selected&&!neighbors.has(n.id)?.28:1)*(cluster&&!highlighted.has(n.id)?.18:1);ctx.fill();ctx.globalAlpha=1;
 if(cluster&&highlighted.has(n.id)){ctx.beginPath();ctx.arc(p.x,p.y,r+4,0,2*Math.PI);ctx.strokeStyle='#cf7b1b';ctx.lineWidth=3;ctx.stroke();}
 if(n.id===selected){ctx.beginPath();ctx.arc(p.x,p.y,r,0,2*Math.PI);ctx.strokeStyle='#132f3a';ctx.lineWidth=3;ctx.stroke();}
 if(el('labels').checked&&(visibleNodes.length<=160||n.id===selected||(el('search').value&&n._matches))){ctx.font=n.id===selected?'bold 12px system-ui':'11px system-ui';const box={x:p.x+r+3,y:p.y-9,w:ctx.measureText(n.label).width,h:14};const collision=labelBoxes.some(b=>box.x<b.x+b.w&&box.x+box.w>b.x&&box.y<b.y+b.h&&box.y+box.h>b.y);if(!collision||n.id===selected||(el('search').value&&n._matches)){ctx.globalAlpha=cluster&&!highlighted.has(n.id)&&n.id!==selected?.3:1;ctx.fillStyle='#203e48';ctx.fillText(n.label,box.x,p.y+4);labelBoxes.push(box);ctx.globalAlpha=1;}}
 }
}
function compactDetails(parent, entries){
 const dl=document.createElement('dl');
 for(const[k,v]of entries){if(v===undefined||v===null||v==='')continue;const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=k.replaceAll('_',' ');dd.textContent=String(v);dl.append(dt,dd);}
 parent.append(dl);
}
function renderMembers(){
 const panel=el('protein-members');panel.hidden=!isDomain||!selected;
 if(panel.hidden)return;
 const n=byId.get(selected),pair=selectedPair,ids=pair?[pair.source,pair.target]:[selected];
 el('member-title').textContent=pair?'Shared proteins':'Proteins';
 el('member-context').textContent=pair?ids.map(id=>byId.get(id).name||byId.get(id).label).join(' + '):'';
 const members=data.members,expected=Number(pair?pair.proteins_with_both:n.proteins_with_domain);
 memberRows=[];el('member-note').textContent='';
 if(!members?.available){
  el('member-status').textContent='Protein list unavailable';
  el('member-note').textContent=members?.reason||'Regenerate with the original protein_architectures.tsv to recover carrier IDs.';
 }else{
  let indices=members.by_domain[ids[0]]||[];
  if(pair){const other=new Set(members.by_domain[ids[1]]||[]);indices=indices.filter(i=>other.has(i));}
  memberRows=indices.map(i=>members.proteins[i]);
  const complete=memberRows.length===expected;
  el('member-status').textContent=`${memberRows.length.toLocaleString()} / ${expected.toLocaleString()} proteins${complete?'':' · Partial list'}`;
  if(!complete)el('member-note').textContent='The HTML protein limit applies. Missing from this preview does not mean missing from the analysis. Use the complete architecture table or network_reports.py --pair for all matches.';
 }
 el('member-notes').hidden=!el('member-note').textContent;
 el('member-search').hidden=!members?.available;el('download-members').hidden=!members?.available;
 el('download-members').textContent='Download TSV';
 el('download-members').title='Download all embedded matches for this domain or pair, including matches hidden by this search.';
 renderMemberList();
}
function renderMemberList(){
 const q=el('member-search').value.trim().toLowerCase(),cluster=activeCluster(),rows=memberRows.filter(r=>Object.values(r).some(v=>String(v).toLowerCase().includes(q)));
 if(cluster)rows.sort((a,b)=>Number(b.cluster===cluster.id)-Number(a.cluster===cluster.id)||a.protein_id.localeCompare(b.protein_id));
 const list=el('member-list');list.replaceChildren();
 for(const row of rows.slice(0,100)){
  const li=document.createElement('li'),name=document.createElement('strong');
  name.className='member-name';name.textContent=row.annotation||row.protein_id;name.title=row.annotation||row.protein_id;li.append(name);
  if(row.annotation){const id=document.createElement('p');id.className='muted member-id';id.textContent=row.protein_id;li.append(id);}
  const badge=document.createElement('span');badge.className='cluster-badge'+(cluster&&row.cluster===cluster.id?' active-cluster':'')+(!row.cluster?' unknown-cluster':'');badge.textContent=row.cluster||'Cluster unavailable';badge.title='Protein cluster';li.append(badge);
  const details=document.createElement('details'),summary=document.createElement('summary');summary.textContent='Details';details.append(summary);
  compactDetails(details,[['protein_id',row.protein_id],['cluster',row.cluster],['domains',row.domain_names],['architecture',row.architecture_names],['NR hit',row.nr_hit],['annotation source',row.annotation_status],['full header',row.annotated_header]]);
  li.append(details);list.append(li);
 }
 const status=el('member-list-status');status.hidden=!q&&rows.length<=100;
 status.textContent=`Showing ${Math.min(100,rows.length)} / ${rows.length} matches. TSV includes all ${memberRows.length} embedded proteins.`;
}
function selectedDetails(){
 const card=el('details'),stats=el('statistics'),panel=el('neighbors');card.replaceChildren();stats.replaceChildren();panel.replaceChildren();
 el('technical-details').hidden=!selected;
 if(!selected||!byId.has(selected)){const p=document.createElement('p');p.className='muted';p.textContent='Click a node to explore.';card.append(p);renderMembers();return;}
 const n=byId.get(selected),name=document.createElement('p');name.className='selection-name';
 const accession=String(n.accession||n.id),match=isDomain&&accession.match(/^(PF[0-9]{5})(?:\.[0-9]+)?$/);
 if(isDomain&&n.domain_source==='CDD'&&/^[0-9]+$/.test(String(n.cdd_pssm_id))){const a=document.createElement('a');a.id='cdd-link';a.href=`https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=${n.cdd_pssm_id}`;a.target='_blank';a.rel='noopener noreferrer';a.textContent=n.label;a.title='Read this CDD model on NCBI';name.append(a);}
 else if(match){const a=document.createElement('a');a.id='pfam-link';a.href=`https://www.ebi.ac.uk/interpro/entry/pfam/${match[1]}/`;a.target='_blank';a.rel='noopener noreferrer';a.textContent=n.label;a.title='Read this Pfam entry on InterPro';name.append(a);}else name.textContent=n.label;
 card.append(name);
 if(!isDomain&&n.annotation){const annotation=document.createElement('p');annotation.className='selection-annotation';annotation.textContent=n.annotation;card.append(annotation);}
 const count=document.createElement('p');count.className='muted';count.textContent=isDomain?`${Number(n.proteins_with_domain).toLocaleString()} ${proteinNoun(n.proteins_with_domain)}`:[n.length_aa?`${n.length_aa} aa`:'',n.cluster].filter(Boolean).join(' · ');if(isDomain&&activeCluster()){const c=activeCluster();count.textContent+=` · ${c.domains[n.id]||0} in ${c.id}`;}card.append(count);
 compactDetails(stats,Object.entries(n).filter(([k])=>!['x','y','size','_matches','label'].includes(k)));
 const attached=visibleEdges.filter(e=>e.source===selected||e.target===selected),h=document.createElement('h2');h.textContent=`Connections (${attached.length})`;panel.append(h);
 let extra=null;if(attached.length>8){extra=document.createElement('details');const summary=document.createElement('summary');summary.textContent=`More connections (${Math.min(100,attached.length)-8})`;extra.append(summary);}
 for(const [index,e] of attached.slice(0,100).entries()){
  const target=e.source===selected?e.target:e.source,p=document.createElement('div'),b=document.createElement('button');p.className='connection';b.className='neighbor-target';b.textContent=byId.get(target).label;b.onclick=()=>selectNode(target);p.append(b);
  const metrics=isDomain?`Jaccard: ${Number(e.jaccard).toFixed(3)}; lift: ${Number(e.lift).toFixed(2)}`:`Identity: ${e.identity_percent}%; coverage A/B: ${(100*e.coverage_a).toFixed(1)}% / ${(100*e.coverage_b).toFixed(1)}%; E-value: ${e.evalue}; weight: ${Number(e.weight).toFixed(2)}`;
  p.title=metrics;
  if(isDomain){const shared=document.createElement('button');shared.className='shared-proteins';shared.textContent=`${e.proteins_with_both} ${proteinNoun(e.proteins_with_both)}`;shared.title='Show proteins containing both domains';shared.onclick=()=>{selectedPair=e;el('member-search').value='';renderMembers();el('protein-members').scrollIntoView({block:'nearest',behavior:'smooth'});};p.append(shared);}
  else{const identity=document.createElement('span');identity.className='muted';identity.textContent=`${e.identity_percent}% identity`;p.append(identity);}
  (index<8?panel:extra).append(p);
  const entry=document.createElement('p');entry.textContent=`${byId.get(target).label}: ${metrics}`;stats.append(entry);
 }
 if(extra)panel.append(extra);
 if(attached.length>100){const more=document.createElement('p');more.className='muted';more.textContent='Showing 100 connections; filter the network to narrow the list.';panel.append(more);}
 if(selectedPair&&!attached.includes(selectedPair))selectedPair=null;
 renderMembers();
}
function selectNode(id){selected=id;selectedPair=null;el('member-search').value='';selectedDetails();draw();}
function update(autoFit=false){
 const q=el('search').value.trim().toLowerCase(),g=el('group').value,minimum=Number(el('threshold').value),pairs=Math.max(1,Number(el('pair').value)||1);
 el('threshold-value').textContent=isDomain?minimum.toFixed(2):String(minimum);
 let matches=new Set();for(const n of data.nodes){n._matches=!q||Object.entries(n).filter(([k])=>!['x','y','_matches'].includes(k)).some(([,v])=>String(v).toLowerCase().includes(q));if(n._matches)matches.add(n.id);}
 const edges=data.edges.filter(e=>e.value>=minimum&&(!isDomain||e.proteins_with_both>=pairs));
 if(q&&el('include-neighbors').checked){const near=new Set(matches);for(const e of edges){if(matches.has(e.source))near.add(e.target);if(matches.has(e.target))near.add(e.source);}matches=near;}
 visibleNodes=data.nodes.filter(n=>(!g||n.group===g)&&matches.has(n.id));let ids=new Set(visibleNodes.map(n=>n.id));visibleEdges=edges.filter(e=>ids.has(e.source)&&ids.has(e.target));
 if(el('hide-isolated').checked){const connected=new Set(visibleEdges.flatMap(e=>[e.source,e.target]));visibleNodes=visibleNodes.filter(n=>connected.has(n.id));}
 if(selected&&!visibleNodes.some(n=>n.id===selected))selected=null;
 el('empty').hidden=visibleNodes.length>0;el('status').textContent=`Visible: ${visibleNodes.length.toLocaleString()} nodes, ${visibleEdges.length.toLocaleString()} edges | Embedded export: ${data.nodes.length.toLocaleString()} nodes, ${data.edges.length.toLocaleString()} edges. View filters do not rerun clustering.`;
 const list=el('node-list');list.replaceChildren();for(const n of visibleNodes.slice(0,100)){const b=document.createElement('button');b.textContent=n.label;b.dataset.nodeId=n.id;if(activeCluster())b.className=Object.hasOwn(activeCluster().domains,n.id)?'cluster-domain-match':'cluster-domain-other';b.onclick=()=>selectNode(n.id);list.append(b);}el('list-status').textContent=`Listed ${Math.min(100,visibleNodes.length)} of ${visibleNodes.length} visible nodes; search to find others.`;
 const legend=el('legend');legend.replaceChildren();for(const g of [...new Set(visibleNodes.map(n=>n.group))].sort().slice(0,40)){const s=document.createElement('span'),dot=document.createElement('i');dot.className='swatch';dot.style.background=color(g);s.append(dot,document.createTextNode(g));legend.append(s);}
 highlightStatus();selectedDetails();if(autoFit)fit();else draw();
}
function pick(x,y){let found=null,distance=Infinity;for(const n of visibleNodes){const p=screen(n),d=Math.hypot(x-p.x,y-p.y),r=Math.max(5,Math.min(24,(n.size||6)*Math.sqrt(camera.scale)))+4;if(d<r&&d<distance){found=n;distance=d;}}return found;}
function coords(e){const r=canvas.getBoundingClientRect();return{x:e.clientX-r.left,y:e.clientY-r.top};}
canvas.addEventListener('pointerdown',e=>{const p=coords(e);drag={...p,node:pick(p.x,p.y),moved:false};canvas.setPointerCapture(e.pointerId);});
canvas.addEventListener('pointermove',e=>{const p=coords(e);if(drag){const dx=p.x-drag.x,dy=p.y-drag.y;if(Math.abs(dx)+Math.abs(dy)>1)drag.moved=true;if(drag.node){drag.node.x+=dx/camera.scale;drag.node.y+=dy/camera.scale;}else{camera.x+=dx;camera.y+=dy;}drag.x=p.x;drag.y=p.y;draw();return;}const n=pick(p.x,p.y),tip=el('hover');tip.style.display=n?'block':'none';canvas.style.cursor=n?'pointer':'grab';if(n){tip.textContent=`${n.label}\n${n.group}`;tip.style.left=Math.max(0,Math.min(width-300,p.x+12))+'px';tip.style.top=Math.max(0,p.y-35)+'px';}});
canvas.addEventListener('pointerup',()=>{if(drag?.node&&!drag.moved)selectNode(drag.node.id);drag=null;});canvas.addEventListener('pointercancel',()=>drag=null);canvas.addEventListener('pointerleave',()=>el('hover').style.display='none');
canvas.addEventListener('wheel',e=>{e.preventDefault();const p=coords(e),old=camera.scale;camera.scale=Math.max(.005,Math.min(12,old*Math.exp(-e.deltaY*.001)));camera.x=p.x-(p.x-camera.x)*camera.scale/old;camera.y=p.y-(p.y-camera.y)*camera.scale/old;draw();},{passive:false});
function relax(){if(relaxing){relaxing=false;return;}relaxing=true;let count=0;el('relax').textContent='Stop layout';const tick=()=>{
 if(!relaxing||count++>=180){relaxing=false;el('relax').textContent='Relax layout';return;}
 const grid=new Map(),cell=70,force=new Map(visibleNodes.map(n=>[n.id,{x:0,y:0}]));for(const n of visibleNodes){const key=`${Math.floor(n.x/cell)},${Math.floor(n.y/cell)}`;if(!grid.has(key))grid.set(key,[]);grid.get(key).push(n);}
 for(const n of visibleNodes){const f=force.get(n.id),gx=Math.floor(n.x/cell),gy=Math.floor(n.y/cell);for(let x=gx-1;x<=gx+1;x++)for(let y=gy-1;y<=gy+1;y++)for(const m of grid.get(`${x},${y}`)||[]){if(n===m)continue;let dx=n.x-m.x,dy=n.y-m.y;if(dx===0&&dy===0){dx=.1;dy=.1;}const d2=dx*dx+dy*dy+10;f.x+=dx*70/d2;f.y+=dy*70/d2;}}
 for(const e of visibleEdges){const a=byId.get(e.source),b=byId.get(e.target),dx=b.x-a.x,dy=b.y-a.y,d=Math.hypot(dx,dy)||1,k=(d-65)*.008;force.get(a.id).x+=dx/d*k;force.get(a.id).y+=dy/d*k;force.get(b.id).x-=dx/d*k;force.get(b.id).y-=dy/d*k;}
 for(const n of visibleNodes){const f=force.get(n.id);n.x+=Math.max(-6,Math.min(6,f.x));n.y+=Math.max(-6,Math.min(6,f.y));}draw();requestAnimationFrame(tick);};requestAnimationFrame(tick);}
function download(content,name,type){const a=document.createElement('a'),url=URL.createObjectURL(new Blob([content],{type}));a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
for(const id of ['search','group','threshold','pair','include-neighbors','hide-isolated'])el(id).addEventListener('input',()=>update(true));el('protein-cluster').addEventListener('input',()=>update(false));el('labels').addEventListener('change',draw);el('fit').onclick=fit;el('relax').onclick=relax;
el('reset').onclick=()=>{relaxing=false;selected=null;el('relax').textContent='Relax layout';el('search').value='';el('group').value='';el('protein-cluster').value='';el('threshold').value='0';el('pair').value='1';el('include-neighbors').checked=true;el('hide-isolated').checked=false;el('labels').checked=true;update(true);};
el('download-data').onclick=()=>download(JSON.stringify({definition:data.definition,nodes:visibleNodes,edges:visibleEdges,highlighted_cluster:activeCluster()?.id||null,highlighted_domains:visibleNodes.filter(n=>highlightedDomains().has(n.id)).map(n=>n.id)},null,2),'visible_network.json','application/json');
el('download-table').onclick=()=>{const keys=[...new Set(visibleNodes.flatMap(n=>Object.keys(n)))].filter(k=>!['x','y','_matches','size'].includes(k));const cell=v=>String(v??'').replace(/[\t\r\n]/g,' ');download([keys.join('\t'),...visibleNodes.map(n=>keys.map(k=>cell(n[k])).join('\t'))].join('\n'),'visible_nodes.tsv','text/tab-separated-values');};
el('member-search').addEventListener('input',renderMemberList);
el('download-members').onclick=()=>{const keys=['protein_id',...(memberRows.some(r=>r.cluster)?['cluster']:[]),'annotation','nr_hit','annotation_status','domain_names','architecture_names','annotated_header'],cell=v=>String(v??'').replace(/[\t\r\n]/g,' ');download([keys.join('\t'),...memberRows.map(r=>keys.map(k=>cell(r[k])).join('\t'))].join('\n')+'\n','domain_carrier_proteins.tsv','text/tab-separated-values');};
el('download-image').onclick=()=>{draw();const a=document.createElement('a');a.href=canvas.toDataURL('image/png');a.download='network.png';a.click();};
new ResizeObserver(()=>{width=plot.clientWidth;height=plot.clientHeight;canvas.width=width*devicePixelRatio;canvas.height=height*devicePixelRatio;fit();}).observe(plot);
window.networkView={data,getVisible:()=>({nodes:visibleNodes,edges:visibleEdges}),selectNode,fit,getHighlight:()=>({cluster:activeCluster()?.id||null,domains:[...highlightedDomains()],visible_domains:visibleNodes.filter(n=>highlightedDomains().has(n.id)).map(n=>n.id),total_domains:activeCluster()?.total_domains||0}),getCamera:()=>({...camera})};update();
})();
