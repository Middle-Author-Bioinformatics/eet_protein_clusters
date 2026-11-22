"""Collection-wide plots and an HTML report for standalone protein FASTA input."""
import html
from pathlib import Path


def draw_prevalence(rows,path,max_domains=30):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    shown=sorted(rows,key=lambda r:(-r['proteins_with_domain'],r['accession']))[:max_domains]
    fig,ax=plt.subplots(figsize=(11,max(3,.32*len(shown)+1.8)))
    if shown:
        labels=[f'{r["name"]} ({r["accession"]})' for r in shown]
        values=[r['percent_proteins'] for r in shown]
        ax.barh(range(len(shown)),values,color='#248c84')
        ax.set_yticks(range(len(shown)),labels,fontsize=8)
        ax.invert_yaxis()
        for i,r in enumerate(shown):
            ax.text(values[i]+.3,i,f'{r["proteins_with_domain"]:,}',va='center',fontsize=8)
        ax.set_xlim(0,max(values)*1.18 if max(values) else 1)
    else:
        ax.text(.5,.5,'No domains detected under the selected thresholds',ha='center',va='center',transform=ax.transAxes)
        ax.set_yticks([])
    ax.set_xlabel('Input proteins carrying the domain (%)')
    ax.set_title(f'Domain prevalence across the protein collection\nShowing {len(shown)} of {len(rows)} detected domain families')
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(path,format='svg',bbox_inches='tight');plt.close(fig)


def _svg(path):
    if not path.exists():return ''
    content=path.read_text()
    return content[content.index('<svg'):]


def render_collection_report(out,summary,domains):
    """All domain families appear in the table, including those omitted from plots."""
    def esc(value):return html.escape(str(value),quote=True)
    def number(value):
        if value is None or value=='NA':return 'Unavailable'
        return f'{int(value):,}'
    cards=[('Input proteins',summary['total_proteins']),('With detected domains',summary.get('proteins_with_domains')),
           ('Distinct domain families',summary['distinct_domains']),('Detected domain copies',summary['domain_copies'])]
    text='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Protein collection domain profile</title>
<style>body{font:15px system-ui;margin:0;background:#f5f8f8;color:#203c40}main{max-width:1300px;margin:auto;padding:28px}h1,h2{color:#155c60}header p{color:#527074}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px}.card,section{background:white;padding:18px;border-radius:8px;margin:18px 0;box-shadow:0 1px 3px #0001}.card{margin:0}.card strong{display:block;font-size:30px;color:#176b69}.card span{font-size:13px}.scroll{overflow:auto;max-height:620px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;text-align:left;border-bottom:1px solid #e1e9e9}th{position:sticky;top:0;background:#e2eeee}input{padding:10px;border:1px solid #bed4d4;border-radius:4px;min-width:280px;margin-bottom:12px}svg{max-width:100%;height:auto}a{color:#086e75}nav a{display:inline-block;margin:6px 16px 6px 0}.muted{color:#61777a}.note{padding:10px;border-left:3px solid #278d85;background:#eef7f5}</style></head><body><main><header><h1>Protein collection domain profile</h1>'''
    text+=f'<p>{esc(summary["dataset"])} · Single protein FASTA · {number(summary["total_proteins"])} input proteins</p></header>'
    text+='<div class="cards">'+''.join(f'<div class="card"><strong>{number(v)}</strong><span>{esc(k)}</span></div>' for k,v in cards)+'</div>'
    text+='<p class="note">Prevalence uses every input protein, including proteins without detected domains. Repeated copies contribute to copy counts; each domain contributes once per protein to co-occurrence. These are descriptive annotations, not proof of EET activity or physical interaction.</p>'
    if summary.get('proteins_with_domains') in [None,'NA']:
        text+='<p class="muted">The number of proteins with any detected domain cannot be recovered from the supplied aggregate tables. It is unavailable in this reconstructed report.</p>'
    files=[('collection_summary.tsv','Collection summary'),('domain_prevalence.tsv','Domain prevalence'),('protein_architectures.tsv','Protein architectures'),('domain_hits.tsv','Domain hits'),('cooccurrence_edges.tsv','Co-occurrence edges'),('domain_network.graphml','GraphML'),('domain_prevalence.svg','Prevalence plot'),('domain_network.svg','Network SVG'),('manifest.json','Manifest')]
    text+='<nav>'+''.join(f'<a download href="{name}">{caption}</a>' for name,caption in files if (out/name).exists())+'</nav>'
    if (out/'domain_network.html').exists():text+='<p><a href="domain_network.html">Open interactive domain network</a> | <a href="domain_network_communities.tsv">Domain names and communities</a></p>'
    plot=_svg(out/'domain_prevalence.svg')
    if plot:text+='<section><h2>Most prevalent domains</h2>'+plot+'</section>'
    text+='<section><h2>Domain breakdown</h2><p>Protein prevalence and copy counts are separate: a protein carrying three copies of a domain contributes one protein and three copies.</p><label for="domains-search">Find a domain </label><input id="domains-search" placeholder="Accession or domain name"><div class="scroll"><table><thead><tr><th>Domain accession</th><th>Name</th><th>Proteins</th><th>Prevalence</th><th>Copies</th><th>Copies / positive protein</th></tr></thead><tbody>'
    for r in sorted(domains,key=lambda r:(-r['proteins_with_domain'],r['accession'])):
        text+=f'<tr><td>{esc(r["accession"])}</td><td>{esc(r["name"])}</td><td>{number(r["proteins_with_domain"])}</td><td>{r["percent_proteins"]:.2f}%</td><td>{number(r["domain_copies"])}</td><td>{r["mean_copies_when_present"]:.2f}</td></tr>'
    if not domains:text+='<tr><td colspan="6">No domains detected under the selected models and thresholds.</td></tr>'
    text+='</tbody></table></div></section>'
    network=_svg(out/'domain_network.svg')
    text+='<section><h2>Within-protein domain co-occurrence</h2><p>An edge joins different domains found in the same protein. Node size represents protein prevalence; line width represents Jaccard similarity. The plot may show a subset for readability; the edge table and GraphML contain the full exported data.</p>'
    text+=network if network else '<p>SVG rendering was not requested. Network data remain available as GraphML and TSV.</p>'
    text+='</section><script>document.getElementById("domains-search").addEventListener("input",e=>document.querySelectorAll("tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(e.target.value.toLowerCase())));</script></main></body></html>'
    (out/'index.html').write_text(text)
