# Interactive protein and domain networks

HTML plots are generated automatically by `cluster_proteins.py`,
`cluster_domains.py`, and `genome_pipeline.py`. Each network HTML includes its
JavaScript, data and initial layout. Open the file directly in a modern browser;
no server, internet connection, CDN, web submission or new runtime dependency is
needed. Keep the repository modules, `network_viewer.html`, and
`network_viewer.js` together when running the scripts. Generated network HTML
files can then be shared individually. Inline JSON/PNG downloads keep working;
links to neighboring TSV/GraphML files require sharing those files as well.

## Existing analyses: regenerate reports without searches

After updating the repository, run inside the `eet-clusters` environment:

```bash
python /home/ark/MAB/bin/eet_protein_clusters/network_reports.py \
  --clustering /path/to/MCL2 \
  --domains /path/to/domain_results
```

Either folder is optional. This reads the stored protein graph, membership,
domain GraphML, `protein_architectures.tsv` and model-name inventory; it does not run DIAMOND, MCL or an HMM
search. For older protein outputs that stored only edge weights, it reconstructs
identity/coverage metadata from checksum-verified saved DIAMOND hits, retaining
MCL membership. If those files/checksums are missing, it reports the needed input
rather than inventing identity values. A report-only manifest records parameters
and the source analysis manifest checksum.

### Show annotations from updated FASTA headers

Replacing FASTA files alone does **not** update an existing HTML file. The saved
architecture table supplies domain carriers and does not contain NR descriptions.
Use the annotated FASTAs as an explicit display-label source, then regenerate:

```bash
python network_reports.py \
  --clustering /path/to/MCL2-test \
  --domains /path/to/domain_results \
  --protein-fastas /path/to/MCL2-test/clusters_nr_annotated
```

There is no need to replace the original cluster folder. If you do replace
`MCL2-test/clusters` and supply `--clustering MCL2-test`, report-only regeneration
automatically reads the current headers there when `--protein-fastas` is omitted.
For domain-only regeneration, supply `--protein-fastas` explicitly.

The first whitespace-delimited protein IDs must stay unchanged and globally
unique. Descriptions appended after those IDs are imported as labels; no graph,
domain hit, cluster membership or sequence is recomputed. Renamed FASTA filenames
are acceptable because matching uses protein IDs. An import with no matching IDs
fails with an actionable error. Partial imports keep unmatched proteins as IDs
and report matched/total counts in the log and report manifest.

In the shared-protein list, the primary label is the imported description, with
the original ID underneath. Both names and IDs are searchable. NR accessions,
full headers, domains, architectures and cluster IDs are available under each
protein's **Details** and in TSV downloads. No-hit/missing-title headers retain
their status and display the original ID. Plain descriptive FASTA headers are
also supported; the annotation source is recorded.

`--protein-fastas` also adds annotations to protein-identity node details, the
node TSV and complete `--pair ... --proteins-out ...` exports. Header imports are
checksum-keyed; paths, checksums and mapping counts are recorded in
`network_report_manifest.json`. Editing a header refreshes the metadata cache.
Original source/architecture tables and protein IDs are never rewritten.

This option imports FASTA headers. The separate LLM summarizer produces a
`cluster_annotations.tsv`/JSONL, not renamed FASTAs; that TSV is not implicitly
treated as a protein-header source.

### Protein cluster IDs and domain highlighting

Every carrier row now shows its protein cluster ID as a visible badge, beside
the protein name/ID. It is also searchable and retained in TSV downloads. Missing
assignments are labeled **Cluster unavailable**, rather than guessed from an ID
or a renamed FASTA filename.

In the pooled `all_domains_network.html`, use **Highlight cluster** to select a
protein cluster. Its domains receive orange outlines; other domains fade while
the graph remains visible. Domain-community colors, node sizes, edges and
Jaccard/lift values retain their original meaning. The existing community filter
is separate from this new protein-cluster selector. Reset clears both controls.

The carrier list continues to include all clusters for the selected domain/pair;
members of the highlighted cluster move to the top and their badges are accented.
Search the list by cluster ID to narrow it. The selected-node card shows both the
global carrier count and the count in the highlighted cluster. The toolbar states
visible versus total domain-family counts for that cluster; filters/export caps
can hide domains. Clusters with no detected domains remain selectable and show
zero, and membership counts include analyzed domain-negative proteins.

**Highlighting means domain presence, not within-cluster co-occurrence.** A cluster
can contain domain A in one protein and domain B in another. Both domains will
be highlighted, while an A–B edge in the pooled graph could come entirely from
a different cluster. The selector does not create or reassign edges. To inspect
within-cluster co-occurrence itself, open that cluster's own `domain_network.html`.

The full saved architecture table supplies the cluster/domain presence index,
independently of `--network-max-proteins`. Even a zero-protein HTML preview can
highlight the correct domains. Only domains embedded in the graph can be drawn;
the status displays the full graph's domain count to expose node truncation.
Per-cluster HTML files restrict the selector to their own analysis scope. Complete
cluster/domain index metadata can increase HTML size for very large collections.

Cluster-folder domain analyses already carry their protein-cluster assignments.
If you originally searched one **unclustered protein FASTA**, attach its later
MCL assignments during report regeneration by supplying `--clustering`:

```bash
python network_reports.py \
  --clustering /path/to/MCL2-test \
  --domains /path/to/single_fasta_domain_results
```

This reads `MCL2-test/membership.tsv` by original protein ID and preserves the
existing domain search, architecture table and graph statistics. You may also
regenerate only the domain views using a saved membership file:

```bash
python network_reports.py --domains /path/to/domain_results \
  --membership /path/to/MCL2-test/membership.tsv
```

`--membership` is for domain-only reports; with `--clustering`, its saved mapping
is used automatically. The TSV requires `protein_id` and `cluster` columns and
unique protein IDs. Duplicate, empty, nonmatching, or conflicting assignments
are rejected. If the domain analysis already has cluster labels, an imported
mapping must agree with them. Partial mappings preserve known assignments and
record unassigned-protein counts. The report manifest records mapping provenance
and its checksum. Complete `--pair ... --cluster ...` exports also support imported
membership for a single-FASTA domain analysis.

Without saved cluster assignments, the selector is disabled with an explanation.
An aggregate GraphML alone cannot reconstruct protein/cluster/domain mappings;
the original `protein_architectures.tsv` is still required. No sequence search or
reclustering is performed by this display feature.

To focus on a particular MCL family:

```bash
python /home/ark/MAB/bin/eet_protein_clusters/network_reports.py \
  --clustering /path/to/MCL2 --domains /path/to/domain_results \
  --cluster cluster_000123 \
  --network-max-nodes 2000 --network-max-edges 10000
```

The protein HTML/preview TSVs are replaced; the underlying MCL membership stays
the same. For old domain outputs, report-only mode creates interactive network
files next to the GraphML. Regenerate `cluster_domains.py --resume` to also refresh
its report index, named architecture lists, pooled network and summary columns;
checksum-matching HMM searches are reused.

## New analyses

```bash
python cluster_proteins.py \
  --input AllMHCGenes.faa --output MCL2 \
  --threads 32 --inflation 2.0 --block-size 0.4 --resume

python cluster_domains.py \
  --input MCL2/clusters --hmm /path/to/Pfam-A.hmm \
  --output domain_results --threads 32 --resume --clustering MCL2

# Add domain names to the protein network's node details:
python network_reports.py --clustering MCL2 --domains domain_results
```

The genome wrapper automatically regenerates the protein network after domain
analysis, attaching its domain names and architectures. Its root report links
both types of interactive network. Standalone protein FASTA domain analysis also
gets a native collection-wide interactive domain network.

## Protein amino-acid similarity network

Open `MCL2/protein_network.html` (also linked from `MCL2/index.html`).
When both folders are supplied to `network_reports.py`, the protein and domain
views link to each other. The genome wrapper adds these links automatically.
For a new domain search, `cluster_domains.py --clustering MCL2` adds the return
link; this option does not rerun clustering or compute identity from domains.

- Each node is a protein; node color is its existing MCL cluster. Isolated proteins
  are retained when included in the preview.
- Each edge is a qualifying protein-pair alignment used in the MCL input graph.
  Details include percent amino-acid identity, coverage of each endpoint, E-value,
  bit score and MCL edge weight.
- Directional/duplicate alignments are collapsed by the same maximum graph weight
  used for clustering. Metadata follows that selected HSP; equal graph weights
  prefer the higher bit score. Coverage A/B follows the exported source/target
  endpoint order even when the selected hit ran in the opposite direction.
- The graph still clusters using the configured weight (normally capped
  -log10 E-value), with the configured identity/coverage filters. Displaying or
  filtering identity does not rerun MCL or change its memberships.

Here “amino-acid identity” means pairwise protein alignment identity. **Genome
AAI** averages identities across shared proteins between two genomes and is a
different analysis; this wrapper does not compute it. A cluster color is not an
AAI threshold or a validated functional label.

`protein_network_nodes.tsv` and `protein_network_edges.tsv` describe the embedded
preview. Complete membership remains in `membership.tsv`; complete qualifying
edges and selected alignment metadata remain in `work/proteins.sqlite`, with
original search hits in the existing alignment output. This avoids creating a
second enormous all-pairs text export by default.

## Domain co-occurrence networks and names

- Per protein cluster: `domain_results/clusters/<cluster>/domain_network.html`.
- Across all input cluster members: `domain_results/all_domains_network.html`.
- Single FASTA mode: `domain_results/domain_network.html`.

Nodes display **HMM domain names plus accessions**, e.g. `SLH (PF00395.29)`.
Names come from the supplied HMM model inventory. They are preserved for Pfam
and custom HMMs; unknown/DUF model names are not replaced with guessed functions.
Domain nodes retain counts, repeated-domain copies and prevalence. Edge details
include names of both domains, shared-protein count, Jaccard similarity and lift.
Click a Pfam node and use **Read about [name]** in its details to open that
family's InterPro/Pfam entry. Model versions are retained in the analysis but
removed from the page URL. Custom HMM accessions do not receive invented links.

### Which proteins carry these domains?

Click a domain node to list its carrier proteins, with IDs, MCL cluster (when
available), and named architectures. Under a connection, click **Show shared
proteins** to list proteins containing **both** domains. Search this list or
download its embedded matches as TSV. The screen renders the first 100 search
matches; the TSV includes all embedded matches for the selected node/pair,
regardless of the carrier-list search. Repeated copies do not duplicate IDs.
Lists reflect the graph's protein collection/cluster, independently of the
currently visible community or Jaccard filters.

The right pane starts with the domain name/Pfam link, protein count, and compact
connection buttons. Click a connection's **N proteins** button to open its shared
list. More than eight connections are collapsed under **More connections**. Raw
node statistics and edge metrics remain under **More details**; architectures
and provenance are under each carrier's **Details**. Browse-node and group lists
are collapsed initially. Complete/partial counts remain visible, with longer
explanations under **About this list**. JSON and TSV downloads retain the full
embedded scientific metadata.

Carrier identities come from the selected `domain_set` column in the original
`protein_architectures.tsv`, not from aggregate edge counts. GraphML by itself
cannot reconstruct protein IDs: those views explicitly label protein lists as
unavailable. Missing or inconsistent architecture counts also produce an
unavailable message. Old architecture tables without domain-name columns remain
usable, with accessions as the fallback. A checksum-keyed SQLite lookup in
`work/network_members.sqlite` is rebuilt automatically when architectures change.

The offline HTML embeds at most **20,000 unique domain-carrying proteins** by
default, selected in protein-ID order among carriers of embedded domains.
`--network-max-proteins` changes this cap; `0` disables embedding identities.
Each list states its embedded/full count and labels partial lists. Complete
architectures remain downloadable from the report folder. To export every
protein sharing a particular pair without any HTML limit or new searches:

```bash
python network_reports.py --domains /path/to/domain_results \
  --pair PF00395 PF01476 --proteins-out SLH_LysM_proteins.tsv

# Restrict to one MCL cluster (only for cluster-folder analyses):
python network_reports.py --domains /path/to/domain_results \
  --pair PF00395 PF01476 --cluster cluster_000123 \
  --proteins-out cluster_000123_SLH_LysM.tsv
```

This export mode skips HTML refresh. Use exact versioned accessions if multiple
versions match; an unversioned accession is accepted when it resolves uniquely.
Single-FASTA results retain their native collection format, without cluster
columns in the exported pair table.

An edge means both different domains occur **inside the same protein**. Repeated
copies contribute to copy counts, but count once per protein for co-occurrence.
Separate proteins carrying A and B do not create an edge, including in the pooled
network. Domain-negative proteins remain in the denominator. These statistics are
descriptive associations, not physical interaction or evidence of EET activity.

Color represents **domain co-occurrence communities** computed with weighted
Louvain (Jaccard weights, resolution 1, seed 17) on the complete exported graph
after `--min-pair-count`. These differ from protein MCL clusters. An isolated
domain gets its own community. Display thresholds do not recompute communities.
Community membership and domain names are exported in
`domain_network_communities.tsv` (or `all_domains_network_communities.tsv`) and
GraphML node attributes. GraphML retains the complete exported network.

Additional named lists are included in normal domain output:

- `protein_architectures.tsv`: `domain_names` and `architecture_names` retain
  repeat order and hit coordinates alongside the accession-based columns.
- `cluster_summary.tsv` or `collection_summary.tsv`: `top_domain_name` and
  the cluster/collection's `domain_names` list.
- `cooccurrence_edges.tsv`: `domain_a_name` and `domain_b_name` alongside accessions.
- Genome wrapper `protein_annotations.tsv`: domain names and named architectures.

## Controls and scale

Search matches protein IDs, original headers, domain names/accessions and attached
annotations. Filter by MCL cluster/domain community; optionally include matching
nodes' neighbors, hide isolates, filter identity/Jaccard or require more shared
proteins. Click nodes for details and visible neighbors. Drag nodes, pan on blank
canvas, scroll to zoom, fit/reset the view, optionally relax the layout, download
visible graph JSON, a named node TSV, or export a PNG.

HTML previews default to **1,000 nodes and 5,000 edges** per network. Every file
states embedded vs complete counts and its selection method. Protein previews
sample members round-robin from the largest MCL clusters, then retain the strongest
eligible edges. Domain previews retain the most prevalent families and strongest
shared-protein edges; communities are computed before this display subset.
Use `--network-max-nodes`, `--network-max-edges` and protein `--network-cluster`
(or report-only `--cluster`) to inspect another subset. Sampling/truncation can
hide links, so displayed topology should be interpreted as a preview.

`--no-plots` skips static SVG rendering while keeping interactive HTML.
`--no-interactive` skips HTML network files, while retaining analysis tables,
GraphML and domain community tables. The bounded viewers prevent embedding a
whole million-protein graph, but many per-cluster HTML files still take disk space.
Report-only regeneration can be used for a focused protein family when needed.

Tests use synthetic models/graphs and real open-source search engines. Actual
Chromium interaction tests run in CI, including offline file loading, search,
threshold/group filtering, reset, selection, zoom, dragging and downloads. Browser
validation requires development dependencies only (`requirements-dev.txt`), not
installation of a browser in the analysis environment.
