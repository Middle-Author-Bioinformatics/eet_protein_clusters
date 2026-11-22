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
domain GraphML and model-name inventory; it does not run DIAMOND, MCL or an HMM
search. For older protein outputs that stored only edge weights, it reconstructs
identity/coverage metadata from checksum-verified saved DIAMOND hits, retaining
MCL membership. If those files/checksums are missing, it reports the needed input
rather than inventing identity values. A report-only manifest records parameters
and the source analysis manifest checksum.

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
  --output domain_results --threads 32 --resume

# Add domain names to the protein network's node details:
python network_reports.py --clustering MCL2 --domains domain_results
```

The genome wrapper automatically regenerates the protein network after domain
analysis, attaching its domain names and architectures. Its root report links
both types of interactive network. Standalone protein FASTA domain analysis also
gets a native collection-wide interactive domain network.

## Protein amino-acid similarity network

Open `MCL2/protein_network.html` (also linked from `MCL2/index.html`).

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
