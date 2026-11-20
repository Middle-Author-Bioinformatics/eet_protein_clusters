# Protein clustering and cluster-specific domain profiles

Two scripts for an already selected collection of candidate EET-related proteins:

1. **cluster_proteins.py**: protein FASTA -> DIAMOND all-versus-all similarity ->
   MCL clustering -> one FASTA per cluster, including singletons.
2. **cluster_domains.py**: cluster FASTAs -> Pfam/custom HMM searches -> per-cluster
   domain prevalence, repeated-domain statistics and within-protein co-occurrence networks.

Keep `common.py` alongside the two scripts. No genome coordinates, earlier EET
wrapper, FeGenie checkout, or web submission is required. No new EET activity
classification is made. These scripts also work for other preselected protein sets.

## 1. Install on TheBelly / Linux

After unzipping this package, enter its directory. The simplest installation is a
separate Conda environment; no sudo is needed:

```bash
conda create -n eet-clusters --override-channels \
  -c conda-forge -c bioconda \
  python=3.11 mcl diamond biopython pyhmmer networkx matplotlib
conda activate eet-clusters

mcl --version
mcxload --version
diamond version
python -c "import Bio, pyhmmer, networkx, matplotlib; print('Python dependencies OK')"
```

Alternatively:

```bash
conda env create -f environment.yml
conda activate eet-clusters
```

MCL's own README recommends package-manager installation and explains that its
GitHub development tree requires additional build machinery and cimfomfa. You do
not need to compile that tree for these scripts. See
https://github.com/micans/mcl#installation-mcl-versions-and-linux-packages.

For an existing environment that already has `mcl`, `mcxload` and `diamond`:

```bash
python -m pip install -r requirements.txt
```

PyHMMER includes the HMMER search engine, so this workflow does not require a
separate `hmmsearch` executable or `hmmpress`. No restricted localization software
is required. MCL and DIAMOND remain independently installed external programs;
no third-party executable or model is distributed in this package.

## 2. Cluster the proteins

```bash
python cluster_proteins.py \
  --input /path/to/eet_candidates.faa \
  --output /path/to/eet_mcl \
  --threads 32 \
  --inflation 2.0 \
  --min-identity 30 \
  --min-coverage 0.5 \
  --resume
```

Defaults:

| Parameter | Default | Meaning |
|---|---:|---|
| DIAMOND sensitivity | very-sensitive | Similarity search suitable for exploring divergent proteins |
| E-value | 1e-5 | Maximum reported/retained alignment E-value |
| Identity | 30% | Minimum identity in the qualifying HSP |
| Coverage | 50% on **both** proteins | Based on aligned sequence coordinate spans, not alignment length including gaps |
| MCL inflation | 2.0 | Starting granularity; compare 1.4, 2 and 4 when assessing stability |
| Edge weight | capped -log10(E) | Greater similarity -> greater weight; zero E-values capped at 200 |
| Repeated/directional edges | maximum weight | One symmetric edge per protein pair |
| Proteins without qualifying edges | retained | Each becomes a singleton FASTA |

DIAMOND is explicitly run with **unlimited reported targets** (`--max-target-seqs 0`)
and one best HSP per pair, avoiding an accidental top-25-neighbor network. MCL
receives an undirected sparse matrix created using `mcxload --stream-mirror`.
The protein-to-cluster output is verified to be a complete partition. A missing
non-isolated node from MCL output is an error, not silently restored as a singleton.

Higher inflation generally produces finer clusters; lower inflation coarser
clusters. Graph construction and similarity thresholds also affect the outcome.
These are sequence-similarity groups, not automatically orthologs, equivalent
functions or validated EET families. Multidomain proteins and repeats can connect
families through shared regions. The bidirectional coverage threshold reduces
short-domain bridging but does not eliminate it. Keep full protein sequences for
this first analysis; do not strip the domains you intend to study later.

To explore more divergent families, lower `--min-identity` (e.g. 25) and/or adjust
coverage deliberately; inspect the resulting architectures. Fragmented proteins
may fail the default coverage threshold. The same-family/function interpretation
requires separate validation.

### Input requirements and scale

- One protein per FASTA record; `.faa.gz` is also supported.
- IDs are the first whitespace-delimited token and must be **globally unique**.
  Namespace repeated locus IDs as `genome|locus` before concatenating genome outputs.
  Duplicate IDs fail visibly rather than merging records.
- Sequences are uppercased and terminal stop symbols removed. Internal stops,
  empty sequences and invalid protein characters are rejected.
- Original IDs and full headers are restored in cluster FASTAs. Identical sequences
  with different IDs remain separate records and contribute separately to prevalence.
  Deduplicate first if that is your intended sampling unit; duplicated strains or
  overrepresented taxa can affect every downstream frequency.
- The graph edge-deduplication table is on disk in SQLite; sequences are accessed
  by file offsets when cluster FASTAs are written. MCL still needs RAM for its sparse
  graph, and all-versus-all searches can create very large temporary/alignment files.
- `--block-size 1` reduces DIAMOND's memory use. It does not put a hard total-memory
  limit on DIAMOND, MCL, or the operating system. Threads are configurable.

### Clustering outputs

```text
eet_mcl/
  clusters/cluster_000001.faa
  clusters/cluster_000002.faa
  ...
  membership.tsv          # original protein IDs and cluster labels
  cluster_summary.tsv     # sizes and length distributions
  id_mapping.tsv          # safe internal IDs -> original IDs/full headers
  manifest.json           # versions, input/script hashes, parameters, commands
  run.log
  work/                    # DIAMOND hits/database, graph, MCL result, raw logs, caches
```

Cluster names are sorted by decreasing size, then first input-record index.
All singletons get their own FASTA. For a very large number of singletons this
can create many small files; the membership table is the compact authoritative map.

`--resume` reuses an identical, checksum-verified DIAMOND result. Changing inflation,
identity/coverage filters or edge weighting rebuilds the graph and cluster FASTAs;
changing the input, search sensitivity, search E-value, block size, DIAMOND version
or thread setting invalidates the search cache. Save independent runs in different
output folders if comparing clusterings. Existing generated `cluster_*.faa` files
are replaced during a recomputation so stale clusters cannot enter script 2.
`--force` explicitly recomputes an existing output folder.

An optional `--hits existing_hits.tsv` bypasses DIAMOND. This must contain the
following **custom 12 columns**, with ORIGINAL protein IDs:

```text
qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen
```

This is **not** the default BLAST 12-column order. For example, generate it with
`diamond blastp --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen`.
Each HSP is filtered independently; multiple HSPs are not pooled into coverage.
The script's own `work/diamond_hits.tsv` uses safe INTERNAL IDs and is reused via
`--resume`, not passed as an original-ID `--hits` file. MCL/mcxload are still needed.

## 3. Download a broad domain database

For your exploratory EET collection, start with **all Pfam families**, rather
than restricting the search to the six surface-associated families from the earlier
wrapper. That can reveal cytochrome, binding, enzyme, membrane and uncharacterized
domains across different clusters.

```bash
mkdir -p databases
wget -O databases/Pfam-A.hmm.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget -O databases/Pfam.version.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
gunzip -c databases/Pfam-A.hmm.gz > databases/Pfam-A.hmm
gunzip -c databases/Pfam.version.gz > databases/Pfam.version.txt
cat databases/Pfam.version.txt
```

The full library is a large download and the uncompressed file is much larger.
Keep the exact HMM file, release text and resulting manifest to reproduce results;
`current_release` may change in future. Model accessions/versions and actual
thresholds are recorded in `model_inventory.tsv`; the HMM SHA256 is in the manifest.
The release-text file is retained alongside your database, not inferred from a filename.

A custom concatenated HMM database can be supplied instead. Gathering thresholds
are used whenever available, with both sequence and domain GA criteria. A model
without GA fails early unless you explicitly supply a fallback such as
`--domain-evalue 1e-5`, applied to sequence and independent-domain E-values.
GA models still use GA when this fallback is supplied. Duplicate accessions or
multiple versions of the same family in one database are rejected to avoid double
counting the same family as separate domains.

## 4. Profile domains and generate networks

```bash
python cluster_domains.py \
  --clusters /path/to/eet_mcl/clusters \
  --hmm databases/Pfam-A.hmm \
  --output /path/to/eet_cluster_domains \
  --threads 32 \
  --resume
```

Every protein is searched. The script searches the **combined collection once**
and aggregates results by cluster, rather than launching a complete Pfam scan for
every small cluster. This shares the same target-database size for all E-values.
It reads HMMs in batches; by default the digital protein sequences are held in RAM.
`--low-memory` streams protein targets from disk, trading memory for repeated I/O.
`--model-batch-size` controls the number of HMMs in memory at a time.

`--resume` reuses a checksum-matching completed HMM search; summaries/plots are
regenerated. The key includes protein content, cluster membership/file content,
model file, fallback threshold, PyHMMER version and search-script hashes. Changing
cluster names/membership invalidates the mapping cache. Interrupted HMM searches
are rerun; only **completed** scans are cached. Logs and an in-progress manifest
remain after a failure; no successful manifest is written for an incomplete run.

### What "domain propensity" means here

For each cluster and domain:

- **Proteins with domain / all proteins in cluster**: prevalence (fraction and percent).
- Total domain copies: repeated domains count individually here.
- Mean copies per cluster member and mean copies among domain-positive members.
- Prevalence outside the cluster and a descriptive ratio against that background.

For example, if 12 of 20 proteins carry SLH, prevalence is **60%**, even if those
12 proteins together contain 36 SLH repeats. The copy count is 36 and copies per
positive protein is 3. Proteins with no detected domain remain in the denominator.
A missing hit is "not detected under these models/thresholds", not proof that the
protein lacks the domain or has no EET function.

The outside-cluster ratio is descriptive, **not** a statistical enrichment test:
sequence clustering, phylogeny and duplicated/overrepresented samples violate a
naive independence assumption. `Inf` indicates no observed outside-cluster copies;
`NA` means there are no proteins outside this cluster.

### Network definition

- Node = distinct domain accession/version. Size = prevalence in that cluster.
- Edge = two different domains found in the **same individual protein** within that cluster.
- Edge table: number of proteins with both, cluster fraction, Jaccard similarity,
  lift, and both conditional probabilities. Figure line width represents Jaccard.
- Repeated copies of the same domain count **once per protein** for an edge.
  A protein carrying A-A-B adds one A-B co-occurrence, not two.
- A domain found only in one cluster member and another found only in a different
  member does not create an edge. A cluster with one detected domain has a node
  and no edges; domain-negative clusters have an empty network and a plot annotation.

These networks describe protein architectures. They do not infer protein-protein
interactions, operons, genomic proximity, surface exposure or EET activity.

By default, strongly overlapping matches from different models are resolved using
a **conservative greedy heuristic**, choosing lower independent E-value, then higher
bit score, and suppressing overlaps covering at least 50% of the shorter hit.
Repeated non-overlapping copies of the same domain are retained. All threshold-
passing raw hits and selection/exclusion reasons remain in `domain_hits.tsv`.
This is **not PfamScan's curated clan/nesting algorithm**; E-values are only a
ranking heuristic and legitimate nested domains can be suppressed. Use
`--overlap-policy all` to retain all matches in summaries for architecture review,
and inspect raw coordinates before interpreting such edges.

Useful options:

```bash
# Reduce weak/single-protein edges in graphs; complete pair table is still retained
python cluster_domains.py -i eet_mcl/clusters --hmm databases/Pfam-A.hmm \
  -o eet_cluster_domains --resume --min-pair-count 3

# Explore all overlapping models and avoid rendering many thousands of SVGs
python cluster_domains.py -i eet_mcl/clusters --hmm databases/Pfam-A.hmm \
  -o domain_profiles_all --overlap-policy all --no-plots --threads 32
```

### Domain output layout

```text
eet_cluster_domains/
  index.html                      # local searchable cluster index
  cluster_summary.tsv             # sizes, fractions annotated, domain/edge counts
  domain_prevalence.tsv           # long-format per-cluster/domain summaries
  protein_architectures.tsv       # every protein, ordered repeated domain hits
  domain_hits.tsv                 # raw threshold-passing hits + overlap decisions
  model_inventory.tsv            # versions, GA/fallback settings
  id_mapping.tsv
  manifest.json
  run.log
  clusters/<cluster>/
    domain_prevalence.tsv
    cooccurrence_edges.tsv
    domain_network.graphml        # open in Cytoscape/Gephi
    domain_network.svg            # scalable figure; omitted with --no-plots
  work/                           # normalized sequences, raw domtblout, SQL/cache
```

Open `index.html` locally. Keep the folder together so its links work. Network SVGs
are standalone vector images. To keep dense figures legible, plots show at most
60 of the most prevalent nodes and 300 strongest qualifying edges by default;
the figure title states shown/total counts. GraphML retains all domain nodes and
all edges passing `--min-pair-count`; the edge TSV retains **all** observed pairs.
TSVs are the complete data, including families absent from the displayed subgraph.
Domain-negative rows are present in `cluster_summary.tsv` and protein architectures;
`domain_prevalence.tsv` is sparse and lists only detected domains.

## Validation and references

See `VALIDATION.md` for tests and exactly which external programs were run.
The scripts do not include an experimental EET positive/negative benchmark.

Primary documentation consulted:

- MCL: https://github.com/micans/mcl
- MCL protein-network protocols: https://micans.org/mcl/man/clmprotocols.html
- DIAMOND output, sensitivity and target-limit options:
  https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
- PyHMMER search/GA API:
  https://pyhmmer.readthedocs.io/en/stable/api/hmmer/profile.html
- Official Pfam data: https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/

For publication, cite MCL (van Dongen 2008) and the protein-family application
(Enright, van Dongen & Ouzounis 2002), plus DIAMOND, HMMER/PyHMMER and Pfam for the
analyses actually used. Keep the manifests and chosen thresholds with the results.
