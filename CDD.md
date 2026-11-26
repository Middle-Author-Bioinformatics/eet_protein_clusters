# CDD domain profiles and networks

Use `download_cdd.py` to install a self-contained NCBI Conserved Domain Database
(CDD) snapshot, then `cluster_cdd.py` to profile one protein FASTA or an existing
folder of cluster FASTAs. Existing MCL memberships are preserved: this does not
recluster the proteins. Keep CDD results in a different folder from Pfam results.

## Additional software

In the existing environment, install **NCBI BLAST+** and **rpsbproc**:

```bash
conda install -n eet-clusters -c conda-forge -c bioconda blast rpsbproc
conda activate eet-clusters
rpsblast -version
blast_formatter -version
rpsbproc -version
```

BLAST+ provides `rpsblast` and `blast_formatter`. The downloader uses Python's
standard library; it does not require wget or curl. Analysis reuses the existing
BioPython, NetworkX, and Matplotlib dependencies. CDD uses PSSMs and RPS-BLAST,
not a new HMM database or DIAMOND. On Debian, `rpsblast+` is recognized as an
alternative executable name. These dependencies remain optional for Pfam runs.
The Conda package manager installs platform-specific runtime libraries; avoid
copying an isolated NCBI executable without its required libraries.

NCBI software is public domain (the rpsbproc Bioconda recipe identifies that
license). CDD incorporates models from multiple contributing databases: retain
their provenance and citations. No models, predictor weights, or NCBI executables
are redistributed by this repository.

## Download once

Use the `/home/ark` data filesystem, with ample space for the archive plus its
expanded database; budgeting at least 20 GB free is prudent. This is a storage
allowance, not a measured peak-use requirement. The default is the complete
`Cdd_LE.tar.gz` search database for little-endian systems (Linux x86-64 and current
Intel/Apple Silicon Macs). Supported executable packages vary by platform.

```bash
python /home/ark/MAB/bin/eet_protein_clusters/download_cdd.py \
  --output /home/ark/databases/CDD --resume
```

For only NCBI-curated models, use a **separate** snapshot directory:

```bash
python /home/ark/MAB/bin/eet_protein_clusters/download_cdd.py \
  --output /home/ark/databases/CDD_NCBI --subset ncbi --resume
```

Setup downloads the preformatted search database, `cddid.tbl`, specific bit-score
thresholds, annotation data, family/superfamily links, `cdd.info`, and README files.
It records source URLs, download times, HTTP metadata, and local SHA256 digests.
The database prefix is detected from the extracted files, not guessed from the
archive's name. Safe extraction rejects links, devices, and traversal paths.

`cdd_manifest.json` points to a completed snapshot beneath `snapshots/`. Successful
`--resume` verifies and reuses its search data **without checking for updates**.
An interrupted setup restarts its downloads; HTTP range resumption is not provided.
`--force` creates a new snapshot and updates the manifest only after success, while
preserving prior complete snapshots. Failed updates leave the prior manifest intact.
The archive is removed after extraction unless `--keep-archive` is given.

Search database and annotation-file checksum failures remain fatal. The two
snapshot-root README copies are documentation: if missing or edited, validation
warns after checking the scientific files, and records expected/actual documentation
hashes and status in the analysis manifest. It does not rewrite the original
snapshot manifest or download a new database. A README warning is not evidence
that the database itself is absent or corrupt. Other missing/changed files are
reported together so a missing snapshot can be distinguished from one missing file.

The release file is checked before and after setup to detect a changing release.
Some annotation files update independently of that release: timestamps and hashes
are retained, but this is not an atomic upstream release or signed authenticity
verification. Keep the whole snapshot and its manifest together for reproducibility.

## Search the existing MCL clusters

```bash
python /home/ark/MAB/bin/eet_protein_clusters/cluster_cdd.py \
  --input "$HOME/Marquette/eet/MCL2-test/clusters" \
  --cdd /home/ark/databases/CDD \
  --output "$HOME/Marquette/eet/MCL2-test/eet_cluster_cdd" \
  --threads 32 --resume
```

Every protein is searched together, then results are summarized by cluster. A single
unclustered `.faa` or `.faa.gz` is also supported, with native collection outputs:

```bash
python /home/ark/MAB/bin/eet_protein_clusters/cluster_cdd.py \
  --input "$HOME/Marquette/eet/GramPosDedupeGenes.faa" \
  --cdd /home/ark/databases/CDD \
  --output "$HOME/Marquette/eet/GramPosDedupeGenes_cdd" \
  --threads 32 --resume
```

This does not represent the single FASTA as an artificial MCL cluster.
The same globally unique first-token protein IDs required for Pfam are required here.

## Search and evidence definitions

The pipeline runs:

1. RPS-BLAST with composition-based correction, SEG off, E-value <= **0.01** by
   default (`--evalue` changes it), and `-outfmt 11` to retain the BLAST archive.
   `-max_target_seqs 1000000` avoids the usual top-hit-only limit at current CDD size;
   repeated HSPs are not restricted to one per model.
2. `blast_formatter` on that archive to retain all reported alignment coordinates,
   PSSM IDs, scores, E-values, and model lengths.
3. `rpsbproc -m rep -t doms` with the snapshot's matching metadata. NCBI's
   representative view resolves redundant model alternatives. Native Specific,
   Non-specific (also the historical Non-Specific spelling), Superfamily,
   Multidom, incompleteness and parent-superfamily
   fields are retained. Specific assignments are made by rpsbproc, with the
   corresponding model-specific bit-score cutoff recorded where supplied.

CDD E-values and bit-score criteria are **not Pfam gathering thresholds**. CDD's
imported Pfam PSSMs may represent a different model release from your HMM database.
An overlapping Pfam hit is not an independent confirmation just because CDD also
reports it. Names, accessions, source databases, release metadata and PSSM IDs
are kept, rather than relabeling CDD accessions as the original Pfam HMM results.

Raw RPS alignments and native rpsbproc output remain under `work/`. Summaries
use representative hits, excluding whole-protein `Multidom` models. An overlapping
parent superfamily is suppressed when its child domain is already selected.
Overlapping occurrences of the **same** model are deduplicated; disjoint repeats
are retained. The overlap/shorter-hit fraction defaults to 0.5 and is configurable
with `--overlap-fraction`. No additional generic overlap suppression between
unrelated models is applied after rpsbproc's representative selection.

`domain_hits.tsv` records selected/rejected representative hits and reasons.
It uses CDD-specific fields, not HMM coordinate or independent-domain E-value
labels. Model coverage is calculated only when a native hit exactly matches a
raw PSSM alignment by protein, PSSM ID and query coordinates. Derived/unmatched
hits have unavailable model coverage; query coverage still uses their query span.
`--min-model-coverage` defaults to zero; a positive cutoff excludes unavailable
model coverage with an explicit reason, rather than inventing a coverage value.
Superfamily entries can legitimately have model length zero in CDD metadata.

Repeated copies contribute to copy counts, but count once per protein for
prevalence and co-occurrence. Proteins without hits remain in all denominators.
Two domains in different proteins in one cluster do not create a co-occurrence
edge. A network describes protein architectures, not physical interactions,
localization or confirmed EET activity. Parent superfamily annotations and domain
models are not equivalent evidence; their hit types remain visible in the tables.
Conserved sites and SPARCLE protein naming are not run by this domain-only pipeline.

## Outputs and interactive reports

Cluster-folder input produces `cluster_summary.tsv`, domain prevalence, raw/selected
hit tables, named protein architectures, `model_inventory.tsv`, per-cluster reports,
and a pooled `all_domains_network.html`. Single-FASTA input produces
`collection_summary.tsv`, collection prevalence and `domain_network.html`, without
artificial cluster columns. Both produce `index.html`, GraphML, tables, and a
manifest with commands, versions, parameters, hashes, and skipped analyses.
The inventory contains the full metadata table, including entries not in a selected
search subset; raw alignments/native hits identify the actual models matched.

Networks reuse the Pfam workflow's same-protein counts, Jaccard, lift, community
grouping, carrier lists, and preview limits. CDD model names are displayed, and
clicking a node links to its **CDD** record by PSSM ID. Use the existing report
refresher to add NR annotations and MCL membership without new searches:

```bash
python /home/ark/MAB/bin/eet_protein_clusters/network_reports.py \
  --domains "$HOME/Marquette/eet/MCL2-test/eet_cluster_cdd" \
  --membership "$HOME/Marquette/eet/MCL2-test/membership.tsv" \
  --protein-fastas "$HOME/Marquette/eet/MCL2-test/clusters_nr_annotated"
```

For a standalone FASTA, use its CDD output folder with the same membership argument
to attach a later clustering. Membership must match by original ID; no cluster is
guessed from a protein's annotation. Cluster highlighting indicates domain presence,
not a recalculation of pooled co-occurrence edges. `--no-plots` skips SVGs;
`--no-interactive` skips the interactive HTML. The search still completes.

Analysis `--resume` reuses checksum-verified completed RPS/rpsbproc output when
input, database snapshot, tools and search parameters match. It regenerates
summaries from those results. Report-only rebuilding needs the source architecture
table and GraphML, not a full database or the search executables. Partial searches
are never marked complete. An existing Pfam output folder is rejected.

Advanced import: `--rpsbproc-output FILE --raw-hits FILE` skips external search
tools but still needs the matching snapshot metadata. The raw table format is
`blast_formatter -outfmt '6 qseqid sseqid qstart qend sstart send evalue bitscore qlen slen'`.
Both files must come from the exact same normalized `work/proteins.faa` and database,
with generated `pNNNNNNNNNN` query IDs in its original order. The native file must
be representative `-m rep -t doms` output without `-f`. IDs, lengths, coordinates,
PSSM/accession matches and complete query assessment are validated; a same-length
sequence substitution in an externally supplied search cannot be detected from
these tables. Preserve the original input checksum when using this advanced route.

## Authoritative references

- [NCBI CDD help](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml)
- [NCBI database and download README](https://ftp.ncbi.nih.gov/pub/mmdb/cdd/README)
- [NCBI rpsbproc documentation](https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/README)
- [NCBI preformatted databases](https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/)
- [BLAST+ manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [Bioconda rpsbproc installation](https://bioconda.github.io/recipes/rpsbproc/README.html)
