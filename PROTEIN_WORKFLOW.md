# Protein FASTA to annotated families, domain matrices, HMMs, and trees

`protein_pipeline.py` combines the existing scripts into one resumable command.
It accepts a protein `.faa` (also `.faa.gz`), retains every protein including
singletons and proteins without NR/domain matches, and preserves the original
protein IDs. IDs must be globally unique. Sequences are normalized to uppercase;
terminal stop characters are removed, as in `cluster_proteins.py`.

## Install

For your existing environment, the **additional programs are MAFFT and FastTree**:

```bash
conda activate eet-clusters
conda install -c conda-forge -c bioconda mafft fasttree
python -m pip install -r /home/ark/MAB/bin/eet_protein_clusters/requirements.txt
```

For a new environment, from this repository:

```bash
conda env create -f environment.yml
conda activate eet-clusters
```

The workflow uses DIAMOND, MCL/mcxload, MAFFT, FastTree, and the Python packages
already used for domain profiling. PyHMMER builds the family profiles, so a
separate `hmmbuild` installation is unnecessary. No localization software or LLM
API is required. The genome/localization scripts remain available separately.

CDD additionally requires:

```bash
conda install -c conda-forge -c bioconda blast rpsbproc
```

MAFFT's core is BSD licensed; FastTree's project supplies GPL-3.0 licensing;
PyHMMER is MIT licensed. These dependencies are installed separately, not
redistributed here. See their authoritative installation/API documentation:

- MAFFT: https://mafft.cbrc.jp/alignment/software/
- MAFFT threading/repeatability: https://mafft.ddbj.nig.ac.jp/alignment/software/multithreading.html
- FastTree: https://morgannprice.github.io/fasttree/ and https://github.com/morgannprice/fasttree
- PyHMMER: https://pyhmmer.readthedocs.io/en/latest/api/plan7/pli.html

NR and Pfam/CDD remain user-supplied reference databases. Keep their release
metadata; use [NR_ANNOTATION.md](NR_ANNOTATION.md), the existing Pfam setup in
[README.md](README.md), or [CDD.md](CDD.md) for database preparation. For a fresh
CDD installation:

```bash
python download_cdd.py --output /home/ark/databases/CDD
```

## Run on your protein collection

```bash
conda activate eet-clusters

python /home/ark/MAB/bin/eet_protein_clusters/protein_pipeline.py \
  --input "$HOME/Marquette/eet/GramPosDedupeGenes.faa" \
  --output "$HOME/Marquette/eet/eet_families_pfam" \
  --nr-db /home/ark/databases/nr.dmnd \
  --hmm /home/ark/databases/HMMs/Pfam-A.hmm \
  --threads 32 \
  --tmpdir /home/ark/eet_scratch \
  --resume
```

Pfam is the default. MCL inflation is **2.0**, the clustering identity cutoff is
**35%**, and coverage is **0.5 on both proteins**, matching your earlier runs.
The standalone `cluster_proteins.py` still defaults to 30% identity; this wrapper
explicitly uses 35%. Override with `--inflation`, `--min-identity`, and
`--min-coverage`. Clustering uses a DIAMOND very-sensitive all-versus-all search
(E-value 1e-5), capped −log10(E-value) edge weights and MCL; isolated proteins are
retained. DIAMOND block size is 0.4 for both searches and is configurable.

The equivalent CDD workflow uses a separate output directory:

```bash
python /home/ark/MAB/bin/eet_protein_clusters/protein_pipeline.py \
  --input "$HOME/Marquette/eet/GramPosDedupeGenes.faa" \
  --output "$HOME/Marquette/eet/eet_families_cdd" \
  --nr-db /home/ark/databases/nr.dmnd \
  --domain-db cdd \
  --cdd /home/ark/databases/CDD \
  --threads 32 \
  --tmpdir /home/ark/eet_scratch \
  --resume
```

Defaults for the NR, Pfam, and CDD paths already match the paths above, but
explicit arguments make a run easier to reproduce. Use `--help` for all options.

## What the workflow does

1. **Clusters sequences** with the existing `cluster_proteins.py`.
2. **Annotates cluster FASTAs** using `annotate_clusters_nr.py`. A single DIAMOND
   search against NR covers all proteins. The best qualifying **reported** hit
   supplies the description; ties and alignment statistics remain in the detailed
   annotation table. Descriptions are not restricted to informative names or
   assumed to be gene symbols. The default report limit is 25 targets
   (`--nr-max-target-seqs 0` requests all). NR thresholds have their own
   `--nr-evalue`, `--nr-min-identity`, `--nr-query-cover`, and `--nr-subject-cover`
   arguments; coverage percentages here are independent of clustering coverage.
3. **Profiles domains** using `cluster_domains.py` (Pfam gathering thresholds
   preferred) or `cluster_cdd.py` (RPS-BLAST plus rpsbproc representative hits).
   Custom Pfam-format HMMs without GA require an explicit `--domain-evalue`.
   CDD uses `--cdd-evalue` (0.01). `--min-model-coverage` defaults to zero.
   The existing overlap/redundancy selection logic is retained, and all raw hits,
   repeated-domain coordinates, scores and search provenance remain available.
4. **Builds the requested tables**, joins exact protein IDs, refreshes the
   annotated interactive networks, and performs the family sequence analyses.

### Per-cluster table

`families/cluster_000001/proteins.tsv` has:

| Column | Contents |
|---|---|
| 1: `protein_id` | Original protein ID |
| 2: `annotation` | NR top-hit product/description; an explicit no-hit/unavailable label when needed |
| 3 onward | One column per selected domain found anywhere in this cluster, headed `Domain_name [accession.version]`; values 1 or 0 |

Repeated copies count **once** in this matrix; copy counts and coordinates remain
in `domains/domain_hits.tsv` and the existing domain summaries. A zero means
not detected under this run's thresholds, rather than proof of biological
absence. A cluster without domain matches has just the two metadata columns.
`domain_columns.tsv` is a machine-readable name/accession key. CDD accessions and
names remain native CDD identifiers, not relabeled Pfam identifiers.

### Alignment and HMM curation

MAFFT runs `--auto --amino --anysymbol --thread N --threadit 0`. Disabling
iterative-refinement threading follows MAFFT's repeatability guidance. Safe IDs
are used internally, and aligned residues and the complete ID set are validated
before restoring original IDs. Full alignments retain all members, including
identical sequences.

HMM training uses a **separate**, auditable alignment:

1. Select one representative of each exact sequence to reduce duplicate bias.
2. Keep columns with at least 50% standard-amino-acid occupancy among current
   unique representatives (`--min-column-occupancy`). Ambiguous residues and gaps
   do not count as occupancy evidence.
3. Exclude representatives covering less than 50% of the retained core
   (`--min-sequence-coverage`). Repeat column/sequence filtering until stable;
   eliminated columns and rows are not reintroduced.
4. Build a profile only with at least **3 unique retained sequences** and
   **20 curated columns** (`--min-hmm-sequences`, `--min-hmm-columns`).

This retains conserved positions useful to HMMs, rather than selecting only
phylogenetically informative sites. It is reproducible automated curation,
not a claim that every family has an optimal biological alignment. Inspect the
full/curated alignments and `curation_sequences.tsv` / `curation_columns.tsv`.
The audit records original 1-based alignment columns, occupancy, duplicate
representatives, fragment exclusions, and exclusion rounds.

PyHMMER builds amino-acid profiles with position-based weighting,
entropy-based effective sequence counts, fast architecture (`symfrac=0.5`), and
seed 17. HMMER may choose fewer match states than curated columns after weighting;
both lengths are recorded. The combined `family_profiles.hmm` contains only
successfully built models. Profiles have **no calibrated GA cutoff** and are
family models, not validated EET detectors. Diverse domain architectures are
flagged for review; a global family HMM may need further subfamily splitting.

### Protein trees

FastTree runs `-lg -gamma -seed 17` on the **full MAFFT alignment**, retaining all
protein IDs, including sequences excluded from HMM training. It produces an
unrooted approximate maximum-likelihood tree: topology/branch optimization uses
the CAT approximation and `-gamma` subsequently rescales lengths under Gamma20.
Local SH-like support values are not bootstrap support. No outgroup is assumed.

Every family gets a Newick file unless `--skip-tree` is requested. Singletons,
two-member families, families with fewer than three distinct sequences, and
alignments without variable standard-AA columns get an explicitly labeled
**unresolved membership representation**, with no invented lengths or supports.
These are not inferred evolutionary trees. `tree.nwk` quotes original protein
IDs; `tree.safe.nwk` and `id_mapping.tsv` support software with restrictive ID
parsers. Exact tip membership is checked after tree generation.

`--skip-hmm` and `--skip-tree` are available; MAFFT alignment and the curation
audit still run. FastTree is optional with `--skip-tree`. Small families retain
tables and alignments even when no HMM can be built.

## Outputs and resuming

| Path | Contents |
|---|---|
| `index.html` | Searchable cluster overview and links to all family products |
| `cluster_summary.tsv` | Sizes, domains, curation counts, HMM/tree status, review flags |
| `clustering/` | Existing MCL results, memberships, original cluster FASTAs, protein identity networks |
| `nr_annotated/` | Same-named annotated FASTAs and detailed `annotations.tsv` |
| `domains/` | Complete Pfam or CDD domain results and co-occurrence networks |
| `families/cluster_*/` | Protein/domain table, alignments, curation audit, HMM if eligible, Newick, family report |
| `family_profiles.hmm` | Combined built family HMMs; absent when none are eligible |
| `logs/` | Logs for each child script, including search errors |
| `manifest.json` | Parameters, input hash, software versions, commands, repository commit, database/child provenance and skipped analyses |

Family `sequence_manifest.json` files additionally record actual MAFFT/FastTree
commands, PyHMMER builder settings, sequence and output checksums, and reasons
for skipped inference. Per-family MAFFT/FastTree diagnostics remain with the
alignment/tree. Failed runs write `manifest.failed.json` without a successful
root manifest. Dependency/model checks occur before expensive searches.

Re-run the same command with `--resume` after an interruption. Existing child
scripts verify their own search caches; family sequence products are reused
only when sequences, settings, tool versions and output hashes match. Updating
NR descriptions or domain thresholds refreshes tables without needlessly
rebuilding identical family alignments. Changing HMM curation settings rebuilds
family sequence products. `--force` recomputes generated results. Use separate
root output directories for Pfam versus CDD; switching databases within a run is
rejected to avoid mixing their results. Do not use an old standalone MCL/domain
output folder as the root pipeline output.

NR uses the existing default stat fingerprint (path/size/time), with DIAMOND
database metadata in its child manifest; `--hash-nr` additionally hashes the full
database at potentially substantial I/O cost. Pfam hashes and model thresholds,
or the checked CDD snapshot metadata, remain in the respective domain manifests.
Protein similarity networks show alignment identity, not genome-wide AAI.
Domain co-occurrence means domains detected in the same protein, not interaction.
No sequence similarity, domain, NR description, profile, or tree here proves EET
activity or surface exposure.

## Validation

`tests/test_protein_pipeline.py` covers exact-ID joins; column order and domain
names; repeated-domain binary coding; no-hit/zero-domain proteins; fragment,
duplicate, ambiguous-residue and low-occupancy curation; alignment residue
validation; singleton/small/degenerate families; punctuation-safe tree ID
round trips; corruption and settings invalidating caches; safe input paths;
failed-stage manifests; and CDD command dispatch.

The small end-to-end integration test executes real DIAMOND, MCL, PyHMMER,
MAFFT and FastTree using synthetic proteins, a tiny local NR-like database, and
test domain HMMs. It also searches a generated HMM to verify it is usable and
checks complete tree tip membership. CDD parsing/selection continues to be tested
using the repository's native rpsbproc fixtures; the new wrapper's CDD dispatch is
tested without downloading a full database. This does not benchmark biological
accuracy against your full NR/Pfam/CDD databases or validate EET function.

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```
