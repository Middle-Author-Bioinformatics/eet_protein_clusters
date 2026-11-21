# Whole-proteome genome workflow

`genome_pipeline.py` runs **all proteins** through gene extraction, DIAMOND
all-versus-all similarity, MCL, the existing per-cluster domain analysis, and
localization. It adds FeGenie-compatible heme/hematite motif counts for every
protein, including domain-negative and motif-negative proteins. It does not
select proteins by motif before clustering or label them confirmed EET proteins.

This wrapper belongs in `eet_protein_clusters` because it orchestrates the two
existing clustering/domain scripts and preserves their outputs. It does not
vendor either FeGenie or MHCscan or require their installation. Their specific
motif definitions and output organization informed the independent motif module.

## Install

```bash
git clone https://github.com/Middle-Author-Bioinformatics/eet_protein_clusters.git
cd eet_protein_clusters
git switch codex/genome-proteome-workflow   # PR branch until merged
conda env create -f environment.yml
conda activate eet-clusters
python genome_pipeline.py --help
```

The environment includes DIAMOND, MCL/mcxload, PyHMMER, Pyrodigal, BioPython,
NetworkX and Matplotlib. Existing environments can instead use
`python -m pip install -r requirements-genomes.txt`, with DIAMOND, MCL and
mcxload separately installed and on PATH. Standalone domain-only users can
continue installing just `requirements.txt`.

PSORTb and the supplementary predictors have separate installations. Keep the
wrapper modules together in the repository; copying only `genome_pipeline.py` is
not sufficient. The wrapper
never downloads them or their weights, accepts licenses, or uploads sequences.
Only bacterial positive/negative envelope models are supported here. For unusual
envelopes, select the appropriate PSORTb model externally and import results;
Gram staining alone can be misleading.

| Tool | Contribution and installation | License/access considerations |
|---|---|---|
| PSORTb 3.0 | Primary whole-proteome localization. Install standalone `psort` following https://psort.org/downloads/INSTALL.html; its official download page also recommends containers: https://psort.org/downloads/. Run `psort --version` before the workflow. On Debian/Ubuntu with the package available, `sudo apt-get update` followed by `sudo apt-get install psortb` supplies `psort`. | Standalone GNU GPL; required dependencies/models must also be installed. Official instructions list Perl/BioPerl, IO::String, BLAST and PFTOOLS. This wrapper requires a working installation, not merely a script copied into PATH. |
| SignalP 6.0 | Supplementary Sec/Tat/pilin and lipoprotein SP classes. Obtain a license/package/models from https://services.healthtech.dtu.dk/services/SignalP-6.0/. Follow the authors' installation instructions: https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md. Install in a **separate** environment as needed; the published package instructions support Python up to 3.10. Expose `signalp6` through an executable launcher or config. | Academic download is restricted to eligible use; commercial users need the applicable DTU license. No package or weights redistributed here. |
| Phobius 1.01 | Specialized joint transmembrane topology/SP prediction; `phobius.pl -short`. Obtain a licensed local installation at https://phobius.sbc.su.se/download. | Academic/nonprofit restrictions; other users must contact the provider. No automatic installation. |
| LipoP 1.0 | Optional independent lipoprotein evidence; native short output import or a configured local command. See https://services.healthtech.dtu.dk/services/LipoP-1.0/. | Subject to DTU download/license terms. The provider calls it outdated and recommends current SignalP for lipoprotein signals. It is trained on Gram-negative proteins; Gram-positive performance has been reported, but treat its provenance separately. It is not required. |

SignalP 6 `LIPO`/`TATLIPO` provides lipoprotein evidence without requiring LipoP.
Phobius supplies full-sequence topology, rather than equating an N-terminal
hydrophobic SP with a retained membrane helix. Other predictors can be run outside
the workflow; do not relabel their output as one of these supported native formats.
A future adapter should preserve their source and class definitions.

For example, a SignalP launcher pointing to a separately installed, licensed
Python 3.10 environment can be an executable file named `signalp6` on PATH:

```bash
#!/usr/bin/env bash
exec /absolute/path/to/signalp-env/bin/signalp6 "$@"
```

Configure installation/model details below; there is no requirement to install
SignalP into the clustering Python environment. For the official PSORTb container, use the provider
[command-line wrapper](https://github.com/brinkmanlab/psortb_commandline_docker)
to generate native long results, then use the import route below. That wrapper
has a different CLI (`-i`, `-r`) from standalone `psort`; do not rename it to
`psort` and assume the default command applies. The provider recommends its
pre-built image because an old build dependency URL no longer works.

PSORTb itself is serial per
genome; `--threads` is passed to DIAMOND, MCL, PyHMMER and SignalP. Genomes are
processed in a deterministic order, one at a time for gene/localization steps.

## Domain database

The existing [database instructions](README.md#3-download-a-broad-domain-database)
also apply here. Download a snapshot once and keep the version file:

```bash
mkdir -p databases
wget -O databases/Pfam-A.hmm.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget -O databases/Pfam.version.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
gunzip -c databases/Pfam-A.hmm.gz > databases/Pfam-A.hmm
gunzip -c databases/Pfam.version.gz > databases/Pfam.version.txt
sha256sum databases/Pfam-A.hmm databases/Pfam.version.txt
```

The `current_release` URL changes over time: retain these exact files to reproduce
a run. `--hmm-metadata databases/Pfam.version.txt` records the release text and
checksum. The HMM checksum, accession versions, model lengths and sequence/domain
thresholds are always recorded in the nested domain manifest/model inventory.
User-supplied concatenated HMMs work. Models with GA use their gathering thresholds;
models without GA require an explicit `--domain-evalue` fallback. All predicted
proteins form one common HMM target database, then results are split by cluster.
Repeated domains are retained and domain prevalence counts distinct proteins.
Co-occurrence means domains in the **same protein**, not different members of a
cluster. The existing documented greedy overlap policy still applies.

## Run a folder of assembled genomes

Place one genome's contigs in each `.fna`, `.fa` or `.fasta` file; `.gz` is accepted.
Keep output outside the input folder. For a uniformly Gram-positive collection:

```bash
python genome_pipeline.py \
  --input /absolute/path/to/genomes \
  --input-type assembly \
  --output /absolute/path/to/genome_analysis \
  --hmm databases/Pfam-A.hmm \
  --hmm-metadata databases/Pfam.version.txt \
  --gram positive --threads 32 --inflation 2.0 \
  --resume
```

The default `--psortb run` requires working `psort` or imported native results.
`--signalp auto` and `--phobius auto` run locally if found or configured, otherwise
record unavailable analyses. LipoP auto imports provided results or runs an explicit
config; it is otherwise skipped. A detected/configured predictor that fails causes
an actionable error, rather than silently making the run appear successful.
`--signalp run --phobius run` requires both supplementary analyses explicitly.

To perform clustering, domains and motifs while localization tools are unavailable:

```bash
python genome_pipeline.py -i /path/to/genomes --input-type assembly \
  -o /path/to/analysis --hmm databases/Pfam-A.hmm --gram negative \
  --psortb skip --signalp skip --phobius skip --lipop skip -t 16 --resume
```

Such a run reports localization as **unavailable**, never as negative evidence.
After installing/configuring predictors, rerun the same command with the desired
predictor flags and `--resume`; the expensive matching stages are reused.

For a mixed collection, create a tab-separated metadata file. `genome` is the
filename stem, removing `.gz` and the final extension (e.g. `isolate_1.fna.gz` ->
`isolate_1`). `gram` selects the envelope model; it is never guessed:

```tsv
genome	gram
isolate_1	positive
isolate_2	negative
```

```bash
python genome_pipeline.py -i /path/to/genomes -o /path/to/analysis \
  --hmm databases/Pfam-A.hmm --genome-metadata genomes.tsv -t 32 --resume
```

### Annotation policy and mappings

- Assembly inputs use Pyrodigal **meta mode**, masked ambiguous regions, open
  contig ends, and models restricted to the requested translation table (11 by
  default; 4 optional). No training randomness is involved. Incomplete genes are
  retained and flagged. This practical default supports fragmented genomes; it
  does not claim to optimize gene prediction for every complete isolate genome.
- GenBank `.gb`, `.gbk`, `.gbff` inputs preserve annotated non-pseudogene CDSs.
  Supplied translations take precedence. Missing translations are computed from
  those CDSs using `codon_start` and `transl_table`; annotation source records that
  distinction. Translation exceptions without supplied translations cause an error.
  No de novo calls are mixed into a GenBank proteome.
- `--input-type auto` recognizes assemblies/GenBank, but duplicate basenames and
  accompanying `.faa` files are rejected as ambiguous. Choose an explicit type
  when a folder contains several annotation representations of the same genomes.
- `--input-type proteins` accepts `.faa/.fa/.fasta` proteomes (also gzipped).
  Coordinates are blank and completeness is unknown; no gene prediction occurs.
- All predictor IDs are deterministic alphanumeric `g00001p00000001` values.
  `inputs/protein_mapping.tsv` maps these to genome, original ID/header, source,
  coordinates, strand, full GenBank location, partial status and translation table.
  IDs repeated across genomes remain distinct. A GenBank ID repeated at different
  loci is preserved, but localization imports using that ambiguous original ID fail.
- Coordinates are 1-based inclusive bounding intervals. The full `location` field
  retains joins/fuzzy boundaries, including origin-spanning annotations. The
  wrapper does not perform neighborhood or operon inference or wrap linear contigs.
  Partial proteins can lack targeting signals, so negative results need caution.

### Predictor configuration, versions, and native imports

`--localization-config predictors.json` can specify `version`, `model_files`,
`argv` and `output` for any supported predictor. Model paths are relative to the
JSON file; argv paths should be absolute or on PATH. Use JSON arrays of arguments;
no shell is invoked. Supported placeholders are `{input}`, `{outdir}`, `{gram}`,
`{gram_flag}` (`-p`/`-n`) and `{threads}`. Output is `stdout` or a relative filename
beneath `{outdir}`. Wrapper stdout/stderr are logged together, so a stdout-producing
launcher must avoid progress text that mimics result rows.

```json
{
  "psortb": {
    "version": "3.0.6",
    "argv": ["/opt/psort/bin/psort", "{gram_flag}", "-o", "long", "{input}"],
    "output": "stdout"
  },
  "signalp": {
    "version": "6.0i (your installed build)",
    "argv": ["/opt/signalp-env/bin/signalp6", "--fastafile", "{input}",
             "--organism", "other", "--output_dir", "{outdir}", "--format", "none",
             "--mode", "fast", "--torch_num_threads", "{threads}", "--write_procs", "1"],
    "output": "prediction_results.txt"
  }
}
```

Add `"model_files": ["/absolute/path/to/actual/model-file"]` to hash the actual
installed model files; list all files that determine prediction, not directories.
No model paths or version strings are invented by the wrapper. It attempts
`--version` for PSORTb/SignalP when not supplied. Other unavailable versions and
missing model checksums are explicitly recorded. For a container/launcher,
record its image digest/version and model files in config. Include the digest in
`--predictor-cache-tag` if weights cannot be accessed directly. **Change that tag
when updating models not covered by recorded model checksums.** Default caches
cannot detect hidden model changes behind an unchanged executable/version.
LipoP command options vary by installation; use its local help to create an argv
array producing native short output and set `output` accordingly, or import it.

Native import files in `--localization-import /path/to/results`:

| Filename | Accepted format |
|---|---|
| `<genome>.psortb.tsv` | PSORTb long TSV `SeqID, Final_Localization, Final_Score`; terse `SeqID, Localization, Score` also accepted |
| `<genome>.signalp.tsv` | SignalP 6 `prediction_results.txt`, including its `# ID` tab header, `Prediction` and `CS Position` |
| `<genome>.phobius.txt` | Phobius `-short`: ID, TM count, SP `Y/0`, topology |
| `<genome>.lipop.txt` | LipoP native short lines: `# ID SpII score=... margin=... cleavage=...` (also SpI/TMH/CYT) |

Results may use the safe IDs in `inputs/proteomes/*.faa` or unambiguous original
IDs within each genome. Imports override local execution for that genome/tool,
unless the tool is explicitly `skip`. Use `auto` for an import-only optional tool
if it may be absent for some genomes. Missing rows remain unavailable and are
counted in warnings, not treated as no-SP/zero-TM. Unknown IDs, conflicting or
duplicate rows cause an error. Inputs/outputs are hashed; each import is copied
under `localization/` and marked imported with its provenance. Supply imported
software versions/model files through the same config for complete provenance.

## Interpretation and motifs

The primary PSORTb localization and score are unchanged. A separate
`localization_interpretation` refines only Unknown/unavailable baselines:

| Evidence | Feature interpretation |
|---|---|
| SignalP LIPO/TATLIPO or LipoP SpII | Putative membrane-attached lipoprotein; membrane identity/exposure unresolved |
| Phobius TM count >0 | Putative membrane protein; topology retained |
| SignalP/Phobius/LipoP signal peptide without attachment evidence | Export signal detected, destination unresolved |
| No supporting evidence | PSORTb Unknown remains Unknown; absent PSORTb remains unavailable |

Confident PSORTb calls are retained even when features conflict, and disagreements
are flagged. Secretion, membrane attachment, cell-wall association and extracellular
exposure are separate columns. SP alone does not establish periplasm vs extracellular
location. A Gram-negative lipoprotein is not assigned inner/outer membrane or an
outward-facing surface from its lipobox alone. Supplementary rules are transparent
heuristics, not a newly validated localization classifier. Ambiguous or incomplete
results should be inspected in native files. Neither motif nor proximity to a
surface-associated domain demonstrates EET.

Motifs were verified against pinned source commits (also in every manifest):

- [FeGenie.py](https://github.com/Middle-Author-Bioinformatics/FeGenie/blob/034f554fb12528b43f7515b5c066de610a6ec6e8/FeGenie.py#L2474): `C..CH`, `C...CH`, `C....CH`, `C.{14}CH`, `C.{15}CH`;
  b-type pattern `G.[HR].C[PLAV]G`; hematite-associated pattern
  `[STC][AVILMFYWH][ST]P[ST]`. `fegenie_heme_c_total` sums all five c-type patterns.
- [MHCscan cytoscan.py](https://github.com/Middle-Author-Bioinformatics/MHCscan/blob/2fe42012d46f3c57281b856331dc51e3777b1872/bin/cytoscan.py#L376):
  `mhcscan_hemes` sums CXXCH + CXXXCH, matching its `hemes` column.
- `cxxch_count` is separately reported as the canonical c-type motif count.
  `canonical_mhc_ge2` in summaries counts proteins with >=2 of these motifs;
  it does **not** filter the clustering input.

Counts use **non-overlapping occurrences within each pattern**, matching both
source implementations. Different patterns are searched independently; overlapping
variant matches may contribute to their separate counts and the broader total.
`motif_hits.tsv` gives every counted hit's regex, sequence and 1-based inclusive
coordinates. These short patterns occur by chance; counts are not measured heme
stoichiometry, validated hematite binding, or EET evidence by themselves. The
hematite pattern is the one implemented by FeGenie, not a universal mineral-binding
diagnostic. FeGenie is AGPL-3.0 and MHCscan GPL-3.0; no source code or licensed
predictor/model is copied into this repository.

## Outputs, resumption and scale

```text
analysis/
  protein_annotations.tsv        # Every protein: mappings, cluster, motif counts,
                                 # domain architecture, source localization, inference
  motif_hits.tsv                 # All counted pattern matches and coordinates
  genome_summary.tsv
  cluster_evidence_summary.tsv   # Counts; includes negative/unavailable members
  localization_distribution.tsv  # PSORTb and interpretation fractions per genome/cluster
  index.html                     # Self-contained searchable preview and data/report links
  manifest.json                  # Parameters, commands, versions, hashes, skipped analyses
  inputs/                        # Combined proteins.faa, proteomes/, mapping, extraction cache
  clustering/                    # clusters/*.faa, membership, original clustering outputs
  domains/                       # Per-cluster domain tables, GraphML, SVG, HTML + raw hits
  localization/                  # Native output, parsed source evidence, commands and caches
  work/                          # Disk-backed annotation join
```

The root report previews the first 1,000 proteins (clearly labeled) to keep it
usable for whole proteomes. TSV downloads contain every protein. The existing
per-cluster domain reports and networks remain linked, including repeats and
negative-member denominators. No unexplained confidence score is added.

`--resume` reuses checksum-matching gene extraction, DIAMOND, HMM and localization
results. MCL/cluster export and reporting are regenerated; clustering parameters
can be changed without rerunning compatible alignments. Membership/HMM/sequence
changes invalidate dependent searches. Native imports are reread on every run.
`--force` recomputes stages. An interrupted run leaves logs and an in-progress
manifest; failed predictor outputs are never marked as completed caches. Always
reuse the same input basenames for stable generated IDs; adding an earlier-sorting
genome can renumber later IDs and correctly invalidate downstream caches.

MCL inflation must be **greater than 1**; a coarse value such as 1.2/1.4 can be
explored. Identity/coverage thresholds also change clusters. All-versus-all whole
proteome matching can be expensive. `--block-size` tunes DIAMOND memory;
`--low-memory` streams HMM targets; `--no-plots` avoids SVG rendering for many
clusters. Full Pfam searches and thousands of per-cluster networks can take substantial
time/disk space. No million-protein performance claim is made.

## Reproducible synthetic demonstration and tests

```bash
python examples/make_genome_demo.py --output /tmp/eet-genome-demo
python genome_pipeline.py -i /tmp/eet-genome-demo/genomes --input-type assembly \
  -o /tmp/eet-genome-demo-analysis --hmm /tmp/eet-genome-demo/models.hmm \
  --gram positive --psortb skip --signalp skip --phobius skip --lipop skip \
  --threads 2 --no-plots --resume
python -m pytest -q
```

The demo generates artificial genomes and HMMs, runs real gene prediction,
clustering and domain searches, and explicitly leaves localization unassessed.
Tests use real Pyrodigal/PyHMMER and, when available, real DIAMOND/MCL. Predictor
adapters use labeled test doubles/native-format fixtures, not real PSORTb,
SignalP, Phobius or LipoP model inference. See [VALIDATION.md](VALIDATION.md).
