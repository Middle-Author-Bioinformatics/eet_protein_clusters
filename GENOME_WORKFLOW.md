# Whole-proteome genome workflow

`genome_pipeline.py` runs all proteins through gene extraction, DIAMOND/MCL
clustering, the existing domain analysis, **SignalP 4 and Phobius**, and per-protein
FeGenie/MHCscan heme and hematite motif counting. PSORTb and LipoP are **skipped by
default**. No protein is excluded because it lacks a motif or domain.

## Install and use the existing TheBelly predictors

```bash
git clone https://github.com/Middle-Author-Bioinformatics/eet_protein_clusters.git
cd eet_protein_clusters
git switch codex/diamond-scratch-directory  # until this PR is merged
conda env create -f environment.yml
conda activate eet-clusters
```

For an existing clone, fetch and switch to the PR branch instead:

```bash
git fetch origin
git switch codex/diamond-scratch-directory
conda activate eet-clusters
```

Keep the repository's Python modules together. The environment supplies DIAMOND,
MCL/mcxload, Pyrodigal, PyHMMER, BioPython, NetworkX and Matplotlib. Alternatively
install `requirements-genomes.txt` with pip and install DIAMOND/MCL separately.

Your `signalp -h` output is the SignalP **4** interface. The wrapper can use its
existing environment without activating it or installing SignalP into eet-clusters:

```bash
conda run --no-capture-output -n signalP signalp -V
```

`--signalp4-env signalP` resolves the environment prefix and uses `conda run` for
both the version probe and prediction. Environment names are case-sensitive. A full
prefix such as `/home/ark/miniconda3/envs/signalP` is also accepted. If `conda` is not
on PATH, add `--conda-executable /home/ark/miniconda3/bin/conda`.

MHCscan's [Phobius implementation](https://github.com/Middle-Author-Bioinformatics/MHCscan/tree/main/bin)
uses `phobius.pl -short`. Its script resolves these files beside itself:
`decodeanhmm.64bit`, `phobius.model`, and `phobius.options`. Point `--phobius-dir` at
the **complete** existing bin folder. The GitHub decoder's executable bit may need
restoring after cloning:

```bash
chmod +x /home/ark/MAB/bin/MHCScan/bin/decodeanhmm.64bit
```

The wrapper invokes Perl directly, checks the bundle, hashes the script/decoder/
model/options, and supplies sequences through stdin. This handles input FASTA
paths containing spaces despite the legacy script's internal shell command.
The installation directory itself must have no spaces or shell metacharacters.
Version banners and warnings go to separate stderr logs, never the prediction TSV.
It does not modify or redistribute the MHCscan predictor bundle.

After preparing a domain database, run:

```bash
python genome_pipeline.py \
  --input /absolute/path/to/genomes --input-type assembly \
  --output /absolute/path/to/genome_analysis \
  --hmm /absolute/path/to/Pfam-A.hmm \
  --gram positive --threads 32 \
  --signalp4-env signalP --signalp run \
  --phobius-dir /home/ark/MAB/bin/MHCScan/bin --phobius run \
  --resume
```

PSORTb and LipoP require no installation and do not run in this command. Use
`--gram negative` for uniformly Gram-negative inputs. For mixed collections,
replace `--gram` with `--genome-metadata genomes.tsv`:

```tsv
genome	gram
isolate_1	positive
isolate_2	negative
```

`genome` is the input filename stem, removing `.gz` and the final extension.
The Gram group is never guessed. SignalP 4 receives `-t gram+` or `-t gram-` per
genome, with `-f short -s best -c 70` (the predictor's default network/truncation
settings made explicit). Do not change the networks to `notm` for a whole proteome.

`--signalp auto` / `--phobius auto` run when installed/configured, otherwise mark
results unavailable. The explicit `run` flags above require the requested tools.
A selected predictor that fails produces an actionable error. `--signalp skip`
and `--phobius skip` allow clustering/domains/motifs to finish without those tools.

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

## Native results and custom installations

`--localization-import /path/to/results` accepts:

| Filename | Format |
|---|---|
| `<genome>.signalp.tsv` | SignalP **4 short** whitespace table (12 columns, native comments) or SignalP **6** `prediction_results.txt` TSV; generation is detected from the content |
| `<genome>.phobius.txt` | Phobius short ID / TM count / SP Y-or-0 / topology |

The `.tsv` suffix is the import convention even though SignalP 4 separates columns
with whitespace. IDs may be the safe IDs in `inputs/proteomes/*.faa` or unambiguous
original IDs within each genome. Imports override execution unless explicitly
skipped. Missing rows remain unavailable; duplicate, unknown or ambiguous IDs
fail validation. No Conda environment is required for an import-only run.

`--signalp-version 4` is the execution default (`signalp` on PATH).
`--signalp-version 6` explicitly selects `signalp6` and its existing adapter;
`--signalp4-env` is only for generation 4. SignalP 6 installations remain optional.

Custom `--localization-config predictors.json` entries support `argv` (an argument
array, no shell), `output` (`stdout` or a relative file beneath `{outdir}`),
`version`, `version_argv`, `model_files`, and `stdin` (boolean). SignalP also accepts
`format`: `signalp4` or `signalp6`. For example:

```json
{
  "signalp": {
    "format": "signalp4",
    "argv": ["/home/ark/miniconda3/bin/conda", "run", "--no-capture-output",
             "-n", "signalP", "signalp", "-t", "{signalp_gram}", "-f", "short", "{input}"],
    "version_argv": ["/home/ark/miniconda3/bin/conda", "run", "--no-capture-output",
                     "-n", "signalP", "signalp", "-V"],
    "output": "stdout"
  }
}
```

Prefer the convenience flags for the installations described above. Do not supply
both a custom argv and its corresponding convenience flag. Other placeholders are
`{input}`, `{outdir}`, `{gram}`, `{gram_flag}` (legacy PSORTb), and `{threads}`.
Model paths are resolved relative to the JSON file. Supply actual installed SignalP
network/model files in `model_files` for full model provenance. The Conda option
tracks the SignalP executable automatically; it cannot detect every hidden weight
file. Add those files or change `--predictor-cache-tag` after updating untracked
models. Exact commands, versions, executable/model hashes, raw outputs, parsed
evidence, missing rows and skipped analyses are preserved in the manifest.

The earlier PSORTb/LipoP adapters are retained for compatibility behind explicit
`--psortb run` / `--lipop run`, but neither is needed by this workflow. Their native
imports/configuration remain accepted. Their historical TSV columns remain present
and say unavailable when skipped; the default HTML preview now focuses on SignalP
and Phobius. No PSORTb Unknown call is invented when it was not run.

## Interpretation

SignalP and Phobius evidence stays separate. SignalP 4 columns include generation,
SP yes/no, predicted cleavage boundary, D-score, D cutoff and networks used. Ymax
is the **first mature residue**, so Ymax=21 is reported as cleavage 20-21. A negative
call has no assigned cleavage boundary, although raw native scores/positions remain
available. The predictor's own Y/N decision is used rather than replacing it with
an invented threshold. Phobius retains TM count, SP flag and full topology.

| Evidence | Transparent feature interpretation |
|---|---|
| Phobius TM count >0 | Putative membrane protein; topology remains available |
| SignalP or Phobius SP without retained TM evidence | Export signal detected, destination unresolved |
| Predictors disagree about SP | Both calls retained and disagreement flagged |
| No SP and no TM detected | `unresolved_no_positive_targeting_evidence`; do not call cytoplasmic from these negatives |
| Missing predictor | Unavailable, not a negative prediction |

Secretion/export, membrane attachment, cell-wall association and extracellular
exposure remain separate. SignalP 4's `SignalP-TM` **network label is not a count
of retained TM helices**; Phobius provides that evidence. Neither an SP alone nor
Phobius's `o` state for a soluble protein establishes extracellular exposure.

SignalP 4 does **not** assess bacterial lipoprotein signals or provide SignalP 6's
LIPO/TATLIPO/PILIN/TAT classification. With the requested SignalP 4 + Phobius setup,
`lipoprotein_assessment` stays unavailable. Generic SP-positive proteins are never
relabeled lipoproteins. See the [SignalP 4.1 documentation](https://services.healthtech.dtu.dk/services/SignalP-4.1/)
for its output and limitations. Existing eligible installations are used as supplied;
no licensed predictor or weights are downloaded or redistributed here.

## Motifs

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
proteome matching can be expensive. `--tmpdir /path/to/large/scratch` is forwarded
to DIAMOND; its default is `$TMPDIR` or `clustering/work`. Put `--output` on a
large filesystem too, because alignments/graphs remain there. The wrapper owns and
cleans only an isolated scratch subdirectory. A failed DIAMOND pass restarts with
`--resume`; changing scratch alone preserves a completed search.
`--block-size 0.4` reduces DIAMOND scratch/memory usage relative to the retained
2.0 default, but can slightly change alignments. See the
[storage troubleshooting instructions](README.md#temporary-storage-and-disk-full-errors).
`--low-memory` streams HMM targets; `--no-plots` avoids SVG rendering for many
clusters. Full Pfam searches and thousands of per-cluster networks can take substantial
time/disk space. No million-protein performance claim is made.

## Reproducible synthetic demo and validation

```bash
python examples/make_genome_demo.py --output /tmp/eet-genome-demo
python genome_pipeline.py -i /tmp/eet-genome-demo/genomes --input-type assembly \
  -o /tmp/eet-genome-demo-analysis --hmm /tmp/eet-genome-demo/models.hmm \
  --gram positive --signalp skip --phobius skip --threads 2 --no-plots --resume
python -m pytest -q
```

The demo generates artificial genomes and HMMs. It uses real gene prediction,
clustering and domain searches; localization is explicitly skipped unless actual
predictors are selected. The MHCscan Phobius bundle was also run on the generated
24-protein collection and through the complete wrapper. SignalP 4's Conda route
was tested using labeled executable test doubles and native-format fixtures;
the user's SignalP installation is not accessible from this execution environment.
See [VALIDATION.md](VALIDATION.md) for exact checks and remaining limits.
