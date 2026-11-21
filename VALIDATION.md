# Validation

## Native collection output (version 1.2.0)

The current suite reports **21 passed, 2 skipped** on the same Linux/Python 3.12
analysis stack. The existing DIAMOND/MCL executable tests remain skipped in this
session. Existing folder-mode domain tests still pass.

The single-FASTA tests now assert a genuinely collection-wide schema: no public
cluster/outside-comparison columns, no generated cluster folder or cluster-summary
file, collection-specific GraphML metadata, a prevalence chart and a domain-level
HTML table. Raw hit/statistic values are compared to equivalent folder analysis
using the appropriate mode-specific columns; repeat counts and negative-protein
denominators agree. Native-to-folder and folder-to-native resume transitions remove
stale output formats. No-plots mode removes both SVGs.

An additional cache test proves that a matching published v1.1.0 cache migrates to
the new output without invoking HMMER. An independent smoke test also executed the
actual published v1.1.0 scripts from git, then resumed using v1.2.0 and verified the
raw domtblout file was unchanged. Reuse validates input/model/search settings and
raw/database checksums. Migration is disabled if the search-core implementation
changes; presentation-only edits do not invalidate new-format search caches.

The four supplied aggregate result files were also reformatted for visual review:
4,984 input proteins, 142 detected domain families, 7,292 domain copies and 156
observed domain pairs. Numerical prevalence and edge statistics were preserved;
GraphML was read back and the prevalence chart was inspected as a raster preview.
No HMM searches were rerun for that reconstruction. Overall numbers of annotated
and unannotated proteins and per-protein architectures cannot be recovered from
aggregate domain/pair counts; those values are explicitly unavailable in the preview.
A real script run has the protein-level information and reports those counts.

## Standalone FASTA extension (version 1.1.0)

At version 1.1.0, the suite reported **20 passed, 2 skipped** on Linux/Python 3.12 with
BioPython 1.88, PyHMMER 0.12.3, NetworkX 3.7 and Matplotlib 3.10.8.
The two existing external DIAMOND/MCL tests skipped because those executables
were unavailable in this session; no clustering implementation was changed.
The domain tests executed the real PyHMMER search engine against controlled HMMs.

Seven new tests cover:

- Identical prevalence, raw hits, architectures and summaries for a single FASTA
  and the equivalent one-cluster folder, including repeated domains and proteins
  without hits in the denominator.
- Standalone GraphML/SVG/co-occurrence exports, explicit input-mode metadata,
  single-collection report wording and unavailable outside-collection statistics.
- Gzipped FASTA input, a path containing spaces, actual subprocess CLI execution,
  completed-search cache reuse and preservation of original protein IDs.
- Removal of stale root SVGs with `--no-plots`, and removal of standalone exports
  when resuming the same output folder with a cluster-folder input.
- A domain-negative single-FASTA collection with a zero-domain summary and empty graph.
- Empty and duplicate-ID FASTAs rejected without a successful completion manifest.
- Missing input paths and output directories containing input files rejected safely.

`--input`, `-i`, and the existing `--clusters` alias are exercised. Full Pfam-A
and a user-supplied biological collection were not run for this extension.

## Original workflow validation (version 1.0.0)

Validated on Linux with Python 3.12, BioPython 1.88, PyHMMER 0.12.3,
NetworkX 3.7, Matplotlib 3.10.8 and pytest 9.2.0.

**15 tests passed**, including execution of the actual external programs:

- DIAMOND **2.2.8**, official Linux release executable.
- MCL and mcxload **22-282**, Debian distribution build with its libtingea dependency.
- PyHMMER's actual HMMER search engine, including sequence and domain GA thresholds.

No mock executable replaced DIAMOND, MCL, mcxload, or the HMM engine in those tests.
The Conda installation command/environment solve was not run in this execution
environment; the tested binaries were provisioned independently for verification.

Run the suite from the unzipped package:

```bash
python -m pip install pytest
python -m pytest tests -q
```

External DIAMOND/MCL tests skip explicitly if those executables are absent from
PATH. Domain tests require the dependencies in requirements.txt. The suite itself
does not download databases or binaries.

## Tested behavior

- Protein IDs normalized safely and original headers retained; duplicate original
  IDs rejected; original IDs that resemble internal IDs mapped without ambiguity.
- Similarity identity, E-value and **both-protein coverage** filters; invalid
  numerical fields/coordinates and unknown IDs rejected.
- Reciprocal and repeated alignment edges deduplicated using maximum weight;
  zero E-values capped rather than becoming infinite weights.
- MCL output verified as a partition; duplicate members and missing non-isolated
  vertices rejected. Isolates and entirely edgeless inputs retained as singletons.
- Real DIAMOND -> mcxload -> MCL execution on synthetic families: two two-member
  families plus one unrelated singleton. Every input protein appeared once.
- DIAMOND cache reuse across a changed MCL inflation setting.
- Repeated domains retained as copies while counting domain presence once per
  protein for prevalence and pair statistics.
- Default overlapping-model suppression, optional retention of all overlaps,
  and minimum model-coverage filtering.
- Real HMM searches against two synthetic HMMs/controlled protein architectures:
  A-A-B, A-only, B-only and domain-negative proteins.
- Denominators include domain-negative proteins. For a three-member cluster,
  A prevalence was 2/3 with three copies and 1.5 copies per positive protein;
  A-B occurred in one protein with Jaccard 0.5 and lift 1.5.
- Domains found in different proteins in the same cluster do **not** create a
  within-protein co-occurrence edge.
- Domain-negative clusters generate zero-domain summaries and empty-network plots.
- Missing-GA models produce an error unless an explicit E-value fallback is supplied.
- Low-memory target streaming, one-model batching, GraphML, SVG, no-plot mode,
  minimum edge counts, architecture tables and zero-hit output were exercised.
- Cached domain searches reused for identical inputs; cluster-name/membership changes
  invalidate the cached mapping and stale generated summaries are removed.

## Scope and limitations

The HMM tests used controlled models generated from synthetic amino-acid sequences,
not an experimental EET validation panel. A full Pfam-A scan was **not** run during
validation, and no user protein collection was supplied. Search/database coverage
and sensitivity in the user's biological data therefore remain unmeasured.

No macOS installation, million-protein workload, real-cohort MCL parameter sweep,
or interactive browser-rendering benchmark was performed. GraphML was read back
with NetworkX, SVG files were generated, and TSV values were asserted. Dense
plots are deliberately limited; complete data remain in the TSVs/GraphML.

Pfam overlap handling is a documented greedy heuristic rather than the curated
PfamScan clan/nested-domain algorithm. Prevalence and network statistics are
observational, not physical interactions or EET activity calls. Real data should
be inspected for partial proteins, sequence redundancy, taxon imbalance, domain
model coverage and clustering sensitivity before biological interpretation.

## Genome workflow (version 1.3.0)

The genome wrapper was validated with **36 passing tests** under Python 3.12,
BioPython 1.88, Pyrodigal 3.7.1, PyHMMER 0.12.3, DIAMOND 2.2.8 and MCL 22-282.
The external-engine tests were run, not skipped. CI now installs DIAMOND/MCL and
Pyrodigal and runs the suite on Python 3.11 and 3.12.

New checks cover:

- Exact FeGenie motif regexes, independent variant totals, canonical CXXCH vs
  MHCscan CXXCH+CXXXCH counts, negative sequences, source-compatible non-overlap
  semantics and 1-based match coordinates.
- Real deterministic Pyrodigal predictions from synthetic bacterial coding
  sequences, restricted translation-table models and reproducible extraction.
- GenBank supplied/derived translations, alternative initiator translation,
  reverse strands, compound/fuzzy locations, exclusion of pseudogenes, stable safe
  IDs and original IDs repeated across genomes or at different annotated loci.
- Protein-only folder inputs preserve IDs and explicitly lack genomic coordinates.
- Explicit Gram metadata, ambiguous input rejection, correct PSORTb `-p`/`-n`
  selection and SignalP organism/mode/thread arguments.
- Native PSORTb long/terse, SignalP 6, Phobius short and LipoP short parsing;
  duplicate/unknown/ambiguous ID rejection; partial imports retain missing calls
  as unavailable; requested missing executables give actionable errors.
- Unknown PSORTb calls gain only feature-supported interpretations. SP-only
  results do not imply extracellular exposure. Confident baseline calls are kept
  when supplementary predictions disagree; disagreements are flagged.
- Predictor command execution, native output validation, checksum caches,
  corrupted-output invalidation, import precedence, and configuration validation.
- A whole-genome **GenBank-to-clustering-to-domain-to-annotation** integration test
  runs actual DIAMOND/MCL/PyHMMER against synthetic A-A-B/domain-negative proteins.
  It verifies repeated domains, all-protein retention, known motif counts, mixed
  Gram classes, partial localization imports, joined mappings and reports.
  A labeled executable test double emits PSORTb-format Unknown calls to test
  orchestration; it does not run PSORTb models. Resume reuses unchanged HMM and
  localization outputs, while changing a genome invalidates dependent caches.

A separate documented CLI demonstration generated two artificial assemblies,
called genes with real Pyrodigal, ran DIAMOND/MCL and the synthetic HMM search,
and produced **24 proteins in 2 clusters**, 6 proteins with >=2 canonical CXXCH
motifs and 6 proteins with the FeGenie hematite pattern. Localization was explicitly
skipped and all 24 proteins were reported as unassessed. This demonstration is
reproducible with `examples/make_genome_demo.py` and the commands in
[GENOME_WORKFLOW.md](GENOME_WORKFLOW.md). The artificial HMMs are not Pfam models
or evidence of EET.

**Not run:** real PSORTb, SignalP 6, Phobius or LipoP model inference; a full Pfam-A
scan; an experimental localization or EET benchmark; a million-protein workload;
or the optional provider containers. Predictor format fixtures and test doubles
validate adapters and rules, not their biological accuracy or the local installation
of licensed tools. No user genomes were provided. SignalP/Phobius/LipoP and their
weights are not redistributed. Supplementary refinements are transparent feature
annotations, not a validated replacement for PSORTb.

## SignalP 4 + MHCscan Phobius update (version 1.3.1)

**45 tests passed locally**, with the real DIAMOND/MCL tests enabled. Nine new
cases cover SignalP 4 positive/negative native short output, D-score/cutoff and
network retention, cleavage coordinates (Ymax=21 -> 20-21), malformed/duplicate
records, Gram-specific arguments, Conda environment discovery/dispatch/version
probes, skipped-tool defaults, unavailable lipoprotein assessment, Phobius bundle
validation, executable permissions and separation of stderr from prediction tables.
A full wrapper test uses a clearly labeled Conda/SignalP executable test double,
imports Phobius output, runs real clustering/HMM searches, and checks resumption
from an output directory containing spaces.

The **actual Phobius 1.01 script, decoder, options and model** already present in
MHCscan were run locally on the 24-protein artificial assembly demo, both directly
and through `genome_pipeline.py --phobius-dir`. All 24 result rows parsed and were
reported assessed. The script's version banner and uninitialized-value warnings
were retained in stderr logs rather than contaminating the machine-readable table.
The synthetic collection had zero predicted TM helices; it is an execution/format
check, not an accuracy benchmark. The bundle was copied only into temporary test
storage; no predictor binary/model is added to this repository.

The real SignalP 4 installation on TheBelly is not accessible here. Its CLI was
verified against the user-provided help and its native output against the DTU
SignalP 4.1 documentation. SignalP model inference was therefore **not** run;
Conda dispatch uses a labeled test double. PSORTb and LipoP are now explicitly
skipped by default and were not run. The earlier version 1.3.0 limitations above
are historical; actual Phobius execution was added in this update. Full Pfam-A
and biological accuracy/scale benchmarks remain unperformed.

## DIAMOND scratch storage update (version 1.3.2)

**53 tests passed locally**, including real DIAMOND 2.2.8 and MCL 22-282. Eight
new storage tests cover explicit scratch paths (including spaces), `$TMPDIR`
fallback, output/work fallback, precedence, cleanup of owned scratch while leaving
unrelated files intact, injected disk-full failure, invalidation of an earlier
completion stamp, restart after failure, reuse of completed searches when only
scratch changes, invalid scratch paths and preservation of actionable errors when
an exhausted output filesystem also prevents writing the failure manifest.

The real DIAMOND test passes `--tmpdir` beneath a custom parent and uses block size
0.4. The genome integration test verifies forwarding the scratch option to the
clustering subprocess. Partial alignments are never treated as a completed cache.
No physical disk was filled during testing; disk-full conditions were injected.
The user's large `AllMHCGenes.faa` dataset and TheBelly filesystems are not
accessible here, so this does not establish peak disk requirements or prove which
mount, quota or inode limit caused that specific run to fail. The established
block-size default 2.0 is retained; the documented recovery command explicitly
selects 0.4 and notes the potential for slightly different alignments.
