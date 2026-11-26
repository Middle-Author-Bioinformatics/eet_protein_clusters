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

## Interactive networks (version 1.4.0)

Local validation: **59 tests passed; 1 actual-browser test skipped**. Real
DIAMOND/MCL and PyHMMER integration checks remain enabled. New analysis checks
verify selected-HSP identity/coverage provenance, reversed coverage endpoints,
retained MCL memberships/isolates, named and repeated domain architectures,
pooled same-protein counts, Louvain communities/isolates, empty networks, display
limits, HTML/JSON escaping, and checksum-verified reconstruction of old graphs
without rerunning DIAMOND. `--no-interactive` cleanup is checked while reusing
cached HMM searches.

CI passed all **60 tests**, including the actual-browser test, on Python 3.11
and 3.12. The browser test is configured in CI on both Python versions. It opens the
self-contained HTML offline in Chromium, tests domain-name/accession search,
neighbor expansion, group and edge-threshold filters, hide-isolated, reset,
selection/details, zoom, canvas dragging, JSON download and PNG export, and checks
for JavaScript errors. CI uploads a screenshot artifact. These checks use a small
synthetic graph; they are not biological accuracy or large-browser benchmarks.

The local Chromium launcher could not start because this execution environment
blocks its process-singleton socket; the alternate headless-shell download was
also unavailable. Therefore local browser rendering is not claimed. CI provides
the actual browser execution result. No new localization model or full Pfam-A
analysis was run for this reporting change.

Protein MCL edge weights and scientific filtering are unchanged. Identity metadata
follows the selected maximum-weight HSP (ties prefer higher bit score). Protein
views do not calculate genome AAI. Domain communities use weighted Louvain on
Jaccard-weighted exported co-occurrence graphs (resolution 1, seed 17); UI filters
change only the view. Bounded HTML exports state their sampling and displayed vs
complete counts. Complete memberships, selected protein-edge metadata and domain
GraphML remain available independently of display limits.

## Pfam links and carrier proteins (version 1.4.1)

Local validation passed **66 tests**, with the actual-browser test skipped because
of the local Chromium launch limitation described above. CI passed **all 67 tests**
on Python 3.11 and 3.12, including the actual offline browser interactions.

Tests check carrier identities against selected same-protein domain sets, repeated
domains, pooled/per-cluster scope, native single-FASTA outputs, domain-negative
proteins, full pair exports, embedding caps including zero, cache invalidation,
duplicate IDs, missing architectures, mismatched domain and pair counts, and
relative links between domain and protein-identity reports. Real PyHMMER searches
exercise automatic embedding and report-only rebuilding without a rescan.
The offline browser test additionally checks version-stripped Pfam URLs, custom
accessions without invented links, node/pair lists, carrier search, TSV downloads
independent of that search, and partial/unavailable states.

The user's aggregate domain GraphML has no protein IDs or protein alignment
identities. Its updated standalone preview can show domain links and counts;
carrier lists and the protein-identity view require the original architecture
and clustering outputs. No identities are inferred from aggregate statistics.

## Standalone NR header annotation

`tests/test_nr_annotation.py` tests same-named output FASTAs, unchanged first-token
IDs, original descriptions, record ordering, unchanged sequence lines including
case/terminal stops/CRLF, gzip inputs, duplicate IDs across files, empty files,
no-hit records, missing database titles, score ranking/ties without a preference
for named functions, invalid alignments/sequences, safe output paths, stale-owned
file cleanup, and checksum-verified resume invalidation for changed inputs,
database fingerprints, search settings, or corrupt raw hits. An injected
disk-full failure checks partial-file/stamp invalidation and owned scratch cleanup.

Two integration tests run actual DIAMOND against small synthetic protein
databases: title extraction, equal-score hits, no-hit and ID-only-title records,
resume, and a short stronger match versus a broader lower-scoring qualifying
match. The latter checks that coverage filtering can retain the broader hit.
No full NR search or benchmark on TheBelly was run; its database and proteins
are unavailable in this workspace. Database content/annotation quality and NR-scale
runtime/space requirements are therefore not established by these tests.

## LLM annotation synthesis

`tests/test_llm_summaries.py` uses mocked HTTP responses: no live OpenAI requests,
API credentials or paid model evaluations are used. It tests header parsing,
species-bracket grouping, copied-hit counts, unknown/empty clusters, bounded
evidence with explicit omissions, dry-run behavior, manual/environment key input
without saved secrets, response parsing, evidence-reference validation, incomplete
responses/refusals, retry/authentication behavior, cache invalidation, partial
failures, and joining NR quality metrics to matching FASTA records. A standalone
`python -S ... --help` invocation checks operation without third-party packages.
These are software and prompt-construction checks, not evidence that the selected
model produces biologically accurate labels. Human curation and comparison with
annotation/domain evidence are required before treating outputs as annotations.

## Annotated network proteins and compact side pane (version 1.4.2)

Tests check header-label parsing, NR/no-hit/missing-title states, gzip and plain
FASTA imports, exact ID matching, duplicate/renamed ID rejection, partial imports,
metadata-cache invalidation, annotation export, automatic use of the current
clustering FASTAs, and report-only regeneration without an HMM rescan or changes
to architecture tables, memberships or alignment identities. HTML JSON escaping
is tested for hostile annotation text. Browser checks exercise annotation search
and TSV export, collapsed statistics/protein details, and expandable connection
lists that initially show eight entries. Actual user NR-annotated FASTAs were not
available; annotation imports are validated with synthetic fixtures.

## Cluster provenance and domain presence highlights (version 1.4.3)

Tests verify complete cluster/domain counts despite zero/limited embedded
proteins, node caps, domain-negative clusters, per-cluster scope, native
single-FASTA membership imports, partial/missing/changed/conflicting/invalid
assignments, exact protein-ID matching, and no-rescan report regeneration.
An integration check uses `--clustering` to attach saved membership to a native
domain analysis. Separate proteins carrying A and B highlight both domains but
do not become a same-protein pair in a cluster-specific export.
Browser checks cover visible cluster badges, selection/reset, separate community
filters, prioritizing highlighted-cluster proteins without removing other rows,
domain presence beyond the protein-preview cap, disabled legacy mapping controls,
and unchanged pooled edges. These tests use synthetic data; the user's actual
membership and architecture files were not available here.

## CDD domain profiling (version 1.5.0)

Local full-suite validation: **134 tests passed, 1 browser test skipped**. Real
DIAMOND 2.2.8, MCL/mcxload 22-282, PyHMMER, and the optional real CDD search test
were enabled. The browser interactions are exercised separately by GitHub CI;
no local Chromium execution is claimed.

The CDD adapter was run end to end with real **RPS-BLAST/BLAST formatter 2.12.0+**
and **rpsbproc 0.5.1**, using NCBI's CDD 3.22 SMART RPS database and matching CDD
annotation files. Positive controls used the actual CDD consensus sequences for
**LysM (smart00257, PSSM 197609)** and **bacterial SH3 (smart00287, PSSM 214600)**.
One artificial protein concatenated LysM-LysM-SH3b with linkers; other proteins
contained LysM only, SH3b only, or an alanine-only negative control. The controls
are database-consensus constructs, not characterized EET proteins.

The real workflow retained all four proteins, selected five domain occurrences,
counted the two LysM copies individually, created exactly one co-occurrence for
the combined protein, and created no edge for separate proteins in the second
cluster. Aligned-model coverage was 1.0 for each selected control hit. The more
generic overlapping SH3 model remained in the raw alignment table but was not
selected by the representative view. A second run reused completed search output
without rewriting its archive. These tests exposed and now cover the current
native `SESSION ... blastp ...` format and zero-length superfamily metadata.

Automated fixtures cover both native session spellings, current/historical
Non-specific hit-type capitalization, original/internal query
mapping, invalid coordinates/lengths/accessions, incomplete or omitted query
blocks, repeated domains versus duplicate alignments, model-specific threshold
retention, derived superfamily IDs, generic fallback annotations, hierarchy and
whole-protein-model exclusion, unavailable coverage, domain-negative denominators,
per-cluster versus pooled counts, native single-FASTA outputs, saved membership
imports and report-only rebuilding. Downloader fixtures check safe extraction,
complete snapshots, checksum invalidation, no-network completed reuse, subset
mismatches, and preservation of a previous snapshot after an injected failed
update. Existing Pfam outputs cannot be overwritten by a CDD resume.

The browser test checks CDD model links in addition to the existing Pfam links,
without changing carrier selection or cluster highlighting. The normal CI suite
does not download CDD. To include the real search control locally, set
`EET_CDD_TEST_DATABASE` to its prepared SMART validation snapshot and
`EET_CDD_TEST_PROTEINS` to the four-control cluster directory, with the three NCBI
executables on PATH. The optional test otherwise skips explicitly.

Not performed: a complete Cdd_LE/NCBI-subset download through the new setup CLI,
a scan of the user's cohort, an EET-function accuracy benchmark, a Conda solve,
or a macOS execution test. A real SMART database and all needed metadata were
downloaded from NCBI for the execution check; the larger downloader's archive and
failure logic use controlled fixtures. CDD imported Pfam models are not counted
as independent biological confirmation of the separate Pfam HMM analysis.

## CDD documentation integrity fix (version 1.5.1)

Snapshot-root `CDD_README.txt` and `rpsbproc_README.txt` copies are informational;
missing/edited copies no longer block an otherwise verified search database.
Checks still cover every scientific file recorded in the original manifest.
Missing/changed database or annotation files remain fatal and are reported
together, so a missing snapshot is not misdiagnosed as only a missing README.
Documentation discrepancies retain expected/actual digests and status in the
analysis manifest; the source snapshot manifest and its hashes are not edited.

Four regression cases exercise missing and edited README copies through the
analysis workflow, their warning/provenance records, unchanged source manifests,
simultaneously missing/corrupted scientific inputs, and strict checking of an
unexpected README-like filename outside the known snapshot-root locations.
Local validation passed **130 tests**, with **9 explicitly skipped** because the
external executables, browser, or optional downloaded-CDD control were unavailable.
The user's database contents were not available here; the reported README absence
was confirmed from the user's diagnostic output, not reproduced from their disk.
