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
