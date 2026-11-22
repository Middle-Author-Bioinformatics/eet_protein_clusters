# Annotate cluster FASTAs with their top NR hits

`annotate_clusters_nr.py` is self-contained: Python **3.9+** and an installed
DIAMOND executable are its only dependencies. It can be copied and run on its
own; it does not import repository modules, BioPython, MCL or localization tools.
The existing `eet-clusters` environment already supplies DIAMOND. If needed:

```bash
conda install -n eet-clusters -c conda-forge -c bioconda diamond
conda activate eet-clusters
diamond version
```

Run from any directory, using the repository's script path:

```bash
mkdir -p /home/ark/eet_scratch

python /home/ark/MAB/bin/eet_protein_clusters/annotate_clusters_nr.py \
  --input "$HOME/Marquette/eet/MCL2-test/clusters" \
  --output "$HOME/Marquette/eet/MCL2-test/clusters_nr_annotated" \
  --db /home/ark/databases/nr.dmnd \
  --threads 32 --tmpdir /home/ark/eet_scratch --resume
```

The script reads top-level `.faa`, `.fa`, `.fasta` and gzipped equivalents. It
combines all records into **one DIAMOND blastp run**, using safe internal query
IDs, and restores the original filenames, file membership, record order, protein
IDs and descriptions. Original sequence lines (including case and terminal stops)
are copied unchanged. Gzipped inputs produce gzipped annotated files of the same
name; compression bytes can differ. A duplicate ID in different input files is
kept separately, and the TSV identifies records by filename and ordinal.

For DIAMOND's search input only, sequences are uppercased and terminal `*` symbols
removed. Empty sequences, internal stops and invalid protein characters fail with
the source file/protein identified. Empty files are retained when other files
contain proteins; an entirely empty protein collection is rejected.

An input header such as:

```text
>GCA_000006785.2_ASM678v2AAK34500.1
```

becomes, **for an illustrative hit**:

```text
>GCA_000006785.2_ASM678v2AAK34500.1 | nr_hit=WP_123456789.1 | nr_annotation=Protein description [Organism name]
```

The original first whitespace-delimited ID stays the same. Any original header
description is retained before the appended NR fields. Hit statistics are in the
TSV, keeping FASTA headers reasonably short. These are top-hit annotations, not
experimentally validated functional assignments. Partial matches and proteins
with different domain architectures need inspection before assigning function.

## Hit definition and thresholds

Defaults: **sensitive**, E-value **<=1e-5**, block size **0.4**, and no additional
identity or coverage filter. Coverage is reported for both query and subject.
To require broader matches, add, for example:

```bash
--query-cover 50 --subject-cover 50 --min-identity 30
```

All three thresholds use **percentages from 0 to 100**, not fractions. Coverage
uses the coordinate span of the best single HSP; disjoint HSPs are not pooled.
DIAMOND receives the filters explicitly, and the parser checks them again.

DIAMOND runs with `--max-target-seqs 25 --max-hsps 1`, retaining a target candidate
pool before the script selects one annotation. The script selects the highest
qualifying reported bit score, then resolves
ties by lower E-value, greater query coverage, greater subject coverage, greater
identity, and lexical subject ID/title. `best_score_ties` records how many
qualifying reported targets shared the highest bit score, within this candidate
limit; it is not a count of every equally scoring NR entry. `--max-target-seqs`
can increase the pool, or `0` requests all reported targets (potentially very
large output). The script deliberately does not use `--top 0`: in a real DIAMOND
test, that option returned no hit with a 90% query-coverage threshold when the
strongest match was short, whereas a target candidate pool retained a lower-scoring
full-length match that met the threshold. There is no preference
for a named function over a hypothetical/uncharacterized annotation. Self/exact
matches in NR are permitted. DIAMOND is a heuristic search; this is its best
reported qualifying alignment, not an exhaustive alignment guarantee.

Annotations come from DIAMOND's `stitle` field, including the organism text when
present. A matching leading subject ID is removed from the display annotation;
the complete original title remains in `annotations.tsv`. No taxonomy database
or accession-to-taxid mapping is needed. If NR was built with ID-only FASTA
headers, the annotation is marked `hit_title_unavailable`; the script cannot
recover a description absent from the database. Rebuild such a database from
reference FASTA with full descriptive headers if descriptions are needed.

Proteins without a qualifying match are retained with:

```text
>original_ID original description | nr_status=no_qualifying_hit
```

## Outputs and reruns

The output folder directly contains the annotated FASTAs with their original
names, plus:

- `annotations.tsv`: every input protein, including unmatched records; source
  file/record/ID, status, hit accession, annotation, full title, identity, query and
  subject coverage, E-value, bit score, alignment coordinates/lengths and ties.
- `cluster_summary.tsv`: counts of proteins, matched/unmatched proteins and hits
  lacking titles per file.
- `manifest.json`: completion/failure status, commands, versions, parameters,
  input/output checksums and database provenance.
- `run.log`, `work/diamond_blastp.log`, and database/version logs.
- `work/nr_hits.tsv`: raw candidate hits with the custom columns listed in
  the manifest; query IDs map through `work/queries.sqlite`.

The input folder is never modified. Use a separate sibling output folder rather
than an input folder or a parent/child of it. Nonempty outputs are accepted only
with `--resume` or `--force` and a manifest identifying this script as their owner.
When input files are removed, a successful rerun removes only stale FASTAs listed
as generated in the prior manifest; unrelated files are retained.

`--resume` reuses a completed, checksum-verified hit file only when the inputs,
database fingerprint, DIAMOND version/executable and search settings match.
Changing scratch alone preserves a completed cache. Interrupted/failed searches
restart; partially written results are not cached. `--force` requests a new
search. A failed manifest marks any remaining older FASTAs as incomplete for the
current run; a successful manifest has `status: complete`.

By default, the NR fingerprint records the resolved file path, size and modification
time, plus DIAMOND `dbinfo` output in its log. This avoids reading an enormous NR
file solely to hash it on every run. Add `--db-version NR_RELEASE_OR_BUILD_DATE`
for a release label, or `--hash-db` to include the full database SHA256 in the
manifest/cache key. Keep the same database snapshot for reproducibility; size and
timestamp alone cannot detect every possible database replacement.

Scratch uses `--tmpdir`, then `$TMPDIR`, otherwise `OUTPUT/work`. Each run owns a
unique scratch subdirectory and removes it on ordinary completion/failure. An
abrupt kill or power loss may leave scratch behind. Put the output and scratch
on filesystems with sufficient space; `--block-size` is configurable. The error
includes the DIAMOND log tail if the search fails.

DIAMOND's authoritative option definitions:
https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
