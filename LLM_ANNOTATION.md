# Draft one annotation per protein cluster using the OpenAI API

`summarize_clusters_llm.py` is a standalone Python **3.9+** script using only the
standard library. No OpenAI SDK or new Conda environment is required. It reads a
folder of annotated protein FASTAs and produces one proposed family label for
each file, with rationale, evidence, conflicts and limitations.

## Is this scientifically reasonable?

This is an **annotation synthesis and triage** step. It can reconcile synonymous
descriptions, choose a suitably broad family label, and expose inconsistent
annotations. The available evidence is still the existing homolog annotations;
an LLM summary does not independently establish a protein's molecular function.
NR names can be incorrect or propagated between related strains, and many copies
of one name/accession are not independent confirmations. Sequence clusters can
also merge proteins with shared domains but different full-protein functions.

The prompt therefore permits `consensus`, `mixed`, and `insufficient_evidence`.
These describe consistency of the evidence, not confidence probabilities. It
does not force a precise annotation onto unknown or heterogeneous clusters, and
forbids assigning EET activity, localization, substrates or reaction specificity
without supplied evidence. Completely generic/unmatched/empty clusters receive
a local unresolved annotation without an API call. A DUF family label is retained
as family information, without pretending its function is known.

Related published work exists:

- **InterPro/Pfam** is the closest precedent: LLMs synthesize protein-family
  descriptions from existing annotation evidence. Its newer workflow also uses
  domain architectures, curated proteins and literature, with provenance and
  curation. This script does not reproduce that richer workflow.
  https://doi.org/10.1093/nar/gkag1114
  https://www.ebi.ac.uk/training/events/llm-generated-summaries-protein-classification-interpro/
- **Hu et al.**, *Evaluation of large language models for discovery of gene set
  function*, evaluated GPT-4 and other LLMs for common-function descriptions.
  These are gene-set analyses rather than homologous protein-family annotations.
  https://doi.org/10.1038/s41592-024-02525-x
- **GeneAgent** uses GPT-4 and database-based verification for gene-set analysis.
  The verification step is relevant; this script has no live biological database
  or literature retrieval and cannot claim the same verification.
  https://doi.org/10.1038/s41592-025-02748-6

For your EET collection, review the output alongside Pfam architectures, alignment
coverage, localization evidence and representative sequence alignments. Before
generalizing labels across all clusters, compare a blinded, manually reviewed
sample with a simple most-common-description baseline, including mixed clusters,
DUFs, low-coverage hits and generic-only clusters. No biological-accuracy benchmark
or paid live-model evaluation has been performed for this script.

## API key and first run

Paste your key near the top of your **local** copy:

```python
API_KEY = "your-openai-api-key"
```

Alternatively leave that constant empty and set `OPENAI_API_KEY` in your shell.
No command-line key is needed. Saved prompts, responses and manifests exclude
the authorization header/key. The script sends annotation summaries, counts,
representative headers/IDs, and optional alignment statistics to OpenAI. It
**does not send amino-acid sequences**. It uses the Responses API with structured
JSON output and `store: false`; that flag does not by itself assert zero data
retention. Normal API usage charges and your account's model access apply.

Start with a local preview (no key, network request or API charge):

```bash
python summarize_clusters_llm.py \
  --input "$HOME/Marquette/eet/MCL2-test/clusters_nr_annotated" \
  --output "$HOME/Marquette/eet/MCL2-test/cluster_llm_annotations" \
  --dry-run --max-clusters 10
```

Inspect `work/*.request.json`, then run a paid pilot on those ten clusters:

```bash
python summarize_clusters_llm.py \
  --input "$HOME/Marquette/eet/MCL2-test/clusters_nr_annotated" \
  --output "$HOME/Marquette/eet/MCL2-test/cluster_llm_annotations" \
  --resume --max-clusters 10
```

Remove `--max-clusters 10` for the whole folder; completed pilot responses are
reused with `--resume`. The default is the pinned model
`gpt-5.4-mini-2026-03-17`, configurable with `--model`. Other models must support
the Responses API and strict structured outputs. `--reasoning-effort` is optional
and must also be supported by your chosen model. The code was checked against
official API documentation; no live request was made without your key.

## Evidence preparation

Supported inputs are top-level `.faa`, `.fa`, `.fasta` and gzipped equivalents.
The previous NR annotator's `nr_hit`, `nr_annotation`, and no-hit fields are read
explicitly. Plain `>ID description` headers also work; their descriptions are
clearly labeled as header-only evidence. First-token IDs without descriptions
are counted as unannotated. Empty files get an unresolved output row.

All headers are counted. Labels differing only in whitespace/case and a final NR
organism-style bracket are grouped; the bracket text is retained separately and
is not treated as validated taxonomy. Qualifiers such as putative, fragment,
partial and predicted are retained. Synonyms are left for the LLM to reconcile;
the program does not merge different functions by guessing.

Each request includes up to 60 distinct label groups, prioritizing informative
labels and their frequencies, plus two example headers per group. Complete
label counts and up to two full examples per label are saved in
`*.evidence.json`. Raw input FASTAs remain the full header source. The bounded
request reports omitted label/protein counts; individual labels/example headers
longer than 1000 characters are clipped and explicitly flagged. Increase
`--max-labels` / `--max-evidence-chars` if necessary. There is no random protein
sampling and omitted proteins remain in all-protein denominators.

If the input folder contains the NR annotator's `annotations.tsv`, alignment
statistics are included automatically: identity and coverage distributions,
low-coverage counts, E-values and bit scores. Table rows are checked against
filename, record ordinal, protein ID, hit accession, status and annotation.
Mismatched tables are rejected. Use `--nr-table PATH` for a matching table stored
elsewhere, or `--no-nr-table` to deliberately use headers alone. The LLM does not
invent missing quality statistics, domains, localization results or literature.

## Results and audit trail

- `cluster_annotations.tsv`: one row per selected file, with one concise
  `annotation_prediction`, evidence status, rationale, supporting/conflicting
  annotation IDs, limitations, proposed checks, model, response ID and token usage.
- `cluster_annotations.jsonl`: the same results with nested arrays/objects intact
  and the path to the corresponding evidence file.
- `work/*.evidence.json`: local complete label counts/statistics and examples.
- `work/*.request.json`: exact request body, including the full editable system
  prompt and bounded evidence; no authorization credentials.
- `work/*.response.json`: returned API JSON, including unresolved/refused or
  incomplete responses when the API returned parseable JSON.
- `manifest.json`: source checksums, model/settings, prompt hash, selected/omitted
  file counts, run status, new-request token usage, and output checksums.

Annotation IDs link the model's support/conflict claims to the supplied label
groups. The program rejects absent/duplicate evidence IDs, generic labels used
as positive support, multiline/overlong labels, missing consensus support,
refusals, incomplete responses and invalid JSON/schema. These checks validate
format and evidence references; they do not prove the biological claim is true
or eliminate prompt injection/hallucinations. Expert review remains necessary.

`--resume` requires matching source content, prompt, evidence, model and request
settings, plus an intact raw-response checksum. Errors are not cached as accepted
annotations. `--force` makes new requests and can incur charges again. Sequential
requests use a default 0.2-second pause; transient errors retry up to three times.
Authentication/model/configuration/quota errors stop promptly. Other exhausted
API errors get explicit error rows and a nonzero exit status so `--resume` can
retry only the failed clusters. Retried/time-out requests can have been processed
by the server; the program cannot guarantee exactly-once billing.

Token usage in the manifest is for new returned responses during that run,
including responses rejected by local validation; cached rows retain their
original response usage. It is not a complete billing ledger. Audits from older
runs can remain in `work`; the current TSV/JSONL lists the active selected files.
Input FASTAs are never rewritten, and no model result is automatically promoted
to a curated annotation or applied to every protein.

Official API references:
https://developers.openai.com/api/docs/guides/structured-outputs
https://developers.openai.com/api/docs/models/gpt-5.4-mini
