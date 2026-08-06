# PROVENANCE — EXP26000892 / cDNA003 HBV enrichment analysis

This repository mixes three kinds of file, and the filenames alone do not distinguish
them. Read this before assuming any script is the canonical version of itself.

| origin | meaning |
|---|---|
| **upstream** | Written for EXP26000559 (cDNA001) by Matt Wolpert. Synced from `s3://tunetx-analysis-…/mattwolpert/2026_07_02_EXP26000559_cDNA001/scripts/`. Present here unmodified except where noted. |
| **adapted** | Derived from an upstream script, renamed, and changed enough that the original would not produce these results. |
| **new** | Written for this analysis. No upstream equivalent. |
| **reimplemented** | Rewritten from a written specification because the original could not be located. **Not** byte-equivalent to the canonical version. |

If the repository was committed as recommended — upstream files in the first commit,
this analysis in the second — then `git diff HEAD~1 -- <file>` is the authoritative
record of what changed, and takes precedence over this document.

---

## Three things to know before using anything here

### 1. `hbv_transcript_classify.py` is NOT Matt's canonical module

The original could not be found on the analysis instance, in the S3 analysis bucket, or
in `__pycache__`. It was rewritten from the classification rules supplied in writing by
the original author. Every documented rule is implemented as stated and the module
carries a 30-case self-test, but it has never been diffed against the original.

It cannot be renamed — `phase3_classify.py` imports it by that exact name. Its `VERSION`
string is `hbv-transcript-classify-REIMPL-v1`, which propagates into every output header.

**Before comparing transcript-class counts against EXP26000465 or EXP26000559, re-run
Phase 3 with the canonical module and confirm the counts are unchanged.**

### 2. Several upstream scripts are present but were never executed

They arrived with the S3 sync of the whole `scripts/` directory. Their presence is not
evidence they were part of this analysis. See the `run?` column below.

### 3. `--kit-name EXP26000892` is deliberate, not a bug

`dorado demux` requires `--kit-name`, and the value must match the `name`/`kit` fields in
`barcode_arrangement_EXP26000892.toml`. It is not an ONT kit name. It was deliberately
chosen so it cannot collide with a built-in kit — see the TOML comments and
nanoporetech/dorado#1548.

---

## File manifest

### Demultiplexing (all new — the custom barcodes are unique to this experiment)

| file | origin | run? | notes |
|---|---|---|---|
| `barcode_arrangement_EXP26000892.toml` | new | ✅ | dorado custom arrangement. Masks verified base-for-base against all 8 ordered oligos (MWolig0127–0134). Scoring thresholds derived from the measured minimum pairwise barcode distance of 12/24 nt, not guessed. |
| `barcode_sequences_EXP26000892.fasta` | new | ✅ | Barcodes named `MW01`–`MW04`, **not** `BC01`–`BC04`: dorado ships built-in barcodes with the latter names and rejects a custom file that redefines them. Sequences are unchanged; only labels differ. |
| `run_dorado_demux.sh` | new | ✅ | Runs strict (`--barcode-both-ends`) and permissive passes. Contains a verification gate that aborts if the custom arrangement was silently ignored. **The gate as written is faulty** — it greps for `MW0[1-4]` in the summary, but dorado 2.x names output directories `barcode01`… and emits a `sequencing_summary.txt` whose barcode column does not contain those names. It fired a false alarm on a successful run. Verify success by counting barcode directories instead. |
| `reshape_demux_output.sh` | new | ✅ | dorado 2.x mirrors MinKNOW's nested layout (`<sample>/<run>/bam_pass/barcodeNN/`); this flattens it into the `fastq_pass/<barcode>/` layout `phase1_align.sh` globs. |
| `demux_qc.sh` | new | ✅ | Single streaming pass over the 46 GB demux summary. **The `barcode_score` analysis in it is broken** — it bins at width 1 assuming a 0–100 score, but dorado 2.x reports a 0–1 float, so that output is meaningless. Length, position and per-barcode counts are correct. |
| `count_barcode_motifs.py` | new | ✅ | Whole-read scan for barcode motifs and internal adapters. Used to establish that unclassified long reads are concatemers rather than lost transcripts. |
| `demux_custom_barcodes.py` | new | ❌ | Independent Python demultiplexer, written before confirming dorado supports custom arrangements. Validated on synthetic reads only; **never run on real data**. Retained as a cross-check implementation. |

### Pipeline setup

| file | origin | run? | notes |
|---|---|---|---|
| `setup_project_EXP26000892.sh` | new | ✅ | Creates the directory layout, symlinks `rundata` and `analysis/refs`, installs the sample sheet, and rewrites the hardcoded `PROJECT_ROOT` in the upstream phase scripts (leaving `.bak` copies). Does **not** cover scripts added after it was written. |
| `samples_EXP26000892.tsv` | new | ✅ | 7-column format the upstream scripts expect, plus a `barcode_seq` column they ignore. |

### Reference

| file | origin | run? | notes |
|---|---|---|---|
| `phase0_build_ref.sh` | upstream | ❌ | The hg38+2×HBV reference already existed at `/data/refs` and was reused. |
| `make_hbv_2x_ref.py` | upstream | ❌ | As above. |

### Alignment and QC

| file | origin | run? | notes |
|---|---|---|---|
| `phase1_align.sh` | upstream | ✅ | Run unmodified apart from the `PROJECT_ROOT` rewrite. |
| `phase1_5_readqc.sh` | upstream | ❌ | Cannot run on this experiment: it reads `$RUNDATA/sequencing_summary_*.txt` (absent after computational demux), hardcodes column positions that differ in dorado's summary, and requires a `poly_tail_length` column that does not exist in this run. Superseded by ↓ |
| `phase1_5_readqc_demux.sh` | adapted | ✅ | Reads dorado's demux summary; locates columns **by name**; drops poly(A) metrics, which have no source (poly(A) estimation was not enabled at basecalling, so neither summary contains `poly_tail_length`). Writes to the same output path as the original. |
| `phase1_5_umi.sh` | upstream | ❌ | Genome-wide pychopper for library-level duplication. Skipped: nothing downstream consumes its output, and the targeted Phase 2 route produces the `hbv.umi.bam` it would have supplied. |

### HBV extraction and UMI deduplication

| file | origin | run? | notes |
|---|---|---|---|
| `phase2_hbv_umi.sh` | upstream | ❌ | Superseded by ↓ |
| `phase2_hbv_umi_v2.sh` | adapted | ✅ | Threads from `nproc` rather than hardcoded 4 (the original assumed tens-to-hundreds of HBV reads; this run has 0.4–1.3 M); per-stage timing; samtools ≥1.16 version gate; `PROJECT_ROOT` as an env override; and an optional homologous-concatemer filter. **The filter was OFF for the EXP26000892 results in this repository** — see below. |
| `filter_concatemers.py` | new | ⚠️ | Removes reads carrying an adapter or barcode motif >200 bp from both ends. Must run **before** pychopper, which trims the primers away. Run once on the 1 ng library for validation, then reverted; **not applied to the published EXP26000892 numbers**. Measured rate: 0.97% of HBV reads but 3.5% of unique molecules (3,800 → 3,667). Enabled from EXP26000993 onward. |
| `phase2_extract_hbv.py` | upstream | ❌ | Alternative HBV extraction route; the `phase2_hbv_umi` path was used instead. |
| `phase2b_dedup_read_ids.sh` | upstream | ❌ | Not used. |

### Transcript classification

| file | origin | run? | notes |
|---|---|---|---|
| `hbv_transcript_classify.py` | **reimplemented** | ✅ | See warning above. Lives in `/data/shared/`, imported by `phase3_classify.py` via `PROJECT_ROOT.parent/"shared"`. |
| `phase3_classify.py` | upstream | ✅ | Run unmodified. Classifies **reads**, not molecules — see below. |
| `phase3b_classify_molecules.py` | new | ✅ | Classifies UMI-deduplicated **molecules**, and reports reads-per-molecule within each class. Written because duplication in this experiment is length-biased, which the upstream read-level approach assumes away. |

### Quantification and reporting

| file | origin | run? | notes |
|---|---|---|---|
| `phase4_quantify.py` | upstream | ✅ | Run **twice**, unmodified, by swapping which parquet sits at `hbv_classified.parquet` — once for read-level, once for molecule-level. Outputs renamed `_reads` / `_molecules` afterward. |
| `phase5_report.py` | upstream | ❌ | Its narrative describes EXP26000559's 8-library poly(A)×PCR matrix. Superseded by ↓ |
| `phase5_report_EXP26000892.py` | adapted | ✅ | Reuses the upstream figure functions; new narrative for this experiment; `--level reads\|molecules`; poly(A) panel removed; adds bounded dose-response and jackpotting figures. |

### Additional analyses (all new)

| file | origin | run? | notes |
|---|---|---|---|
| `cross_library_umi_check.py` | new | ✅ | Shared UMI-and-position analysis across libraries. Measured <1% cross-library carry-over. Validated against synthetic data with a planted 25% contamination rate (recovered 23.6%). |
| `umi_saturation_curve.sh` | new | ✅ | Subsamples to 10/25/50/75/100% and deduplicates at each depth, with both `directional` and `unique`. Established that **no library reached saturation** — all molecule counts are lower bounds. |
| `library_qc_funnel.py` | new | ❌ | Consolidated per-library funnel. Written but not run. |
| `test_midstrand.sh` | new | ❌ | A/B test of `midstrand_flank_score`. Written but not run — the question it addressed was answered another way. |

---

## Execution order actually used

```
run_dorado_demux.sh            (strict pass; ~8.5 h at 4 threads on 131 M reads)
reshape_demux_output.sh
demux_qc.sh
count_barcode_motifs.py        (unclassified + one classified library, as control)
setup_project_EXP26000892.sh
phase1_align.sh                (one library at a time)
phase1_5_readqc_demux.sh
phase2_hbv_umi_v2.sh           (CONCAT_FILTER=0 for these results)
phase3_classify.py             → read-level
phase3b_classify_molecules.py  → molecule-level
cross_library_umi_check.py
umi_saturation_curve.sh
phase4_quantify.py             ×2, parquet swapped between runs
phase5_report_EXP26000892.py   → two reports
```

## Known limitations of these results

- **Molecule counts are lower bounds.** No library reached saturation; deeper sequencing would recover more.
- **UMI counts are inflated at low input.** PCR jackpotting is severe below 10 ng — in the 0.1 ng library one molecule holds 440,452 reads and the top 1% of molecules hold 98.7% of all reads. Sequencing errors in the 28-nt UMI of such a molecule generate thousands of spurious singletons that deduplication does not fully collapse. Complexity is therefore reported as a range, bounded below by distinct alignment coordinates.
- **MAPQ is not a quality filter on this reference.** A read contained within one copy of the 2× HBV contig matches both copies equally and receives MAPQ 0 by construction. `hbv_reads_mapq20` effectively counts reads that uniquely span the copy junction, biasing toward long transcripts.
- **Homologous concatemers were not removed** from these results (0.97% of reads, 3.5% of molecules, measured on the 1 ng library).
- **All libraries entered hybridisation capture at equal mass** (375 ng each). This measures whether a low-input library still contains HBV, not how much HBV can be recovered from a given input.
- **The transcript classifier is a reimplementation.** See above.
