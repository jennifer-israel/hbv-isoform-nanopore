# Pipeline overview: sequencing run to HBV molecule counts

What happens between a finished nanopore run and a per-library table of unique HBV molecules,
and why each step is there. No commands — see the demux runbook for those, and the QC metrics
document for the checks that run alongside.

The endpoint is a read funnel like this, one row per library:

| bin | sample | single-end | both-ends | kept | HBV | % HBV | filtered | molecules | dup |
|---|---|---|---|---|---|---|---|---|---|
| bc08 | YecgNT | 1,543,835 | 931,051 | 60.3% | 262,919 | 28.24% | 254,678 | **192,779** | 24.3% |

Every column is a stage below.

---

## What makes this assay unusual

Three things drive most of the pipeline's shape.

**The barcodes are custom and internal.** Libraries are barcoded by PCR with oligos carried
inside the cDNA construct, not by ONT native barcoding. MinKNOW therefore does no
demultiplexing — every read lands in `bam_pass/` unclassified — and basecalling must run with
`--no-trim` so the adapters and barcodes survive to be read later. If that flag is missed, the
run cannot be demultiplexed at all and has to be re-basecalled from POD5.

**HBV is a circular 3.2 kb genome.** Transcripts routinely wrap past the origin, so a linear
reference would break them into two alignments. The reference doubles HBV end-to-end to
6,364 bp, letting wrapping transcripts align contiguously.

That choice has consequences that shape the rest of the pipeline. A read which doesn't span the
junction matches both copies equally and gets **MAPQ 0 by construction**. This is why alignment
happens in **two passes against two references**, long reads against the 2× and the leftover
short fraction against a 1×.

**HBV RNA is a tiny fraction of the input.** Libraries are enriched by hybridisation capture
against HBV probes. That works well, but probe capacity is finite, so pooling a high-viral-load
library alongside a low-load one lets the former consume the capture and suppress the latter.

---

## Stage 1 — Run to S3

The sequencer writes basecalled reads as BAM (`bam_pass/`, `bam_fail/`) and the same reads
again as FASTQ (`fastq_pass/`, `fastq_fail/`). Only `bam_pass` is needed: it's the same data as
`fastq_pass` and the BAM retains tags, so pulling both roughly doubles the transfer for nothing.

## Stage 2 — Demultiplexing

`dorado demux` with a custom barcode arrangement: a TOML describing the construct's flanking
sequences and a FASTA of the run's barcode oligos.

The important choice here is **requiring the barcode at both ends**. A read carrying agreeing
barcodes at both ends is confidently assigned. A read where the two ends disagree is a fusion
of molecules from two different libraries, and it's discarded.

Running demux **both ways** — single-end and both-ends — gives the funnel's first three columns
and, with them, a direct measurement of how often a barcode is detected at all.

What both-ends demultiplexing does *not* catch: a fusion of two molecules from the **same**
library. Both its outer barcodes agree, so it looks like a perfectly valid read. Those are
removed later, at stage 5.

## Stage 3 — Reshape and project layout

Dorado writes a nested directory tree; the reshape step flattens it and renames barcodes to the
`custom_bcNN` convention the downstream scripts expect. A sample sheet maps each barcode to a
library, its input amount, PCR cycle count and poly(A) route.

Two things belong here rather than later. Barcode-to-library identity should be **confirmed
against the reads**, because dorado names its output directories positionally and nothing
downstream can detect a transposed sample sheet. And the `sample_key` chosen now names every
output directory from this point on, so renaming it later orphans completed work.

## Stage 4 — First alignment pass

minimap2 in spliced-alignment mode against a composite reference: **hg38 plus the doubled HBV
genome**. Both are needed in one index — aligning to HBV alone would recruit host reads with
partial similarity, and the host component gives the denominator for on-target percentage.

This pass handles reads long enough to be placed spliced. What it can't place is set aside for
the second pass at stage 6 rather than discarded.

Two consequences worth knowing. Reads from **chimeric mouse models** won't map, since no mouse
genome is present, which depresses the overall mapping rate without saying anything about HBV
recovery. And **HBV as a share of *mapped* reads is a misleading metric** when a large fraction
of the sample can't map at all — HBV over *total* reads is the honest on-target figure.

A separate read-QC pass summarises read length and quality per library. Read length turns out to
be the variable that most often limits what can be asked of a dataset: transcript classification
needs reads long enough to span a transcript, and fragmented input can't provide that no matter
how good the capture.

## Stage 5 — HBV extraction, concatemer removal, UMI deduplication

Four steps in one pass, and the order matters.

**Extract** reads aligning to the HBV contig. Everything after this operates on hundreds of
thousands of reads rather than the whole run.

**Remove homologous concatemers** — the same-library fusions that survived demultiplexing. They
are found by looking for an adapter or barcode motif more than 200 bp from *both* read ends,
which cannot happen within an intact molecule. Left in, they cause two problems: a fused pair
carries two UMIs, so deduplication keeps one and silently loses the other; and on the doubled
reference a fused pair can align across the junction and be scored as a genuine
genome-spanning transcript.

**Orient and extract the UMI** with pychopper, which also trims the primers. This is why
concatemer removal has to come first — after trimming, the adapter sequence is gone whether or
not the read was a fusion, so the evidence no longer exists.

**Deduplicate** on the UMI plus alignment position. This is the step that converts read counts
into molecule counts. Because PCR amplification is exponential, a molecule amplified early can
account for a large share of a library's reads, so reads are a measure of amplification luck
and molecules are the quantity of interest.

## Stage 6 — Second alignment pass: the short fraction

Stage 4 discards reads it can't place in spliced mode, and in fragmented input — serum and
plasma especially — that can be most of the data. Those reads are real RNA, so they get a
second pass with different settings.

**Against a single-copy HBV reference, not the doubled one.** This is the key difference, and
it follows from why the reference was doubled in the first place. A short fragment cannot span
the copy junction, so on the 2× reference it matches **both copies equally well** — it
multi-maps, MAPQ collapses to 0, and it becomes indistinguishable from a genuinely ambiguous
alignment. On a single-copy reference there is one locus and the fragment maps uniquely.
Doubling buys nothing for a molecule too short to wrap, and costs the ability to trust the
alignment.

**And in short-read mode rather than spliced mode**, since a 60–150 bp fragment has no introns
to find and the spliced preset's gap penalties are tuned for a different problem.

The 1× reference is derived from the 2× by extracting its first 3,182 bp, so the two are
guaranteed to be the same sequence rather than separately sourced.

The short fraction then goes through the same UMI extraction and deduplication as the long
fraction, giving its own molecule count.

## Stage 7 — Combining the two passes

The two passes measure the same libraries at different fragment lengths, and they combine
differently depending on the question.

**For abundance, add them.** Total unique HBV molecules is the long-fraction count plus the
short-fraction count. Omitting the short fraction understates recovery by a variable amount —
it has run from a couple of percent of molecules in intact cellular RNA to a large share in
degraded serum.

**For composition, use the long fraction alone.** Every transcript class is gated on a minimum
span, and the smallest gate is several hundred bases. Short-fraction molecules satisfy none of
them, so folding them into a classification inflates the unclassified bucket without adding
information. Report the short-fraction count alongside rather than inside.

**To compare positions between the passes, coordinates must be folded.** The long fraction is
in 2× coordinates spanning 6,364 bp; the short fraction is in 1× coordinates spanning 3,182 bp.
A long-fraction position maps into the short fraction's frame as `((pos − 1) mod 3182) + 1`.
Comparing them without folding silently compares two different coordinate systems — positions
past 3,182 in one frame don't exist in the other.

**The short fraction is also informative on its own.** Those fragments were never efficiently
captured, being shorter than the probes, so their HBV content relative to the long fraction
indicates what the capture is and isn't reaching, and their positional distribution shows
whether the short material is random degradation or has structure.

## Stage 8 — Downstream

With molecules in hand: transcript classification assigns each molecule to an HBV transcript
class from its 5′ end, span and splice structure; quantification builds per-library recovery and
composition tables; a report assembles the figures. Classification is the step most often
blocked by read length, since every transcript class is gated on a minimum span.

---

## The read funnel

| column | what it is | derivation |
|---|---|---|
| **bin** | barcode | demux output directory |
| **sample** | library | sample sheet |
| **single-end** | reads assigned needing a barcode at one end | demux without `--barcode-both-ends` |
| **both-ends** | reads assigned needing agreeing barcodes at both ends | demux with `--barcode-both-ends` |
| **kept** | both-ends ÷ single-end | how often the second end was also readable |
| **HBV** | reads aligning to the HBV contig | stage 4, first pass |
| **% HBV** | HBV ÷ both-ends | on-target fraction |
| **filtered** | HBV reads surviving concatemer removal | stage 5 |
| **molecules** | distinct molecules after UMI deduplication | stage 5 |
| **dup** | 1 − molecules ÷ filtered | duplication rate |

Every column here describes the **first pass**. Where the short fraction matters, the table
extends with the second pass and a combined total:

| additional column | what it is |
|---|---|
| **short HBV** | reads from the unmapped fraction aligning to the 1× reference |
| **short molecules** | distinct molecules after deduplicating those |
| **total molecules** | long + short — the abundance figure |
| **short %** | short molecules ÷ total molecules |

`short %` is the number that says whether the second pass was worth running for a given sample
type. It's small for intact cellular RNA and large for degraded serum.

### Reading across a row

**kept** is an estimate of single-end barcode detection probability. If detection at one end
succeeds with probability *p*, single-end assignment scales with *p* and both-ends with *p*², so
their ratio is *p* — measured, rather than assumed. It varies with read length and quality,
since a short or degraded read is less likely to present both ends intact.

**% HBV** is capture performance, and it is only comparable between libraries that shared a
capture pool. A library pooled with a much larger one is competing for finite probe capacity, so
a low value may describe the pool rather than the sample.

**filtered vs HBV** is the concatemer loss. Detection needs a motif more than 200 bp from both
ends, so short-read libraries report a low rate whether or not fusions occurred — read it
alongside read length.

**dup** rises with sequencing depth for a fixed library, so it is not comparable between
libraries sequenced to different depths. Two libraries with the same duplication rate at very
different depths do not have the same complexity.

**molecules** is the endpoint, and it is a floor unless the library is saturated. If additional
sequencing would still be finding new molecules, the count reflects how much was sequenced as
much as what was in the sample.

### What a row can tell you at a glance

A high **kept** with a low **% HBV** means the library sequenced cleanly but capture didn't
enrich it — either genuinely low viral load, or suppression by a larger library in the pool.

A low **kept** points at read length or end quality rather than at HBV.

A high **dup** with few **molecules** means the library was amplified past its complexity: more
sequencing will add reads without adding information.

A large gap between **HBV** and **filtered** means fusion events are consuming the library, which
points at ligation conditions and adapter-to-insert ratio.

---

## Where things go wrong

| symptom | usually means |
|---|---|
| nothing demultiplexes | adapters were trimmed at basecalling; needs re-basecalling from POD5 |
| one library far below the others | pool imbalance, or that library's input |
| high mapping rate, low % HBV | capture didn't enrich — check pool composition before the sample |
| low mapping rate overall | reference is missing a genome present in the sample (e.g. mouse) |
| % HBV looks fine, molecules very low | amplification, not capture — check duplication and input |
| classification mostly unclassified | read length, upstream of anything the pipeline can fix |
| short fraction is most of the library | fragmented input; the first pass alone will understate recovery |
| short-fraction alignments all MAPQ 0 | second pass ran against the 2× reference instead of the 1× |
| coverage profiles disagree between passes | 2× and 1× coordinates compared without folding |

Two causes recur: **too few starting molecules relative to PCR cycles**, and **pool imbalance**.
Neither is fixable after sequencing, which is why both are worth checking early.

---

## Companion documents

- **demux runbook** — commands for stages 1 through 4
- **QC metrics** — per-library rates for concatemers, jackpotting and barcode hopping, with the
  gates that stop bad data propagating
