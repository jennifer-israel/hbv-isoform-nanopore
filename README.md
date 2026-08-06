# hbv-isoform-nanopore
HBV transcript isoform detection and quantification from long-read nanopore RNA-seq

Analysis pipeline for detecting and quantifying HBV transcript isoforms from long-read nanopore sequencing of cDNA libraries. Covers computational demultiplexing of custom PCR barcodes, alignment to a doubled HBV reference so genome-wrapping transcripts align contiguously, UMI-based deduplication, and transcript classification by transcription start site and splice structure. Applies to both hybridisation-enriched and un-enriched libraries; enrichment status is recorded per experiment in config/.

Derived from the EXP26000559 (cDNA001) pipeline by Matt Wolpert; see PROVENANCE.md for what is upstream, what was adapted, and what is new.
