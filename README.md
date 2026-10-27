# RNase R processes ribosomal RNA precursors in *Haloferax volcanii*

**Code accompanying the manuscript** by Michael Jüttner, Felix Grünberger, Rosalie Kieffer, Viola Goebel, Dina Grohmann, Lisa-Katharina Maier, Sébastien Ferreira-Cerca and Anita Marchfelder.  
*Frontiers in Microbiology* (submitted 2026; manuscript ID **1973540**). Publication details and citation should be updated once the article is published.

This repository contains a **minimal selection of main-figure computational analyses** (not all exploratory analyses or supplementary figures). The scripts were reduced from the authors' original R code to focus on main experimental and long-read sequencing results.

## Analysis scripts

| Manuscript figure | Script | Output |
|---|---|---|
| Fig. 2 | `scripts/Fig2_growth.R` | Growth-curve summary |
| Fig. 3 | `scripts/Fig3_differential_expression.R` | DESeq2 differential expression and arCOG-annotated MA plot |
| Fig. 4B–E | `scripts/Fig4_transcript_ends.R` | Primary 5′/3′ end overlaps and terminal read distributions |
| Fig. 5B–C | `scripts/Fig5_transcript_classification.R` | Full-length, 5′-/3′-processed and degraded read classifications |
| Fig. 8B–C | `scripts/Fig8_rRNA_precursors.R` | 16S/23S intermediate classification, DESeq2 and RNA-type proportions |

Figure 1 (protein-domain/structural comparison), Figures 6 and 10 (models), and Northern blot Figures 7 and 9 are **not produced by these R scripts**; their illustrations/images are not included. The source-generated supplementary analyses were excluded intentionally.

## Data

PCR-cDNA sequencing data: **ENA [PRJEB87392](https://www.ebi.ac.uk/ena/browser/view/PRJEB87392)**. The code assumes processed reference annotations, BAM alignments, peak/coverage tables, sample metadata, and growth measurements are supplied under `data/`. These files are **not included in the repository**.

Required processed input files include:

- `data/genome/hvo.gff`, `data/genome/hvo.fasta`, `data/genome/hvo_annotation_arcog.tsv`
- `data/samples/ont_run_info.xlsx`, `data/growth/Neue Wachstumskurven.xlsx`
- `data/Rdata/tss_peaks_all_longread.tsv`, `data/Rdata/tts_peaks_all_longread_cov5.tsv`
- `data/Rdata/total_frame_anno_g.tsv`, `data/Rdata/counts_all_tab.tsv`
- `data/mapped_fl_dedup/` (UMI-deduplicated BAMs for Fig. 3)
- `data/mapped_fl_only_umi_removed/`, `data/mapped_fl_only_umi_removed_circ16S/`, `data/mapped_fl_only_umi_removed_circ23S/` (BAMs for Fig. 8)
- `data/bamCoverage/` (strand-specific end5/end3 CPM bedGraph files for Fig. 4)

The *Methods* section of the manuscript describes upstream read processing (Dorado, pychopper, Cutadapt, Minimap2, SAMtools, UMI-tools and end-position calling). **This code package starts from the mapped/processed data** and is not a complete raw-POD5-to-figure workflow.

## Running the scripts

Run each R script from the repository root in an R session with the necessary CRAN/Bioconductor dependencies installed. Output PDFs are written to `figures/`. The shared R functions and imports are in `scripts/helpers.R`. Scripts may require substantial memory to process whole-genome coverage tracks.

> **Author-review status (October 2026):** the scripts were trimmed without access to input datasets, and the revised versions have not been rerun or validated against the final assembled manuscript figures. Fig. 4E in particular uses an `intragenic` 3′-end class in its source code although the caption refers to primary termini. This must be checked before publication. This README does not claim full computational reproducibility without the missing inputs.

## Code availability

Repository: `https://github.com/felixgrunberger/RNaseR_Haloferax` (the URL specified in the manuscript’s Code availability statement). The present code is a minimal draft for author review; results must be checked against the final figures before the repository is made public.
