# RNase R processes ribosomal RNA precursors in *Haloferax volcanii*

Michael Jüttner, Felix Grünberger, Rosalie Kieffer, Viola Goebel, Dina Grohmann, Lisa-Katharina Maier, Sébastien Ferreira-Cerca and Anita Marchfelder

*Frontiers in Microbiology* · Manuscript ID 1973540

## About this repository

This repository contains the Nanopore PCR-cDNA processing workflow and R scripts used to analyse transcript abundance, RNA ends, degradation intermediates and rRNA maturation in wild-type and Δ*rnr* *Haloferax volcanii*.

The code covers the principal sequencing analyses and the generation of selected main-manuscript figures.

## Nanopore sequencing and mapping

PCR-cDNA libraries were prepared using the Oxford Nanopore SQK-PCB114.24 kit, with a custom 3′ adapter strategy to preserve transcript-end information, and sequenced on R10.4.1 flow cells.

The processing steps are collected in [`pipeline/nanopore_processing.sh`](pipeline/nanopore_processing.sh):

| Step | Processing |
|---|---|
| `basecall` | Dorado 0.8.0; SUP model `sup@v5.0.0` |
| `demux` | Dorado demultiplexing with `SQK-PCB114-24` |
| `trim` | Pychopper full-length read orientation with custom primers; Cutadapt linker and 5′-end trimming |
| `map` | Minimap2 alignment to the *H. volcanii* reference genome |
| `dedup` | UMI correction, PCR duplicate removal and UMI primer trimming |
| `remap` | Alignment of deduplicated reads to the reference genome |
| `circular` | Mapping to permuted 16S and 23S rRNA references |
| `coverage` | Strand-specific 5′ and 3′ end coverage using deepTools |
| `peaks` | Enriched end-position detection using `termseq_peaks` |

The pipeline uses Minimap2 `-ax map-ont --MD`, SAMtools for BAM processing, UMI-tools for PCR duplicate removal and deepTools `bamCoverage` with CPM normalization and single-base offsets. The `peaks` step operates on `barcode20`–`barcode23`, representing the rRNA-depleted wild-type and Δ*rnr* libraries used for transcript-end analyses.

### Running the workflow

Set the input directory and run each processing stage in sequence:

```bash
export DATA_DIR=/path/to/data
export POD5_DIR=/path/to/pod5
export PRIMER_FASTA=/path/to/pychopper_custom.fasta
export SAMPLE_INFO=pipeline/sample_info.csv
export THREADS=16
export DORADO_DEVICE=cuda:0,1

bash pipeline/nanopore_processing.sh basecall
bash pipeline/nanopore_processing.sh demux
bash pipeline/nanopore_processing.sh trim
bash pipeline/nanopore_processing.sh map
bash pipeline/nanopore_processing.sh dedup
bash pipeline/nanopore_processing.sh remap
bash pipeline/nanopore_processing.sh circular
bash pipeline/nanopore_processing.sh coverage
bash pipeline/nanopore_processing.sh peaks
```

The workflow expects `data/genome/hvo.fasta` and a custom Pychopper primer FASTA. For circular rRNA mapping, the reference names are read from `pipeline/sample_info.csv`; reference FASTA files are stored relative to `data/genome/`. The corresponding files (`circ/circ16S_hvo.fasta` and `circ/circ23S_hvo.fasta`) contain permuted rRNA templates. These references and the experimental input data are not included.

The end-position annotation and filtering used for the final analysis require the processed peak tables described below.

## Analysis scripts

| Figure | Script | Analysis |
|---|---|---|
| Figure 2 | [`Fig2_growth.R`](scripts/Fig2_growth.R) | Growth under different culture conditions |
| Figure 3 | [`Fig3_differential_expression.R`](scripts/Fig3_differential_expression.R) | Differential gene expression and arCOG annotation |
| Figure 4B–E | [`Fig4_transcript_ends.R`](scripts/Fig4_transcript_ends.R) | Transcript-end overlap and positional coverage profiles |
| Figure 5B–C | [`Fig5_transcript_classification.R`](scripts/Fig5_transcript_classification.R) | Full-length and processed mRNA classes |
| Figure 8B–C | [`Fig8_rRNA_precursors.R`](scripts/Fig8_rRNA_precursors.R) | rRNA precursor classification, differential abundance and RNA composition |

Common R imports and annotation functions are defined in [`scripts/helpers.R`](scripts/helpers.R). Scripts are run from the repository root and write figures to `figures/`.

The input directory contains genome annotations (`genome/`), sample metadata (`samples/`), processed count and peak tables (`Rdata/`), deduplicated BAM alignments and strand-specific coverage files (`bamCoverage/`). Growth measurements are read from `data/growth/`.

## Data availability

Nanopore PCR-cDNA sequencing data are available from the European Nucleotide Archive under accession [PRJEB87392](https://www.ebi.ac.uk/ena/browser/view/PRJEB87392).

## Code availability

https://github.com/felixgrunberger/RNaseR_Haloferax
