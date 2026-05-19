# Nextflow + nf-core/smrnaseq Workflow for Ultra-Short RNA Fragment (USRF) Discovery and Small RNA-seq Analysis

Pipeline development and downstream analysis framework for small RNA sequencing with emphasis on ultra-short RNA fragment (USRF) detection, residual read characterization, and iterative filtering of canonical small RNA species.

---

# Overview

This repository contains workflow scripts, quality-control analysis, filtering utilities, and downstream exploratory analyses for small RNA-seq datasets processed using:

- Nextflow
- nf-core/smrnaseq
- Bowtie / Bowtie2
- FastQC / MultiQC
- fastp
- mirtrace
- samtools

The project focuses specifically on preserving and interrogating ultra-short RNA fragments (approximately 8–20 nt) that are typically discarded in conventional small RNA-seq preprocessing pipelines.

The workflow was developed on SLURM-based HPC infrastructure and emphasizes reproducibility, modular filtering, and scalable execution across immune-cell sequencing datasets.

---

# Scientific Motivation

Conventional small RNA-seq workflows are generally optimized for canonical miRNAs and related annotated small RNA species. Reads shorter than mature miRNAs are often aggressively filtered during preprocessing or ignored during downstream analysis.

This project explores whether residual ultra-short sequencing reads that remain after canonical filtering steps may contain reproducible biological signal.

Core questions include:

- Which RNA fragments persist after mature miRNA and hairpin filtering?
- Are residual ultra-short reads reproducible across samples and immune cell types?
- Do specific sequence classes remain enriched after iterative filtering?
- Can these residual reads represent biologically meaningful RNA species rather than sequencing noise or degradation artifacts?

---

# Workflow Architecture

The pipeline combines standard nf-core/smrnaseq processing with custom residual-read extraction and downstream filtering steps.

## High-level workflow

```text
Raw small RNA FASTQ
        ↓
fastp preprocessing
        ↓
FastQC / MultiQC QC
        ↓
nf-core/smrnaseq alignment
        ↓
miRNA + hairpin quantification
        ↓
Residual unmapped read extraction
        ↓
Iterative filtering against:
    - mature miRNA
    - hairpin miRNA
    - genome
    - tRNA databases
    - additional RNA references
        ↓
Ultra-short RNA fragment analysis
````

---

# Repository Structure

```text
.
├── 01.nfcore_smrnaseq/
│   ├── fastp/
│   ├── fastqc/
│   ├── genome_quant/
│   ├── mirna_quant/
│   ├── mirtrace/
│   ├── multiqc/
│   └── pipeline_info/
│
├── 03.shortRNA_filtering/
│   ├── logs/
│   ├── qc/
│   └── filtering utilities
│
├── launch/
│   ├── sample sheets
│   └── launch metadata
│
├── logs/
│   └── SLURM workflow logs
│
├── *.ipynb
│   └── downstream QC and exploratory analyses
│
├── run_smrnaseq_keep8.sh
├── run_smrnaseq_kraken.sh
├── residual_run.sh
└── filter_mature_hairpin_filteredResiduals.sh
```

---

# Key Features

## Ultra-short read preservation

Unlike standard small RNA-seq preprocessing workflows, this pipeline intentionally preserves reads below typical miRNA size thresholds.

Custom workflow variants retain reads as short as approximately 8 nt for downstream investigation.

---

## Residual read analysis

The workflow extracts reads remaining after:

* mature miRNA alignment
* hairpin alignment
* genome alignment

Residual fractions are then iteratively filtered against additional RNA databases to isolate candidate USRF populations.

---

## HPC-oriented workflow design

The pipeline was developed and tested in a SLURM HPC environment using:

* Nextflow
* Conda/Mamba environments
* nf-core modular workflows

The repository includes workflow launch scripts and execution logs demonstrating reproducible execution on cluster infrastructure.

---

## Comprehensive QC reporting

The project incorporates multiple QC layers including:

* FastQC
* MultiQC
* fastp reports
* mirtrace summaries
* alignment statistics
* residual filtering summaries

Notebook-based exploratory QC analyses are included for visualization and downstream interpretation.

---

# Included Data

This repository intentionally excludes:

* raw FASTQ sequencing files
* large intermediate workflow outputs
* Nextflow work directories
* cached Conda environments
* large reference databases and genome indices

The repository instead preserves:

* workflow logic
* launch scripts
* QC summaries
* filtering utilities
* downstream analyses
* reproducible execution metadata

---

# Example Workflow Execution

## Example nf-core/smrnaseq execution

```bash
nextflow run nf-core/smrnaseq \
    -profile singularity \
    --input sample_index.csv \
    --protocol illumina \
    --aligner bowtie \
    --genome GRCh38 \
    --save_trimmed \
    --minlength 8
```

## Example custom launch

```bash
bash run_smrnaseq_keep8.sh
```

---

# Technologies

## Workflow / Infrastructure

* Nextflow
* nf-core/smrnaseq
* SLURM
* Conda / Mamba

## Bioinformatics Tools

* Bowtie
* Bowtie2
* samtools
* fastp
* FastQC
* MultiQC
* mirtrace

## Downstream Analysis

* Python
* Jupyter
* pandas
* matplotlib

---

# Potential Applications

Potential downstream applications of the framework include:

* discovery of non-canonical short RNA species
* characterization of degradation-derived RNA populations
* immune-cell-specific short RNA profiling
* residual small RNA analysis
* exploratory biomarker discovery
* development of alternative small RNA filtering paradigms

---

# Notes

This repository is intended primarily as:

1. A reproducible workflow framework for small RNA-seq processing
2. A development environment for USRF analysis methodology
3. A reference implementation for preserving and analyzing ultra-short sequencing reads

The repository is under active development and exploratory analysis strategies may evolve over time.

---

# Author

Joshua Rivera
Dana-Farber Cancer Institute / Harvard Medical School

Computational biology, single-cell analysis, and scalable bioinformatics workflow development.

```
```
