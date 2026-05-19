````markdown id="9t6x1m"
# Nextflow + nf-core/smrnaseq Pipeline for Ultra-Short RNA Fragment (USRF) Analysis

A scalable small RNA-seq workflow and downstream analysis framework for detecting, filtering, and characterizing ultra-short RNA fragments (USRFs) using Nextflow and nf-core/smrnaseq.

---

## Overview

This repository contains workflow scripts, QC analyses, filtering utilities, and exploratory notebooks for small RNA sequencing datasets processed with a modified nf-core/smrnaseq pipeline designed to preserve ultra-short reads.

While conventional small RNA-seq pipelines are primarily optimized for canonical miRNAs, this project focuses on the residual short RNA population that remains after iterative filtering of known small RNA classes.

The workflow was developed in an HPC environment using SLURM, Nextflow, and Conda/Mamba-based execution.

---

## Goals

The primary goals of this project are:

- Preserve and analyze reads in the ~8–20 nt range
- Characterize residual RNA populations after canonical filtering
- Explore reproducibility of ultra-short fragments across immune cell types
- Develop modular filtering workflows for iterative residual analysis
- Build reproducible Nextflow-based small RNA-seq infrastructure

---

## Core Technologies

### Workflow Infrastructure

- Nextflow
- nf-core/smrnaseq
- SLURM
- Conda / Mamba

### Bioinformatics Tools

- fastp
- FastQC
- MultiQC
- mirtrace
- Bowtie / Bowtie2
- samtools

### Downstream Analysis

- Python
- Jupyter notebooks
- pandas
- matplotlib

---

## Workflow Summary

```text id="0qyn5o"
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

## Repository Structure

```text id="72rjfa"
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
│   └── workflow launch metadata
│
├── logs/
│   └── SLURM execution logs
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

## Key Features

### Ultra-short read preservation

The workflow intentionally preserves reads below conventional miRNA size thresholds to enable exploratory analysis of ultra-short RNA species.

Custom pipeline configurations retain reads as short as approximately 8 nucleotides.

---

### Residual read extraction

The pipeline isolates reads remaining after:

* mature miRNA alignment
* hairpin alignment
* genome alignment

Residual populations are then subjected to iterative filtering against additional RNA references.

---

### Reproducible HPC execution

The workflow was developed for scalable execution on SLURM-based HPC infrastructure using:

* Nextflow orchestration
* Conda/Mamba environments
* nf-core modular workflows

Execution metadata and launch scripts are preserved for reproducibility.

---

### Integrated QC reporting

Included QC layers include:

* fastp reports
* FastQC summaries
* MultiQC aggregation
* mirtrace outputs
* alignment statistics
* residual filtering summaries

---

## Example Workflow Execution

### nf-core/smrnaseq

```bash id="oz9nqx"
nextflow run nf-core/smrnaseq \
    -profile singularity \
    --input sample_index.csv \
    --protocol illumina \
    --aligner bowtie \
    --genome GRCh38 \
    --save_trimmed \
    --minlength 8
```

### Custom USRF workflow

```bash id="c6e6uv"
bash run_smrnaseq_keep8.sh
```

---

## Included vs Excluded Data

This repository includes:

* workflow scripts
* launch files
* QC summaries
* filtering utilities
* exploratory notebooks
* execution metadata

This repository excludes:

* raw FASTQ files
* large intermediate outputs
* Nextflow work directories
* cached Conda environments
* genome indices and large reference databases

---

## Potential Applications

Potential downstream applications include:

* non-canonical short RNA discovery
* residual RNA profiling
* immune-cell-specific short RNA analysis
* degradation-pattern analysis
* exploratory biomarker discovery
* alternative small RNA filtering strategies

---

## Notes

This repository represents an active exploratory workflow framework for small RNA-seq analysis and ultra-short RNA fragment characterization.

The emphasis is on reproducible workflow engineering, iterative residual filtering, and scalable exploratory analysis rather than finalized biological conclusions.

---

## Author

Joshua Rivera
Dana-Farber Cancer Institute / Harvard Medical School

Computational biology, bioinformatics workflow engineering, and single-cell/multi-omic analysis.

```
```
