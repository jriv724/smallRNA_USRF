#!/usr/bin/env bash
set -euo pipefail

BASE="/samurlab1/Joshua/smallRNA_USRF"

mkdir -p \
  "$BASE/03.shortRNA_filtering/00.unmapped_reads" \
  "$BASE/03.shortRNA_filtering/01.databases" \
  "$BASE/03.shortRNA_filtering/02.indices" \
  "$BASE/03.shortRNA_filtering/03.filtered" \
  "$BASE/03.shortRNA_filtering/logs"

# -----------------------------
# 1) Copy non-miRNA residual reads per sample
# -----------------------------
for f in "$BASE"/02.residual_reads/*_mature_hairpin.unmapped.fastq.gz; do
  sample=$(basename "$f" _mature_hairpin.unmapped.fastq.gz)
  cp -v "$f" \
    "$BASE/03.shortRNA_filtering/00.unmapped_reads/${sample}.unmapped_reads.fastq.gz"
done

# -----------------------------
# 2) Download known RNA references
# -----------------------------
cd "$BASE/03.shortRNA_filtering/01.databases"

# tRNA: GtRNAdb hg38
wget -c http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa \
  -O hg38_gtrnadb_tRNAs.fa

# RNAcentral (large, catch-all)
wget -c ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz \
  -O rnacentral_active.fasta.gz

# extract human subset (quick filter)
zcat rnacentral_active.fasta.gz | awk '
  /^>/ {keep = ($0 ~ /9606|Homo sapiens/)}
  keep {print}
' > rnacentral_human.fa