#!/usr/bin/env bash
set -euo pipefail

BASE="/samurlab1/Joshua/smallRNA_USRF"
FILT="$BASE/03.shortRNA_filtering"

IN="$FILT/00.unmapped_reads"
DB="$FILT/01.databases"
IDX="$FILT/02.indices"
OUT="$FILT/03.filtered"
LOG="$FILT/logs"
QC="$FILT/qc"

THREADS=8

mkdir -p "$IDX" "$OUT" "$LOG" "$QC"

command -v bowtie-build >/dev/null || { echo "ERROR: bowtie-build not found"; exit 1; }
command -v bowtie >/dev/null || { echo "ERROR: bowtie not found"; exit 1; }

mkdir -p "$IDX/tRNA" "$IDX/RNAcentral_human"

if [[ ! -s "$IDX/tRNA/tRNA.1.ebwt" ]]; then
  bowtie-build "$DB/hg38_gtrnadb_tRNAs.fa" "$IDX/tRNA/tRNA" \
    > "$LOG/build_tRNA.log" 2>&1
fi

if [[ ! -s "$IDX/RNAcentral_human/RNAcentral_human.1.ebwt" ]]; then
  bowtie-build "$DB/rnacentral_human.fa" "$IDX/RNAcentral_human/RNAcentral_human" \
    > "$LOG/build_RNAcentral_human.log" 2>&1
fi

echo -e "sample\tstage\treads" > "$QC/filter_counts.tsv"

for fq in "$IN"/*.unmapped_reads.fastq.gz; do
  sample=$(basename "$fq" .unmapped_reads.fastq.gz)
  echo "Processing $sample"

  SDIR="$OUT/$sample"
  mkdir -p "$SDIR"

  input_n=$(zcat "$fq" | awk 'END{print NR/4}')
  echo -e "$sample\tinput\t$input_n" >> "$QC/filter_counts.tsv"

  # -----------------------------
  # tRNA filter
  # -----------------------------
  zcat "$fq" | bowtie -q -v 0 -k 1 --best -p "$THREADS" \
    --un "$SDIR/${sample}.not_tRNA.fastq" \
    --al "$SDIR/${sample}.tRNA_matched.fastq" \
    "$IDX/tRNA/tRNA" \
    - \
    "$SDIR/${sample}.tRNA.sam" \
    > "$LOG/${sample}.tRNA.stdout.log" \
    2> "$LOG/${sample}.tRNA.stderr.log"

  gzip -f "$SDIR/${sample}.not_tRNA.fastq" "$SDIR/${sample}.tRNA_matched.fastq"

  not_tRNA_n=$(zcat "$SDIR/${sample}.not_tRNA.fastq.gz" | awk 'END{print NR/4}')
  tRNA_n=$(zcat "$SDIR/${sample}.tRNA_matched.fastq.gz" | awk 'END{print NR/4}')

  echo -e "$sample\ttRNA_matched\t$tRNA_n" >> "$QC/filter_counts.tsv"
  echo -e "$sample\tnot_tRNA\t$not_tRNA_n" >> "$QC/filter_counts.tsv"

  # -----------------------------
  # RNAcentral human filter
  # -----------------------------
  zcat "$SDIR/${sample}.not_tRNA.fastq.gz" | bowtie -q -v 0 -k 1 --best -p "$THREADS" \
    --un "$SDIR/${sample}.not_RNAcentral.fastq" \
    --al "$SDIR/${sample}.RNAcentral_matched.fastq" \
    "$IDX/RNAcentral_human/RNAcentral_human" \
    - \
    "$SDIR/${sample}.RNAcentral.sam" \
    > "$LOG/${sample}.RNAcentral.stdout.log" \
    2> "$LOG/${sample}.RNAcentral.stderr.log"

  gzip -f "$SDIR/${sample}.not_RNAcentral.fastq" "$SDIR/${sample}.RNAcentral_matched.fastq"

  final_n=$(zcat "$SDIR/${sample}.not_RNAcentral.fastq.gz" | awk 'END{print NR/4}')
  rnacentral_n=$(zcat "$SDIR/${sample}.RNAcentral_matched.fastq.gz" | awk 'END{print NR/4}')

  echo -e "$sample\tRNAcentral_matched\t$rnacentral_n" >> "$QC/filter_counts.tsv"
  echo -e "$sample\tcandidate_remaining\t$final_n" >> "$QC/filter_counts.tsv"

  ln -sf "$SDIR/${sample}.not_RNAcentral.fastq.gz" \
    "$SDIR/${sample}.candidate_shortRNA.fastq.gz"
done

echo "Done."
echo "Counts: $QC/filter_counts.tsv"
echo "Filtered outputs: $OUT"