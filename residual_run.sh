PROJECT_DIR="/samurlab1/Joshua/smallRNA_USRF"
RUN_DIR="${PROJECT_DIR}/01.nfcore_smrnaseq"
RESID_DIR="${PROJECT_DIR}/02.residual_reads"

mkdir -p "${RESID_DIR}"

find "${RUN_DIR}" -type f \( \
    -name "*unmapped*.fastq.gz" -o \
    -name "*trimmed*.fastq.gz" -o \
    -path "*/fastp/*" -name "*.fastq.gz" \
\) | while read -r f; do
    ln -sf "$f" "${RESID_DIR}/$(basename "$f")"
done