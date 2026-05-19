#!/usr/bin/env bash
#SBATCH --job-name=smrnaseq_usrf
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/samurlab1/Joshua/smallRNA_USRF/logs/%x_%j.out
#SBATCH --error=/samurlab1/Joshua/smallRNA_USRF/logs/%x_%j.err

####
#### this script includes the corrected CATS 3' primer sequence and corrected min max length  8 -> 10 and 100 -> 20 for fastp. 
#### It also includes a check to make sure the R1 and R2 files match up in the sample sheet, and that all files exist. It pre-builds the conda environments with mamba to avoid issues with the nf-core pipeline trying to build them with micromamba. It also links the residual read files to a separate directory for easier access.
####

set -euo pipefail

PROJECT_DIR="/samurlab1/Joshua/smallRNA_USRF"
RAW_SHEET="${PROJECT_DIR}/sample_index.csv"
RUN_DIR="${PROJECT_DIR}/01.nfcore_smrnaseq"
WORK_DIR="${PROJECT_DIR}/work_smrnaseq"
LAUNCH_DIR="${PROJECT_DIR}/launch"
LOG_DIR="${PROJECT_DIR}/logs"
TMP_DIR="${PROJECT_DIR}/tmp"

GENOME="GRCh38"
MIRTRACE_SPECIES="hsa"
PIPELINE_VER="2.4.1"

MINLEN=10
MAXLEN=20
THREE_PRIME_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

CONDA_BASE="/samurlab1/Joshua/joshMiniforge3"
CONDA_ENV_NAME="nextflow"
CONDA_ENV_PREFIX="${CONDA_BASE}/envs/${CONDA_ENV_NAME}"

CONDA_PYTHON="${CONDA_ENV_PREFIX}/bin/python"
CONDA_NEXTFLOW="${CONDA_ENV_PREFIX}/bin/nextflow"
CONDA_JAVA="${CONDA_ENV_PREFIX}/bin/java"
CONDA_MAMBA="${CONDA_ENV_PREFIX}/bin/mamba"

export PATH="${CONDA_ENV_PREFIX}/bin:${PATH}"
export NXF_HOME="${PROJECT_DIR}/.nextflow"
export NXF_WORK="${WORK_DIR}"
export NXF_CONDA_CACHEDIR="${PROJECT_DIR}/.conda_cache"
export TMPDIR="${TMP_DIR}"
export TEMP="${TMP_DIR}"
export TMP="${TMP_DIR}"
export MAMBA_NO_BANNER=1
export PYTHONNOUSERSITE=1

PROFILE="mamba,illumina"

mkdir -p \
  "${PROJECT_DIR}" \
  "${RUN_DIR}" \
  "${WORK_DIR}" \
  "${LAUNCH_DIR}" \
  "${LOG_DIR}" \
  "${TMP_DIR}" \
  "${NXF_HOME}" \
  "${NXF_CONDA_CACHEDIR}"

echo "=== environment check ==="

[[ -x "${CONDA_NEXTFLOW}" ]] || { echo "ERROR: nextflow not found at ${CONDA_NEXTFLOW}"; exit 1; }
[[ -x "${CONDA_JAVA}" ]]     || { echo "ERROR: java not found at ${CONDA_JAVA}"; exit 1; }
[[ -x "${CONDA_PYTHON}" ]]   || { echo "ERROR: python not found at ${CONDA_PYTHON}"; exit 1; }
[[ -x "${CONDA_MAMBA}" ]]    || { echo "ERROR: mamba not found at ${CONDA_MAMBA}"; exit 1; }

echo "nextflow: ${CONDA_NEXTFLOW}"
echo "java:     ${CONDA_JAVA}"
echo "python:   ${CONDA_PYTHON}"
echo "mamba:    ${CONDA_MAMBA}"

"${CONDA_NEXTFLOW}" -version
"${CONDA_JAVA}" -version
"${CONDA_PYTHON}" --version
"${CONDA_MAMBA}" --version

"${CONDA_PYTHON}" - <<'PY'
mods = ["pandas"]
missing = []
for m in mods:
    try:
        __import__(m)
    except Exception:
        missing.append(m)
if missing:
    raise SystemExit(f"ERROR: missing Python modules in launcher env: {missing}")
print("Launcher Python dependencies look OK")
PY

echo "Using profile: ${PROFILE}"
echo "Using 3-prime adapter: ${THREE_PRIME_ADAPTER}"
echo "Using fastp length bounds: ${MINLEN}-${MAXLEN} nt"

SE_SHEET="${LAUNCH_DIR}/sample_index.smrnaseq.single_end.csv"
PAIR_AUDIT="${LAUNCH_DIR}/sample_index.pair_audit.tsv"

export RAW_SHEET
export SE_SHEET
export PAIR_AUDIT

"${CONDA_PYTHON}" - <<'PY'
import os
import pandas as pd

raw_sheet = os.environ["RAW_SHEET"]
se_sheet = os.environ["SE_SHEET"]
audit = os.environ["PAIR_AUDIT"]

df = pd.read_csv(raw_sheet)

required = {"sample", "fastq_1"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"Missing required columns: {sorted(missing)}")

if "fastq_2" not in df.columns:
    df["fastq_2"] = ""

rows = []
for _, r in df.iterrows():
    sample = str(r["sample"])
    r1 = str(r["fastq_1"]) if pd.notna(r["fastq_1"]) else ""
    r2 = str(r["fastq_2"]) if pd.notna(r["fastq_2"]) else ""
    rows.append({
        "sample": sample,
        "fastq_1": r1,
        "fastq_2": r2,
        "r1_exists": os.path.exists(r1),
        "r2_exists": os.path.exists(r2) if r2 else False,
        "r2_matches_r1_name": (
            os.path.basename(r2).replace("_R2", "_R1") == os.path.basename(r1)
        ) if r2 else False,
    })

audit_df = pd.DataFrame(rows)
audit_df.to_csv(audit, sep="\t", index=False)

if not audit_df["r1_exists"].all():
    bad = audit_df.loc[~audit_df["r1_exists"], ["sample", "fastq_1"]]
    raise SystemExit("One or more R1 files do not exist:\n" + bad.to_string(index=False))

se_df = df[["sample", "fastq_1"]].copy()
se_df.to_csv(se_sheet, index=False)

print(f"Wrote samplesheet: {se_sheet}")
print(f"Wrote pair audit:  {audit}")
PY

echo "=== pulling nf-core/smrnaseq assets ==="
"${CONDA_NEXTFLOW}" pull nf-core/smrnaseq -r "${PIPELINE_VER}"

rebuild_env_from_yml() {
    local env_name="$1"
    local env_yml="$2"
    local env_prefix="$3"
    local check_exe="$4"

    echo "=== checking ${env_name} environment file ==="
    [[ -f "${env_yml}" ]] || { echo "ERROR: ${env_name} environment file not found at ${env_yml}"; exit 1; }

    echo "${env_name} env file: ${env_yml}"
    echo "${env_name} env dir:  ${env_prefix}"

    if [[ -d "${env_prefix}" && ! -x "${env_prefix}/bin/${check_exe}" ]]; then
        echo "Removing partial ${env_name} env: ${env_prefix}"
        rm -rf "${env_prefix}"
    fi

    if [[ ! -x "${env_prefix}/bin/${check_exe}" ]]; then
        echo "=== prebuilding ${env_name} env with exact mamba binary ==="
        "${CONDA_MAMBA}" env create \
            --yes \
            --prefix "${env_prefix}" \
            --file "${env_yml}"
    else
        echo "${env_name} env already present."
    fi

    [[ -x "${env_prefix}/bin/${check_exe}" ]] || {
        echo "ERROR: ${env_name} executable ${check_exe} missing after env build"
        exit 1
    }
}

rebuild_env_from_spec() {
    local env_name="$1"
    local env_prefix="$2"
    local check_exe="$3"
    shift 3
    local specs=("$@")

    echo "=== checking ${env_name} env ==="
    echo "${env_name} env dir: ${env_prefix}"

    if [[ -d "${env_prefix}" && ! -x "${env_prefix}/bin/${check_exe}" ]]; then
        echo "Removing partial ${env_name} env: ${env_prefix}"
        rm -rf "${env_prefix}"
    fi

    if [[ ! -x "${env_prefix}/bin/${check_exe}" ]]; then
        echo "=== prebuilding ${env_name} env with exact mamba binary ==="
        "${CONDA_MAMBA}" create \
            --yes \
            --prefix "${env_prefix}" \
            "${specs[@]}"
    else
        echo "${env_name} env already present."
    fi

    [[ -x "${env_prefix}/bin/${check_exe}" ]] || {
        echo "ERROR: ${env_name} executable ${check_exe} missing after env build"
        exit 1
    }
}

FASTQC_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/nf-core/fastqc/environment.yml"
FASTQC_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-e88d3b603ce8089926a3a3fee5b86d02"
rebuild_env_from_yml "FastQC" "${FASTQC_ENV_YML}" "${FASTQC_ENV_PREFIX}" "fastqc"
"${FASTQC_ENV_PREFIX}/bin/fastqc" --version || true

FASTX_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-e3041154361d4096fb2d21a1924036f7"
rebuild_env_from_spec \
    "FASTX Toolkit" \
    "${FASTX_ENV_PREFIX}" \
    "fastx_collapser" \
    "bioconda::fastx_toolkit=0.0.14"

SEQCLUSTER_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-2ec791992f67c095c8990f681e82a3c9"
rebuild_env_from_spec \
    "seqcluster collapse" \
    "${SEQCLUSTER_ENV_PREFIX}" \
    "seqcluster" \
    "python=3.11" \
    "setuptools<81" \
    "bioconda::seqcluster=1.2.9"

SEQKIT_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-3d528e9269ca03d0ab9096382b8a8b88"
rebuild_env_from_spec \
    "seqkit fq2fa" \
    "${SEQKIT_ENV_PREFIX}" \
    "seqkit" \
    "bioconda::seqkit"

MIRDEEP_MAPPER_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/nf-core/mirdeep2/mapper/environment.yml"
MIRDEEP_MAPPER_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-0f4632abd08c7b22f16beba6b7b27ccb"
rebuild_env_from_yml \
    "mirdeep2 mapper" \
    "${MIRDEEP_MAPPER_ENV_YML}" \
    "${MIRDEEP_MAPPER_ENV_PREFIX}" \
    "mapper.pl"

SAMTOOLS_SORT_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/nf-core/samtools/sort/environment.yml"
SAMTOOLS_SORT_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-fa23c14bfe20d281e05f8173f5eedb7c"
rebuild_env_from_yml \
    "samtools sort" \
    "${SAMTOOLS_SORT_ENV_YML}" \
    "${SAMTOOLS_SORT_ENV_PREFIX}" \
    "samtools"
"${SAMTOOLS_SORT_ENV_PREFIX}/bin/samtools" --version | head -n 1 || true

MIRTOP_GFF_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-13b8b84647f7fbd2832bc10a1fda7eea"
rebuild_env_from_spec \
    "mirtop gff" \
    "${MIRTOP_GFF_ENV_PREFIX}" \
    "mirtop" \
    "python=3.11" \
    "setuptools<81" \
    "bioconda::mirtop=0.4.30" \
    "bioconda::samtools=1.21"

"${MIRTOP_GFF_ENV_PREFIX}/bin/mirtop" --version || true
"${MIRTOP_GFF_ENV_PREFIX}/bin/samtools" --version | head -n 1 || true

PIVOT_LONGER_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/local/pivot/longer/environment.yml"
PIVOT_LONGER_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-8431c62a66b876a27f12f64049a69ae4"
rebuild_env_from_yml \
    "pivot longer" \
    "${PIVOT_LONGER_ENV_YML}" \
    "${PIVOT_LONGER_ENV_PREFIX}" \
    "R"

PIVOT_WIDER_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/local/pivot/wider/environment.yml"
PIVOT_WIDER_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-3db4a9c331614be839da48ed4ceac3e8"
rebuild_env_from_yml \
    "pivot wider" \
    "${PIVOT_WIDER_ENV_YML}" \
    "${PIVOT_WIDER_ENV_PREFIX}" \
    "R"

EDGER_QC_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/local/edger_qc/environment.yml"
EDGER_QC_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-f0767c44ad966901614d792f2c3ac61d"
rebuild_env_from_yml \
    "edger qc" \
    "${EDGER_QC_ENV_YML}" \
    "${EDGER_QC_ENV_PREFIX}" \
    "R"

DATATABLE_MERGE_ENV_YML="${NXF_HOME}/assets/nf-core/smrnaseq/modules/local/datatable_merge/environment.yml"
DATATABLE_MERGE_ENV_PREFIX="${NXF_CONDA_CACHEDIR}/env-c7043b656647942706ba45e4c88a9dbc"
rebuild_env_from_yml \
    "datatable merge" \
    "${DATATABLE_MERGE_ENV_YML}" \
    "${DATATABLE_MERGE_ENV_PREFIX}" \
    "R"

cd "${LAUNCH_DIR}"

echo "=== launching pipeline ==="
set +e
"${CONDA_NEXTFLOW}" run nf-core/smrnaseq \
    -r "${PIPELINE_VER}" \
    -profile "${PROFILE}" \
    --input "${SE_SHEET}" \
    --genome "${GENOME}" \
    --mirtrace_species "${MIRTRACE_SPECIES}" \
    --outdir "${RUN_DIR}" \
    --save_intermediates \
    --save_trimmed_fail \
    --three_prime_adapter "${THREE_PRIME_ADAPTER}" \
    --fastp_min_length "${MINLEN}" \
    --fastp_max_length "${MAXLEN}" \
    -work-dir "${WORK_DIR}" \
    -resume
NF_EXIT=$?
set -e

RESID_DIR="${PROJECT_DIR}/02.residual_reads"
mkdir -p "${RESID_DIR}"

find "${RUN_DIR}" -type f \( \
    -name "*unmapped*.fastq.gz" -o \
    -name "*trimmed*.fastq.gz" -o \
    -path "*/fastp/*" -name "*.fastq.gz" \
\) | while read -r f; do
    ln -sf "$f" "${RESID_DIR}/$(basename "$f")"
done

echo "Done."
echo "Single-end samplesheet used: ${SE_SHEET}"
echo "Pair audit:                  ${PAIR_AUDIT}"
echo "nf-core output:              ${RUN_DIR}"
echo "Residual read links:         ${RESID_DIR}"
echo "Nextflow exit code:          ${NF_EXIT}"

exit "${NF_EXIT}"