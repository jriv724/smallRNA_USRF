#!/usr/bin/env bash
#SBATCH --job-name=smrnaseq_usrf_patched
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/samurlab1/Joshua/smallRNA_USRF/logs_isolated/%x_%j.out
#SBATCH --error=/samurlab1/Joshua/smallRNA_USRF/logs_isolated/%x_%j.err

set -euo pipefail

#############################################
# Project paths
#############################################
PROJECT_DIR="/samurlab1/Joshua/smallRNA_USRF"
RAW_SHEET="${PROJECT_DIR}/sample_index.csv"

RUN_TAG="${RUN_TAG:-smrnaseq_$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${PROJECT_DIR}/isolated_runs/${RUN_TAG}"

RUN_DIR="${RUN_ROOT}/01.nfcore_smrnaseq"
WORK_DIR="${RUN_ROOT}/work_smrnaseq"
LAUNCH_DIR="${RUN_ROOT}/launch"
LOG_DIR="${RUN_ROOT}/logs"
TMP_DIR="${RUN_ROOT}/tmp"

GENOME="GRCh38"
MIRTRACE_SPECIES="hsa"
PIPELINE_VER="2.4.1"
MINLEN=10
MAXLEN=20

#############################################
# Local pipeline clone + prebuilt envs
#############################################
PIPELINE_DIR="${PROJECT_DIR}/pipelines/nf-core-smrnaseq-${PIPELINE_VER}"
ENV_ROOT="${PROJECT_DIR}/prebuilt_envs/nf-core-smrnaseq-${PIPELINE_VER}"
PATCH_MANIFEST="${LAUNCH_DIR}/conda_patch_manifest.tsv"
PATCHED_CONFIG="${LAUNCH_DIR}/patched_conda.config"
PATCH_REPORT="${LAUNCH_DIR}/conda_patch_report.txt"

#############################################
# Exact known-good tool paths
#############################################
CONDA_BASE="/samurlab1/Joshua/joshMiniforge3"
CONDA_ENV_NAME="nextflow"
CONDA_ENV_PREFIX="${CONDA_BASE}/envs/${CONDA_ENV_NAME}"

CONDA_PYTHON="${CONDA_ENV_PREFIX}/bin/python"
CONDA_NEXTFLOW="${CONDA_ENV_PREFIX}/bin/nextflow"
CONDA_JAVA="${CONDA_ENV_PREFIX}/bin/java"
CONDA_MAMBA="${CONDA_ENV_PREFIX}/bin/mamba"
CONDA_GIT="$(command -v git || true)"

#############################################
# Nextflow / conda caches
#############################################
export PATH="${CONDA_ENV_PREFIX}/bin:${PATH}"
export NXF_HOME="${PROJECT_DIR}/.nextflow"
export NXF_WORK="${WORK_DIR}"
export NXF_CONDA_CACHEDIR="${PROJECT_DIR}/.conda_cache"
export TMPDIR="${TMP_DIR}"
export TEMP="${TMP_DIR}"
export TMP="${TMP_DIR}"
export MAMBA_NO_BANNER=1
export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX="${PROJECT_DIR}/.mamba_root"
export CONDA_PKGS_DIRS="${PROJECT_DIR}/.conda_pkgs"
export NXF_ANSI_LOG=false

PROFILE="illumina"

mkdir -p \
  "${PROJECT_DIR}" \
  "${RUN_DIR}" \
  "${WORK_DIR}" \
  "${LAUNCH_DIR}" \
  "${LOG_DIR}" \
  "${TMP_DIR}" \
  "${NXF_HOME}" \
  "${NXF_CONDA_CACHEDIR}" \
  "${ENV_ROOT}" \
  "${MAMBA_ROOT_PREFIX}" \
  "${CONDA_PKGS_DIRS}" \
  "${PROJECT_DIR}/pipelines"

#############################################
# Environment sanity checks
#############################################
echo "=== environment check ==="

[[ -x "${CONDA_NEXTFLOW}" ]] || { echo "ERROR: nextflow not found at ${CONDA_NEXTFLOW}"; exit 1; }
[[ -x "${CONDA_JAVA}" ]]     || { echo "ERROR: java not found at ${CONDA_JAVA}"; exit 1; }
[[ -x "${CONDA_PYTHON}" ]]   || { echo "ERROR: python not found at ${CONDA_PYTHON}"; exit 1; }
[[ -x "${CONDA_MAMBA}" ]]    || { echo "ERROR: mamba not found at ${CONDA_MAMBA}"; exit 1; }
[[ -n "${CONDA_GIT}" ]]      || { echo "ERROR: git not found in PATH"; exit 1; }

echo "nextflow: ${CONDA_NEXTFLOW}"
echo "java:     ${CONDA_JAVA}"
echo "python:   ${CONDA_PYTHON}"
echo "mamba:    ${CONDA_MAMBA}"
echo "git:      ${CONDA_GIT}"

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
echo "Isolated run root: ${RUN_ROOT}"

action_log() {
  echo
  echo "=== $* ==="
}

#############################################
# Build R1-only samplesheet
#############################################
action_log "building R1-only samplesheet"
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

#############################################
# Clone exact pipeline version locally
#############################################
action_log "cloning nf-core/smrnaseq locally"
if [[ ! -d "${PIPELINE_DIR}/.git" ]]; then
  "${CONDA_GIT}" clone --branch "${PIPELINE_VER}" --depth 1 https://github.com/nf-core/smrnaseq.git "${PIPELINE_DIR}"
else
  echo "Pipeline clone already exists: ${PIPELINE_DIR}"
  (
    cd "${PIPELINE_DIR}"
    "${CONDA_GIT}" fetch --tags --depth 1 origin "refs/tags/${PIPELINE_VER}:refs/tags/${PIPELINE_VER}" || true
    "${CONDA_GIT}" checkout -f "${PIPELINE_VER}"
  )
fi

#############################################
# Patch conda directives -> existing env dirs
#############################################
action_log "inventorying and patching conda directives"
export PIPELINE_DIR ENV_ROOT PATCH_MANIFEST PATCHED_CONFIG PATCH_REPORT

"${CONDA_PYTHON}" - <<'PY'
import csv
import hashlib
import os
import re
from pathlib import Path

pipeline_dir = Path(os.environ["PIPELINE_DIR"]).resolve()
env_root = Path(os.environ["ENV_ROOT"]).resolve()
manifest_path = Path(os.environ["PATCH_MANIFEST"]).resolve()
config_path = Path(os.environ["PATCHED_CONFIG"]).resolve()
report_path = Path(os.environ["PATCH_REPORT"]).resolve()

env_root.mkdir(parents=True, exist_ok=True)

# Match: optional whitespace + conda + quoted string, single-line only.
conda_re = re.compile(r'^(?P<indent>\s*)conda\s+(?P<quote>["\'])(?P<value>.*?)(?P=quote)\s*$')

def sanitize_rel_path(p: Path) -> str:
    rel = p.relative_to(pipeline_dir)
    return str(rel).replace(os.sep, "__").replace(".yml", "").replace(".yaml", "")

def spec_env_name(spec: str) -> str:
    h = hashlib.sha1(spec.encode()).hexdigest()[:16]
    return f"spec__{h}"

def resolve_env_target(nf_file: Path, value: str):
    raw = value.strip()
    if "environment.yml" in raw or "environment.yaml" in raw:
        if raw.startswith("${moduleDir}/"):
            rel = raw.split("${moduleDir}/", 1)[1]
            yml_path = (nf_file.parent / rel).resolve()
        elif raw.startswith("./") or raw.startswith("../") or raw.endswith(".yml") or raw.endswith(".yaml"):
            yml_path = (nf_file.parent / raw).resolve()
        else:
            yml_path = (pipeline_dir / raw).resolve()
        if not yml_path.exists():
            raise FileNotFoundError(f"Could not resolve environment file '{raw}' from {nf_file}")
        env_prefix = env_root / sanitize_rel_path(yml_path)
        return {
            "kind": "yaml",
            "source": str(yml_path),
            "env_prefix": str(env_prefix),
            "conda_value": str(env_prefix),
        }
    else:
        env_prefix = env_root / spec_env_name(raw)
        return {
            "kind": "spec",
            "source": raw,
            "env_prefix": str(env_prefix),
            "conda_value": str(env_prefix),
        }

manifest_rows = []
report_lines = []
seen_keys = set()
patched_files = 0
patched_directives = 0

for nf_file in sorted(pipeline_dir.rglob("*.nf")):
    text = nf_file.read_text()
    new_lines = []
    changed = False
    for line in text.splitlines(True):
        line_no_nl = line.rstrip("\n")
        m = conda_re.match(line_no_nl)
        if not m:
            new_lines.append(line)
            continue

        target = resolve_env_target(nf_file, m.group("value"))
        indent = m.group("indent")
        quote = m.group("quote")
        replacement = f"{indent}conda {quote}{target['conda_value']}{quote}"
        if line.endswith("\n"):
            replacement += "\n"

        new_lines.append(replacement)
        changed = True
        patched_directives += 1
        key = (target["kind"], target["source"], target["env_prefix"])
        if key not in seen_keys:
            seen_keys.add(key)
            manifest_rows.append(target)
        report_lines.append(f"PATCHED\t{nf_file.relative_to(pipeline_dir)}\t{m.group('value')}\t{target['env_prefix']}")

    if changed:
        nf_file.write_text("".join(new_lines))
        patched_files += 1

# Also scan configs for any visible conda directives and report them.
for cfg_file in sorted(pipeline_dir.rglob("*.config")):
    for idx, line in enumerate(cfg_file.read_text().splitlines(), start=1):
        if "conda" in line:
            report_lines.append(f"CONFIG_SCAN\t{cfg_file.relative_to(pipeline_dir)}:{idx}\t{line}")

with manifest_path.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["kind", "source", "env_prefix", "conda_value"], delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

config_path.write_text(
    """
conda {
  enabled = true
  useMamba = true
  cacheDir = '/samurlab1/Joshua/smallRNA_USRF/.conda_cache'
  createTimeout = '2h'
}

docker.enabled = false
singularity.enabled = false
apptainer.enabled = false
podman.enabled = false
charliecloud.enabled = false
""".strip() + "\n"
)

report_header = [
    f"Patched files: {patched_files}",
    f"Patched conda directives: {patched_directives}",
    f"Unique env targets: {len(manifest_rows)}",
    "",
]
report_path.write_text("\n".join(report_header + report_lines) + "\n")

print(f"Patched files: {patched_files}")
print(f"Patched conda directives: {patched_directives}")
print(f"Unique env targets: {len(manifest_rows)}")
print(f"Manifest: {manifest_path}")
print(f"Config:   {config_path}")
print(f"Report:   {report_path}")
PY

#############################################
# Build all prebuilt envs from patched manifest
#############################################
action_log "prebuilding all required conda envs"

build_yaml_env() {
  local yml="$1"
  local prefix="$2"

  if [[ -x "${prefix}/bin/python" || -x "${prefix}/bin/R" || -x "${prefix}/bin/bash" ]]; then
    echo "exists: ${prefix}"
    return 0
  fi

  rm -rf "${prefix}"
  mkdir -p "$(dirname "${prefix}")"

  echo "building YAML env -> ${prefix}"
  "${CONDA_MAMBA}" env create --yes --prefix "${prefix}" --file "${yml}"
}

build_spec_env() {
  local spec="$1"
  local prefix="$2"

  if [[ -x "${prefix}/bin/python" || -x "${prefix}/bin/R" || -x "${prefix}/bin/bash" ]]; then
    echo "exists: ${prefix}"
    return 0
  fi

  rm -rf "${prefix}"
  mkdir -p "$(dirname "${prefix}")"

  echo "building SPEC env -> ${prefix}"
  read -r -a spec_parts <<< "${spec}"
  "${CONDA_MAMBA}" create --yes --prefix "${prefix}" "${spec_parts[@]}"
}

while IFS=$'\t' read -r kind source env_prefix conda_value; do
  [[ "${kind}" == "kind" ]] && continue
  if [[ "${kind}" == "yaml" ]]; then
    build_yaml_env "${source}" "${env_prefix}"
  elif [[ "${kind}" == "spec" ]]; then
    build_spec_env "${source}" "${env_prefix}"
  else
    echo "ERROR: unknown manifest kind '${kind}'"
    exit 1
  fi
done < "${PATCH_MANIFEST}"

#############################################
# Quick validation of prebuilt env count
#############################################
action_log "validating prebuilt env inventory"
REQ_COUNT=$(tail -n +2 "${PATCH_MANIFEST}" | wc -l | awk '{print $1}')
HAVE_COUNT=$(find "${ENV_ROOT}" -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1}')
echo "Required env targets: ${REQ_COUNT}"
echo "Present env dirs:     ${HAVE_COUNT}"
[[ "${HAVE_COUNT}" -ge 1 ]] || { echo "ERROR: no prebuilt envs found"; exit 1; }

#############################################
# Run patched local pipeline
#############################################
action_log "launching patched local pipeline"
cd "${LAUNCH_DIR}"

"${CONDA_NEXTFLOW}" run "${PIPELINE_DIR}" \
  -r "${PIPELINE_VER}" \
  -profile "${PROFILE}" \
  -c "${PATCHED_CONFIG}" \
  --input "${SE_SHEET}" \
  --genome "${GENOME}" \
  --mirtrace_species "${MIRTRACE_SPECIES}" \
  --outdir "${RUN_DIR}" \
  --save_intermediates \
  --save_trimmed_fail \
  --fastp_min_length "${MINLEN}" \
  --fastp_max_length "${MAXLEN}" \
  -work-dir "${WORK_DIR}" \
  -resume

#############################################
# Collect residual candidate FASTQs
#############################################
action_log "collecting residual candidate FASTQs"
RESID_DIR="${RUN_ROOT}/02.residual_reads"
mkdir -p "${RESID_DIR}"

find "${RUN_DIR}" -type f \( \
    -name "*unmapped*.fastq.gz" -o \
    -name "*trimmed*.fastq.gz" -o \
    -path "*/fastp/*" -name "*.fastq.gz" \
\) | while read -r f; do
    ln -sf "$f" "${RESID_DIR}/$(basename "$f")"
done

#############################################
# Summary
#############################################
action_log "done"
echo "Single-end samplesheet used: ${SE_SHEET}"
echo "Pair audit:                  ${PAIR_AUDIT}"
echo "Patched manifest:            ${PATCH_MANIFEST}"
echo "Patched conda config:        ${PATCHED_CONFIG}"
echo "Patch report:                ${PATCH_REPORT}"
echo "Pipeline clone:              ${PIPELINE_DIR}"
echo "Prebuilt env root:           ${ENV_ROOT}"
echo "nf-core output:              ${RUN_DIR}"
echo "Residual read links:         ${RESID_DIR}"
