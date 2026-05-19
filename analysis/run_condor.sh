#!/bin/bash
set -euo pipefail

DATASET="${1:?Missing dataset argument}"

WORKDIR="${_CONDOR_SCRATCH_DIR:-$PWD}"
cd "${WORKDIR}"

export HOME="${WORKDIR}"
export XDG_CACHE_HOME="${WORKDIR}"
export XDG_CONFIG_HOME="${WORKDIR}"
export MPLCONFIGDIR="${WORKDIR}"
export NUMBA_CACHE_DIR="${WORKDIR}"
export PIP_CACHE_DIR="${WORKDIR}"

echo "=== Job started on $(hostname) ==="
date
echo "WORKDIR=${WORKDIR}"
echo "DATASET=${DATASET}"

export X509_USER_PROXY="${WORKDIR}/x509up_u147757"

export XRD_NETWORKSTACK=IPv4
export XRD_REQUESTTIMEOUT=120
export XRD_REDIRECTLIMIT=10

# --- conda env ---
tar -xzf py38.tgz
source bin/activate
conda-unpack

# --- analysis code ---
tar -xzf analysis.tgz

# 여기 중요 👇
cd analysis
set -x
python3 run.py \
  -p stop_new2024 \
  -m KNU_2024_v4 \
  -w 1 \
  -d "${DATASET}"
set +x

echo "=== After run.py: pwd=$(pwd) ==="
echo "=== List hists/stop_new2024 ==="
ls -lah hists/stop_new2024 || true

echo "=== Find futures files ==="
find . -maxdepth 6 -name "*.futures" -ls || true

# futures 위치도 경로 맞춰야 함
FUT="hists/stop_new2024/${DATASET}.futures"
if [ ! -f "${FUT}" ]; then
  echo "ERROR: expected futures not found: ${FUT}" >&2
  echo "Dumping tree (maxdepth 4)..." >&2
  find . -maxdepth 4 -type d -print >&2
  exit 2
fi

cp -v "${FUT}" "${WORKDIR}/${DATASET}.futures"

echo "=== Job finished ==="
date