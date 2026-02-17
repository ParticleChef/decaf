#!/bin/bash
set -euo pipefail

DATASET="${1:?Missing dataset argument}"

WORKDIR="${_CONDOR_SCRATCH_DIR:-$PWD}"
cd "${WORKDIR}"

# ----- stable scratch env -----
export HOME="${WORKDIR}"
export XDG_CACHE_HOME="${WORKDIR}"
export XDG_CONFIG_HOME="${WORKDIR}"
export MPLCONFIGDIR="${WORKDIR}"
export NUMBA_CACHE_DIR="${WORKDIR}"
export PIP_CACHE_DIR="${WORKDIR}"
export PYTHONUNBUFFERED=1

# ----- proxy shipped via transfer_input_files -----
export X509_USER_PROXY="${WORKDIR}/x509up_u147757"

# ----- XRootD -----
export XRD_NETWORKSTACK=IPv4
export XRD_REQUESTTIMEOUT=120
export XRD_REDIRECTLIMIT=10

echo "=== START $(date) ==="
echo "HOST=$(hostname)"
echo "WORKDIR=${WORKDIR}"
echo "DATASET=${DATASET}"
echo "ls(workdir) before unpack:"
ls -lah

# ----- unpack conda env locally (no 'source activate') -----
echo "=== Unpack py38.tgz ==="
tar -xzf py38.tgz

export PATH="${WORKDIR}/bin:${PATH}"
export LD_LIBRARY_PATH="${WORKDIR}/lib:${LD_LIBRARY_PATH:-}"

PY="${WORKDIR}/bin/python3"
[ -x "${PY}" ] || PY="${WORKDIR}/bin/python"
echo "Using PY=${PY}"
"${PY}" -V

# conda-pack일 때만 존재. 없어도 OK.
if [ -x "${WORKDIR}/bin/conda-unpack" ]; then
  echo "=== conda-unpack ==="
  "${WORKDIR}/bin/conda-unpack" || true
fi

# ----- unpack analysis code -----
echo "=== Unpack analysis.tgz ==="
tar -xzf analysis.tgz
test -d analysis
cd analysis
echo "PWD=$(pwd)"
ls -lah | head

# ----- run -----
echo "=== Run ==="
"${PY}" -u run.py \
  -p stop_new2024 \
  -m KNU_2024_v4 \
  -w 1 \
  -d "${DATASET}"

# ----- stage-out: copy futures to a fixed name at WORKDIR root -----
FUT="hists/stop_new2024/${DATASET}.futures"
echo "=== Locate futures ==="
ls -lah "hists/stop_new2024" || true

if [ ! -f "${FUT}" ]; then
  echo "ERROR: futures not found: ${FUT}" >&2
  echo "Dumping candidates:" >&2
  find . -maxdepth 6 -name "*.futures" -ls >&2 || true
  exit 2
fi

cp -v "${FUT}" "${WORKDIR}/out.futures"

echo "=== END $(date) ==="
