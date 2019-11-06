#!/bin/bash
#
#SBATCH --job-name=reconstruction
#SBATCH --ntasks=1 --exclusive
#SBATCH --array=1-19
#SBATCH --export=NONE

# This script has to be submitted from the WGMS subdir of the initialization tree.

set -e

export WORKDIR="/work/${SLURM_JOB_USER}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export S_WORKDIR="$WORKDIR/wgms"
export OUTDIR="$HOME/initialization/out"

mkdir -p "$S_WORKDIR"
mkdir -p "$OUTDIR"

export TEMP_BIAS=0

echo "Workdir for this run: $S_WORKDIR"
echo "$SLURM_ARRAY_TASK_ID"

# I should be one value from I = 0-3
export I=0
export REGION="$SLURM_ARRAY_TASK_ID"

export OGGM_DOWNLOAD_CACHE="/home/data/download"
export OGGM_DOWNLOAD_CACHE_RO=1

singularity exec docker://oggm/oggm:20191104 bash -s <<EOF
  set -e
  pwd
  python3 -m venv --system-site-packages "$WORKDIR/venv"
  source "$WORKDIR/venv/bin/activate"
  pip install --upgrade pip setuptools
  pip install ..

  exec python3 ./wgms_experiments.py
EOF

echo "Processing DONE"

cp -r "${S_WORKDIR}/." "${OUTDIR}"
rm -rf "${WORKDIR}"

echo "DONE"
