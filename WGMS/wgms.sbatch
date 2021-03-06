#!/bin/bash
#
#SBATCH --job-name=reconstruction
#SBATCH --partition=low
#SBATCH --exclusive
#SBATCH --array=1-19

# Abort whenever a single step fails. Without this, bash will just continue on errors.
#set -e

# Load the required environment modules for OGGM
module load python/3.6.4 oggm-binary-deps/1

# Activate our local OGGM virtualenv
source /home/users/julia/python3_env/bin/activate

# On every node, when slurm starts a job, it will make sure the directory
# /work/username exists and is writable by the jobs user.
# We create a sub-directory there for this job to store its runtime data at.
S_WORKDIR="/work/$SLURM_JOB_USER/wgms/temp_0"
OUTDIR="/home/users/julia/initialization/out"

TEMP_BIAS = 0

mkdir -p "$S_WORKDIR"
mkdir -p "$OUTDIR"

export S_WORKDIR
export OUTDIR

echo "Workdir for this run: $S_WORKDIR"
echo "$SLURM_ARRAY_TASK_ID"

# I should be one value from I = 0-3
export I=0
export REGION = $SLURM_ARRAY_TASK_ID

export OGGM_DOWNLOAD_CACHE="/home/data/download"
export OGGM_DOWNLOAD_CACHE_RO=1


# Run the actual job. The srun invocation starts it as individual step for slurm.
srun -n 1 -c "${SLURM_JOB_CPUS_PER_NODE}" python3 ./wgms_experiments.py
# Print a final message so you can actually see it being done in the output log.
echo "DONE"

# Once a slurm job is done, slurm will clean up the /work directory on that node from any leftovers $
# So copy any result data you need from there back to your home dir!
# $SLURM_SUBMIT_DIR points to the directory from where the job was initially commited.

# Copy any neccesary result data.
cp -r "${S_WORKDIR}" "${OUTDIR}"


