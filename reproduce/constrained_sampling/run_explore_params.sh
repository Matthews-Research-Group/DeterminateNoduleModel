#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -J bash_sample
#SBATCH -o log/job_%A_%a.out    # Save stdout (unique for each job)
#SBATCH -e log/job_%A_%a.err
#SBATCH --array=1-10                  # Submit 10 job instances
#SBATCH --ntasks=1                       # Number of tasks per job
#SBATCH --cpus-per-task=5               # Adjust CPU allocation
#SBATCH --mem=50g                         # Adjust memory allocation
#SBATCH -N 1
#SBATCH --mail-user=your_email
#SBATCH --mail-type=ALL
#SBATCH -D /your_desired_path

# ----------------Load Modules--------------------
module load MATLAB/2020b-IGB-gcc-8.2.0
# ----------------Commands------------------------
export TMPDIR=$PWD/tmp_$SLURM_JOB_ID
mkdir -p $TMPDIR

echo "SLURM Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "TMPDIR: $TMPDIR"

matlab -nodisplay -nosplash -nodesktop -r \
"cluster = parcluster('local'); \
cluster.JobStorageLocation = getenv('TMPDIR'); \
parpool(cluster, str2double(getenv('SLURM_CPUS_PER_TASK'))); \
explore_params($SLURM_ARRAY_TASK_ID); exit;"


