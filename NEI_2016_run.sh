#!/bin/bash
# Bash script to submit standard run to queue

#number of nodes to use. can set to >1 for MPI jobs, otherwise leave at 1 to make sure all allocated CPUs are on same node.
#SBATCH -N 1

#number of CPUs to utilize (this job is one task)
#SBATCH -c 24

#memory per node. job will crash if limit is exceeded. Default is 1G per allocated core
#SBATCH --mem=60G

#time limit after which the job will be killed. default time limit is 60 minutes. Specified as HH:MM:SS or D-HH:MM
#SBATCH -t 03-00:00:00

#partition on which to run the job.
#SBATCH -p edr

#output from the job
#SBATCH -o my_run-%j.out
#SBATCH -e my_run-%j.err

module avail
#load the GC bash
conda init bash
conda activate conda_env
echo "running"
python GriddedEpa2GCfinal.py 02

#configure stack size and remove limit
ulimit -s unlimited
export OMP_STACKSIZE=400m

#number of CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

exit $?
