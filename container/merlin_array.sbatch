#!/bin/bash
# One cohort subject per array task on Merlin. Submit on a merlin login node:
#   sbatch --array=1-N container/merlin_array.sbatch
# (N must match the count staged by merlin_stage.sh; -n uses TASK_MAX.)
# Aggregate when done:
#   head -1 out/batch_summary_task1.csv > out/summary.csv
#   tail -q -n +2 out/batch_summary_task*.csv >> out/summary.csv
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 0:30:00
#SBATCH -o /home-nfs/schwat1/pet-fdg/out/task-%A_%a.out
R=/home-nfs/schwat1/pet-fdg
singularity exec -B $R/bids:/bids:ro -B $R/out:/work -B $R/templates:/templates \
  $R/pet-fdg.sif pet-fdg batch /bids -n ${SLURM_ARRAY_TASK_MAX} \
  --line ${SLURM_ARRAY_TASK_ID} \
  --work /work/batch --summary /work/batch_summary_task${SLURM_ARRAY_TASK_ID}.csv
