#!/bin/bash

#SBATCH --job-name=getTransSNP_%a
#SBATCH -p mignot
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=../SlurmOutputs/getTransSNP/getTransSNP_%a.out
#SBATCH --error=../SlurmOutputs/getTransSNP/getTransSNP_%a.err
#SBATCH --nodes=1
#SBATCH -t 3:00:00

cd ..
python3 getTransSNP.py -ChrIndex $SLURM_ARRAY_TASK_ID
