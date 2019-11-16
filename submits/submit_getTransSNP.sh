#!/bin/bash

#SBATCH --job-name=getTransSNP_%a
#SBATCH -p owners
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=20G
#SBATCH --output=../SlurmOutputs/getTransSNP/getTransSNP_%a.out
#SBATCH --error=../SlurmOutputs/getTransSNP/getTransSNP_%a.err
#SBATCH --nodes=1
#SBATCH -t 03:00:00

cd ..
python3 getTransSNP.py -ChrIndex $SLURM_ARRAY_TASK_ID