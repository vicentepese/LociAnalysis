#!/bin/bash

#SBATCH --job-name=getSNP_%a
#SBATCH -p owners
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=../SlurmOutputs/getSNP/getSNP_CHR%a.out
#SBATCH --error=../SlurmOutputs/getSNP/getSNP_CHR%a.err
#SBATCH --nodes=1
#SBATCH -t 01:00:00