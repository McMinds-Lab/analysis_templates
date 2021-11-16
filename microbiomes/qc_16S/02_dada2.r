#!/bin/bash
#SBATCH --job-name=02_dada2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mail-user=salexander4@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/02_dada2.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20
#SBATCH --time=01:00:00
