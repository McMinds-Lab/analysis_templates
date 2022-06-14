#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=7-00:00:00
#SBATCH --mem=187G
#SBATCH --job-name=zeroinfl
#SBATCH --output=file.path
#SBATCH --cpus-per-task=24

outdir="file.path"

mkdir -p ${outdir}/logs

module purge
module load hub.apps/anaconda3
source activate dekupl

dekupl-run --configfile file.path/config.json -j24 --resources ram=187000 -p
