
in_kmers=$1
sampledat=$2
threshold=$3
formula=$4
keycolumn=$5
outdir=$6
n_processes=$7
maxmem=$8

subdir=${outdir}/02_id_diffs

mkdir ${subdir}

## get samplewise offsets and latent factors
cat <<EOF > ${subdir}/02a_id_diffs_newwave.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --mem=${maxmem}
#SBATCH --time=7-00:00:00
#SBATCH --job-name=02_id_diffs
#SBATCH --ntasks=${n_processes}
#SBATCH --cpus-per-task=1
#SBATCH --output=${subdir}/id_diffs_newwave.log

module purge
module load hub.apps/R/4.1.1

Rscript 02_id_diffs_newwave.r ${in_kmers} ${sampledat} ${threshold} ${formula} ${keycolumn} ${outdir} ${n_cores}

EOF

sbatch ${subdir}/02_id_diffs_newwave.sbatch
