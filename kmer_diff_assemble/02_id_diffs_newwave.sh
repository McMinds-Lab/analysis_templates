
in_kmers=$1
sampledat=$2
threshold=$3
formula=$4
keycolumn=$5
outdir=$6
n_nodes=$7 ## will not actually work with more than 1 node at the moment
autorun=$8

subdir=${outdir}/02_id_diffs

mkdir ${subdir}

## get samplewise offsets and latent factors
cat <<EOF > ${subdir}/02_id_diffs_newwave.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --job-name=02_id_diffs
#SBATCH --nodes=${n_nodes}
#SBATCH --exclusive
#SBATCH --output=${subdir}/id_diffs_newwave.log

nodenames=\$(scontrol show hostname \$SLURM_NODELIST | tr '\n' ' ')

module purge
module load hub.apps/R/4.1.1

Rscript 02_id_diffs_newwave.r ${in_kmers} ${sampledat} ${threshold} ${formula} ${keycolumn} ${outdir} "\${nodenames}"

EOF

if $autorun; then
   sbatch ${subdir}/02_id_diffs_newwave.sbatch
fi
