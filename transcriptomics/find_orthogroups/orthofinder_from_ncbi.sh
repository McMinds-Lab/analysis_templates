source scripts/immune_allometry/de_novo/local.env

mkdir -p ${outdir}/04_orthofinder/logs

cat <<EOF > ${outdir}/04_orthofinder/04_orthofinder.sbatch
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --time=6-20:00:00
#SBATCH --mem=${maxram}
#SBATCH --ntasks=${nprocesses}
#SBATCH --mail-user=${email}
#SBATCH --job-name=04_orthofinder
#SBATCH --output=${outdir}/04_orthofinder/logs/04_orthofinder.out
#SBATCH --error=${outdir}/04_orthofinder/logs/04_orthofinder.err
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}

module purge
module load hub.apps/anaconda3/2020.11
conda activate orthofinder

orthofinder -t ${nprocesses} -d -M msa -A mafft -T iqtree \
  -f ${outdir}/03_corset/corset_longest/ \
  -o ${outdir}/04_orthofinder/of_out

EOF

if $autorun; then
    sbatch ${outdir}/04_orthofinder/04_orthofinder.sbatch
fi
