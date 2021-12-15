# get local variables
source local.env

mkdir -p ${outdir}/02_dada2

mkdir -p ${outdir}/02_dada2/data

cat <<EOF > ${outdir}/02_dada2/02_dada2.sbatch
#!/bin/bash
#SBATCH --job-name=02_dada2
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/02_dada2/02_dada2.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=20G
#SBATCH --time=01:00:00

module load hub.apps/anaconda3
source activate dada2
Rscript 02_dada2.r

EOF

if $autorun; then
    sbatch ${outdir}/02_dada2/02_dada2.sbatch
fi
