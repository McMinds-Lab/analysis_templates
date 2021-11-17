# get local variables
source local.env

mkdir -p ${outdir}/02_dada2

cat <<EOF > ${outdir}/02_dada2/02_dada2.sbatch
#!/bin/bash
#SBATCH --job-name=02_dada2
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/02_dada2/02_dada2.log
#SBATCH --ntasks=${nthreads}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20
#SBATCH --time=01:00:00

module load apps/R
Rscript 02_dada2.r

EOF

if $autorun; then
    sbatch ${outdir}/01_init_QC/01_init_QC.sbatch
fi
