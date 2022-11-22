indir=$1
outdir=$2
ntasks=$3

mkdir -p ${outdir}/02_dada2/
echo "bash $0 $@" > ${outdir}/02_dada2/this_command.sh
cp 02_dada2.r ${outdir}/02_dada2/

cat <<EOF > ${outdir}/02_dada2/02_dada2.sbatch
#!/bin/bash
#SBATCH --job-name=02_dada2
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=${outdir}/02_dada2/02_dada2.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=20G
#SBATCH --time=3-00:00:00

module load hub.apps/R
Rscript ${outdir}/02_dada2/02_dada2.r ${nthreads} ${indir} ${outdir}/02_dada2/

EOF

if $autorun; then
    sbatch ${outdir}/02_dada2/02_dada2.sbatch
fi
