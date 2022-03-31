# get local variables
source local.env

indir=$1
outdir=$2

mkdir -p ${outdir}/02_dada2/
echo "bash $0 $@" > ${outdir}/02_dada2/this_command.sh
cp 02_dada2.r ${outdir}/02_dada2/

cat <<EOF > ${outdir}/02_dada2/02_dada2.sbatch
#!/bin/bash
#SBATCH --job-name=02_dada2
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/02_dada2/02_dada2.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=${maxram}
#SBATCH --time=3-00:00:00

module load hub.apps/R
Rscript ${outdir}/02_dada2/02_dada2.r ${nthreads} ${indir} ${outdir}/02_dada2/

EOF

if $autorun; then
    sbatch ${outdir}/02_dada2/02_dada2.sbatch
fi
