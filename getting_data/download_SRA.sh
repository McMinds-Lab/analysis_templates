## two positional arguments specifying 1) a file with single NCBI Run ID per line (e.g. SRR15960006) and 2) the output directory
readarray -t samples < $1
outdir=$2

# get local variables
source local.env

mkdir -p ${outdir}
cp $1 ${outdir}/sample_list.txt

cat <<EOF > ${outdir}/download_SRA.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=${maxram}
#SBATCH --ntasks=${nthreads}
#SBATCH --job-name=download_SRA
#SBATCH --output=${outdir}/download_SRA.log
#SBATCH --array=0-$((${#samples[@]}-1))%10

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

module purge
module load apps/sratoolkit/2.10.7

fasterq-dump -e ${nthreads} \${sample} -O ${outdir}/\${sample}
gzip -c ${outdir}/\${sample}/\${sample}_1.fastq > ${outdir}/\${sample}/\${sample}_1.fastq.gz
rm ${outdir}/\${sample}/\${sample}_1.fastq
gzip -c ${outdir}/\${sample}/\${sample}_2.fastq > ${outdir}/\${sample}/\${sample}_2.fastq.gz
rm ${outdir}/\${sample}/\${sample}_2.fastq

EOF

if $autorun; then

sbatch ${outdir}/download_SRA.sbatch

fi

