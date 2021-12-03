# get local variables
source local.env

## two positional arguments specifying 1) an SRA project ID starting with SRP and 2) the output directory
srp=$1
outdir=$2

mkdir -p ${outdir}/logs

module purge
module load hub.apps/anaconda3
source activate entrez-direct

esearch -db sra -query ${srp} | efetch -format runinfo > ${outdir}/runInfo.csv

samples=($(grep SRR ${outdir}/runInfo.csv | cut -d ',' -f 1))

cat <<EOF > ${outdir}/download_SRA.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=${maxram}
#SBATCH --job-name=download_SRA
#SBATCH --output=${outdir}/logs/download_SRA_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%${nthreads}

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

subdir=\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 13)

mkdir -p ${outdir}/\${subdir}/\${sample}

module purge
module load hub.apps/anaconda3
conda activate entrez-direct

biosample=\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 26)
esearch -db biosample -query \${biosample} | efetch > ${outdir}/\${subdir}/\${sample}/metadata.txt

conda deactivate
module purge
module load apps/sratoolkit/2.10.7

fasterq-dump -e 1 \${sample} -t ${outdir}/\${subdir}/\${sample} -O ${outdir}/\${subdir}/\${sample}

gzip --best -c ${outdir}/\${subdir}/\${sample}/\${sample}_1.fastq > ${outdir}/\${subdir}/\${sample}/\${sample}_1.fastq.gz &
gzip --best -c ${outdir}/\${subdir}/\${sample}/\${sample}_2.fastq > ${outdir}/\${subdir}/\${sample}/\${sample}_2.fastq.gz

rm ${outdir}/\${subdir}/\${sample}/\${sample}_1.fastq
rm ${outdir}/\${subdir}/\${sample}/\${sample}_2.fastq

EOF

if $autorun; then

sbatch ${outdir}/download_SRA.sbatch

fi

