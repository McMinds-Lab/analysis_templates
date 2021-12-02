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
#SBATCH --mem=10G
#SBATCH --job-name=download_SRA
#SBATCH --output=${outdir}/download_SRA.log

module purge
module load apps/sratoolkit/2.10.7

for sample in ${samples[@]}; do

fasterq-dump --split-3 \${sample} -O ${outdir}/\${sample}
gzip -c ${outdir}/\${sample}/\${sample}_1.fastq > ${outdir}/\${sample}/\${sample}_1.fastq.gz
rm ${outdir}/\${sample}/\${sample}_1.fastq
gzip -c ${outdir}/\${sample}/\${sample}_2.fastq > ${outdir}/\${sample}/\${sample}_2.fastq.gz
rm ${outdir}/\${sample}/\${sample}_2.fastq

done

EOF

if $autorun; then

sbatch ${outdir}/download_SRA.sbatch

fi

