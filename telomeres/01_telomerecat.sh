# get local variables
source local.env

## three positional arguments specifying 1) the directory of SRA downloads (with runInfo.csv in it), 2) the output analysis directory, and 3) either the string NONE or the path to a reference genome
indir=$1
outdir=$2
reference=$3

mkdir -p ${outdir}/01_telomerecat/logs

if [ ${reference} != NONE ]; then

cat <<EOF > ${outdir}/01_telomerecat/01_telomerecat_build_index.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --mem=${maxram}
#SBATCH --output=${outdir}/01_telomerecat/logs/01_telomerecat_index.log
module purge
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
source activate bowtie2

zcat ${reference} > ${outdir}/01_telomerecat/reference.fasta
bowtie2-build --threads ${nthreads} ${outdir}/01_telomerecat/reference.fasta ${outdir}/01_telomerecat/reference
EOF

jid=$(sbatch ${outdir}/01_telomerecat/01_telomerecat_build_index.sbatch | cut -d ' ' -f4)
dependency="#SBATCH --dependency=afterok:${jid}"
else
dependency=''
fi

samples=($(grep SRR ${indir}/runInfo.csv | grep 'WGA\|WGS' | cut -d ',' -f 1))

cat <<EOF > ${outdir}/01_telomerecat/01_telomerecat.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=${maxram}
#SBATCH --job-name=01_telomerecat
#SBATCH --output=${outdir}/01_telomerecat/logs/01_telomerecat_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%10
${dependency}

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

mkdir ${outdir}/01_telomerecat/\${sample}/

#trim adapters
module purge
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
source activate bbtools

in1=(${indir}/*/\${sample}/\${sample}_1.fastq.gz)
in2=(${indir}/*/\${sample}/\${sample}_2.fastq.gz)

bbduk.sh \
  in1=\${in1} \
  in2=\${in2} \
  out1=${outdir}/01_telomerecat/\${sample}/trimmed_R1.fastq \
  out2=${outdir}/01_telomerecat/\${sample}/trimmed_R2.fastq \
  ref=adapters,artifacts,phix \
  ktrim=l k=23 mink=11 hdist=2 qtrim=r minlength=40 forcetrimright=70 tpe maq=10 trimq=10 maxgc=0.6 mingc=0.4
  
echo \$((\$(wc -l ${outdir}/01_telomerecat/\${sample}/trimmed_R1.fastq | cut -d' ' -f1) / 4)) > ${outdir}/01_telomerecat/\${sample}_nreads.txt

if [ ${reference} = NONE ]; then

##convert fastqs to unaligned bams
module purge
module load hub.apps/anaconda3
source activate picard

echo picard
picard FastqToSam \
  -F1 ${outdir}/01_telomerecat/\${sample}/trimmed_R1.fastq \
  -F2 ${outdir}/01_telomerecat/\${sample}/trimmed_R2.fastq \
  -O ${outdir}/01_telomerecat/\${sample}/unaligned_read_pairs.bam \
  -SM \${sample} \
  -SO unsorted

##telomerecat requires a specific header format for some reason
conda deactivate
module purge
module load apps/samtools/1.3.1
samtools view -H ${outdir}/01_telomerecat/\${sample}/unaligned_read_pairs.bam > ${outdir}/01_telomerecat/\${sample}/myheader.sam
printf "@SQ\tSN:fake_contig\tLN:1\n" >> ${outdir}/01_telomerecat/\${sample}/myheader.sam
samtools reheader ${outdir}/01_telomerecat/\${sample}/myheader.sam ${outdir}/01_telomerecat/\${sample}/unaligned_read_pairs.bam > ${outdir}/01_telomerecat/\${sample}/\${sample}.bam
rm ${outdir}/01_telomerecat/\${sample}/unaligned_read_pairs.bam ${outdir}/01_telomerecat/\${sample}/myheader.sam

else

## align reads to reference
module purge
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
source activate bowtie2

bowtie2 --threads ${nthreads} \
  -x ${outdir}/01_telomerecat/reference \
  -1 ${outdir}/01_telomerecat/\${sample}/trimmed_R1.fastq \
  -2 ${outdir}/01_telomerecat/\${sample}/trimmed_R2.fastq |
samtools view -b - > ${outdir}/01_telomerecat/\${sample}/\${sample}.bam

fi

## run telomerecat
module purge
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
source activate telomerecat
telomerecat bam2length --output ${outdir}/01_telomerecat/\${sample}.csv ${outdir}/01_telomerecat/\${sample}/\${sample}.bam

rm -rf ${outdir}/01_telomerecat/\${sample}

EOF

if $autorun; then
    sbatch ${outdir}/01_telomerecat/01_telomerecat.sbatch
fi
