## two positional arguments specifying the directory of SRA downloads (with sample_list.txt in it) and the output analysis directory
readarray -t samples < $1/sample_list.txt
indir=$1
outdir=$2

# get local variables
source local.env

mkdir -p ${outdir}/01_telomerecat/logs

cat <<EOF > ${outdir}/01_telomerecat/01_telomerecat.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=10G
#SBATCH --job-name=01_telomerecat
#SBATCH --output=${outdir}/01_telomerecat/logs/01_telomerecat_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%10

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

##convert fastqs to bams
module purge 
module load hub.apps/anaconda3
source activate picard

echo picard

picard FastqToSam \
  -F1 ${outdir}/\${sample}/\${sample}_1.fastq.gz \
  -F2 ${outdir}/\${sample}/\${sample}_2.fastq.gz \
  -O ${outdir}/\${sample}/unaligned_read_pairs.bam \
  -SM \${sample} \
  -SO unsorted

##telomerecat requires a specific header format for some reason
conda deactivate
module purge 
module load apps/samtools/1.3.1
samtools view -H ${outdir}/\${sample}/unaligned_read_pairs.bam > ${outdir}/\${sample}/myheader.sam
printf "@SQ\tSN:fake_contig\tLN:1\n" >> ${outdir}/\${sample}/myheader.sam
samtools reheader ${outdir}/\${sample}/myheader.sam ${outdir}/\${sample}/unaligned_read_pairs.bam > ${outdir}/\${sample}/unaligned_read_pairs.bam_tmp
mv ${outdir}/\${sample}/unaligned_read_pairs.bam_tmp ${outdir}/\${sample}/unaligned_read_pairs.bam

## run telomerecat
module purge 
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
source activate telomerecat
telomerecat bam2length --output ${outdir}/\${sample}.csv ${outdir}/\${sample}/unaligned_read_pairs.bam

EOF

if $autorun; then
    sbatch ${outdir}/01_telomerecat/01_telomerecat.sbatch
fi
