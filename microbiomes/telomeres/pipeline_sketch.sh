## define variables
source local.sh

## download sequence data (below is command for single sample)
## https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/
module purge 
module load apps/sratoolkit/2.10.7
fasterq-dump --split-3 ${sample} -O ${output_dir}/${sample}

## convert fastqs to bams
module purge 
module load hub.apps/anaconda3
conda activate picard
picard FastqToSam -F1 ${output_dir}/${sample}/${sample}_1.fastq -F2 ${output_dir}/${sample}/${sample}_2.fastq -O ${output_dir}/${sample}/unaligned_read_pairs.bam -SM ${sample} -SO unsorted

## telomerecat requires a specific header format for some reason
module load apps/samtools/1.3.1
samtools view -H ${output_dir}/${sample}/unaligned_read_pairs.bam > ${output_dir}/${sample}/myheader.sam
printf "@SQ\tSN:fake_contig\tLN:1\n" >> ${output_dir}/${sample}/myheader.sam 
samtools reheader ${output_dir}/${sample}/myheader.sam ${output_dir}/${sample}/unaligned_read_pairs.bam > ${output_dir}/${sample}/unaligned_read_pairs.bam_tmp
mv ${output_dir}/${sample}/unaligned_read_pairs.bam_tmp ${output_dir}/${sample}/unaligned_read_pairs.bam

## run telomerecat
conda deactivate
conda activate telomerecat
telomerecat bam2length --output ${output_dir}/${sample}.csv ${output_dir}/${sample}/unaligned_read_pairs.bam
