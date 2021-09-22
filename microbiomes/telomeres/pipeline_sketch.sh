## define variables
source local.sh

## download sequence data (below is command for single sample)
## https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/
module purge 
module load apps/sratoolkit/2.10.7
fasterq-dump --split-3 ${sample} -O ${output_dir}/${sample}

## assemble genome. IRL this should be done with concatenated set of all samples from study, not once-per-sample. Otherwise of course use a reference genome if available
module purge 
module load apps/spades/3.11.1
spades.py --pe1-1 ${output_dir}/${sample}/${sample}_1.fastq --pe1-2 ${output_dir}/${sample}/${sample}_2.fastq -o ${output_dir}/${sample}/spades

## align to genome
module purge 
module load hub.apps/anaconda3
conda activate bowtie2
module load apps/samtools/1.3.1
bowtie2-build ${output_dir}/${sample}/spades/scaffolds.fasta ${output_dir}/${sample}/spades/genome
bowtie2 -x ${output_dir}/${sample}/spades/genome -1 ${output_dir}/${sample}/${sample}5_1.fastq -2 ${output_dir}/${sample}/${sample}5_2.fastq | samtools view -bS - > ${output_dir}/${sample}/aligned_reads.bam

## run telomerecat
conda deactivate
conda activate telomerecat
telomerecat bam2length --output ${output_dir}/${sample}.csv ${output_dir}/${sample}/aligned_reads.bam
