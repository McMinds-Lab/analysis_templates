source local.sh

module purge 
module load apps/sratoolkit/2.10.7

## download sequence data (below is command for single sample)
## https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/
fasterq-dump --split-3 SRR6793695 -O ${output_dir}/SRR6793695

module purge 
module load apps/spades/3.11.1
spades.py --pe1-1 ${output_dir}/SRR6793695/SRR6793695_1.fastq --pe1-2 ${output_dir}/SRR6793695/SRR6793695_2.fastq -o ${output_dir}/SRR6793695/spades

module purge 
module load hub.apps/anaconda3
conda activate bowtie2
module load apps/samtools/1.3.1
bowtie2-build ${output_dir}/SRR6793695/spades/scaffolds.fasta ${output_dir}/SRR6793695/spades/genome
bowtie2 -x ${output_dir}/SRR6793695/spades/genome -1 ${output_dir}/SRR6793695/SRR6793695_1.fastq -2 ${output_dir}/SRR6793695/SRR6793695_2.fastq | samtools view -bS - > ${output_dir}/SRR6793695/aligned_reads.bam

conda deactivate
conda activate telomerecat
telomerecat bam2length --output ${output_dir}/SRR6793695.csv ${output_dir}/SRR6793695/aligned_reads.bam
