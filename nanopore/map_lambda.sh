# https://help.nanoporetech.com/en/articles/6632022-what-is-the-lambda-dna-control-lmd

conda activate ncbi_datasets

datasets download virus genome accession NC_001416.1 --filename ~/outputs/lambda/lambda_genome.zip

unzip /Users/Ryan/outputs/lambda/lambda_genome.zip -d /Users/Ryan/outputs/lambda/lambda_genome
rm /Users/Ryan/outputs/lambda/lambda_genome.zip

## map reads to lambda reference
conda activate minimap2
mkdir /Users/Ryan/outputs/sym1
minimap2 -ax map-ont /Users/Ryan/outputs/lambda/lambda_genome/ncbi_dataset/data/genomic.fna /Users/Ryan/outputs/nano/pass/*.fastq.gz > /Users/Ryan/outputs/sym1/lambdamap.sam

## find unmapped reads
conda activate samtools
samtools view -S -f4 /Users/Ryan/outputs/sym1/lambdamap.sam | cut -f1 | sort | uniq > /Users/Ryan/outputs/sym1/notlambda_reads.txt

conda activate seqtk
seqtk subseq <(zcat </Users/Ryan/outputs/nano/pass/*.fastq.gz) /Users/Ryan/outputs/sym1/notlambda_reads.txt > /Users/Ryan/outputs/sym1/notlambda_reads.fastq

## assemble and or annotate remainder
conda activate flye
flye -m 1000 --meta --nano-raw /Users/Ryan/outputs/sym1/notlambda_reads.fastq -o /Users/Ryan/outputs/sym1/flye_out

## annotate
conda activate mmseqs2
mmseqs databases UniProtKB/Swiss-Prot /Users/Ryan/outputs/sym1/swissprot /Users/Ryan/outputs/sym1/tmp
mmseqs easy-taxonomy /Users/Ryan/outputs/sym1/flye_out/assembly.fasta /Users/Ryan/outputs/sym1/swissprot /Users/Ryan/outputs/sym1/mmseqtax /Users/Ryan/outputs/sym1/tmp

## consider aligning to this and then looking at remainder: https://www.ncbi.nlm.nih.gov/nuccore/HG966617.1
