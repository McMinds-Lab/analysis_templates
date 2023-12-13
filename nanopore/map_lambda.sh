# https://help.nanoporetech.com/en/articles/6632022-what-is-the-lambda-dna-control-lmd

conda activate ncbi_datasets

mkdir ~/outputs/lambda/

datasets download virus genome accession NC_001416.1 --filename ~/outputs/lambda/lambda_genome.zip

unzip ~/outputs/lambda/lambda_genome.zip -d ~/outputs/lambda/lambda_genome
rm ~/outputs/lambda/lambda_genome.zip

## map reads to lambda reference
conda activate samtools
samtools fastq ~/outputs/nano/duplex.bam | gzip --best > ~/outputs/nano/duplex.fastq.gz
conda activate minimap2
minimap2 -ax map-ont ~/outputs/lambda/lambda_genome/ncbi_dataset/data/genomic.fna ~/outputs/nano/duplex.fastq.gz > ~/outputs/nano/lambdamap.sam

## find unmapped reads
conda activate samtools
samtools view -S -f4 ~/outputs/nano/lambdamap.sam | cut -f1 | sort | uniq > ~/outputs/nano/notlambda_reads.txt

conda activate seqtk
seqtk subseq <(zcat <~/outputs/nano/duplex.fastq.gz) ~/outputs/nano/notlambda_reads.txt > ~/outputs/nano/notlambda_reads.fastq

## assemble and or annotate remainder
conda activate flye
flye -m 1000 --meta --nano-raw ~/outputs/nano/notlambda_reads.fastq -o ~/outputs/nano/flye_out

## annotate
conda activate mmseqs2
mmseqs databases --compressed 1 GTDB ~/outputs/nano/GTDB ~/outputs/nano/tmp
mmseqs easy-taxonomy --compressed 1 ~/outputs/nano/flye_out/assembly.fasta ~/outputs/nano/GTDB ~/outputs/nano/mmseqtax ~/outputs/nano/tmp

## consider aligning to this and then looking at remainder: https://www.ncbi.nlm.nih.gov/nuccore/HG966617.1
