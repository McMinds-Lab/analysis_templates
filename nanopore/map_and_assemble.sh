# https://help.nanoporetech.com/en/articles/6632022-what-is-the-lambda-dna-control-lmd

conda activate ncbi_datasets

datasets download genome accession GCF_014640905.1 --filename ~/outputs/nano/pyruvatibacter_genome.zip

unzip ~/outputs/nano/pyruvatibacter_genome.zip -d ~/outputs/nano/pyruvatibacter_genome
rm ~/outputs/nano/pyruvatibacter_genome.zip

## map reads to reference
conda activate minimap2
minimap2 -ax map-ont ~/outputs/nano/pyruvatibacter_genome/ncbi_dataset/data/*/*.fna ~/outputs/nano/duplex.fastq.gz > ~/outputs/nano/pyruvatibacter_map.sam

## find mapped reads
conda activate samtools
samtools view -S -F4 ~/outputs/nano/pyruvatibacter_map.sam | cut -f1 | sort | uniq > ~/outputs/nano/pyruvatibacter_reads.txt

conda activate seqtk
seqtk subseq <(zcat <~/outputs/nano/duplex.fastq.gz) ~/outputs/nano/pyruvatibacter_reads.txt > ~/outputs/nano/pyruvatibacter_reads.fastq

## assemble and or annotate remainder
conda activate flye
flye -m 1000 --nano-hq ~/outputs/nano/pyruvatibacter_reads.fastq -o ~/outputs/nano/flye_out_pyruvatibacter

circlator all --data_type nanopore-corrected --assemble_not_only_assembler --assemble_not_careful ~/outputs/nano/flye_out_pyruvatibacter/assembly.fasta ~/outputs/nano/pyruvatibacter_reads.fastq ~/outputs/nano/flye_out_pyruvatibacter/circlator_out
