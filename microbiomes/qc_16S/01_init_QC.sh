#!/bin/bash
#SBATCH --job-name=01_init_qc
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mail-user=salexander4@usf.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/01_init_qc.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20
#SBATCH --time=01:00:00

outdir=$WORK/output/dada2_tutorial
in_fwd=
in_rev=
barcodes_fwd=barcodes_fwd.fasta
barcodes_rev=barcodes_rev.fasta
nthreads=1

mkdir -p ${outdir}/1_init_QC/demultiplexed/
mkdir -p ${outdir}/1_init_QC/merged

# trim indices and primers from sequences, demultiplex, and discard any sequences that don't contain both full barcodes and primers
module purge
module load hub.apps/anaconda3
source activate cutadapt
cutadapt \
  --no-indels \
  --discard-untrimmed \
  --pair-filter=any \
  -g file:${barcodes_fwd} \
  -G file:${barcodes_rev} \
  -o ${outdir}/trimmed/{name1}-{name2}_R1.fastq \
  -p ${outdir}/trimmed/{name1}-{name2}_R2.fastq \
  ${in_fwd} \
  ${in_rev}

for file in ${outdir}/1_init_QC/demultiplexed/*_R1.fastq; do

  # trim "R1" from filenames to get Sample IDs that match mapping file
  filename=$(basename $file)
  sampleid=${filename/_*/} ## need to somehow translate the filenames to sample IDs (this doesn't actually do that)

  # merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)

  conda deactivate
  module purge
  module load apps/vsearch
  vsearch \
    --fastq_mergepairs ${file} \
    --reverse ${file/R1/R2} \
    --fastq_allowmergestagger \
    --fasta_width 0 \
    --threads ${nthreads} \
    --fastqout ${outdir}/1_init_QC/merged/${sampleid}.fastq

done

