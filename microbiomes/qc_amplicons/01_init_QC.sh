# get local variables
source local.env

## five positional arguments specifying 1) input fwd fastq, 2) input rev fastq, 3) forward primer sequence, 4) reverse primer sequence, 5) input fwd barcodes file, 6) input rev barcodes file, and 7) the output directory
in_fwd=$1
in_rev=$2
primer_fwd=$3
primer_rev=$4
barcodes_fwd=$5
barcodes_rev=$6
outdir=$7

mkdir -p ${outdir}/01_init_QC

cat <<EOF > ${outdir}/01_init_QC/01_init_QC.sbatch
#!/bin/bash
#SBATCH --job-name=01_init_QC
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/01_init_QC/01_init_QC.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=20G
#SBATCH --time=24:00:00

mkdir -p ${outdir}/01_init_QC/demultiplexed/
mkdir -p ${outdir}/01_init_QC/merged

# trim indices and primers from sequences, demultiplex, and discard any sequences that don't contain both full barcodes and primers
module purge
module load hub.apps/anaconda3
source activate cutadapt

# find all pairs with both primers in forward orientation
cutadapt \
  --cores=${nthreads} \
  --action=none \
  -g "${primer_fwd};min_overlap=${#primer_fwd}" \
  -G "${primer_rev};min_overlap=${#primer_rev}" \
  --output ${outdir}/01_init_QC/oriented_R1.fastq.gz \
  --paired-output ${outdir}/01_init_QC/oriented_R2.fastq.gz \
  --untrimmed-output ${outdir}/01_init_QC/disoriented_R1.fastq.gz \
  --untrimmed-paired-output ${outdir}/01_init_QC/disoriented_R2.fastq.gz \
  ${in_fwd} \
  ${in_rev}

# find all pairs with both primers in reverse orientation
cutadapt \
  --cores=${nthreads} \
  --action=none \
  -g "${primer_fwd};min_overlap=${#primer_fwd}" \
  -G "${primer_rev};min_overlap=${#primer_rev}" \
  --output ${outdir}/01_init_QC/reoriented_R1.fastq.gz \
  --paired-output ${outdir}/01_init_QC/reoriented_R2.fastq.gz \
  --discard-untrimmed \
  ${outdir}/01_init_QC/disoriented_R2.fastq.gz \
  ${outdir}/01_init_QC/disoriented_R1.fastq.gz

# demultiplex, allowing a single error in each 8bp index
cutadapt \
  -e 0.15 \
  --overlap 8 \
  --no-indels \
  -g file:${barcodes_fwd} \
  -G file:${barcodes_rev} \
  --output ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}_R1.fastq.gz \
  --paired-output ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}_R2.fastq.gz \
  <(zcat ${outdir}/01_init_QC/oriented_R1.fastq.gz ${outdir}/01_init_QC/reoriented_R1.fastq.gz) \
  <(zcat ${outdir}/01_init_QC/oriented_R2.fastq.gz ${outdir}/01_init_QC/reoriented_R2.fastq.gz)

rm ${outdir}/01_init_QC/demultiplexed/*unknown*

for file in ${outdir}/01_init_QC/demultiplexed/*_R1.fastq.gz; do

  # trim "R1" from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/_R1*/}
  
  # trim primers
  source activate cutadapt
  cutadapt \
    --cores=${nthreads} \
    -g ${primer_fwd} \
    -G ${primer_rev} \
    --output ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R1.fastq.gz \
    --paired-output ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R2.fastq.gz \
    \${file} \
    \${file/R1/R2}
  
  # merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)
  source activate vsearch
  vsearch \
    --threads ${nthreads}\
    --fastq_mergepairs ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R1.fastq.gz \
    --reverse ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R2.fastq.gz \
    --fastq_allowmergestagger \
    --fastq_maxdiffs 50 \
    --fasta_width 0 \
    --fastqout_notmerged_fwd ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
    --fastqout_notmerged_rev ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq \
    --fastqout - | gzip --best > ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz

  # add concatenated fwd and reverse reads that could not be merged, assuming that they were not merged because they had a gap between reads
  vsearch \
    --fastx_join ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
    --reverse ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq \
    --join_padgap '' \
    --fastqout - | gzip --best >> ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
  
done

EOF

if $autorun; then

sbatch ${outdir}/01_init_QC/01_init_QC.sbatch

fi
