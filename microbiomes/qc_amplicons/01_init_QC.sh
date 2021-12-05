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
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00

mkdir -p ${outdir}/01_init_QC/demultiplexed/
mkdir -p ${outdir}/01_init_QC/merged

# trim indices and primers from sequences, demultiplex, and discard any sequences that don't contain both full barcodes and primers
module purge
module load hub.apps/anaconda3
source activate cutadapt

cutadapt \
  --action=none \
  -g "${primer_fwd};min_overlap=${#primer_fwd}" \
  -G "${primer_rev};min_overlap=${#primer_rev}" \
  --output ${outdir}/01_init_QC/oriented_R1.fastq \
  --paired-output ${outdir}/01_init_QC/oriented_R2.fastq \
  --untrimmed-output ${outdir}/01_init_QC/disoriented_R1.fastq \
  --untrimmed-paired-output ${outdir}/01_init_QC/disoriented_R2.fastq \
  ${in_fwd} \
  ${in_rev}

cutadapt \
  --action=none \
  -g "${primer_fwd};min_overlap=${#primer_fwd}" \
  -G "${primer_rev};min_overlap=${#primer_rev}" \
  --output ${outdir}/01_init_QC/reoriented_R1.fastq \
  --paired-output ${outdir}/01_init_QC/reoriented_R2.fastq \
  ${outdir}/01_init_QC/disoriented_R2.fastq \
  ${outdir}/01_init_QC/disoriented_R1.fastq

cutadapt \
  -e 0.15 \
  --overlap 8 \
  --no-indels \
  -g file:${barcodes_fwd} \
  -G file:${barcodes_rev} \
  -o ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}_R1.fastq \
  -p ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}_R2.fastq \
  <(cat ${outdir}/01_init_QC/oriented_R1.fastq ${outdir}/01_init_QC/reoriented_R1.fastq) \
  <(cat ${outdir}/01_init_QC/oriented_R2.fastq ${outdir}/01_init_QC/reoriented_R2.fastq)

rm ${outdir}/01_init_QC/demultiplexed/*unknown*

for file in ${outdir}/01_init_QC/demultiplexed/*_R1.fastq; do

  # trim "R1" from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/_R1*/}
  
  cutadapt \
    -g ${primer_fwd} \
    -G ${primer_rev} \
    --output ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R1.fastq \
    --paired-output ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R2.fastq \
    \${file} \
    \${file/R1/R2}
  
  # merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)
  module load apps/vsearch
  vsearch \
    --fastq_mergepairs ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R1.fastq \
    --reverse ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R2.fastq \
    --fastq_allowmergestagger \
    --fasta_width 0 \
    --threads ${nthreads} \
    --fastqout - | gzip --best > ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
    
done

EOF

if $autorun; then

sbatch ${outdir}/01_init_QC/01_init_QC.sbatch

fi
