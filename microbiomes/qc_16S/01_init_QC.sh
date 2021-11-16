# get local variables
source local.env

mkdir -p ${outdir}/01_init_QC

cat <<EOF > ${outdir}/01_init_QC/01_init_QC.sbatch
#!/bin/bash
#SBATCH --job-name=01_init_QC
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/01_init_QC/01_init_QC.log
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20
#SBATCH --time=01:00:00

# get local variables
source local.env

mkdir -p ${outdir}/01_init_QC/demultiplexed/
mkdir -p ${outdir}/01_init_QC/merged

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
  -o ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}_R1.fastq \
  -p ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}_R2.fastq \
  ${in_fwd} \
  ${in_rev}
# double check that reads are oriented consistently (does the above need to be re-run with the forward and reverse indices and or primers switched?)

for file in ${outdir}/01_init_QC/demultiplexed/*_R1.fastq; do

  # trim "R1" from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/-*/} ## double check that all files have matching name1 and name2 or else this could overwrite a good file with a bad one

  # merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)

  conda deactivate
  module purge
  module load apps/vsearch
  vsearch \
    --fastq_mergepairs \${file} \
    --reverse \${file/R1/R2} \
    --fastq_allowmergestagger \
    --fasta_width 0 \
    --threads ${nthreads} \
    --fastqout ${outdir}/01_init_QC/merged/\${sampleid}.fastq

done

EOF

if $autorun; then
    sbatch ${outdir}/01_init_QC/01_init_QC.sbatch
fi
