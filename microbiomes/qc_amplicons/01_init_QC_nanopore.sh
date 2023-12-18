## six positional arguments specifying 1) input basecalls bam (with simplex basecalling), 2) forward primer sequence, 3) reverse primer sequence, 4) orienting sequence next to barcode on forward end, 5) orienting sequence next to barcode on reverse end, 6) input linked barcodes file, and 7) the output directory
in_bam=$1
primer_fwd=$2
primer_rev=$3
twostep_fwd=$4
twostep_rev=$5
barcodes_file=$6
outdir=$7

mkdir -p ${outdir}/01_init_QC
echo "bash $0 $@" > ${outdir}/01_init_QC/this_command.sh

cat <<EOF > ${outdir}/01_init_QC/01_init_QC.sbatch
#!/bin/bash
#SBATCH --job-name=01_init_QC
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=${outdir}/01_init_QC/01_init_QC.log
#SBATCH --ntasks=24
#SBATCH --mem=20G
#SBATCH --time=24:00:00

mkdir -p ${outdir}/01_init_QC/demultiplexed/
mkdir -p ${outdir}/01_init_QC/trimmed/
module purge
module load hub.apps/anaconda3

# trim indices and primers from sequences, demultiplex, and discard any sequences that don't contain both full barcodes and primers

source activate samtools-1.19
# extract reads from dorado basecalling file
samtools fastq ${in_bam} | gzip > ${outdir}/01_init_QC/reads.fastq.gz

# find reverse complement of reverse primer
primer_rev_rc=$(echo ${primer_rev} | tr ACGTRYSWKMBVDHacgtryswkmbvdh TGCAYRSWMKVBHDtgcayrswmkvbhd | rev)
twostep_rev_rc=$(echo ${twostep_rev} | tr ACGTRYSWKMBVDHacgtryswkmbvdh TGCAYRSWMKVBHDtgcayrswmkvbhd | rev)

# find seq lengths
overlap_fwd=$((${#twostep_fwd}+8))
overlap_rev=$((${#twostep_rev}+8))

source activate cutadapt-4.6
# find all pairs with both adaptors plus at least enough bp to see whole index, reorient them, then trim everything before the index, then demultiplex
# revcomp option exists but has weird interaction with trimming so doing it manually with pipe
cutadapt \
  --cores=24 \
  --revcomp \
  --action=none \
  -g ${twostep_fwd}...\${twostep_rev_rc} \
  ${outdir}/01_init_QC/reads.fastq.gz |
cutadapt \
  --cores=24 \
  --action=retain \
  --discard-untrimmed \
  -g "N{8}${twostep_fwd};min_overlap=$((${#twostep_fwd}+8))...\${twostep_rev_rc}N{8};min_overlap=$((${#twostep_rev}+8))" \
  - |
cutadapt \
  --cores=24 \
  --no-indels \
  -g file:${barcodes_file} \
  --output ${outdir}/01_init_QC/demultiplexed/{name}.fastq.gz \
  -

rm ${outdir}/01_init_QC/demultiplexed/*unknown*

for file in ${outdir}/01_init_QC/demultiplexed/*.fastq.gz; do

  # trim extension from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/.fastq.gz/}
  
  # trim primers
  cutadapt \
    --cores=24 \
    --discard-untrimmed \
    -g ${primer_fwd}...\${primer_rev_rc} \
    --output ${outdir}/01_init_QC/trimmed/\${sampleid}.fastq.gz \
    \${file}
  
done

EOF

sbatch ${outdir}/01_init_QC/01_init_QC.sbatch
