## nine positional arguments specifying 1) input basecalls bam (with simplex basecalling), 2) the maximum length an amplicon can be (will trim to this length to avoid chimeric reads) 3) forward primer sequence, 4) reverse primer sequence, 5) orienting sequence next to barcode on forward end, 6) orienting sequence next to barcode on reverse end, 7) input linked barcodes file, 8) the porechop-formatted adapters file, and 9) the output directory
in_bam=$1
maxlen=$2
primer_fwd=$3
primer_rev=$4
twostep_fwd=$5
twostep_rev=$6
barcodes_file=$7
adapters_file=$8
outdir=$9

mkdir -p ${outdir}/01_init_QC
echo "bash $0 $@" > ${outdir}/01_init_QC/this_command.sh

cat <<EOF > ${outdir}/01_init_QC/01_init_QC.sbatch
#!/bin/bash
#SBATCH --job-name=01_init_QC
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=${outdir}/01_init_QC/01_init_QC.log
#SBATCH --ntasks=24
#SBATCH --mem=175G
#SBATCH --time=48:00:00

mkdir -p ${outdir}/01_init_QC/demultiplexed/
mkdir -p ${outdir}/01_init_QC/trimmed/
module purge
module load hub.apps/anaconda3

# trim indices and primers from sequences, demultiplex, and discard any sequences that don't contain both full barcodes and primers

source activate samtools-1.19
# extract reads from dorado basecalling file
samtools fastq -@ 24 ${in_bam} | pigz -p 24 > ${outdir}/01_init_QC/reads.fastq.gz

# find adapters in the middle of the reads and split chimeras (also trims adapters from ends)
source activate porechop_abi
porechop_abi \
  --threads 24 \
  --min_split_read_size 50 \
  --extra_end_trim 0 \
  --middle_threshold 75.0 \
  --extra_middle_trim_good_side 0 \
  --extra_middle_trim_bad_side 0 \
  --discard_database \
  --custom_adapters ${adapters_file} \
  --input ${outdir}/01_init_QC/reads.fastq.gz \
  --output ${outdir}/01_init_QC/reads_split.fastq.gz

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
  -g "N{8}${twostep_fwd};max_error_rate=0.25...\${twostep_rev_rc}N{8};max_error_rate=0.25" \
  ${outdir}/01_init_QC/reads_split.fastq.gz |
cutadapt \
  --cores=24 \
  --length=${maxlen} \
  - |
cutadapt \
  --cores=24 \
  -e 0.25 \
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
