## six positional arguments specifying 1) input basecalls bam (with duplex basecalling), 2) forward primer sequence, 3) reverse primer sequence, 4) input linked barcodes file, and 5) the output directory
in_bam=$1
primer_fwd=$2
primer_rev=$3
barcodes_file=$4
outdir=$5

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
# extract duplex and simplex reads from dorado basecalling file (excluding simplex parents of duplex reads, dx:-1)
samtools view -b -d dx:1 ${in_bam} | samtools fastq - | gzip --best > ${outdir}/01_init_QC/reads.fastq.gz
samtools view -b -d dx:0 ${in_bam} | samtools fastq - | gzip --best >> ${outdir}/01_init_QC/reads.fastq.gz

# find reverse complement of reverse primer
primer_rev_rc=$(echo ${primer_rev} | tr ACGTRYSWKMBVDHacgtryswkmbvdh TGCAYRSWMKVBHDtgcayrswmkvbhd | rev)

source activate cutadapt-4.6
# find all pairs with both primers plus at least 8 extra bp in forward orientation, trim everything before the adapter
cutadapt \
  --cores=24 \
  --revcomp \
  --action=retain \
  --discard-untrimmed \
  -g "N{8}${primer_fwd};min_overlap=$((${#primer_fwd}+8))...\${primer_rev_rc}N{8};min_overlap=$((${#primer_rev}+8))" \
  --output ${outdir}/01_init_QC/oriented.fastq.gz \
  ${outdir}/01_init_QC/reads.fastq.gz

# demultiplex, allowing a single error in each 8bp index
cutadapt \
  --cores=24 \
  -e 1 \
  --overlap 8 \
  -g file:${barcodes_file} \
  --output ${outdir}/01_init_QC/demultiplexed/{name}.fastq.gz \
  ${outdir}/01_init_QC/oriented.fastq.gz

rm ${outdir}/01_init_QC/demultiplexed/*unknown*

for file in ${outdir}/01_init_QC/demultiplexed/*.fastq.gz; do

  # trim extension from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/.fastq.gz/}
  
  # trim primers
  source activate cutadapt-4.6
  cutadapt \
    --cores=24 \
    -g ${primer_fwd}...\${primer_rev_rc} \
    --output ${outdir}/01_init_QC/trimmed/\${sampleid}.fastq.gz \
    \${file}
  
done

## may want to add logic to parse duplex/simplex reads such that:
## 1. if a duplex read exists that is shorter than the simplex template, stitch the overhang onto the duplex read and use only that result
## 2. if no duplex read exists, check that adjacent reads aren't simply reverse complements that didn't make the duplex command's thresholds (eg if it's less than 50% of the template, the second strand probably won't be incorporated as a duplex read, but we don't want to count it as a separate read; we would want to ignore it)

EOF

if $autorun; then

sbatch ${outdir}/01_init_QC/01_init_QC.sbatch

fi
