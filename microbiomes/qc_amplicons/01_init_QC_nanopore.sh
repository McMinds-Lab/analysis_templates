## six positional arguments specifying 1) input basecalls bam (with duplex basecalling), 2) forward primer sequence, 3) reverse primer sequence, 4) input fwd barcodes file, 5) input rev barcodes file, and 6) the output directory
in_bam=$1
primer_fwd=$2
primer_rev=$3
barcodes_fwd=$4
barcodes_rev=$5
outdir=$6

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

source activate SAMTOOLS
# extract duplex and simplex reads from dorado basecalling file (excluding simplex parents of duplex reads, dx:-1)
samtools view -d dx:1 ${in_bam} | samtools fastq -@ 24 -c 9 > ${outdir}/01_init_QC/reads.fastq.gz
samtools view -d dx:0 ${in_bam} | samtools fastq -@ 24 -c 9 >> ${outdir}/01_init_QC/reads.fastq.gz

source activate cutadapt-4.6

# find all pairs with both primers plus at least 8 extra bp in forward orientation, trim everything before the adapter
cutadapt \
  --cores=${nthreads} \
  --revcomp \
  --action=retain \
  -g "N{8}${primer_fwd};min_overlap=$((${#primer_fwd}+8))" \
  -G "N{8}${primer_rev};min_overlap=$((${#primer_rev}+8))" \
  --output ${outdir}/01_init_QC/oriented.fastq.gz \
  ${outdir}/01_init_QC/reads.fastq.gz

# demultiplex, allowing a single error in each 8bp index
cutadapt \
  --cores=${nthreads} \
  -e 1 \
  --overlap 8 \
  -g file:${barcodes_fwd} \
  -G file:${barcodes_rev} \
  --output ${outdir}/01_init_QC/demultiplexed/{name1}-{name2}.fastq.gz \
  ${outdir}/01_init_QC/oriented.fastq.gz

rm ${outdir}/01_init_QC/demultiplexed/*unknown*

for file in ${outdir}/01_init_QC/demultiplexed/*.fastq.gz; do

  # trim extension from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/.fastq.gz/}
  
  # trim primers
  source activate cutadapt-4.6
  cutadapt \
    --cores=${nthreads} \
    -g ${primer_fwd} \
    -G ${primer_rev} \
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
