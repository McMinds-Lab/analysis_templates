# get local variables
source local.env

## two positional arguments specifying 1) the directory containing fastqs in samplewise subfolders a la download_SRA.sh, and 2) the output analysis directory
indir=$1
outdir=$2

mkdir -p ${outdir}/01_BBDuk/logs
mkdir ${outdir}/01_BBDuk/trimmed

samples=($(grep SRR ${indir}/runInfo.csv | grep 'WGA\|WGS\|RNA-Seq' | cut -d ',' -f 1))

cat <<EOF > ${outdir}/01_BBDuk/01_BBDuk.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=${maxram}
#SBATCH --job-name=01_BBDuk
#SBATCH --output=${outdir}/01_BBDuk/logs/01_BBDuk_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%12

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]} ## each array job has a different sample

module purge
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
source activate bbtools

in1=(${indir}/*/\${sample}/\${sample}_1.fastq.gz)
in2=(${indir}/*/\${sample}/\${sample}_2.fastq.gz)

out1=${outdir}/01_BBDuk/trimmed/\${sample}_1.fastq
out2=${outdir}/01_BBDuk/trimmed/\${sample}_2.fastq

## qtrim = trim the 3' ends of reads based on quality scores
## ktrim = trim both 3' ends of reads based on matches to sequencing adapters and artifacts
## k = use 23-mers to identify adapters and artifacts
## mink = don't look for adapter matches below this size
## hdist = allow two mismatches for adapter detection
## minlength = default length of a single kmer downstream in dekupl; if a read is trimmed shorter than this just discard it
## trimq = trim reads once they reach quality scores of 20 (for de-kupl I think it may pay to be stringent here; maybe even more than 20)
## tbo = trim read overhangs if they completely overlap
## tpe = if kmer trimming happens, trim paired reads to same length
## ecco = perform error-correction using pair overlap
bbduk.sh \
  in1=\${in1} \
  in2=\${in2} \
  out1=\${out1} \
  out2=\${out2} \
  ref=adapters,artifacts \
  qtrim=r \
  ktrim=r \
  k=23 \
  mink=11 \
  hdist=2 \
  minlength=31 \
  trimq=20 \
  ftl=10 \
  tbo \
  tpe \
  ecco

gzip -c \${out1} > \${out1/.fastq/.fastq.gz} &
gzip -c \${out2} > \${out2/.fastq/.fastq.gz}

rm \${out1}
rm \${out2}

EOF

if $autorun; then
    sbatch ${outdir}/01_BBDuk/01_BBDuk.sbatch
fi
