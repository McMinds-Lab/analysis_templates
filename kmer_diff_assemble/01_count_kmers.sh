## inspired in part by https://github.com/Transipedia/dekupl-joinCounts/tree/146604e1debee6759b3f7eeb604d031c5d4113b3

## two positional arguments specifying 1) the directory containing fastqs in samplewise subfolders a la download_SRA.sh, and 2) the output analysis directory
# may be good idea to add input of sample sheet which could be used to select the specific files from the input? also add QC prior to kmer counting?
indir=$1
outdir=$2
n_processes=$3
n_threads=$4

subdir=${outdir}/01_jellyfish

mkdir -p ${subdir}/logs
mkdir ${subdir}/counts
mkdir ${subdir}/trimmed

samples=($(grep SRR ${indir}/runInfo.csv | grep 'WGA\|WGS\|RNA-Seq' | cut -d ',' -f 1))

cat <<EOF > ${subdir}/01a_jellyfish.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=160G
#SBATCH --ntasks=${n_threads}
#SBATCH --job-name=01a_jellyfish
#SBATCH --output=${subdir}/logs/01_jellyfish_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%${n_processes}

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]} ## each array job has a different sample

in1=(${indir}/*/\${sample}/\${sample}_1.fastq.gz)
in2=(${indir}/*/\${sample}/\${sample}_2.fastq.gz)

## set up named pipes instead of temporary intermediate files
mkdir -p ${subdir}/temp
pipe1=${subdir}/temp/\${sample}_p1.fastq
pipe2=${subdir}/temp/\${sample}_p2.fastq
pipe3=${subdir}/temp/\${sample}_p3.tsv

mkfifo \${pipe1} \${pipe2} \${pipe3}

## conda seems to need extra help loading this package...
module purge
module load apps/jellyfish/2.2.6
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
source activate bbtools

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
  threads=${n_threads} \
  overwrite=true \
  in1=\${in1} \
  in2=\${in2} \
  out1=\${pipe1} \
  out2=\${pipe2} \
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
  ecco &

## although BBDuk is documented to gzip its output, it hasn't been consistent for me in the past
pigz < \${pipe1} > ${subdir}/trimmed/\${sample}_1.fastq.gz &
pigz < \${pipe2} > ${subdir}/trimmed/\${sample}_2.fastq.gz

#paired-end, unstranded data. theoretically I think in future I would like to first merge reads, then quality-control them (incl. adapter and quality trimming), then feed three files to jellyfish for each sample (merged, unmerged R1, unmerged R2). Downstream would need to be able to handle that new format
jellyfish count \
  -t ${n_threads} \
  -m 31 \
  -s 10000 \
  -o \${pipe3} \
  -C \
  <(zcat ${subdir}/trimmed/\${sample}_1.fastq.gz) \
  <(zcat ${subdir}/trimmed/\${sample}_2.fastq.gz) &

## sort doesn't actually parallelize for pipes unless you make the buffer big (S5G)
jellyfish dump \
  -c \${pipe3} |
  sort -S5G --parallel ${n_threads} -k 1 |
  pigz > ${subdir}/counts/\${sample}.tsv.gz

rm -f \${pipe1} \${pipe2} \${pipe3}

EOF

cat <<EOF > ${subdir}/01b_merge.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=160G
#SBATCH --ntasks=${n_threads}
#SBATCH --job-name=01b_merge
#SBATCH --output=${subdir}/logs/01b_merge.log

samples=(${samples[@]})
prefix=${subdir}/counts/

files=(\${samples[@]/#/\${prefix}})
files=(\${files[@]/%/'.tsv.gz'})

multi_join() {

  # unzip first arg if it's the top level of the recursion (e.g. its a gzipped file not a raw pipe)
  if [ -z \${in_recursion} ]; then
    f1=<(zcat \$1)
    export in_recursion=1
  else
    f1=\$1
  fi
  f2=\$2
  shift 2

  # if the number of remaining files is greater than two, continue the recursion; else join these two and end
  if [ \$# -gt 0 ]; then
    join -o auto -a 1 -a 2 -e 0 \${f1} <(zcat \${f2}) | multi_join - \$@
  else
    join -o auto -a 1 -a 2 -e 0 \${f1} <(zcat \${f2})
  fi

}

## concatenate a header with sample names to a body formed by joining the samplewise jellyfish outputs
cat \
  <(printf '%s\t' 'kmer' \${samples[@]} | sed 's/\t$/\n/g') \
  <(multi_join \${files[@]} | tr ' ' '\t') |
  pigz > ${subdir}/counts_matrix.tsv.gz

EOF

## run the samplewise commands as a slurm array, capturing the jobid, and run the merging code after waiting for the entire array to finish
arrayID=$(sbatch ${subdir}/01a_jellyfish.sbatch | cut -d' ' -f 4)
echo 'Submitted batch job' ${arrayID}
sbatch --dependency=afterany:${arrayID} ${subdir}/01b_merge.sbatch
