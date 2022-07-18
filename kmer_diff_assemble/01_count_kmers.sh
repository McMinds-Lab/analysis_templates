## inspired in part by https://github.com/Transipedia/dekupl-joinCounts/tree/146604e1debee6759b3f7eeb604d031c5d4113b3

## two positional arguments specifying 1) the directory containing fastqs in samplewise subfolders a la download_SRA.sh, and 2) the output analysis directory
# may be good idea to add input of sample sheet which could be used to select the specific files from the input? also add QC prior to kmer counting?
indir=$1
outdir=$2
n_processes=$3
n_threads=$4

mkdir -p ${outdir}/01_jellyfish/logs
mkdir ${outdir}/01_jellyfish/counts
mkdir ${outdir}/01_jellyfish/trimmed

samples=($(grep SRR ${indir}/runInfo.csv | grep 'WGA\|WGS\|RNA-Seq' | cut -d ',' -f 1))

cat <<EOF > ${outdir}/01_jellyfish/01a_jellyfish.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=160G
#SBATCH --job-name=01a_jellyfish
#SBATCH --output=${outdir}/01_jellyfish/logs/01_jellyfish_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%${n_processes}

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]} ## each array job has a different sample

in1=(${indir}/*/\${sample}/\${sample}_1.fastq.gz)
in2=(${indir}/*/\${sample}/\${sample}_2.fastq.gz)

mkdir ${outdir}/01_jellyfish/temp
pipe1_1=${outdir}/01_jellyfish/temp/\${sample}_p1gz
pipe1_2=${outdir}/01_jellyfish/temp/\${sample}_p1jf
pipe2_1=${outdir}/01_jellyfish/temp/\${sample}_p2gz
pipe2_2=${outdir}/01_jellyfish/temp/\${sample}_p2jf
pipe3=${outdir}/01_jellyfish/temp/\${sample}_p3

mkfifo \$pipe1_1 \$pipe1_2 \$pipe2_1 \$pipe2_2 \$pipe3

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
  overwrite=true \
  in1=\${in1} \
  in2=\${in2} \
  out1=\$pipe1_1 \
  out2=\$pipe2_1 \
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

cat \$pipe1_1 | tee \$pipe1_2 | gzip > ${outdir}/01_jellyfish/trimmed/\${sample}_1.fastq.gz &
cat \$pipe2_1 | tee \$pipe2_2 | gzip > ${outdir}/01_jellyfish/trimmed/\${sample}_2.fastq.gz &

#paired-end, unstranded data. theoretically I think I would like to first merge reads, then quality-control them (incl. adapter and quality trimming), then feed three files to jellyfish for each sample (merged, unmerged R1, unmerged R2). Downstream would need to be able to handle that new format
jellyfish count -t ${n_threads} -m 31 -s 10000 -o \$pipe3 -C \$pipe1_2 \$pipe2_2 &
jellyfish dump -c \$pipe3 | sort --parallel ${n_threads} -k 1 | gzip > ${outdir}/01_jellyfish/counts/\${sample}.tsv.gz

rm -rf ${outdir}/01_jellyfish/temp

EOF

cat <<EOF > ${outdir}/01_jellyfish/01b_merge.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=160G
#SBATCH --job-name=01b_merge
#SBATCH --output=${outdir}/01_jellyfish/logs/01b_merge.log

samples=(${samples[@]})
prefix=${outdir}/01_jellyfish/counts/

files=(\${samples[@]/#/\${prefix}})
files=(\${files[@]/%/'.tsv.gz'})

join_rec() {

  # first arg should be one when function is initially called; internally it is set to zero
  # unzip first arg if it's not being piped
  if [ \$1 -eq 1 ]; then
    f1=<(zcat \$2)
  else
    f1=\$2
  fi
  f2=\$3
  shift 3

  # if the number of remaining files is greater than two, continue the recursion; else join these two and end
  if [ \$# -gt 0 ]; then
    join -o auto -a 1 -a 2 -e 0 "\$f1" <(zcat "\$f2") | join_rec 0 - "\$@"
  else
    join -o auto -a 1 -a 2 -e 0 "\$f1" <(zcat "\$f2")
  fi

} #https://stackoverflow.com/questions/10726471/join-multiple-files

cat <(printf "\$(printf '%s\t' 'kmer' "\${samples[@]}")\n") <(join_rec 1 \${files[@]}) | tr ' ' '\t') | gzip > ${outdir}/01_jellyfish/counts_matrix.tsv.gz

EOF

arrayID=$(sbatch ${outdir}/01_jellyfish/01a_jellyfish.sbatch | cut -d' ' -f 4)
sbatch --dependency=afterany:$arrayID ${outdir}/01_jellyfish/01b_merge.sbatch
