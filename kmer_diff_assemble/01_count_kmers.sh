## inspired in part by https://github.com/Transipedia/dekupl-joinCounts/tree/146604e1debee6759b3f7eeb604d031c5d4113b3

## two positional arguments specifying 1) the directory containing fastqs in samplewise subfolders a la download_SRA.sh, and 2) the output analysis directory
# may be good idea to add input of sample sheet which could be used to select the specific files from the input? also add QC prior to kmer counting?
indir=$1
outdir=$2
n_processes=$3
n_threads=$4

mkdir -p ${outdir}/01_jellyfish/logs
mkdir ${outdir}/01_jellyfish/counts

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

out=${outdir}/01_jellyfish/counts/\${sample}

module purge
module load apps/jellyfish/2.2.6

#paired-end, unstranded data
jellyfish count -t ${n_threads} -m 31 -s 10000 -o \${out}.jf -C <(zcat \${in1}) <(zcat \${in2})
jellyfish dump -c \${out}.jf | sort --parallel ${n_threads} -k 1 | gzip > \${out}.tsv.gz

rm \${out}.jf

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

files=(\${samples[@]/#/\'${outdir}/01_jellyfish/counts/\'})
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
