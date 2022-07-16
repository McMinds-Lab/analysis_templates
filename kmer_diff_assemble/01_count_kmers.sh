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

out=${outdir}/01_jellyfish/counts/\${sample}.tsv

module purge
module load apps/jellyfish/2.2.6

#paired-end, unstranded data
jellyfish count -t ${n_threads} -m 31 -s 10000 -o \${out/.tsv/.jf} -C <(zcat \${in1}) <(zcat \${in2})
jellyfish dump -c \${out/.tsv/.jf} | sort --parallel ${n_threads} -k 1 > \${out}

rm \${out/.tsv/.jf}

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

join -o auto -a 1 -a 2 -e '0' ${outdir}/01_jellyfish/counts/\${samples[0]}.tsv ${outdir}/01_jellyfish/counts/\${samples[1]}.tsv > ${outdir}/01_jellyfish/temp_1.tsv

for i in \$(seq 2 \$((\${#samples[@]}-1))); do

join -o auto -a 1 -a 2 -e '0' ${outdir}/01_jellyfish/temp_\$(($i-1)).tsv ${outdir}/01_jellyfish/counts/\${samples[\$i]}.tsv > ${outdir}/01_jellyfish/temp_\${i}.tsv
rm ${outdir}/01_jellyfish/temp_\$(($i-1)).tsv

done

cat <(printf "\$(printf '%s\t' 'kmer' "\${samples[@]}")\n" ${outdir}/01_jellyfish/temp_\$((\${#samples[@]}-1)).tsv > ${outdir}/01_jellyfish/counts_matrix.tsv
rm ${outdir}/01_jellyfish/temp_\$((\${#samples[@]}-1)).tsv

EOF

arrayID=$(sbatch ${outdir}/01_jellyfish/01a_jellyfish.sbatch | cut -d' ' -f 4)
sbatch --dependency=afterany:$arrayID ${outdir}/01_jellyfish/01b_merge.sbatch
