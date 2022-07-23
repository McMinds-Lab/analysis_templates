## calculate normalization factors just like (or even using tools from) de-kupl
## run rscript with custom statistics to identify significantly different kmers

in_kmers=$1
outdir=$2
n_processes=$3
n_threads=$4
maxmem=$5
nbatches=$6
threshold=$7
sampledat=$8
formula=$9
keycolumn=${10}

subdir=${outdir}/02_id_diffs

mkdir -p ${subdir}/logs
mkdir ${subdir}/temp

## get samplewise offsets and latent factors
cat <<EOF > ${subdir}/02a_id_diffs_zinbwave.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --mem=${maxmem}
#SBATCH --time=2-00:00:00
#SBATCH --job-name=02a_id_diffs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${n_threads}
#SBATCH --output=${subdir}/logs/02a_id_diffs_zinbwave.log

pipe1=${subdir}/temp/p1
mkfifo \${pipe1}

## split input into chunks before feeding to R.
#first shuffle lines. skip header for shuffling, then pipe result to both a single file that adds back the header, and to split, which adds the header to each result
#if i want to do any filtering it might be good to do before this?
nlines=\$(zcat ${in_kmers} | wc -l)
paste <(shuf -i2-\${nlines}) <(zcat ${in_kmers} | tail -n +2) |
  sort -S20G --parallel ${n_threads} -T ${subdir}/temp --compress-program pigz -k 1 -n |
  cut -f2- |
  tee \${pipe1} |
  cat <(zcat ${in_kmers} | head -1) - |
  pigz > ${in_kmers/.tsv.gz/_shuffled.tsv.gz} &

split -d \
  -l \$((\$nlines / $nbatches)) \
  --filter='cat <(zcat ${in_kmers} | head -1) - | gzip > \$FILE' \
  --additional-suffix='.gz' \
  \${pipe1} \
  ${subdir}/temp/chunk_

module purge
module load hub.apps/R/4.1.1

Rscript 02a_id_diffs_zinbwave.r ${subdir}/temp/chunk_ ${sampledat} ${threshold} ${formula} ${keycolumn} ${outdir}

EOF

## get significance for each gene independently
cat <<EOF > ${subdir}/02b_id_diffs_zeroinfl.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --mem=${maxmem}
#SBATCH --time=7-00:00:00
#SBATCH --job-name=02b_id_diffs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${n_threads}
#SBATCH --output=${subdir}/logs/02b_id_diffs_zeroinfl_%a.log
#SBATCH --array=0-$((${#nbatches[@]}))%${n_processes}

module purge
module load hub.apps/R/4.1.1

Rscript 02b_id_diffs_zeroinfl.r ${outdir} \${SLURM_ARRAY_TASK_ID}

EOF

## collate results and clean up
cat <<EOF > ${subdir}/02c_id_diffs_collate.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --job-name=02c_id_diffs_collate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=${subdir}/logs/02c_id_diffs_collate.log

cat ${outdir}/02_id_diffs/temp/results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_results.tsv.gz &
cat ${outdir}/02_id_diffs/temp/sig_results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_sig_results.tsv.gz &
cat ${outdir}/02_id_diffs/temp/pos_sig_results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_pos_sig_results.tsv.gz &
cat ${outdir}/02_id_diffs/temp/neg_sig_results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_neg_sig_results.tsv.gz &

wait

rm -rf ${outdir}/02_id_diffs/temp/

EOF

arrayID1=$(sbatch ${subdir}/02a_id_diffs_zinbwave.sbatch | cut -d' ' -f 4)
echo 'Submitted batch job' ${arrayID1}
arrayID2=$(sbatch --dependency=afterany:${arrayID1} ${subdir}/02b_id_diffs_zeroinfl.sbatch | cut -d' ' -f 4)
echo 'Submitted batch job' ${arrayID2}
sbatch --dependency=afterany:${arrayID2} ${subdir}/02c_id_diffs_collate.sbatch
