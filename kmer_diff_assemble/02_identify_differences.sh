## calculate normalization factors just like (or even using tools from) de-kupl
## run rscript with custom statistics to identify significantly different kmers

in_kmers=$1
outdir=$2
n_processes=$3
n_threads=$4
maxmem=$5
splitsize=$6
threshold=$7
sampledat=$8
formula=$9
keycolumn=$10

subdir=${outdir}/02_id_diffs

mkdir -p ${subdir}/logs
mkdir ${subdir}/temp

cat <<EOF > ${subdir}/02_id_diffs.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --mem-per-cpu=1G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=02_id_diffs
#SBATCH --ntasks=${n_processes}
#SBATCH --cpus-per-task=${n_threads}

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
  -l ${splitsize} \
  --filter='cat <(zcat ${in_kmers} | head -1) - | gzip > $FILE' \
  --additional-suffix='.gz' \
  \${pipe1} \
  ${subdir}/temp/chunk_

## temporary use of dekupl env because i somehow got Rmpi installed there and can't figure it out in the main hub R module
module purge
module load hub.apps/anaconda3
source activate dekupl
module load compilers/gcc/4.8.1 mpi/openmpi/1.6.1

mpirun -np ${n_processes} Rscript 02_identify_differences.r ${subdir}/temp/chunk_ ${sampledat} ${threshold} ${formula} ${keycolumn} ${outdir}

cat ${outdir}/02_id_diffs/temp/results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_results.tsv.gz
cat ${outdir}/02_id_diffs/temp/sig_results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_sig_results.tsv.gz
cat ${outdir}/02_id_diffs/temp/pos_sig_results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_pos_sig_results.tsv.gz
cat ${outdir}/02_id_diffs/temp/neg_sig_results_chunk_*.tsv.gz > ${outdir}/02_id_diffs/all_neg_sig_results.tsv.gz

rm -rf ${outdir}/02_id_diffs/temp/

EOF

sbatch ${subdir}/02_id_diffs.sbatch
