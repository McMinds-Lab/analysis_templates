## identify reads that contain significant kmers
indir=$1
in_kmers=$2
outdir=$3
sample_regex_1=$4
sample_regex_2=$5
sample_regex_m=$6

in1=${outdir}/in_1.fastq.gz
in2=${outdir}/in_2.fastq.gz
inm=${outdir}/in_m.fastq.gz

cat $(find ${indir}/${sample_regex_1}) > ${in1}
cat $(find ${indir}/${sample_regex_2}) > ${in2}
cat $(find ${indir}/${sample_regex_m}) > ${inm}

zcat ${in_kmers} | tail -n +2 | awk '$5 > 0' | cut -f 1 > ${outdir}/kmers.txt

## search for both kmer itself and its reverse complement, in both reads. save matching lines with line numbers in a single file
tr ACGTacgt TGCAtgca < ${in_kmers} | rev > ${outdir}/kmers_rc.txt

## put this group of tasks in a background subshell so the merged filtration can start right away
(
## identify lines that have matches in either forward or reverse reads
zcat ${in1} ${in2} |
  grep -F -n -B 1 -A 2 -f ${outdir}/kmers.txt -f ${outdir}/kmers_rc.txt > ${outdir}/keepers_unmerged.txt

## use line nnumbers to filter unmerged reads
## save line numbers of matches into variable idx; "-" delimiter for lines before and after match; ":" delimiter for line with match itself.
awk -F '[-:]' 'NR==FNR {idx[$1]; next} FNR in idx' ${outdir}/keepers_unmerged.txt <(zcat ${in1}) | pigz > ${outdir}/filt_1.fastq.gz &
awk -F '[-:]' 'NR==FNR {idx[$1]; next} FNR in idx' ${outdir}/keepers_unmerged.txt <(zcat ${in2}) | pigz > ${outdir}/filt_2.fastq.gz
) &

## filter merged reads directly
zcat ${inm} |
  grep -F -B 1 -A 2 -f ${outdir}/kmers.txt -f ${outdir}/kmers_rc.txt |
  pigz > ${outdir}/filt_m.fastq.gz

## wait for all the background tasks to finish
wait

rm -f ${in1} ${in2} ${inm}

## assemble reads with megahit. same for RNA or DNA? better newer tools now?
module purge
module load apps/megahit/1.2.9/
megahit \
  --k-step 10 \
  -1 ${outdir}/filt_1.fastq.gz \
  -2 ${outdir}/filt_2.fastq.gz \
  -r ${outdir}/filt_m.fastq.gz \
  -o ${outdir}/assembled_genome
