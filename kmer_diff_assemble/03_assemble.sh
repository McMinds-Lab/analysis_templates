## identify reads that contain significant kmers
indir=$1
in_kmers=$2
outdir=$3
sample_regex_1=$4
sample_regex_2=$5

in1=${outdir}/in_1.fastq.gz
in2=${outdir}/in_2.fastq.gz

cat $(find ${indir}/${sample_regex_1}) > ${in1}
cat $(find ${indir}/${sample_regex_2}) > ${in2}

zcat ${in_kmers} | tail -n +2 | awk '$5 > 0' | cut -f 1 > ${outdir}/kmers.txt

## search for both kmer itself and its reverse complement, in both reads. save matching lines with line numbers in a single file
tr ACGTacgt TGCAtgca < ${in_kmers} | rev > ${outdir}/kmers_rc.txt
zcat ${in1} | grep -F -n -B 1 -A 2 -f ${outdir}/kmers.txt -f ${outdir}/kmers_rc.txt > ${outdir}/keepers.txt
zcat ${in2} | grep -F -n -B 1 -A 2 -f ${outdir}/kmers.txt -f ${outdir}/kmers_rc.txt >> ${outdir}/keepers.txt

## save line numbers of matches into variable idx; "-" delimiter for lines before and after match; ":" delimiter for line with match itself.
awk -F '[-:]' 'NR==FNR {idx[$1]; next} FNR in idx' ${outdir}/keepers.txt <(zcat ${in1}) | gzip > "${outdir}/filt_1.fastq.gz"
awk -F '[-:]' 'NR==FNR {idx[$1]; next} FNR in idx' ${outdir}/keepers.txt <(zcat ${in2}) | gzip > "${outdir}/filt_2.fastq.gz"

## assemble reads with megahit. same for RNA or DNA? better newer tools now?
module purge
module load apps/megahit/1.2.9/
megahit \
  --k-step 10 \
  -1 ${outdir}/filt_1.fastq.gz \
  -2 ${outdir}/filt_2.fastq.gz \
  -o ${outdir}/assembled_genome
