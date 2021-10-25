source ~/.bashrc

outdir=/raid1/home/micro/mcmindsr/ryan/20170228_swarm_plus_MED/
nthreads=9

mkdir ${outdir}/swarm/

## dereplicate sequences
vsearch --derep_fulllength ${outdir}/seqs.fna --output ${outdir}/seqs_derep.fasta --uc seqs_derep.uc --relabel_sha1 --relabel_keep --sizeout --threads ${nthreads} --fasta_width 0

## run swarm, with fastidious option (merge rare and similar OTUs) and with breaking
swarm-2.1.12-linux-x86_64 \
-d 1 -z -f -t ${nthreads} \
-w ${outdir}/swarm/swarm_rep_set.fasta \
-l ${outdir}/swarm/swarm.log \
-s ${outdir}/swarm/swarm_statistics.txt \
-o ${outdir}/swarm/swarms.txt \
-u ${outdir}/swarm/swarms.uc \
< ${outdir}/seqs_derep.fasta

## sort swarm reps
vsearch --sortbysize ${outdir}/swarm/swarm_rep_set.fasta --fasta_width 0 --output ${outdir}/swarm/swarm_rep_set_sorted.fasta

## (c) detect and remove chimeras using uchime denovo
vsearch --uchime_denovo ${outdir}/swarm/swarm_rep_set_sorted.fasta --fasta_width 0 --nonchimeras ${outdir}/swarm/swarm_rep_set_sorted_no_chimeras.fasta --threads ${nthreads}

## get all seqs from each OTU above a certain size and that isn't a chimera (code modified from https://github.com/torognes/swarm )
INPUT_SWARM="${outdir}/swarm/swarms.txt"
INPUT_STATS="${outdir}/swarm/swarm_statistics.txt"
INPUT_FASTA="${outdir}/seqs_derep.fasta"
MIN_READS=100
OUTPUT_FOLDER="${outdir}/swarms_fastas"
AMPLICONS="${outdir}/swarms_fasta_temp.txt"
TEMPLATE="${outdir}/swarm/swarm_rep_set_sorted_no_chimeras.fasta"

mkdir "${OUTPUT_FOLDER}"

sort -k2,2nr "${INPUT_STATS}" | \
awk -v MIN_READS="${MIN_READS}" 'NR==FNR {dict[$1]=1} NR!=FNR && $2 >= MIN_READS && $3 in dict {print $3}' \
RS='>' FS=';' "${TEMPLATE}" \
RS='\\n' FS='\\s' - | \
while read SEED ; do
	grep -m 1 "^${SEED}" "${INPUT_SWARM}" | \
	tr " " "\n" | sed -e 's/^/>/' > "${AMPLICONS}"
	grep -A 1 -F -f "${AMPLICONS}" "${INPUT_FASTA}" | \
	sed -e '/^--$/d' > "${OUTPUT_FOLDER}/${SEED}.fasta"
done
rm "${AMPLICONS}"



## (d) remove extremely short reads, which are most likely host mitochondria
vsearch --fastx_filter ${outdir}/swarm/swarm_rep_set_sorted_no_chimeras.fasta --fastq_minlen 220 --fastaout ${outdir}/swarm/swarm_rep_set_sorted_no_chimeras_no_shortreads.fna --fasta_width 0 --threads ${nthreads}
