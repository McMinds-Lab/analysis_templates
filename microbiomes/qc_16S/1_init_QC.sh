source ~/.bashrc

outdir=/raid1/home/micro/mcmindsr/ryan/20170228_swarm_plus_MED2/
indir=/raid1/home/micro/mcmindsr/labhome/web_downloads/GCMP/data/16S/raw_data/fastqs/
nthreads=1
fwdprimer=GTGYCAGCMGCCGCGGTAA
revprimer=GGACTACNVGGGTWTCTAAT

mkdir -p ${outdir}/trimmed/
mkdir -p ${outdir}/filtered/

for file in ${indir}/*_R1*fastq.gz; do

    # trim "R1" from filenames to get Sample IDs that match mapping file
    filename=$(basename $file)
    sampleid=${filename/_*/}

	## trim primers from sequences, and discard any sequences that don't contain both full primers exactly at the beginnings of the reads.
    # (a) merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)
    ## (b) filter merged reads based on 'expected errors' and a minimum read length. This function determines the likelihood that each read has an error based on all the (merged) quality scores. If even one error is expected to occur within the read (including even the presence of a single 'N'), it throws out the entire read. No trimming occurs, which is theoretically better for downstream clustering and denoising applications. This command also changes the fasta headers in the output so they fit the format ">sampleID_arbitraryreadnumber"
	cutadapt --discard-untrimmed --pair-filter=any -g "^${fwdprimer}" -G "^${revprimer}" -o ${outdir}/trimmed/${sampleid}_R1.fastq -p ${outdir}/trimmed/${sampleid}_R2.fastq ${file} ${file/R1/R2}

	vsearch --fastq_mergepairs ${outdir}/trimmed/${sampleid}_R1.fastq --reverse ${outdir}/trimmed/${sampleid}_R2.fastq --fastq_allowmergestagger --fasta_width 0 --threads ${nthreads} --fastqout - |

	vsearch --fastq_filter - --fastq_maxns 1 --fastq_maxee 1 --fastq_ascii 33 --fasta_width 0 --relabel ${sampleid}_ --threads ${nthreads} --fastaout ${outdir}/filtered/${sampleid}.fasta

	rm ${outdir}/trimmed/${sampleid}_R1.fastq ${outdir}/trimmed/${sampleid}_R2.fastq

done

cat ${outdir}/filtered/*.fasta > ${outdir}/seqs.fna

rm -r ${outdir}/trimmed/ ${outdir}/filtered/