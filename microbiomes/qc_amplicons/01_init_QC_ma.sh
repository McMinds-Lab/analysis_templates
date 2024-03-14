## five positional arguments specifying: 
## 1) directory with input fastqs
## 2) forward primer sequence
## 3) reverse primer sequence
## 4) output directory
## 5) number of threads
## 6) lower quality score threshold for iterative trimming and merging
## 7) upper quality score threshold for iterative trimming and merging
## 8) increment for iterative trimming and merging
## 9) quality threshold to keep and concatenate unmerged read pairs (as fraction of bases expected to be wrong in a read pair, including 0 to mean do not concatenate or 1 to mean concatenate all, or 'fwd' to mean just keep the fwd read) 

indir=$1
primer_fwd=$2
primer_rev=$3
outdir=$4
nthreads=$5
thresh_low=$6
thresh_hi=$7
thresh_inc=$8
maxee_rate=$9

mkdir -p ${outdir}/01_init_QC
echo "bash $0 $@" > ${outdir}/01_init_QC/this_command.sh

cat <<EOF > ${outdir}/01_init_QC/01_init_QC.sbatch
#!/bin/bash
#SBATCH --job-name=01_init_QC
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --output=${outdir}/01_init_QC/01_init_QC.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=20G
#SBATCH --time=24:00:00

mkdir -p ${outdir}/01_init_QC/trimmed
mkdir -p ${outdir}/01_init_QC/merged

# trim primers from sequences and discard any sequences that don't contain both full primers
module purge
module load hub.apps/anaconda3

for file in ${indir}/*_R1_001.fastq.gz; do

  # trim "R1" from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/_L001_R1*/}
  
  # trim primers and keep only the reads with the primers at exactly the beginning
  source activate cutadapt-4.6
  cutadapt \
    -e 0.25 \
    --cores=${nthreads} \
    -g ^${primer_fwd} \
    -G ^${primer_rev} \
    --output ${outdir}/01_init_QC/trimmed/\${sampleid}_R1.fastq.gz \
    --paired-output ${outdir}/01_init_QC/trimmed/\${sampleid}_R2.fastq.gz \
    --discard-untrimmed \
    \${file} \
    \${file/R1/R2}
  
  # merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)
  source activate vsearch
  vsearch \
    --threads ${nthreads} \
    --fastq_mergepairs ${outdir}/01_init_QC/trimmed/\${sampleid}_R1.fastq.gz \
    --reverse ${outdir}/01_init_QC/trimmed/\${sampleid}_R2.fastq.gz \
    --fastq_allowmergestagger \
    --fastq_maxdiffs 1000 \
    --fasta_width 0 \
    --fastq_maxns 0 \
    --fastqout_notmerged_fwd ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
    --fastqout_notmerged_rev ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq \
    --fastqout - | gzip --best > ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
  
  cp ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq ${outdir}/01_init_QC/merged/\${sampleid}_latest_R1.fastq
  cp ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq ${outdir}/01_init_QC/merged/\${sampleid}_latest_R2.fastq
  
  for thresh in {${thresh_low}..${thresh_hi}..${thresh_inc}}; do
  
    echo "Retrying merge after quality trimming with threshold \${thresh}\n"
    
    # quality trim reads that did not merge. 
    vsearch \
      --fastx_filter ${outdir}/01_init_QC/merged/\${sampleid}_latest_R1.fastq \
      --reverse ${outdir}/01_init_QC/merged/\${sampleid}_latest_R2.fastq \
      --fastq_truncqual \${thresh} \
      --fastqout ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_Q\${thresh}t_R1.fastq \
      --fastqout_rev ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_Q\${thresh}t_R2.fastq
  
    # retry merging reads that were quality trimmed
    vsearch \
      --threads ${nthreads}\
      --fastq_mergepairs ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_Q\${thresh}t_R1.fastq \
      --reverse ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_Q\${thresh}t_R2.fastq \
      --fastq_allowmergestagger \
      --fastq_maxdiffs 1000 \
      --fasta_width 0 \
      --fastq_maxns 0 \
      --fastqout_notmerged_fwd ${outdir}/01_init_QC/merged/\${sampleid}_tmp_R1.fastq \
      --fastqout_notmerged_rev ${outdir}/01_init_QC/merged/\${sampleid}_tmp_R2.fastq \
      --fastqout - | gzip --best >> ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
      
      mv -f ${outdir}/01_init_QC/merged/\${sampleid}_tmp_R1.fastq ${outdir}/01_init_QC/merged/\${sampleid}_latest_R1.fastq
      mv -f ${outdir}/01_init_QC/merged/\${sampleid}_tmp_R2.fastq ${outdir}/01_init_QC/merged/\${sampleid}_latest_R2.fastq
      
  done
  
  if [ "${maxee_rate}" = 0 ]; then
  
    echo 'discarding unmerged reads\n'

  else 
  
    zcat ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz | awk 'NR%4==1 {print \$1}' | cut -c 2- > ${outdir}/01_init_QC/merged/\${sampleid}_merged_labels.txt
  
    if [ "${maxee_rate}" = "fwd" ]; then
  
      echo 'adding unmerged forward reads\n'
        
      source activate vsearch
      vsearch \
        --fastx_getseqs ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
        --label_substr_match \
        --labels ${outdir}/01_init_QC/merged/\${sampleid}_merged_labels.txt \
        --notmatchedfq ${outdir}/01_init_QC/merged/\${sampleid}_tmp.fastq.gz
      
      source activate cutadapt-4.6
      cutadapt \
        -e 0.5 \
        --minimum-length 100
        --cores=${nthreads} \
        -a $(tr ACGTRYSWKMBVDHacgtryswkmbvdh TGCAYRSWMKVBHDtgcayrswmkvbhd <<< ${primer_rev} | rev) \
        ${outdir}/01_init_QC/merged/\${sampleid}_tmp.fastq.gz |
      gzip --best >> ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
  
      rm ${outdir}/01_init_QC/merged/\${sampleid}_tmp.fastq.gz

    else

      echo "adding concatenated reads if they didnt merge but are good quality\n"
    
      vsearch \
        --fastx_getseqs ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
        --label_substr_match \
        --labels ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_labels.txt \
        --notmatchedfq ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_final_R1.fastq
    
      vsearch \
        --fastx_getseqs ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq \
        --label_substr_match \
        --labels ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_labels.txt \
        --notmatchedfq ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_final_R2.fastq
    
      vsearch \
        --fastx_filter ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_final_R1.fastq \
        --reverse ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_final_R2.fastq \
        --fastq_maxee_rate ${maxee_rate} \
        --fastqout ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_filt_R1.fastq \
        --fastqout_rev ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_filt_R2.fastq
    
      vsearch \
        --fastq_join ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_filt_R1.fastq \
        --reverse ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_filt_R2.fastq \
        --join_padgap '' \
        --fastqout - | gzip --best >> ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
  
    fi
  
  fi
  
done

rm ${outdir}/01_init_QC/merged/*unmerged*
rm ${outdir}/01_init_QC/merged/*latest*

EOF

if $autorun; then

sbatch ${outdir}/01_init_QC/01_init_QC.sbatch

fi
