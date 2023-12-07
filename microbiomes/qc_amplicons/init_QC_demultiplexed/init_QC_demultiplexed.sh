## five positional arguments specifying 1) input fwd fastq, 2) input rev fastq, 3) the output directory, 4) either a path SLURM parameters file or the word 'false' (to use all default slurm parameters, you can specify a path to an empty script. if false, script will be run with 'bash' rather than 'sbatch'), and 5) optional path to a script that will put cutadapt and vsearch in your PATH
in_fwd=$1
in_rev=$2
outdir=$3
slurm_params=$4

## input error checking
if [ $# -lt 4 ]; then
  echo "First seven positional arguments are required"
fi

## if tools are not automatically in path, put them there
if [ $# -ge 5 ]; then
  prepare_path=$5
  source ${prepare_path}
fi

## create a copy of the command used to run this script
mkdir -p ${outdir}/01_init_QC
echo "bash $0 $@" > ${outdir}/01_init_QC/this_command.sh

## if slurm parameters are set, create a file with an sbatch header
if [ ${slurm_params} != "false" ]; then

  scriptname=${outdir}/01_init_QC/01_init_QC.sbatch

  # get slurm parameters
  source ${slurm_params}
  
  cat <<__EOF | sed -e 's/^  //' > ${scriptname}
  #!/bin/bash
  #SBATCH --job-name=01_init_QC
  #SBATCH --partition=${partition:-rra}
  #SBATCH --qos=${qos:-rra}
  #SBATCH --output=${outdir}/01_init_QC/01_init_QC.log
  #SBATCH --ntasks=${nthreads:-1}
  #SBATCH --mem=${mem:-20G}
  #SBATCH --time=${time:-24:00:00}

__EOF

else

  scriptname=${outdir}/01_init_QC/01_init_QC.sh

fi

## add commands to the sbatch or shell script
cat <<EOF >> ${scriptname}

mkdir -p ${outdir}/01_init_QC/merged

for file in ${in_fwd}; do

  # trim "R1" from filenames to get Sample IDs that match mapping file
  filename=\$(basename \$file)
  sampleid=\${filename/_R1*/}
  
  # trim primers
  cutadapt \
    --cores=${nthreads} \
    -g ${primer_fwd} \
    -G ${primer_rev} \
    --output ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R1.fastq.gz \
    --paired-output ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R2.fastq.gz \
    \${file} \
    \${file/R1/R2}
  
  # merge paired-end reads such that short reads, where the read is longer than the insertion (such as mitochondria), are not discarded, and nucleotides are trimmed that extend past the beginning of the paired read (which are just adaptor sequences)
  vsearch \
    --threads ${nthreads}\
    --fastq_mergepairs ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R1.fastq.gz \
    --reverse ${outdir}/01_init_QC/demultiplexed/\${sampleid}_trimmed_R2.fastq.gz \
    --fastq_allowmergestagger \
    --fastq_maxdiffs 100 \
    --fasta_width 0 \
    --fastq_maxns 0 \
    --fastqout_notmerged_fwd ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
    --fastqout_notmerged_rev ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq \
    --fastqout - | gzip --best > ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz

  # add concatenated fwd and reverse reads that could not be merged, assuming that they were not merged because they had a gap between reads
  vsearch \
    --fastq_join ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R1.fastq \
    --reverse ${outdir}/01_init_QC/merged/\${sampleid}_unmerged_R2.fastq \
    --join_padgap '' \
    --fastqout - |
  vsearch \
    --fastx_filter - \
    --fastq_maxns 0 \
    --fastqout - |
  gzip --best >> ${outdir}/01_init_QC/merged/\${sampleid}.fastq.gz
  
done

rm ${outdir}/01_init_QC/merged/*unmerged*

EOF

## if slurm parameters were specified, then use slurm, otherwise run as bash script
if [ ${slurm_params} != "false" ]; then

  sbatch ${scriptname}

else

  bash ${scriptname}

fi
