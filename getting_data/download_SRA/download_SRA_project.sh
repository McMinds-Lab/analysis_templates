## four or five positional arguments specifying 1) an SRA project ID starting with SRP, 2) the output directory, 3) a comma-separated list of acceptable sample types (e.g. 'WGA,AMPLICON') or the word 'all' for no filtering, 4) either a path SLURM parameters file or the word 'true' or the word 'false' (If true, default slurm parameters will be used; if false, script will simply download samples with for loop and won't submit a slurm job), and 5) optional path to a script to set your PATH variable to include the entrez-direct and sratools packages
## make sure you have 'prepare_path.sh' and 'slurm_params.sh' in the same directory as this script, if you need them (example files are included as 'prepare_path.sh.txt' and 'slurm_params.sh.txt')

srp=$1
outdir=$2
accept=$3
slurm_params=$4

## input error checking
if [ $# -lt 4 ]; then
  echo "First four positional arguments are required"
fi

## if tools are not automatically in path, put them there
if [ $# -ge 5 ]; then
  prepare_path=$5
  source ${prepare_path}
fi

## this tells the script where to find this file and download_SRA_sample.sh
scriptdir=$( cd "$(dirname "$0")" ; pwd -P )

## prepare output directory
mkdir -p ${outdir}
cp ${scriptdir}/download_SRA_sample.sh ${outdir}

## download a list of samples, and some of their metadata, associated with the specified project
esearch -db sra -query ${srp} | efetch -format runinfo > ${outdir}/runInfo.csv

## extract sample identifiers
samples=($(awk -v header='LibraryStrategy' -v accept=$(sed 's/^/^/;s/$/$/;s/,/$|^/g' <<< ${accept}) 'BEGIN {FS=","; column=0} NR == 1 {for (i=1;i<=NF;i++) {if ($i==header) {column=i}} } NR > 1 && $column ~ accept {print $1}' ${outdir}/runInfo.csv))

## require the explicit input of 'false' to enable for-loop
if [ ${slurm_params} != "false" ]; then

  ## create a logs directory
  mkdir -p ${outdir}/logs

  if [ ${slurm_params} != "true" ]; then
    # get slurm parameters
    source ${slurm_params}
  fi

  ## create a slrum submission script (strip indentation with sed before saving)
  cat <<__EOF | sed -e 's/^  //' > ${outdir}/download_SRA_samples.sbatch
  #!/bin/bash
  #SBATCH --qos=${qos:-rra}
  #SBATCH --partition=${partition:-rra}
  #SBATCH --time=${time:-6-00:00:00}
  #SBATCH --mem=${mem:-1G}
  #SBATCH --job-name=${job_name:-download_SRA_project}
  #SBATCH --output=${outdir}/logs/download_SRA_project_%a.log
  #SBATCH --array=0-$((${#samples[@]}-1))%${narray:-10}

  ## use task ID to assign each batch job a single unique sample
  samples=(${samples[@]})
  sample=\${samples[\$SLURM_ARRAY_TASK_ID]}
  biosample=\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 26)
  subdir=${outdir}/\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 13)

  ## make sure the web requests of all array jobs aren't submitted at once
  if [ \$SLURM_ARRAY_TASK_ID -lt ${narray:-10} ]; then
    sleep \$((\$SLURM_ARRAY_TASK_ID * 10))
  fi

  bash ${outdir}/download_SRA_sample.sh \${sample} \${biosample} \${subdir} ${prepare_path}

__EOF

  ## submit array job to slurm
  sbatch ${outdir}/download_SRA_samples.sbatch

else

  for sample in ${samples}; do
  
    biosample=$(grep ${sample} ${outdir}/runInfo.csv | cut -d ',' -f 26)
    subdir=${outdir}/$(grep ${sample} ${outdir}/runInfo.csv | cut -d ',' -f 13)

    bash ${outdir}/download_SRA_sample.sh ${sample} ${biosample} ${subdir} ${prepare_path}

  done

fi
