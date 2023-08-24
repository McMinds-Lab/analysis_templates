## three or four positional arguments specifying 1) an SRA project ID starting with SRP, 2) the output directory, 3) either a path SLURM parameters file or the word 'false' (If false script will simply download samples with for loop and won't submit a slurm job), and 4) optional path to a script to set your PATH variable to include the entrez-direct and sratools packages
## make sure you have 'prepare_path.sh' and 'slurm_params.sh' in the same directory as this script, if you need them (example files are included as 'prepare_path.sh.txt' and 'slurm_params.sh.txt')

srp=$1
outdir=$2
slurm_params=$3

## input error checking
if [ $# -lt 3 ]; then
  echo "First three positional arguments are required"
fi

## if tools are not automatically in path, put them there
if [ $# -eq 4 ]; then
  prepare_path=$4
  source ${prepare_path}
fi

## this tells the script where to find this file and download_SRA_sample.sh
scriptdir=$(dirname "$(realpath "$0")")

## prepare output directory
mkdir -p ${outdir}

## download a list of samples, and some of their metadata, associated with the specified project
esearch -db sra -query ${srp} | efetch -format runinfo > ${outdir}/runInfo.csv

## extract sample identifiers
samples=($(grep SRR ${outdir}/runInfo.csv | cut -d ',' -f 1))

## require the explicit input of 'false' to enable for-loop
if [ ${slurm_params} != "false" ]; then

## create a logs directory
mkdir -p ${outdir}/logs

# get slurm parameters
source ${slurm_params}

cat <<EOF > ${outdir}/download_SRA_samples.sbatch
#!/bin/bash
#SBATCH --qos=${qos:-rra}
#SBATCH --partition=${partition:-rra}
#SBATCH --time=${time:-6-00:00:00}
#SBATCH --mem=${mem:-1G}
#SBATCH --job-name=${job_name:-download_SRA_project}
#SBATCH --output=${outdir}/logs/download_SRA_project_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%${nthreads:-1}

## use task ID to assign each batch job a single unique sample
samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}
biosample=\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 26)
subdir=${outdir}/\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 13)

## make sure the web requests of all array jobs aren't submitted at once
if [ \$SLURM_ARRAY_TASK_ID -lt ${nthreads:-1} ]; then
  sleep \$((\$SLURM_ARRAY_TASK_ID * 10))
fi

bash ${scriptdir}/download_SRA_sample.sh \${sample} \${biosample} \${subdir} ${prepare_path}

EOF

sbatch ${outdir}/download_SRA_samples.sbatch

else

  for sample in ${samples}; do
  
    biosample=$(grep ${sample} ${outdir}/runInfo.csv | cut -d ',' -f 26)
    subdir=${outdir}/$(grep ${sample} ${outdir}/runInfo.csv | cut -d ',' -f 13)

    bash ${scriptdir}/download_SRA_sample.sh ${sample} ${biosample} ${subdir} ${prepare_path}

  done

fi
