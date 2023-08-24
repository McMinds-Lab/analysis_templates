# get local variables (I should reformat this to be system-agnostic, such that this step just sets the PATH variable to enable everything else)
source local.env

## two positional arguments specifying 1) an SRA project ID starting with SRP and 2) the output directory
srp=$1
outdir=$2

mkdir -p ${outdir}/logs

## put the entrez-direct package in the PATH variable
module purge
module load hub.apps/anaconda3/2020.11
source activate entrez-direct

## download a list of samples, and some of their metadata, associated with the specified project
esearch -db sra -query ${srp} | efetch -format runinfo > ${outdir}/runInfo.csv

## extract sample identifiers
samples=($(grep SRR ${outdir}/runInfo.csv | cut -d ',' -f 1))

## between here and 'EOF' create a new script that can be run on the server as a batch job for parallel downloads of all the samples
cat <<EOF > ${outdir}/download_SRA.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=${maxram}
#SBATCH --job-name=download_SRA
#SBATCH --output=${outdir}/logs/download_SRA_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%${nthreads}

## use task ID to assign each batch job a single unique sample
samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

subdir=\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 13)

mkdir -p ${outdir}/\${subdir}/\${sample}

module purge
module load hub.apps/anaconda3/2020.11
source activate entrez-direct

## make sure the web requests of all array jobs aren't submitted at once
if [ \$SLURM_ARRAY_TASK_ID -lt ${nthreads} ]; then
  sleep \$((\$SLURM_ARRAY_TASK_ID * 10))
fi

## download sample metadata
biosample=\$(grep \${sample} ${outdir}/runInfo.csv | cut -d ',' -f 26)
esearch -db biosample -query \${biosample} | efetch > ${outdir}/\${subdir}/\${sample}/metadata.txt
esearch -db biosample -query \${biosample} | efetch -format xml > ${outdir}/\${subdir}/\${sample}/metadata.xml

conda deactivate
module purge
module load apps/sratoolkit/2.10.7

## download raw sequences from sample
fasterq-dump -e 1 \${sample} -t ${outdir}/\${subdir}/\${sample} -O ${outdir}/\${subdir}/\${sample}

## compress sequence file
for file in ${outdir}/\${subdir}/\${sample}/*.fastq; do
  gzip --best \${file} &
done
wait

EOF

if $autorun; then

sbatch ${outdir}/download_SRA.sbatch

fi

