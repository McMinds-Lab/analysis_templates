# get local variables
source local.env

## three positional arguments specifying 1) the directory of SRA downloads (with runInfo.csv in it), 2) the output analysis directory, and 3) either the string NONE or the path to a reference genome
indir=$1
outdir=$2
reference=$3

mkdir -p ${outdir}/00_browniecorrector/logs

samples=($(grep SRR ${indir}/runInfo.csv | grep 'WGA\|WGS' | cut -d ',' -f 1))

cat <<EOF > ${outdir}/00_browniecorrector/00_browniecorrector.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=${maxram}
#SBATCH --job-name=00_browniecorrector
#SBATCH --output=${outdir}/00_browniecorrector/logs/00_browniecorrector_%a.log
#SBATCH --array=0-$((${#samples[@]}-1))%10

samples=(${samples[@]})
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

mkdir ${outdir}/00_browniecorrector/\${sample}/

#trim adapters
module purge
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
source activate bbtools

in1=(${indir}/*/\${sample}/\${sample}_1.fastq.gz)
in2=(${indir}/*/\${sample}/\${sample}_2.fastq.gz)

bbduk.sh \
  in1=\${in1} \
  in2=\${in2} \
  out1=${outdir}/00_browniecorrector/\${sample}/trimmed_R1.fastq \
  out2=${outdir}/00_browniecorrector/\${sample}/trimmed_R2.fastq \
  ref=adapters,artifacts,phix \
  ktrim=l k=23 mink=11 hdist=2 qtrim=r minlength=40 tpe trimq=10
  
reformat.sh \
  in1=${outdir}/00_browniecorrector/\${sample}/trimmed_R1.fastq \
  in2=${outdir}/00_browniecorrector/\${sample}/trimmed_R2.fastq \
  out=${outdir}/00_browniecorrector/\${sample}/interleaved.fastq

rm ${outdir}/00_browniecorrector/\${sample}/trimmed_R1.fastq ${outdir}/00_browniecorrector/\${sample}/trimmed_R2.fastq

conda deactivate
source activate browniecorrector

coverage=\$(bc <<< "scale=2;\$(tail -n +2 ${outdir}/00_browniecorrector/\${sample}/interleaved.fastq | awk '1 == NR % 4 {L+=length(\$0)} END {print(L)}')/$(zcat ${reference} | awk '0 == NR % 2 {L+=length($0)} END {print(L)}')")

mcminds_browniecorrector.sh \
  ${outdir}/00_browniecorrector/\${sample}/interleaved.fastq \
  \${coverage} \
  ${outdir}/00_browniecorrector/\${sample}/brownieout \
  TTAGGGTTAGG \
  /shares/omicshub/apps/anaconda3/envs/browniecorrector/browniecorrector/bashScripts/

rm ${outdir}/00_browniecorrector/\${sample}/interleaved.fastq
rm ${outdir}/00_browniecorrector/\${sample}/brownieout/normal.fastq
rm -r ${outdir}/00_browniecorrector/\${sample}/brownieout/clusters
gzip --best -c ${outdir}/00_browniecorrector/\${sample}/brownieout/brownie.fastq > ${outdir}/00_browniecorrector/\${sample}/brownieout/brownie.fastq.gz
rm ${outdir}/00_browniecorrector/\${sample}/brownieout/brownie.fastq
  
EOF

if $autorun; then
    sbatch ${outdir}/00_browniecorrector/00_browniecorrector.sbatch
fi
