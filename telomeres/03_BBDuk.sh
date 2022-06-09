# get local variables
source local.env

## three positional arguments specifying 1) the directory of SRA downloads (with runInfo.csv in it), 2) the output analysis directory, and 3) either the string NONE or the path to a reference genome
indir=$1
outdir=$2
reference=$3

mkdir -p ${outdir}/00_browniecorrector/logs

samples=(${indir})

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
in1=(${indir}/*_1.fastq.gz)
in2=(${indir}/*_2.fastq.gz)
bbduk.sh \
  in1=\${in1} \
  in2=\${in2} \
  out1=${outdir}/00_browniecorrector/${sample}_trimmed_R1.fastq
  out2=${outdir}/00_browniecorrector/${sample}_trimmed_R2.fastq
  ref=adapters,artifacts \
  qtrim=t 
  k=0
  minlf=0 
  hdist=2 
  minlength=X  
  trimq=6
  ftl=10
  ftr
EOF

if $autorun; then
    sbatch ${outdir}/00_browniecorrector/00_browniecorrector.sbatch
fi