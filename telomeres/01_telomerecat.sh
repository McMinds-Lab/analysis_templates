## two positional arguments specifying the directory of SRA downloads (with sample_list.txt in it) and the output analysis directory
readarray -t samples < $1/sample_list.txt
indir=$1
outdir=$2

# get local variables
source local.env

mkdir -p ${outdir}/01_telomerecat

cat <<EOF > ${outdir}/01_telomerecat/01_telomerecat.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=10G
#SBATCH --job-name=01_telomerecat
#SBATCH --output=${outdir}/01_telomerecat/01_telomerecat.log

for sample in ${samples[@]}; do

##convert fastqs to bams
module purge 
module load hub.apps/anaconda3
source activate picard

echo picard

picard FastqToSam \
  -F1 ${output_dir}/\${sample}/\${sample}_1.fastq.gz \
  -F2 ${output_dir}/\${sample}/\${sample}_2.fastq.gz \
  -O ${output_dir}/\${sample}/unaligned_read_pairs.bam \
  -SM \${sample} \
  -SO unsorted

##telomerecat requires a specific header format for some reason
conda deactivate
module purge 
module load apps/samtools/1.3.1
samtools view -H ${output_dir}/${sample}/unaligned_read_pairs.bam > ${output_dir}/${sample}/myheader.sam
printf "@SQ\tSN:fake_contig\tLN:1\n" >> ${output_dir}/${sample}/myheader.sam 
samtools reheader ${output_dir}/${sample}/myheader.sam ${output_dir}/${sample}/unaligned_read_pairs.bam > ${output_dir}/${sample}/unaligned_read_pairs.bam_tmp
mv ${output_dir}/${sample}/unaligned_read_pairs.bam_tmp ${output_dir}/${sample}/unaligned_read_pairs.bam

## run telomerecat
module purge 
module load hub.apps/anaconda3
source /shares/omicshub/apps/anaconda3/etc/profile.d/conda.sh
source activate telomerecat
telomerecat bam2length --output ${output_dir}/${sample}.csv ${output_dir}/${sample}/unaligned_read_pairs.bam
done
