## two positional arguments specifying 1) a filepath to a list of accession numbers, and 2) the output directory
accession_list=$1
outdir=$2

mkdir -p ${outdir}
cp ${accession_list} ${outdir}/accession_list.txt

cat <<EOF > ${outdir}/download_ref_genomes.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=20G
#SBATCH --job-name=download_ref_genomes
#SBATCH --output=${outdir}/download_ref_genomes.log

module purge
module load hub.apps/anaconda3
source activate ncbi_datasets

cd ${outdir}

datasets download genome accession --dehydrated --no-progressbar --include genome,gtf,gbff,gff3,rna,cds,protein,seq-report --inputfile accession_list.txt

unzip ncbi_dataset.zip
rm ncbi_dataset.zip

datasets rehydrate --gzip --directory ncbi_dataset

EOF

if $autorun; then

sbatch ${outdir}/download_ref_genomes.sbatch

fi

