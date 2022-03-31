# get local variables
source local.env

seqpath=$1
outdir=$2
taxpath=$3
dbpath=$4

mkdir -p ${outdir}/BLCA/
echo "bash $0 $@" > ${outdir}/BLCA/this_command.sh
cp clean_BLCA.r ${outdir}/BLCA/

cat <<EOF > ${outdir}/BLCA/BLCA.sbatch
#!/bin/bash
#SBATCH --job-name=BLCA
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/BLCA/BLCA.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=${maxram}
#SBATCH --time=3-00:00:00

module load hub.apps/anaconda3
source activate BLCA

gunzip -c ${seqpath} > ${outdir}/BLCA/seqs.fasta

#2.blca_main.py from rmcminds' fork
2.blca_main.py -i ${outdir}/BLCA/seqs.fasta -r ${taxpath} -q ${dbpath} -o ${outdir}/BLCA/blca_taxa.txt -p ${nthreads} -k Superkingdom,Kingdom,Subkingdom,Superphylum,Phylum,Subphylum,Superclass,Class,Subclass,Infraclass,Superorder,Order,Suborder,Infraorder,Parvorder,Superfamily,Family,Subfamily,Tribe,Subtribe,Genus,Subgenus,Species,Subspecies

module purge
module load hub.apps/R
Rscript ${outdir}/BLCA/clean_BLCA.r ${outdir}/BLCA/blca_taxa.txt ${outdir}/BLCA/seqs.fasta

rm ${outdir}/BLCA/seqs.fasta

EOF

if $autorun; then
    sbatch ${outdir}/BLCA/BLCA.sbatch
fi
