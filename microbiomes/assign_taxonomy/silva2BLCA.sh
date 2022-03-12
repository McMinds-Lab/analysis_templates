## requires python2 with pycogent installed
source local.env

outdir=$1
taxdump=$2
acc2tax=$3
silvapath=$4

if [ ${taxdump} == download ]; then
  taxdump_path=${outdir}/03_ref/ncbi_taxonomy
else
  taxdump_path=${taxdump}
fi
if [ ${acc2tax} == download ]; then
  acc2tax_path=${outdir}/03_ref/ncbi_taxonomy/nucl_gb.accession2taxid
else
  acc2tax_path=${acc2tax}
fi

mkdir -p ${outdir}/03_ref
echo "bash $0 $@" > ${outdir}/03_ref/this_command.sh
cp 03_silva2BLCA.r ${outdir}/03_ref/

cat <<EOF > ${outdir}/03_ref/03_ref.sbatch
#!/bin/bash
#SBATCH --job-name=03_ref
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/03_ref/03_ref.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=${maxram}
#SBATCH --time=3-00:00:00

cd ${outdir}
mkdir -p ${outdir}/03_ref/ncbi_taxonomy
cd ${outdir}/03_ref/ncbi_taxonomy

if [ ${taxdump} == download ]; then
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-01-01.zip
unzip taxdmp_2022-01-01.zip
fi
if [ ${acc2tax} == download ]; then
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz
fi

mkdir ${outdir}/03_ref/silva2BLCA
cd ${outdir}/03_ref/silva2BLCA

if [ ${silvapath} == download ]; then
  wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz
  zcat SILVA_138.1_SSURef_tax_silva.fasta.gz | cut -d'.' -f1 > SILVA_138.1_SSURef_tax_silva_clean.fasta
else
  zcat ${silvapath} | cut -d'.' -f1 > SILVA_138.1_SSURef_tax_silva_clean.fasta
fi

cd ${outdir}/03_ref
module load hub.apps/anaconda3
source activate python2
python2 $(pwd)/entrez_qiime/entrez_qiime.py \
  -i ${outdir}/03_ref/silva2BLCA/SILVA_138.1_SSURef_tax_silva_clean.fasta \
  -n ${taxdump_path} \
  -a ${acc2tax_path} \
  -r superkingdom,kingdom,subkingdom,superphylum,phylum,subphylum,superclass,class,subclass,infraclass,superorder,order,suborder,infraorder,parvorder,superfamily,family,subfamily,tribe,subtribe,genus,subgenus,species,subspecies

module load hub.apps/R
Rscript 03_silva2BLCA.r

source activate BLCA
gunzip -c ${outdir}/03_ref/silva2BLCA/SILVA_138.1_SSURef_tax_silva_clean_euks_BLCA.fasta.gz | makeblastdb -in - -parse_seqids -blastdb_version 5 -title SILVA_138.1_euks_BLCA_db -dbtype nucl -out ${outdir}/03_ref/silva2BLCA/SILVA_138.1_euks_BLCA

gzip --best silva2BLCA/SILVA_138.1_SSURef_tax_silva_clean.fasta

EOF

if $autorun; then
    sbatch ${outdir}/03_ref/03_ref.sbatch
fi
