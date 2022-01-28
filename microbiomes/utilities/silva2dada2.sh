## requires python2 with cogent installed

mkdir ~/data/ref/ncbi_taxonomy
cd ~/data/ref/ncbi_taxonomy

curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -o taxdump.tar.gz
tar -zxvf taxdump.tar.gz -C taxdump
curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz -o nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz


mkdir ~/data/ref/silva/silva2dada2
cd ~/data/ref/silva/silva2dada2

zcat < ~/data/ref/silva/SILVA_138.1_SSURef_tax_silva.fasta.gz | cut -d'.' -f1 > SILVA_138.1_SSURef_tax_silva_clean.fasta

python2 ~/scripts/my_entrez_qiime.py \
  -i SILVA_138.1_SSURef_tax_silva_clean.fasta \
  -n ~/data/ref/ncbi_taxonomy/taxdump \
  -a ~/data/ref/ncbi_taxonomy/nucl_gb.accession2taxid \
  -r superkingdom,kingdom,subkingdom,superphylum,phylum,subphylum,superclass,class,subclass,infraclass,superorder,order,suborder,superfamily,family,subfamily,tribe,subtribe,genus,subgenus,species,subspecies

Rscript ~/scripts/analysis_templates/microbiomes/utilities/silva2dada2.r

gzip --best ~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean.fasta

cd ~/data/ref/
tar cvzf ncbi_taxonomy.tar.gz ncbi_taxonomy


