## this script uses the Rmd at https://zenodo.org/record/3986799#.Ye8aV_XMLxg as a guide for how to do-it-yourself

ref_dir <- '~/data/ref'
base_url <- 'https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/'
ref99 <- 'SILVA_138_SSURef_NR99_tax_silva.fasta.gz'
refFull <- 'SILVA_138_SSURef_tax_silva.fasta.gz'
tax <- 'tax_slv_ssu_138.txt.gz'
  
if(!file.exists(file.path(ref_dir, ref99))) {
download.file(file.path(base_url, ref99), file.path(ref_dir, ref99))
}
if(!file.exists(file.path(ref_dir, refFull))) {
download.file(file.path(base_url, refFull), file.path(ref_dir, refFull))
}
if(!file.exists(file.path(ref_dir, tax))) {
download.file(file.path(base_url, 'taxonomy', tax), file.path(ref_dir, tax))
}

source('dada2_from_silva_functions.r')

makeTaxonomyFasta_SilvaNR_18S(file.path(ref_dir, ref99), file.path(ref_dir, tax), file.path(ref_dir, paste0(sub('.fasta.gz','',ref99),'_dada2.fasta.gz')))
makeTaxonomyFasta_SilvaNR_18S(file.path(ref_dir, refFull), file.path(ref_dir, tax), file.path(ref_dir, paste0(sub('.fasta.gz','',refFull),'_dada2.fasta.gz')))

