ranks <- c('Superkingdom','Kingdom','Subkingdom','Superphylum','Phylum','Subphylum','Superclass','Class','Subclass','Infraclass','Superorder','Order','Suborder','Infraorder','Parvorder','Superfamily','Family','Subfamily','Tribe','Subtribe','Genus','Subgenus','Species','Subspecies')
dna <- Biostrings::readDNAStringSet('~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_dada2_euks.fasta.gz')

trainingSet <- DECIPHER::LearnTaxa(dna,paste0('Root;',names(dna)))

save(trainingSet, file='~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_dada2_euks_IDTAXA.RData')

