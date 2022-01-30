ranks <- c('Superkingdom','Kingdom','Subkingdom','Superphylum','Phylum','Subphylum','Superclass','Class','Subclass','Infraclass','Superorder','Order','Suborder','Superfamily','Family','Subfamily','Tribe','Subtribe','Genus','Subgenus','Species','Subspecies')
abbr <- c('_spk','_k','_sbk','_spp','_p','_sbp','_spc','_c','_sbc','_ic','_spo','_o','_sbo','_spf','_f','_sbf','_t','_sbt','_g','_sbg','_s','_sbs')

tax <- read.table('~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_accession_taxonomy.txt.txt', sep='\t',quote='')
tax <- tax[c(grep('Eukaryota',tax[,2]),sample(grep('Bacteria',tax[,2]),100),sample(grep('Archaea',tax[,2]),100)),]
tax[,2] <- sapply(strsplit(tax[,2],';'), function(x) {
  for(y in length(x):2) {
    if(grepl('uncultured|unidentified|metagenome|environmental|NA|artificial|fungal|endophyte|eukaryote|Incertae|-like|_nr|_flagellate', x[[y]])) {
      x <- x[-y]
    } else {
      break
    } 
  }
  if(length(x) == 0) {
    x <- 'Unknown'
  } else if(length(x) == 1) {
    if(x == 'NA') {
      x <- 'Unknown'
    }
  } else {
    for(y in 2:length(x)) {
      if(x[[y]] == 'NA') {
        x[[y]] <- paste0(sub('.*:','',x[[y-1]]),abbr[[y]])
      }
      x[[y]] <- paste(ranks[[y]],x[[y]],sep=':')
    }
  }
  x[[1]] <- paste('Superkingdom',x[[1]],sep=':')
  paste0(rev(x), collapse=';')
})

write.table(tax, file='~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_accession_taxonomy_BLCA.txt', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

dna <- Biostrings::DNAStringSet(gsub('U','T',Biostrings::readBStringSet('~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean.fasta.gz')))
dna <- dna[tax[,1]]
Biostrings::writeXStringSet(dna,'~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_euks_BLCA.fasta.gz',compress=TRUE,format='fasta',width=10000)
