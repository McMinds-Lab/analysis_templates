abbr <- c('_spk','_k','_sbk','_spp','_p','_sbp','_spc','_c','_sbc','_ic','_spo','_o','_sbo','_spf','_f','_sbf','_t','_sbt','_g','_sbg','_s','_sbs')
dna <- Biostrings::DNAStringSet(gsub('U','T',Biostrings::readBStringSet('~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean.fasta')))
tax <- read.table('~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_accession_taxonomy.txt.txt', sep='\t',quote='')
tax[,2] <- sapply(strsplit(tax[,2],';'), function(x) {
  for(y in length(x):2) {
    if(grepl('uncultured|unidentified|metagenome|environmental|NA', x[[y]])) {
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
    for(z in 2:length(x)) {
      if(x[[z]] == 'NA') {
        x[[z]] <- paste0(x[[z-1]],abbr[[z]])
      }
    }
  }
  paste0(paste(x,collapse=';'),';')
})
names(dna) <- tax[match(names(dna),tax[,1]),2]
Biostrings::writeXStringSet(dna,'~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_dada2.fasta.gz',compress=TRUE)
