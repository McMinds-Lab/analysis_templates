abbr <- c('_spk','_k','_sbk','_spp','_p','_sbp','_spc','_c','_sbc','_ic','_spo','_o','_sbo','_io','_po','_spf','_f','_sbf','_t','_sbt','_g','_sbg','_s','_sbs')
dna <- Biostrings::DNAStringSet(gsub('U','T',Biostrings::readBStringSet('~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean.fasta.gz')))
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
dna <- dna[!names(dna) %in% c('Unknown;','Bacteria;','Archaea;','Eukaryota;')]
euk_inds <- grep('^Eukaryota;',names(dna))
arch_inds <- grep('^Archaea;',names(dna))
bact_inds <- grep('^Bacteria;',names(dna))
dna <- dna[c(euk_inds, sample(arch_inds, min(100,length(arch_inds))), sample(bact_inds, min(100,length(bact_inds)))),]
Biostrings::writeXStringSet(dna,'~/data/ref/silva/silva2dada2/SILVA_138.1_SSURef_tax_silva_clean_dada2_euks.fasta.gz',compress=TRUE,format='fasta',width=10000)
