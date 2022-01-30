args <- commandArgs(TRUE)

taxid_in <- args[[1]]
dna <- Biostrings::readDNAStringSet(args[[2]])

ranks <- c('Superkingdom','Kingdom','Subkingdom','Superphylum','Phylum','Subphylum','Superclass','Class','Subclass','Infraclass','Superorder','Order','Suborder','Superfamily','Family','Subfamily','Tribe','Subtribe','Genus','Subgenus','Species','Subspecies')
taxid_raw <- read.table(taxid_in, sep='\t')
taxid_raw[,1] <- as.character(dna)[taxid_raw[,1]]
taxid_raw <- cbind(taxid_raw[,1],t(simplify2array(sapply(strsplit(taxid_raw[,2],';'), function(x) if(length(x)==1) c(x,rep(NA,43)) else x))))
taxid <- matrix(NA,nrow=nrow(taxid_raw),ncol=length(ranks),dimnames=list(taxid_raw[,1],ranks))
for(i in 1:nrow(taxid)) {
  for(j in ranks) {
    col <- grep(j,taxid_raw[i,])
    if(length(col>0)) {
      col <- col[[1]]
      if(as.numeric(taxid_raw[i,col+1]) > 50) {
        taxid[i,j] <- sub(paste0(j,':'),'',taxid_raw[i,col])
      } else {
        break
      }
    }
  }
  taxid[i,taxid[i,]=='NA'] <- NA
  for(j in 22:2) {
    if(!is.na(taxid[i,j])) {
      if(grepl(paste0(taxid[i,j-1],'_'),taxid[i,j])) {
        taxid[i,j] <- NA
      } else {
        break
      }
    }
  }
}

write.table(taxid, file=sub('.txt','_clean.txt',taxid_in), sep='\t', quote=F)
