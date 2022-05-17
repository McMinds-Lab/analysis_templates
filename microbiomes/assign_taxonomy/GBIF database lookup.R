library("taxize")
library("purrr")
library('jsonlite')
library('data.table')
library("expss")
library("reshape")
taxid <- read.table(file.path, 
                    sep = "\t" , 
                    header = T,
                    row.names = 1)
#taxid <- head(taxid, n=10)
allids <- character(0)

start <- proc.time()

taxid[is.na(taxid$Kingdom), "Kingdom"] <- "Eukaryota"

for (i in 1:nrow(taxid)) {
  for (j in ncol(taxid):1) {
    if(is.na(taxid[i, j])){
      next
    }
    
    tempid <- gnr_resolve((taxid[i,j]),
                          best_match_only=TRUE, 
                          fields ='all',
                          preferred_data_sources=11)
    
    if((length(tempid)==0 & j>1)){
      next
    }
    else if(length(tempid)==0 & j==1){
      allids <- c(allids, NA)
    }
    else if(length(tempid)>0){
      allids <- c(allids, tempid$taxon_id)
      break
    }
  }
}

delta <- proc.time() - start
delta
names(allids) <- rownames(taxid)

tax <- classification(allids, db = 'gbif', return_id = TRUE)
names(tax) <- rownames(taxid)

tax_list <- do.call(rbind, tax)
tax_list[is.na(tax_list$rank), "rank"] <- "domain"
tax_list[is.na(tax_list$rank), "name"] <- "Eukaryote"
tax_list$rank <- factor(tax_list$rank, levels=c('domain',
                                                "kingdom",
                                                'phylum',
                                                'class',
                                                'order',
                                                'family',
                                                'genus',
                                                'species'))

tax_list$asv <- rownames(tax_list)
rownames(tax_list) <- NULL
tax_list$asv <- unlist(map(strsplit(tax_list$asv,'.',fixed=TRUE), 1))
tax_fin <- cast(tax_list,
                formula= asv ~rank,
                value = "name" )
tax_fin$domain <- "Eukaryote"
write.table(tax_fin, file = file.path,
            sep="\t",
            row.names=T)

