## this is modified from the function makeTaxonomyFasta_SilvaNR in dada2 v1.20.0
makeTaxonomyFasta_SilvaNR_18S <- function (fin, ftax, fout, include.species = FALSE, compress = TRUE) {
  ##
  library("dada2", quietly=TRUE)
  library("Biostrings", quietly=TRUE)
  library("ShortRead", quietly=TRUE)
  #
  xset <- DNAStringSet(readRNAStringSet(fin, format = "fasta"))
  taxl <- names(xset)
  names(taxl) <- sapply(strsplit(names(xset), "\\s"), `[`, 
                        1)
  if (any(duplicated(names(taxl)))) 
    stop("Duplicated sequence IDs detected.")
  names(xset) <- names(taxl)
  taxl <- gsub("^[A-Za-z0-9.]+\\s", "", taxl)
  taxl <- gsub(";YM;", ";", taxl)
  taxa <- strsplit(taxl, ";")
  silva.taxa <- read.table(ftax, sep = "\t", col.names = c("Taxon", 
                                                           "V2", "Level", "V4", "V5"), stringsAsFactors = FALSE)
  silva.taxa <- silva.taxa[, c("Taxon", "Level")]
  kingdom <- sapply(strsplit(taxl, ";"), `[`, 1)
  taxl.ba <- taxl[kingdom %in% 'Eukaryota']
  taxa.ba <- taxa[names(taxl.ba)]
  for(i in 1:length(taxa.ba)) {
    spl <- strsplit(taxa.ba[[i]][length(taxa.ba[[i]])], ' ')[[1]]
    if(length(spl) > 2) {spl <- spl[1:2]}
    if(length(spl) == 2) {
      if(spl[[2]] %in% c('sp.','sp','cf.','cf')) {
        spl <- spl[-2]
      } else {
        spl[[2]] <- paste(spl,collapse='_')
      }
    }
    taxa.ba[[i]] <- c(taxa.ba[[i]][-length(taxa.ba[[i]])], spl)
  }
  taxa.ba.mat <- unname(t(do.call(cbind, lapply(taxa.ba, ts))))
  rownames(taxa.ba.mat) <- names(taxl.ba)
  taxa.ba.mat.string <- matrix("UNDEF", nrow = nrow(taxa.ba.mat), 
                               ncol = ncol(taxa.ba.mat))
  rownames(taxa.ba.mat.string) <- names(taxl.ba)
  taxa.ba.mat.string[, 1] <- paste0(taxa.ba.mat[, 1], ";")
  for (col in seq(2, ncol(taxa.ba.mat))) {
    taxa.ba.mat.string[, col] <- paste0(taxa.ba.mat.string[, 
                                                           col - 1], taxa.ba.mat[, col], ";")
  }
  if (any(taxa.ba.mat.string == "UNDEF")) 
    stop("Taxon string matrix was not fully initialized.")
  taxa.ba.mat[grep(paste(c("uncultured",'unidentified','artificial','fungal','endophyte','eukaryote','metagenome','environmental','Incertae','-like','_gen.'),collapse='|'), taxa.ba.mat, ignore.case=TRUE)] <- NA
  set.seed(100)
  N_EUK <- 100
  euk.keep <- sample(names(taxl)[kingdom %in% c("Bacteria", "Archaea")], 
                     N_EUK)
  taxa.euk.mat <- matrix("", nrow = N_EUK, ncol = ncol(taxa.ba.mat))
  rownames(taxa.euk.mat) <- euk.keep
  taxa.euk.mat[, 1] <- "Bacteria_Archaea"
  taxa.euk.mat[, 2:ncol(taxa.euk.mat)] <- NA
  taxa.mat.final <- rbind(taxa.ba.mat, taxa.euk.mat)
  taxa.string.final <- apply(taxa.mat.final, 1, function(x) {
    tst <- paste(x, collapse = ";")
    tst <- paste0(tst, ";")
    tst <- gsub("NA;", "", tst)
    tst
  })
  if (any(is.na(names(taxa.string.final)))) 
    stop("NA names in the final set of taxon strings.")
  if (!all(names(taxa.string.final) %in% names(xset))) 
    stop("Some names of the final set of taxon strings don't match sequence names.")
  xset.out <- xset[names(taxa.string.final)]
  cat(length(xset.out), "reference sequences were output.\n")
  print(table(taxa.mat.final[, 1], useNA = "ifany"))
  writeFasta(ShortRead(unname(xset.out), BStringSet(taxa.string.final)), 
             fout, width = 20000L, compress = compress)
}
