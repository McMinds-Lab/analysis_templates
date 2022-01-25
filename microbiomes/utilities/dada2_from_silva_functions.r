## this is modified from the function makeTaxonomyFasta_SilvaNR in dada2 v1.20.0
makeTaxonomyFasta_SilvaNR_18S <- function (fin, ftax, fout, include.species = FALSE, compress = TRUE) {
  ##
  if(!isNamespaceLoaded('dada2')) {
    loadNamespace("dada2")
    attachNamespace("dada2")
  }
  if(!isNamespaceLoaded('Biostrings')) {
    loadNamespace("Biostrings")
    attachNamespace("Biostrings")
  }
  if(!isNamespaceLoaded('ShortRead')) {
    loadNamespace("ShortRead")
    attachNamespace("ShortRead")
  }
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
  taxa.ba.mat[grep(paste(c("Uncultured", "uncultured",'unidentified','artificial','fungal','endophyte','eukaryote'),collapse='|'), taxa.ba.mat)] <- NA
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


my_assignTaxonomy <- function (seqs, refFasta, minBoot = 50, tryRC = FALSE, outputBootstraps = FALSE, 
             taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", 
                           "Genus", "Species"), multithread = FALSE, verbose = FALSE) {
  ##
  loadNamespace("dada2")
  #
  MIN_REF_LEN <- 20
  MIN_TAX_LEN <- 50
  seqs <- getSequences(seqs)
  if (min(nchar(seqs)) < MIN_TAX_LEN) {
    warning("Some sequences were shorter than ", MIN_TAX_LEN, 
            " nts and will not receive a taxonomic classification.")
  }
  refsr <- readFasta(refFasta)
  lens <- width(sread(refsr))
  if (any(lens < MIN_REF_LEN)) {
    refsr <- refsr[lens >= MIN_REF_LEN]
    warning(paste0("Some reference sequences were too short (<", 
                   MIN_REF_LEN, "nts) and were excluded."))
  }
  refs <- as.character(sread(refsr))
  tax <- as.character(id(refsr))
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x))
  UNITE <- FALSE
  if (all(grepl("FU\\|re[pf]s", tax[1:10]))) {
    UNITE <- TRUE
    cat("UNITE fungal taxonomic reference detected.\n")
    tax <- sapply(strsplit(tax, "\\|"), `[`, 5)
    tax <- gsub("[pcofg]__unidentified;", "_DADA2_UNSPECIFIED;", 
                tax)
    tax <- gsub(";s__(\\w+)_", ";s__", tax)
    tax <- gsub(";s__sp$", ";_DADA2_UNSPECIFIED", tax)
  }
  if (!grepl(";", tax[[1]])) {
    if (length(unlist(strsplit(tax[[1]], "\\s"))) == 3) {
      stop("Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies).")
    }
    else {
      stop("Incorrect reference file format for assignTaxonomy.")
    }
  }
  tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- max(tax.depth)
  for (i in seq(length(tax))) {
    if (tax.depth[[i]] < td) {
      ## looks liek the raw silva fasta needs a semicolon pasted to the end of each taxonomy string if not the final
      tax[[i]] <- paste0(tax[[i]],';')
      #
      for (j in seq(td - tax.depth[[i]])) {
        tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
      }
    }
  }
  genus.unq <- unique(tax)
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol = td, 
                    byrow = TRUE)
  tax.df <- as.data.frame(tax.mat)
  for (i in seq(ncol(tax.df))) {
    tax.df[, i] <- factor(tax.df[, i])
    tax.df[, i] <- as.integer(tax.df[, i])
  }
  tax.mat.int <- as.matrix(tax.df)
  if (is.logical(multithread)) {
    if (multithread == TRUE) {
      RcppParallel::setThreadOptions(numThreads = "auto")
    }
    else {
      RcppParallel::setThreadOptions(numThreads = 1)
    }
  }
  else if (is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
  }
  else {
    warning("Invalid multithread parameter. Running as a single thread.")
    RcppParallel::setThreadOptions(numThreads = 1)
  }
  assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus, 
                                   tax.mat.int, tryRC, verbose)
  bestHit <- genus.unq[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i, 
  ] >= minBoot])
  tax.out <- matrix(NA_character_, nrow = length(seqs), ncol = td)
  for (i in seq(length(seqs))) {
    if (length(taxes[[i]]) > 0) {
      tax.out[i, 1:length(taxes[[i]])] <- taxes[[i]]
    }
  }
  rownames(tax.out) <- seqs
  colnames(tax.out) <- taxLevels[1:ncol(tax.out)]
  tax.out[tax.out == "_DADA2_UNSPECIFIED"] <- NA_character_
  if (outputBootstraps) {
    boots.out <- matrix(boots, nrow = length(seqs), ncol = td)
    rownames(boots.out) <- seqs
    colnames(boots.out) <- taxLevels[1:ncol(boots.out)]
    list(tax = tax.out, boot = boots.out)
  }
  else {
    tax.out
  }
}
