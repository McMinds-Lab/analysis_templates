print(sessionInfo())
cat(paste('dada2 version:', packageVersion('dada2'), '\n'))

args <- commandArgs(TRUE)

nthreads <- as.numeric(args[[1]])
indir <- args[[2]]
outdir <- args[[3]]
taxref <- args[[4]]

merged_reads <- sort(list.files(indir, pattern=".fastq.gz", full.names=TRUE))
sample.names <- sub('.fastq.gz','',basename(merged_reads))

png(file.path(outdir, "quality_profile.png"), height=600, width=600)
dada2::plotQualityProfile(merged_reads[1:10])
dev.off()

#learning errors is a key component, must be done with reads that have no N's
err_merged_reads <- dada2::learnErrors(merged_reads, multithread=nthreads)

pdf(file.path(outdir, "dada_err.pdf"))
dada2::plotErrors(err_merged_reads, nominalQ=TRUE)
dev.off()

#after dada2 has learned errors, user can denoise
derepped <- dada2::derepFastq(merged_reads)
dada_merged <- dada2::dada(derepped, err=err_merged_reads, pool='pseudo', multithread=nthreads)

#Construct ASV
seqtab <- dada2::makeSequenceTable(dada_merged)

write.table(seqtab, file.path(outdir, 'asv.tsv'), sep='\t')

taxid <- dada2::assignTaxonomy(seqtab, taxref, multithread=nthreads, taxLevels = c('Superkingdom','Kingdom','Subkingdom','Superphylum','Phylum','Subphylum','Superclass','Class','Subclass','Infraclass','Superorder','Order','Suborder','Superfamily','Family','Subfamily','Tribe','Subtribe','Genus','Subgenus','Species','Subspecies'))

write.table(taxid, file.path(outdir, 'taxid.tsv'), sep='\t')

save.image(file.path(outdir, 'ASVs.RData'))

