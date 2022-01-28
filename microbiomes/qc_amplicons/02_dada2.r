library(dada2)

args <- commandArgs(TRUE)

nthreads <- as.numeric(args[[1]])
indir <- args[[2]]
outdir <- args[[3]]

merged_reads <- sort(list.files(indir, pattern=".fastq.gz", full.names=TRUE))
sample.names <- sub('.fastq.gz','',basename(merged_reads))

png(file.path(outdir, "quality_profile.png"), height=600, width=600)
plotQualityProfile(merged_reads[1:10])
dev.off()

#learning errors is a key component, must be done with reads that have no N's
err_merged_reads <- learnErrors(merged_reads, multithread=nthreads)

pdf(file.path(outdir, "dada_err.pdf"))
plotErrors(err_merged_reads, nominalQ=TRUE)
dev.off()

#after dada2 has learned errors, user can denoise
derepped <- derepFastq(merged_reads)
dada_merged <- dada(derepped, err=err_merged_reads, pool='pseudo', multithread=nthreads)

#Construct ASV
seqtab <- makeSequenceTable(dada_merged)

write.table(seqtab, file.path(outdir, 'asv.tsv'), sep='\t')

save.image(file.path(outdir, 'ASVs.RData'))

