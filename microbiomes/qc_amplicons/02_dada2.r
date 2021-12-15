library(dada2)

args <- commandArgs(TRUE)

nthreads <- as.numeric(args[[1]])
indir <- args[[2]]
outdir <- args[[3]]

merged_reads <- sort(list.files(indir, pattern=".fastq.gz", full.names=TRUE))
sample.names <- sub('.fastq.gz','',basename(merged_reads))

filt_dir <- file.path(outdir, "filtered")
dir.create(filt_dir, recursive=TRUE)

png(file.path(outdir, "quality_profile.png"), height=600, width=600)
plotQualityProfile(merged_reads[1:10])
dev.off()

merged_filt <- file.path(filt_dir, paste0(sample.names, "_M_filt.fastq.gz"))

names(merged_filt) <- sample.names

filtered_out <- filterAndTrim(
  fwd = merged_reads,
  filt = merged_filt,
  truncQ = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = nthreads
)

#learning errors is a key component, must be done with filtered reads
err_merged_reads <- learnErrors(merged_filt, multithread=nthreads)

pdf(file.path(outdir, "dada_err.pdf"))
plotErrors(err_merged_reads, nominalQ=TRUE)
dev.off()

#after dada2 has learned errors, user can denoise
derepped <- derepFastq(merged_filt)
dada_merged <- dada(derepped, err=err_merged_reads, multithread=nthreads)

#Construct ASV
seqtab <- makeSequenceTable(dada_merged)

write.table(seqtab, file.path(outdir, 'asv.tsv'), sep='\t')

save.image(file.path(outdir, 'ASVs.RData'))

