print(sessionInfo())
cat(paste('dada2 version:', packageVersion('dada2'), '\n'))

args <- commandArgs(TRUE)

nthreads <- as.numeric(args[[1]])
indir <- args[[2]]
outdir <- args[[3]]

merged_reads <- sort(list.files(indir, pattern = ".fastq.gz", full.names = TRUE))
sample.names <- sub('.fastq.gz', '', basename(merged_reads))

filt_reads <- file.path(outdir, "filtered", basename(merged_reads))
# filter out reads that are probably just primer dimer
filt <- dada2::filterAndTrim(merged_reads, filt_reads, truncQ = 0, minLen = 50, compress = TRUE, multithread = nthreads, qualityType='FastqQuality')

filt_reads_exist <- filt_reads[file.exists(filt_reads)]
stopifnot(length(filt_reads_exist) > 0)

# note that the q-scores appear to be inflated relative to Nanoplot results. This is because dada2 takes mean of raw q-scores, while nanoplot takes mean of probabilities (e.g. arithmetic vs geometric means)
png(file.path(outdir, "quality_profile.png"), height=600, width=600)
dada2::plotQualityProfile(filt_reads_exist[1:min(10, length(filt_reads_exist))])
dev.off()

# only reads as many samples as needed to get to 1e8 bases (can be changed), with samples randomized
err_merged_reads <- dada2::learnErrors(filt_reads_exist, multithread = nthreads, randomize = TRUE, MAX_CONSIST = 30, pool = 'pseudo', qualityType='FastqQuality')

pdf(file.path(outdir, "dada_err.pdf"))
dada2::plotErrors(err_merged_reads, nominalQ=TRUE)
dev.off()

#denoise
derepped <- dada2::derepFastq(filt_reads_exist, qualityType='FastqQuality')
dada_merged <- dada2::dada(derepped, err = err_merged_reads, pool = 'pseudo', multithread = nthreads)

#Construct ASV table
seqtab <- dada2::makeSequenceTable(dada_merged)
rownames(seqtab) <- sample.names[file.exists(filt_reads)]

## write fasta with ASV representative seqs
dna <- Biostrings::DNAStringSet(dada2::getSequences(seqtab))
names(dna) <- sprintf(paste0('ASV%0',floor(log10(length(dna))) + 1,'d'),1:length(dna))
Biostrings::writeXStringSet(dna, file.path(outdir, 'ASVs.fasta.gz'), compress = TRUE, format = 'fasta', width = 10000)

## rename asvs and write ASV table
colnames(seqtab) <- names(dna)
write.table(seqtab, file.path(outdir, 'asv.tsv'), sep='\t')

## save entire environment
save.image(file.path(outdir, 'ASVs.RData'))

