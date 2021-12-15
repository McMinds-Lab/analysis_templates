
library(dada2)
source('local_r.env')

merged_reads <- sort(list.files(indir, pattern=".fastq", full.names=T))
sample.names <- sapply(strsplit(basename(merged_reads), "_"), `[`, 1)

png(file.path(outdir, "quality_profile.png"), height=600, width=600)
plotQualityProfile(merged_reads[1:10])
dev.off()
save.image(file.path(outdir, 'quality_profile.png')

filtmerged <- file.path(indir, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filtmerged) <- sample.names
out <- filterAndTrim(merged_reads, filtmerged, maxN=0, maxEE=2, truncQ=0, compress=T, multithread=T)

#provide number for multithreads in local.env

#learning errors is a key component, must be done with filtered reads

err_merged_reads <- learnErrors(merged_reads, multithread=T)
#learnErrors(filtmerged
#There is a bult-in function that allows the user to plot errors, just allows user to visualize the estimated error rates

pdf(file.path(outdir, 'dada_err.pdf')
plotErrors(err_merged_reads, nominalQ=T)
dev.off()

#after dada2 has learned errors, user can denoise
derepFastq #output will be for dada

dada_merged <- dada(filtmerged, err=err_merged_reads, multithread=T)

#Construct ASV

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

write.table(seqtab, file.path(outdir, 'asv.tsv'), sep='\t')

save.image(file.path(outdir, 'ASV.png') #change to above ".RData"

