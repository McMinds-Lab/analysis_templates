
library(dada2)
source('local_r.env')

merged_reads <- sort(list.files(indir, pattern=".fastq", full.names=T))
sample.names <- sapply(strsplit(basename(merged_reads), "_"), `[`, 1)

#png(file.path(outdir, "quality_profile.png"), height=600, width=600)
#plotQualityProfile(merged_reads[1:10])
#dev.off()
#save.image(file.path(outdir, 'quality_profile.png')

filt_dir <- file.path(indir, "filtered")
if (!dir.exists(filt_dir)) dir.create(filt_dir)

merged_filt <- file.path(filt_dir, "filtered", paste0(sample.names, "_M_filt.fastq"))

names(merged_filt) <- sample.names

filtered_out <- filterAndTrim(
fwd = merged_reads,
filt = merged_filt,
truncQ = 0,
rm.phix = T,
compress = T,
multithread = T
)

#provide number for multithreads in local.env

#learning errors is a key component, must be done with filtered reads

err_merged_reads <- learnErrors(merged_filt, multithread=T)

pdf(file.path(outdir, "dada_err.pdf"))
plotErrors(err_merged_reads, nominalQ=T)
dev.off()

#after dada2 has learned errors, user can denoise
derepped <- derepFastq(merged_filt)

dada_merged <- dada(derepped, err=err_merged_reads, multithread=T)

#Construct ASV

seqtab <- makeSequenceTable(dada_merged)
dim(seqtab)

write.table(seqtab, file.path(outdir, 'asv.tsv'), sep='\t')

save.image(file.path(outdir, 'ASV.png')) #change to above ".RData"

