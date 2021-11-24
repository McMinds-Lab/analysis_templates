library(dada2)

#some people run packageVersion(dada2) after to see current version, however it's not helpful

path <- "path1"

#need to change "path1" to where the data can be found

forward_reads <- sort(list.files(path, pattern="R1.fastq", full.names=T))

reverse_reads <- sort(list.files(path, pattern="R2.fastq", full.names=T))

#I believe that dada2 might want us to run its own filter and trim function, however we might get away without using it since prevous cutting will be done

#The following is to assign the filtered reads a file

filtered_forward_reads <- file.path(path, "filtered", paste0(sample.names, "_forward_filtered.fastq.gz"))

filtered_reverse_reads <- file.path(path, "filtered", paste0(sample.names, "_reverse_filtered.fastq.gz"))

out <- filteredAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=T, compress=T, multithread=T)

#learning errors is a key component, must be done with filtered reads

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=T)

err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=T)

#There is a bult-in function that allows the user to plot errors, just allows user to visualize the estimated error rates
#If decided to run this plot-

plotErrors(err_forward_reads, nominalQ=T)

plotErrors(err_reverse_reads, nominalQ=T)

#after dada2 has learned errors, user can denoise

dada_forward <- dada(filtered_forward_reads, err=err_forward_reads, multithread=T)

dada_reverse <- dada(Filtered_reverse_reads, err-err_reverse_reads, multithread=T)

#Inspect the denoising results, if wanted by "dada_forward[[1]]

#After denoising, we can merge the paired reads

mergers <- mergePairs(dada_forward, filtered_forward_reads, dada_reverse, filtered_reverse_reads, verbose=T)

#Again, if wanted to inspect mergers run head(mergers[[1]])

#Construct ASV

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


