## positional arguments:
## 1) filepath to metadata info (tsv with sample name in first column; fwd index name in 4th, and rev index in 5th)
## 2) filepath to Schloss_lab_indices.txt
## 3) comma-separated list of fwd primer names to filter on (in quotes) (keep anything that has at least one of these, or keep everything if empty string)
## 4) comma-separated list of rev primer names to filter on (in quotes) (keep anything that has at least one of these, or keep everything if empty string)
## 5) output path

args <- commandArgs(TRUE)

meta <- read.table(args[[1]], sep='\t', header=TRUE)
idx_tab <- read.table(args[[2]], sep=' ', header=FALSE)
fwd_filter <- strsplit(args[[3]],',')[[1]]
rev_filter <- strsplit(args[[4]],',')[[1]]
out_fp <- args[[5]]

meta_filt <- meta
if(length(fwd_filter) > 0) {
  meta_filt <- meta_filt[meta_filt[,2] %in% fwd_filter,]
}

if(length(rev_filter) > 0) {
  meta_filt <- meta_filt[meta_filt[,3] %in% rev_filter,]
}

res <- sapply(meta_filt[,1], \(x) {
  
  i <- match(x, meta_filt[,1])
  fwd <- idx_tab[match(meta_filt[i,4], idx_tab[,1]),2]
  rev_revcomp <- paste(rev(strsplit(chartr("ACGTRYSWKMBVDHacgtryswkmbvdh", "TGCAYRSWMKVBHDtgcayrswmkvbhd", idx_tab[match(meta_filt[i,5], idx_tab[,1]),2]), NULL)[[1]]), collapse = "")
  
  paste0('^', fwd, '...', rev_revcomp, '$')
  
}, simplify = FALSE)

outvec <- unlist(rbind(paste('>',names(res), sep=''), res))

file_conn <- file(out_fp)
writeLines(outvec, file_conn)
close(file_conn)
