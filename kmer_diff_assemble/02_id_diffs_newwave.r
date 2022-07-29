library(parallel)

args <- commandArgs(TRUE)

countsfile   <- args[[1]]
sampledat <- args[[2]]
threshold <- as.numeric(args[[3]])
formula <- as.formula(args[[4]])
keycolumn <- as.numeric(args[[5]])
outdir <- args[[6]]
nodenames <- strsplit(args[[7]], ' ')

nodenames_expanded <- as.vector(sapply(nodenames, \(x) rep(x,20)))
cl <- makePSOCKcluster(nodenames_expanded)

## read in sample metadata
conditions <- read.table(sampledat, header=TRUE, row.names=1)

## set chunk size for reading in counts
block_rows <- 1e8

## find number of kmers and samples in counts file
cat('finding number of kmers in input')
inconnect <- gzfile(countsfile, 'r')
incolnames <- read.table(inconnect, nrows=1)[-1]
nsamples <- length(incolnames)
nkmers <- 0
while(TRUE) {
  addedkmers <- length(readLines(inconnect, block_rows))
  if(addedkmers) {
    nkmers <- nkmers + addedkmers
  } else {
    break
  }
}
close(inconnect)
mode(nkmers) <- 'integer'
mode(nsamples) <- 'integer'
cat('done')
##

## convert counts file to delayed array
cat('create delayed array container')
dir.create(file.path(outdir,'02_id_diffs', 'hdf5_files'), recursive = TRUE)
HDF5Array::setHDF5DumpDir(file.path(outdir, '02_id_diffs', 'hdf5_files'))
DelayedArray::setAutoRealizationBackend("HDF5Array")
sink <- DelayedArray::AutoRealizationSink(c(nkmers, nsamples), type='integer')
sink_grid <- DelayedArray::RegularArrayGrid(dim(sink), spacings=c(block_rows, nsamples))
cat('container created')

cat('filling container')
inconnect <- gzfile(countsfile, 'r')
for (bid in seq_along(sink_grid)) {
  viewport <- sink_grid[[bid]]
  block <- as.matrix(read.table(inconnect, header=(bid==1), row.names=1, nrows=block_rows))
  mode(block) <- 'integer'
  sink <- write_block(sink, viewport, block)
}
close(sink)
close(inconnect)
cat('container filled; converting container')
counts <- as(sink, "DelayedArray")
cat('converted')
##

## make sure counts and conditions match
cat('filtering samples by conditions file')
filtsamplenames <- incolnames[incolnames %in% rownames(conditions)]
counts <- counts[,incolnames %in% rownames(conditions)]

cat('filtering kmers by prevalence')
counts_grid <- DelayedArray::RegularArrayGrid(dim(counts), spacings=c(1e6, length(filtsamplenames)))
incounts_keep <- unlist(clusterApply(cl, counts_grid, \(viewport) apply(read_block(counts, viewport), 1, sum(x>0)>2)))
counts <- counts[incounts_keep, ]
cat('done')

conditions <- conditions[filtsamplenames,, drop=FALSE]
##

cat('fitting model')
nfit <- NewWave::newFit(counts, 
                        X = model.matrix(formula, data=conditions), 
                        K = 2,
                        children = nodenames_expanded,
                        n_gene_par = 1000,
                        commondispersion = FALSE) ## using character vector for children is undocumented but looking into code it might work
cat('done')

## beta are coefficients of X (samplewise model)
## gamma are coefficients of V (genewise model - intercept gives samplewise 'size factors' that could be used as offsets in other linear models)
## alpha are coefficients of W (which are samplewise latent variables)
W <- NewWave::newW(nfit)
offsets <- as.vector(NewWave::newGamma(zfit))

##if using zinbwave parameters directly, use this code
beta <- t(NewWave::newBeta(zfit))

## for a 2-fold increase in odds, the corresponding logit is log(2). try this for a cutoff?
sig_increase <- rownames(counts)[beta[,keycolumn] > log(threshold)]
sig_decrease <- rownames(counts)[beta[,keycolumn] < -log(threshold)]
##

save.image(file.path(outdir,'02_id_diffs','newwave.RData'))

