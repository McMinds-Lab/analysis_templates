
args <- commandArgs(TRUE)

countsfile   <- args[[1]]
sampledat <- args[[2]]
threshold <- as.numeric(args[[3]])
formula <- as.formula(args[[4]])
keycolumn <- as.numeric(args[[5]])
outdir <- args[[6]]
n_cores <- as.numeric(args[[7]])

## read in sample metadata
conditions <- read.table(sampledat, header=TRUE, row.names=1)

## set chunk size for reading in counts
block_rows <- 1e6

## find number of kmers and samples in counts file
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
##

## convert counts file to delayed array
DelayedArray::setAutoRealizationBackend("HDF5Array")
sink <- DelayedArray::AutoRealizationSink(c(nkmers, nsamples), dimnames=list(paste0('line_',1:nkmers),incolnames), type='integer')
sink_grid <- DelayedArray::RegularArrayGrid(dim(sink), spacings=c(block_rows, nsamples))

inconnect <- gzfile(countsfile, 'r')
for (bid in seq_along(sink_grid)) {
  viewport <- sink_grid[[bid]]
  block <- as.matrix(read.table(inconnect, header=(bid==1), row.names=1, nrows=block_rows))
  mode(block) <- 'integer'
  sink <- write_block(sink, viewport, block)
}
close(sink)
close(inconnect)
counts <- as(sink, "DelayedArray")
##

## make sure counts and conditions match
counts <- counts[,colnames(counts) %in% rownames(conditions)]
counts <- counts[apply(counts,1, \(x) sum(x>0)>0), ]

conditions <- conditions[colnames(counts),, drop=FALSE]
##

nfit <- NewWave::newFit(counts, 
                        X = model.matrix(formula, data=conditions), 
                        K = 2,
                        children=n_cores,
                        n_gene_par = 1000,
                        commondispersion = FALSE)

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

