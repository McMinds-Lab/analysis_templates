library(parallel)

args <- commandArgs(TRUE)

countsfile   <- args[[1]]
sampledat <- args[[2]]
threshold <- as.numeric(args[[3]])
formula <- as.formula(args[[4]])
keycolumn <- as.numeric(args[[5]])
outdir <- args[[6]]
input_children <- strsplit(args[[7]], ' ')[[1]]

if(any(is.na(suppressWarnings(as.numeric(input_children))))) {
  cat('interpreting "children" option to be node names, not number of processes\n')
  firstnodesize <- as.numeric(system(paste0('sinfo --Node --nodes ', input_children[[1]],' --format %c | tail -1'), internal=TRUE))
  children <- rep(input_children[[1]], ceiling(firstnodesize / 2))
  if(length(input_children) > 1) {
    children <- c(children, as.vector(sapply(input_children[-1], \(x) {
      nodesize <- as.numeric(system(paste0('sinfo --Node --nodes ', x,' --format %c | tail -1'), internal=TRUE))
      rep(x, max(1, nodesize - 1))
    })))
  }
  n_children <- length(children)
} else {
  cat('interpreting "children" option to be number of processes on single node, not node names\n')
  children <- as.numeric(input_children)
  n_children <- children
}

## read in sample metadata
conditions <- read.table(sampledat, header=TRUE, row.names=1)

inconnect <- gzfile(countsfile, 'r')
incolnames <- unlist(read.table(inconnect, nrows=1)[-1])
close(inconnect)

if(file.exists(file.path(outdir,'02_id_diffs', 'hdf5_files', 'auto00001.h5'))) {
  
  counts <- HDF5Array::HDF5Array(file.path(outdir,'02_id_diffs', 'hdf5_files', 'auto00001.h5'), 'HDF5ArrayAUTO00001', as.sparse=TRUE, type='integer')
  
} else {
  
  ## find number of kmers and samples in counts file
  cat('finding number of kmers in input\n')
  nkmers <- as.numeric(scan(text = system(paste0('zcat ', countsfile, ' | wc -l'), intern = TRUE), what = '')[[1]]) - 1
  nsamples <- length(incolnames)
  mode(nkmers) <- 'integer'
  mode(nsamples) <- 'integer'
  cat('done\n')
  ##
  
  ## convert counts file to delayed array
  cat('create delayed array container\n')
  block_rows <- ceiling(nkmers / (n_children * 10)) ## make sure there are at least as many blocks as tasks so the workers don't have to read in too much data from mismatching blocks
  dir.create(file.path(outdir, '02_id_diffs', 'hdf5_files'), recursive = TRUE)
  HDF5Array::setHDF5DumpDir(file.path(outdir, '02_id_diffs', 'hdf5_files'))
  DelayedArray::setAutoRealizationBackend("HDF5Array")
  sink <- DelayedArray::AutoRealizationSink(c(nkmers, nsamples), type='integer', as.sparse = TRUE)
  sink_grid <- DelayedArray::RegularArrayGrid(dim(sink), spacings = c(block_rows, nsamples))
  cat('container created\n')
  
  cat('filling container\n')
  inconnect <- gzfile(countsfile, 'r')
  for (bid in seq_along(sink_grid)) {
    viewport <- sink_grid[[bid]]
    block <- as.matrix(read.table(inconnect, header=(bid==1), row.names=1, nrows=block_rows))
    mode(block) <- 'integer'
    sink <- DelayedArray::write_block(sink, viewport, block)
  }
  close(sink)
  close(inconnect)
  cat('container filled\n')
  cat('converting container\n')
  counts <- as(sink, "DelayedArray")
  cat('converted\n')

}

cat('size of input matrix:\n')
print(dim(counts))
##

## make sure counts and conditions match
cat('filtering samples by conditions file\n')
filtsamplenames <- incolnames[incolnames %in% rownames(conditions)]
cat('remaining samples:\n')
print(filtsamplenames)

counts <- counts[,incolnames %in% rownames(conditions)]

cat('filtering kmers by prevalence\n')
cl <- makePSOCKcluster(children)
clusterExport(cl, 'counts')
counts_grid <- DelayedArray::RegularArrayGrid(dim(counts), spacings=c(1e6, length(filtsamplenames)))
cat('  determining which kmers to keep\n')
incounts_keep <- unlist(clusterApplyLB(cl, counts_grid, \(viewport) apply(DelayedArray::read_block(counts, viewport), 1, \(x) sum(x>0)>2)))
stopCluster(cl)
cat('  filtering\n')
counts <- counts[incounts_keep, ]
cat('done filtering kmers\n')


cat('new dimensions:\n')
print(dim(counts))

conditions <- conditions[filtsamplenames,, drop=FALSE]
##

cat('fitting model\n')
nfit <- NewWave::newFit(counts, 
                        X = model.matrix(formula, data=conditions), 
                        K = 2,
                        children = children,
                        n_gene_par = 1000,
                        commondispersion = FALSE,
                        verbose = TRUE) ## using character vector for children is undocumented but my github repo has edits to potentially make it work
cat('done\n')

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

