## in doMPI batches:
## zinbwave to get latent vars and sample intercepts ('size factors'). (maybe only do this on a subset of batches, and summarize the results so they're fixed for all batches of the below tests. mean size factor; procrustes average of latent vars? take first sample as reference, do procrustes rotation of all other samples to match it, take the elementwise averages across all samples, and finally perform PCA to rotate the result again to have orthogonal columns?)
## pscl::zeroinfl with model "~ 1 + offset(size_factor) + study + host_species + latent_vars + is_diseased + has_sctld" to get prevalence significance. still not sure how best to use 'disease exposed' or 'apparently healthy' type samples

library("foreach")
library("MASS")
library("pscl")
library("doMPI")
## also requires RhpcBLASctl

args <- commandArgs(TRUE)

prefix   <- args[[1]]
sampledat <- args[[2]]
threshold <- as.numeric(args[[3]])
formula <- as.formula(args[[4]])
keycolumn <- as.numeric(args[[5]])
outdir <- args[[6]]

binary_only <- TRUE ## stick with this for now; leaving other code in case I want to refer to it in future

# register MPI so clones don't process until the foreach loop
cl <- startMPIcluster() 
registerDoMPI(cl)
opts <- list(preschedule=FALSE)

#maybe worth running with multiple chunks for estimates of W and offsets
files <- Sys.glob(paste0(prefix,'*'))
ref_subset <- sample(files,1)

counts <- as.matrix(read.table(gzfile(ref_subset), header=TRUE, row.names=1))
conditions <- read.table(sampledat, header=TRUE, row.names=1)

counts <- counts[,colnames(counts) %in% rownames(conditions)]
counts <- counts[apply(counts,1, \(x) sum(x>0)>1), ]

conditions <- conditions[colnames(counts),, drop=FALSE]

zfit <- zinbwave::zinbFit(counts, 
                          X = model.matrix(formula, data=conditions), 
                          K = 2, 
                          commondispersion = FALSE)

## beta are coefficients of X (samplewise model)
## gamma are coefficients of V (genewise model - intercept gives samplewise 'size factors' that could be used as offsets in other linear models)
## alpha are coefficients of W (which are samplewise latent variables)
## mu are abundance parameters
## pi are prevalence parameters
W <- zinbwave::getW(zfit)
abundance_offsets <- as.vector(zinbwave::getGamma_mu(zfit))
prevalence_offsets <- as.vector(zinbwave::getGamma_pi(zfit))

##if using zinbwave parameters directly, use this code
beta_mu <- t(zinbwave::getBeta_mu(zfit))
beta_pi <- -t(zinbwave::getBeta_pi(zfit)) ## just like in pscl::zeroinfl, the zero inflation coefficients seem to have signs swapped from their intuitive directions, for some reason

## for a 2-fold increase in odds, the corresponding logit is log(2). try this for a cutoff?
sig_increase <- (beta_pi[,keycolumn] > log(threshold))
sig_decrease <- (beta_pi[,keycolumn] < -log(threshold))
##

invisible(foreach(i=1:length(files), .options.multicore=opts) %dopar% {
  
  RhpcBLASctl::omp_set_num_threads(1)
  RhpcBLASctl::blas_set_num_threads(1)
  
  ## fix the above parameters and optimize zinbwave model for all genes using these values
  ## need to rewrite zinbwave optimizer to do this and use MPI to parallelize across multiple nodes like i did with de-kupl
  
  counts <- as.matrix(read.table(gzfile(files[[i]]), header=TRUE, row.names=1))
  counts <- counts[,colnames(counts) %in% rownames(conditions)]
  counts <- counts[apply(counts,1, \(x) sum(x>0)>1), ]
  
  res <- apply(counts, 1, function(jcounts) {
    
    tryCatch({
    
      if(sd(jcounts) > 0 & !any(is.infinite(jcounts))) {
      
        if(0 %in% jcounts) {
    
          if(binary_only) {
            
            # from library pscl
            full <- zeroinfl(jcounts ~ 1 + offset(abundance_offsets) | colData$condition + offset(prevalence_offsets), dist='negbin', model = FALSE, y = FALSE)
            redz <- zeroinfl(jcounts ~ 1 + offset(abundance_offsets) | 1 + offset(prevalence_offsets), dist='negbin', model = FALSE, y = FALSE)
          
            pz <- pchisq(2 * (logLik(full) - logLik(redz)), df=1, lower.tail=FALSE)

            return(c(pz, -full$coefficients$zero[[2]]))
          
          } else {
            
            # from library pscl
            full <- zeroinfl(jcounts ~ colData$condition + offset(abundance_offsets) | colData$condition + offset(prevalence_offsets), dist='negbin', model = FALSE, y = FALSE)
            redz <- zeroinfl(jcounts ~ colData$condition + offset(abundance_offsets) | 1 + offset(prevalence_offsets), dist='negbin', model = FALSE, y = FALSE)
            red2 <- zeroinfl(jcounts ~ 1 + offset(abundance_offsets) | 1 + offset(prevalence_offsets), dist='negbin', model = FALSE, y = FALSE)
            redc <- zeroinfl(jcounts ~ 1 + offset(abundance_offsets) | colData$condition + offset(prevalence_offsets), dist='negbin', model = FALSE, y = FALSE)
            
            pz <- pchisq(2 * (logLik(full) - logLik(redz)), df=1, lower.tail=FALSE)
            p2 <- pchisq(2 * (logLik(full) - logLik(red2)), df=2, lower.tail=FALSE)
            pc <- pchisq(2 * (logLik(full) - logLik(redc)), df=1, lower.tail=FALSE)
          
            if(pc < pvalue_threshold & pz >= pvalue_threshold) {
              return(c(min(p2,pc), full$coefficients$count[[2]]))
            } else {
              return(c(min(p2,pz), -full$coefficients$zero[[2]])) ## for some reason the 'zero' model coefficient is 'backwards' intuitively
            } ## if diff abund but not diff prev, set coefficient to actual log fc so sign is appropriate
            
          }
    
        } else {
          
          if(binary_only) {
            
            return(c(1,0))
            
          } else {
    
            # from library MASS
            full <- glm.nb(jcounts ~ colData$condition + offset(offsets), model = FALSE, y = FALSE)
            red <- glm.nb(jcounts ~ 1 + offset(offsets), model = FALSE, y = FALSE)

            return(c(pchisq(2 * (logLik(full) - logLik(red)), df=1, lower.tail=FALSE), 
                     full$coefficients[[2]]))
            
          }
    
        }
        
      } else {
        
        return(c(1,0))
        
      }
      
    }, error = \(x) return(c(1,0)))
    
  })
  
  outdf <- data.frame(kmer=rownames(counts), 
                      pvalue=res[1,],
                      diffLogOdds=res[2,],
                      as.data.frame(sapply(1:ncol(counts), function(x) counts[,x] / exp(abundance_offsets[[x]]))))
  write.table(outdf,
              file  = gzfile(file.path(outdir, paste0('results_chunk_', i, 'tsv.gz'))),
              sep   = "\t",
              quote = FALSE,
              col.names = FALSE,
              row.names = FALSE)
  
})

closeCluster(cl)
mpi.finalize()
  
