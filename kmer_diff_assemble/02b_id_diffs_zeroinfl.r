## in doMPI batches:
## zinbwave to get latent vars and sample intercepts ('size factors'). (maybe only do this on a subset of batches, and summarize the results so they're fixed for all batches of the below tests. mean size factor; procrustes average of latent vars? take first sample as reference, do procrustes rotation of all other samples to match it, take the elementwise averages across all samples, and finally perform PCA to rotate the result again to have orthogonal columns?)
## pscl::zeroinfl with model "~ 1 + offset(size_factor) + study + host_species + latent_vars + is_diseased + has_sctld" to get prevalence significance. still not sure how best to use 'disease exposed' or 'apparently healthy' type samples

library("MASS")
library("pscl")

args <- commandArgs(TRUE)

outdir   <- args[[1]]
taskID <- as.numeric(args[[2]])

load(file.path(outdir,'02_id_diffs','zinbwave.RData'))

## fix the above parameters and optimize zinbwave model for all genes using these values
## need to rewrite zinbwave optimizer to do this and use MPI to parallelize across multiple nodes like i did with de-kupl
## for now just use the offsets in pscl::zeroinfl models

counts <- as.matrix(read.table(gzfile(files[[taskID + 1]]), header=TRUE, row.names=1))
counts <- counts[,colnames(counts) %in% rownames(conditions)]
counts <- counts[apply(counts,1, \(x) sum(x>0)>0), ]

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
            file  = gzfile(file.path(outdir, '02_id_diffs', 'temp', paste0('results_chunk_', taskID, '.tsv.gz'))),
            sep   = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

outfilt <- outdf[outdf$pvalue < threshold,]
write.table(outfilt,
            file  = gzfile(file.path(outdir, '02_id_diffs', 'temp', paste0('sig_results_chunk_', taskID, '.tsv.gz'))),
            sep   = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
pos <- outfilt[outfilt$diffLogOdds > 0,]
write.table(pos,
            file  = gzfile(file.path(outdir, '02_id_diffs', 'temp', paste0('pos_sig_results_chunk_', taskID, '.tsv.gz'))),
            sep   = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)  
neg <- outfilt[outfilt$diffLogOdds < 0,]
write.table(neg,
            file  = gzfile(file.path(outdir, '02_id_diffs', 'temp', paste0('neg_sig_results_chunk_', taskID, '.tsv.gz'))),
            sep   = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)  
