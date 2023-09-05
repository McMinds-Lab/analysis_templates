print(sessionInfo())
cat(paste('DESeq2 version:', packageVersion('DESeq2'), '\n'))
cat(paste('cmdstanr version:', packageVersion('cmdstanr'), '\n'))
cat(paste('cmdstan version:', cmdstanr::cmdstan_version(), '\n'))

## get user-edited environmental variables outdir, counts_fp, etc
newargs <- commandArgs(TRUE)

load(newargs[[1]])
taxid_fp <- newargs[[2]]
samplenames <- read.table(newargs[[3]],sep='\t',header=T)
metadat <- read.table(newargs[[4]],sep='\t',header=T)
id_conversion <- read.table(newargs[[5]],sep='\t',header=T)

outdir <- newargs[[6]]
K_s <- as.numeric(newargs[[7]])
nchains <- as.numeric(newargs[[8]])
nthreads <- as.numeric(newargs[[9]])
opencl <- as.logical(newargs[[10]])
opencl_device <- as.numeric(newargs[[11]])
model_dir <- newargs[[12]]
algorithm <- newargs[[13]]
##

rownames(seqtab) <- sapply(rownames(seqtab), \(x) samplenames$Sample[samplenames$Tag == sub('.fastq.gz','',sub('-',':',x))])
seqtab_F <- seqtab[grep('UFL',rownames(seqtab)),]
seqtab_M <- t(sapply(unique(sub('_.*','',rownames(seqtab_F))), \(x) apply(seqtab_F[grep(x,rownames(seqtab_F)),], 2, sum)))

metadat$id_argaly <- id_conversion$id_argaly[match(metadat$ID..Alison, id_conversion$id_alison)]
m2 <- metadat[match(rownames(seqtab_M),metadat$id_argaly),]

taxid <- read.table(taxid_fp,sep='\t',row.names = 1, header=TRUE)

## set up output directory
dir.create(file.path(outdir, 'zip_glm'))
##

counts <- t(seqtab_M)
counts <- counts[apply(counts,1,sd) > 0,]

## merge all taxa with a single occurrence under a given taxonomic level
taxid2 <- taxid[rownames(counts),]
taxid2[is.na(taxid2)] <- 'NA'
taxid2 <- cbind(taxid2, ASV=rownames(taxid2))
for(col in 2:ncol(taxid2)) {taxid2[,col] <- apply(taxid2[,(col-1):col],1,paste,collapse=';')}
newtax <- list()
newtax[[ncol(taxid2)]] <- as.list(setNames(rownames(taxid2),rownames(taxid2)))
newcounts <- list()
newcounts[[ncol(taxid2)]] <- counts
for(level in (ncol(taxid2)-1):1) {
  newtax[[level]] <- list()
  newcounts[[level]] <- matrix(NA,0,ncol(counts))
  for(tax in unique(taxid2[,level])) {
    asvs <- rownames(taxid2)[taxid2[,level] == tax]
    asvs <- asvs[asvs %in% rownames(newcounts[[level+1]])]
    if(length(asvs) > 1) {
      doMerge <- apply(newcounts[[level+1]][asvs,] > 0, 1, sum) <= 1
      for(i in which(doMerge)) {
        uTax <- taxid2[asvs[doMerge],level+1,drop=F]
        if(sum(uTax[,1] == uTax[asvs[i],]) > 1) {
          doMerge[i] <- FALSE
        }
      }
      if(sum(doMerge) == 1) {
        if(length(asvs) == 2) {
          doMerge[!doMerge] <- TRUE
        } else {
          doMerge[doMerge] <- FALSE
        }
      }
      if(any(!doMerge)) {
        newtax[[level]][asvs[!doMerge]] <- newtax[[level+1]][asvs[!doMerge]]
        newcounts[[level]] <- rbind(newcounts[[level]], newcounts[[level+1]][asvs[!doMerge],])
        rownames(newcounts[[level]])[(nrow(newcounts[[level]])-sum(!doMerge)+1):nrow(newcounts[[level]])] <- asvs[!doMerge]
      }
      if(sum(doMerge) > 1) {
        newname <- sample(asvs[doMerge],1)
        newtax[[level]][[newname]] <- unlist(newtax[[level+1]][asvs[doMerge]])
        newcounts[[level]] <- rbind(newcounts[[level]], apply(newcounts[[level+1]][asvs[doMerge],], 2, sum))
        rownames(newcounts[[level]])[nrow(newcounts[[level]])] <- newname
      }
    } else if(length(asvs) == 1) {
      newtax[[level]][[asvs]] <- newtax[[level+1]][[asvs]]
      newcounts[[level]] <- rbind(newcounts[[level]], newcounts[[level+1]][asvs,])
      rownames(newcounts[[level]])[nrow(newcounts[[level]])] <- asvs
    }
  }
}
counts <- newcounts[[1]]

rownames(m2) <- m2$id_argaly
m2 <- m2[colnames(counts)[colnames(counts) %in% rownames(m2)],]
m2$Location <- as.factor(m2$Location)

NS <- ncol(counts)
NF <- nrow(counts)

X_s <- cbind(Intercept=1, model.matrix(~ 0 + Location, m2)) ## standardize continuous variables before placing in model
X_s[,-1] <- apply(X_s[,-1], 2, \(x) x - mean(x))
rownames(X_s) <- rownames(m2)
X_s <- X_s[colnames(counts),]

idx_s <- c(1, rep(2, nlevels(m2$Location)))
NSB <- max(idx_s)
NB_s <- length(idx_s)

## work on this section to make sure columns are unique and so that taxa are filtered out if they apply to ALL samples as well as none
int_taxid <- cbind('Intercept', unique(taxid))
estimables <- lapply(2:(ncol(int_taxid)-1), \(x) {
  y <- unique(int_taxid[,1:x])
  keepers <- apply(y, 1, \(z) sum(apply(y[,1:(x-1),drop=F],1, \(r) identical(z[1:(x-1)],r))) > 1 & !is.na(z[[x]]))
  return(unname(y[keepers,x]))
})

X_f <- cbind(Intercept=1, do.call(cbind, sapply(1:length(estimables), \(x) sapply(estimables[[x]], \(y) as.numeric(taxid[,x] == y)))))
rownames(X_f) <- rownames(taxid)
X_f <- X_f[rownames(counts),]
X_f[is.na(X_f)] <- 0
mode(X_f) <- 'numeric'
X_f_nonUnique <- colSums(X_f)>1 & !duplicated(t(X_f))
X_f <- X_f[,X_f_nonUnique]
X_f[,-1] <- apply(X_f[,-1], 2, \(x) x-mean(x))
##

## order the samples and features such that the ones that seem to define latent dims are first to be constrained in the cholesky factor
counts_rlog <- DESeq2::rlog(counts)
counts_rlog_centered <- t(apply(apply(counts_rlog, 2, \(x) x-mean(x)), 1, \(x) x-mean(x)))
counts_rlog_standard <- t(apply((apply(apply(counts_rlog_centered, 1, \(x) x/sd(x)), 1, \(x) x/sd(x))), 1, \(x) x/sd(x)))
counts_rlog_residuals <- counts_rlog_standard - X_f %*% MASS::ginv(X_f) %*% counts_rlog_standard %*% t(X_s %*% MASS::ginv(X_s))
pr <- prcomp(counts_rlog_residuals)$rotation
sampleOrder <- vector('numeric')
sampsLeft <- rownames(pr)
for(i in 1:(ncol(pr)-1)) {
  winner <- which.max(apply(pr[sampsLeft,i:ncol(pr)],1, \(x) abs(x[1]) / sqrt(sum(x^2))))
  sampleOrder <- c(sampleOrder, sampsLeft[winner])
  sampsLeft <- sampsLeft[-winner]
}
sampleOrder <- c(sampleOrder, sampsLeft)
counts <- counts[,sampleOrder]
X_s <- X_s[colnames(counts),]
##

idx_f = c(1, rep(2, ncol(X_f)-1))
NFB   = max(idx_f)
NB_f  = length(idx_f)

countsbin <- t(as.matrix(counts))
countsbin[countsbin > 0] <- 1

prior_scale_p <- sqrt(exp(mean(log(apply(countsbin,2,var)[apply(countsbin,2,var) > 0]))))
counts_clr <- apply(log(counts), 2, \(x) x-mean(x[!is.infinite(x)]))
logfscales <- log(apply(counts_clr,1,\(x) var(x[!is.infinite(x)])))
prior_scale_a <- sqrt(exp(mean(logfscales[!is.na(logfscales)])))

standat <- list(NS            = NS,
                NB_s          = NB_s,
                NSB           = NSB,
                idx_s         = idx_s,
                X_s           = X_s,
                NF            = NF,
                NB_f          = NB_f,
                NFB           = NFB,
                idx_f         = idx_f,
                X_f           = t(X_f),
                count         = counts,
                prior_scale_a = prior_scale_a,
                prior_scale_p = prior_scale_p,
                K_s           = K_s)

inits <- list(global_scale_prevalence = 0.1,
              global_scale_abundance  = 0.1,
              sd_prevalence_norm      = rep(0.1, (NSB+2)*(NFB+1)-2),
              sd_abundance_norm       = rep(0.1, (NSB+2)*NFB-1),
              beta_prevalence_i       = matrix(0,NB_s+K_s,NB_f),
              beta_prevalence_s       = matrix(0,NB_s+K_s,NF),
              beta_prevalence_f       = matrix(0,NS,NB_f),
              beta_abundance_i        = matrix(0,NB_s+K_s,NB_f-1),
              beta_abundance_s        = matrix(0,NB_s+K_s,NF),
              beta_abundance_f        = matrix(0,NS,NB_f-1),
              residuals               = matrix(0,NS,NF),
              multinomial_nuisance    = apply(counts,2, \(x) weighted.mean(c(log(mean(x)), mean(log(x[x>0]))),c(sum(x==0),sum(x>0)))),
              L_s                     = diag(1,NS,K_s))

relabund <- apply(counts, 2, \(x) x / sum(x))

save.image(file.path(outdir, 'zip_glm', 'zip_glm_setup.RData'))

cmdstanr::write_stan_json(standat, file.path(outdir, 'zip_glm', 'zip_test_data.json'))
cmdstanr::write_stan_json(inits, file.path(outdir, 'zip_glm', 'zip_test_inits.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], file.path(model_dir,'zip_glm')))

sampling_commands <- list(hmc = paste('./zip_glm',
                                      paste0('data file=',path.expand(file.path(outdir, 'zip_glm', 'zip_test_data.json'))),
                                      paste0('init=',path.expand(file.path(outdir, 'zip_glm', 'zip_test_inits.json'))),
                                      'output',
                                      paste0('file=',path.expand(file.path(outdir, 'zip_glm', 'zip_test_data_samples.csv'))),
                                      paste0('refresh=', 1),
                                      'method=sample',
                                      paste0('num_chains=',nchains),
                                      'algorithm=hmc',
                                      'stepsize=0.1',
                                      'engine=nuts',
                                      'max_depth=12',
                                      'adapt t0=10',
                                      'delta=0.8',
                                      'kappa=0.75',
                                      'num_warmup=1000',
                                      'num_samples=1000',
                                      paste0('num_threads=',nthreads),
                                      (paste0('opencl platform=0 device=', opencl_device))[opencl],
                                      sep=' '),
                          advi = paste('./zip_glm',
                                       paste0('data file=',path.expand(file.path(outdir, 'zip_glm', 'zip_test_data.json'))),
                                       'init=0',
                                       'output',
                                       paste0('file=',path.expand(file.path(outdir, 'zip_glm', 'zip_test_data_samples.csv'))),
                                       paste0('refresh=', 100),
                                       'method=variational algorithm=meanfield',
                                       #'grad_samples=1',
                                       #'elbo_samples=1',
                                       'iter=20000',
                                       'eta=0.1',
                                       'adapt engaged=0',
                                       'tol_rel_obj=0.01',
                                       #'eval_elbo=1',
                                       'output_samples=1000',
                                       (paste0('opencl platform=0 device=', opencl_device))[opencl],
                                       sep=' '))

setwd(model_dir)
print(sampling_commands[[algorithm]])
print(date())
system(sampling_commands[[algorithm]])

#stan.fit.var <- cmdstanr::read_cmdstan_csv(Sys.glob(path.expand(file.path(outdir,'zip_glm','zip_test_data_samples_*.csv'))), format = 'draws_array', variables=c('global_scale_prevalence','global_scale_abundance','sd_prevalence_norm','sd_abundance_norm','sd_resid_s','sd_resid_f','L_s'))

#summary(stan.fit.var$post_warmup_sampler_diagnostics)
#summary(stan.fit.var$post_warmup_draws[,,grep('L_s',dimnames(stan.fit.var$post_warmup_draws)[[3]])])

#sd_prevalence_norm <- matrix(c(100,as.data.frame(summary(stan.fit.var$post_warmup_draws[,,grep('sd_pre',dimnames(stan.fit.var$post_warmup_draws)[[3]])]))$mean,1),ncol=NFB+1)
#sd_prevalence_norm^2 / (sum(sd_prevalence_norm^2) - sd_prevalence_norm[1,1]^2)
#sd_abundance_norm <- matrix(c(as.data.frame(summary(stan.fit.var$post_warmup_draws[,,grep('sd_abu',dimnames(stan.fit.var$post_warmup_draws)[[3]])]))$mean,1),ncol=NFB)
#sd_abundance_norm^2 / sum(sd_abundance_norm^2)

L_s <- cbind(apply(stan.fit.var$post_warmup_draws[,1,paste0('L_s[',1:NS,',1]')], 3, mean), apply(stan.fit.var$post_warmup_draws[,1,paste0('L_s[',1:NS,',2]')],3,mean))
pcs <- princomp(L_s)
#plot(pcs$scores, xlab = "PCA1", ylab = "PCA2",axes = TRUE, main = "First samplewise latent variables", col=as.factor(m2$env.features), pch=16)
#points(pcs$scores[1,1],pcs$scores[1,2],pch=24)
#points(pcs$scores[2,1],pcs$scores[2,2],pch=25)

#taxid[do.call(order, as.data.frame(taxid)),]
#  sdnames_prev <- paste0(as.vector(outer(c('Overall_prevalence','Location','Latent','Samples'),c('Richness','Taxonomy','ASVs'),paste,sep='.'))[-c(1,4*3)],'.p')
#  sdnames_abund <- paste0(as.vector(outer(c('Overall_abundance','Location','Latent','Samples'),c('Taxonomy','ASVs'),paste,sep='.'))[-(4*2)],'.a')
#  hi <- stan.fit.var$post_warmup_draws[,1,grep('sd.*_norm',dimnames(stan.fit.var$post_warmup_draws)[[3]]), drop=TRUE]
#  colnames(hi) <- c(sdnames_prev,sdnames_abund)
#  par(mar=c(10, 4, 4, 2) + 0.1); boxplot(hi,las=2)
