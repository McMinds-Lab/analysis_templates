print(sessionInfo())
cat(paste('MASS version:', packageVersion('MASS'), '\n'))
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

rownames(seqtab) <- sapply(rownames(seqtab), function(x) samplenames$Sample[samplenames$Tag == sub('.fastq.gz','',sub('-',':',x))])

seqtab_F <- seqtab[grep('UFL',rownames(seqtab)),]

seqtab_M <- t(sapply(unique(sub('_.*','',rownames(seqtab_F))), function(x) apply(seqtab_F[grep(x,rownames(seqtab_F)),], 2, sum)))

metadat$id_argaly <- id_conversion$id_argaly[match(metadat$ID..Alison, id_conversion$id_alison)]

m2 <- metadat[match(rownames(seqtab_M),metadat$id_argaly),]

taxid <- read.table(taxid_fp,sep='\t',row.names = 1, header=TRUE)

## set up output directory
dir.create(file.path(outdir, 'zip_glm'))
##

counts <- t(seqtab_M)
counts <- counts[apply(counts,1,sd) > 0,]
counts_rlog <- DESeq2::rlog(counts)
sampleOrder <- order(apply(abs(cov2cor(crossprod(t(apply(counts_rlog,1,function(x) x-mean(x)))))),2,mean), decreasing=TRUE)
featureOrder <- order(apply(abs(cov2cor(tcrossprod(apply(counts_rlog,2,function(x) x-mean(x))))),2,mean), decreasing=TRUE)
counts <- counts[featureOrder,sampleOrder]

relabund <- apply(counts, 2, function(x) x / sum(x))

countsbin <- t(as.matrix(counts))
countsbin[countsbin > 0] <- 1

countsmod <- t(as.matrix(counts))
countsmod[countsmod==0] <- 0.5
logcountsmod <- log(countsmod)

rownames(m2) <- m2$id_argaly
m2 <- m2[colnames(counts)[colnames(counts) %in% rownames(m2)],]
m2$Location <- as.factor(m2$Location)

NS    = ncol(counts)
NF    = nrow(counts)

X_s <- cbind(Intercept=1, model.matrix(~ 0 + Location, m2)) ## standardize continuous variables before placing in model
X_s[,-1] <- apply(X_s[,-1], 2, function(x) x - mean(x))
rownames(X_s) <- rownames(m2)
X_s <- X_s[colnames(counts),]

idx_s <- c(1, rep(2, nlevels(m2$Location)))
NSB <- max(idx_s)
NB_s <- length(idx_s)

## work on this section to make sure columns are unique and so that taxa are filtered out if they apply to ALL samples as well as none
hi <- cbind('Intercept',unique(taxid))

estimables <- lapply(2:(ncol(hi)-1), function(x) {
  y <- unique(hi[,1:x])
  keepers <- apply(y, 1, function(z) sum(apply(y[,1:(x-1),drop=F],1, function(r) identical(z[1:(x-1)],r))) > 1 & !is.na(z[[x]]))
  return(unname(y[keepers,x]))
})

X_f <- cbind(Intercept=1, do.call(cbind, sapply(1:length(estimables), function(x) sapply(estimables[[x]], function(y) as.numeric(taxid[,x] == y)))))
rownames(X_f) <- rownames(taxid)
X_f <- X_f[rownames(counts),]
X_f[is.na(X_f)] <- 0
mode(X_f) <- 'numeric'
X_f_nonUnique <- colSums(X_f)>1 & !duplicated(t(X_f))
X_f <- X_f[,X_f_nonUnique]
X_f[,-1] <- apply(X_f[,-1], 2, function(x) x-mean(x))
##

idx_f = c(1, rep(2,ncol(X_f)-1))
NFB   = max(idx_f)
NB_f  = length(idx_f)

prior_scale_p <- sqrt(exp(mean(log(apply(countsbin,2,var)[apply(countsbin,2,var) > 0]))))
prior_scale_a <- sqrt(exp(mean(log(apply(counts,2,function(x) var(log(x[x>0])))))))

qr.xs <- qr(X_s[,-1])
rank_X_s <- qr.xs$rank
R_s <- qr.R(qr.xs)
R_inv_s <- MASS::ginv(R_s)[,1:rank_X_s]

if(rank_X_s < ncol(X_s[,-1])) {
  I_diff_AA_reduced_s <- diag(1,nrow(R_inv_s)) - R_inv_s %*% R_s[1:rank_X_s,]
  qr.is <- qr(I_diff_AA_reduced_s)
  top_rows <- sort(order(abs(diag(qr.R(qr.is))),decreasing=T)[1:(NB_s-rank_X_s)])
  Q_I_s <- qr.Q(qr.is)[, ]
  
  X_sR_inv <- cbind(R_inv_s,Q_I_s)
} else {
  X_sR_inv <- R_inv_s
}

qr.xf <- qr(X_f[,-1])
rank_X_f <- qr.xf$rank
R_f <- qr.R(qr.xf)
R_inv_f <- MASS::ginv(R_f)[,1:rank_X_f]

if(rank_X_f < ncol(X_f[,-1])) {
  I_diff_AA_reduced_f <- diag(1,nrow(R_inv_f)) - R_inv_f %*% R_f[1:rank_X_f,]
  qr.if <- qr(I_diff_AA_reduced_f)
  top_rows <- sort(order(abs(diag(qr.R(qr.if))),decreasing=T)[1:(NB_f-rank_X_f)])
  Q_I_s <- qr.Q(qr.if)[,top_rows]
  
  X_fR_inv <- cbind(R_inv_f,Q_I_f)
} else {
  X_fR_inv <- R_inv_f 
}

standat <- list(NS            = NS,
                NB_s          = NB_s,
                NSB           = NSB,
                idx_s         = idx_s,
                X_s           = X_s,
                X_sR_inv      = X_sR_inv,
                NF            = NF,
                NB_f          = NB_f,
                NFB           = NFB,
                idx_f         = idx_f,
                X_f           = t(X_f),
                X_fR_inv      = t(X_fR_inv),
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
              multinomial_nuisance    = apply(counts,2,function(x) log(mean(x))),
              L_s                     = diag(1,NS,K_s))

save.image(file.path(outdir, 'zip_glm', 'zip_glm_setup.RData'))

cmdstanr::write_stan_json(standat, file.path(outdir, 'zip_glm', 'zip_test_data.json'))
cmdstanr::write_stan_json(inits, file.path(outdir, 'zip_glm', 'zip_test_inits.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], file.path(model_dir,'zip_glm_qr')))

sampling_commands <- list(hmc = paste('./zip_glm_qr',
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

#stan.fit.var <- cmdstanr::read_cmdstan_csv(Sys.glob(path.expand(file.path(outdir,'zip_glm','zip_test_data_samples_*.csv'))),
#                                           format = 'draws_array')

#summary(stan.fit.var$post_warmup_sampler_diagnostics)
#plot(apply(stan.fit.var$post_warmup_draws[,1,paste0('L_s[',1:NS,',1]')], 3, mean), apply(stan.fit.var$post_warmup_draws[,1,paste0('L_s[',1:NS,',2]')],3,mean), xlab = "PCA1", ylab = "PCA2",axes = TRUE, main = "First samplewise latent variables", col=as.factor(m2$env.features), pch=16)
#taxid[do.call(order, as.data.frame(taxid)),]
