## get user-edited environmental variables output_prefix, counts_fp, etc
model_dir <- getwd()
source(file.path(model_dir, 'zip_glm.env'))
K_s <- 5
K_f <- 5
##

## set up output directory
dir.create(output_prefix)
##

counts_raw <- read.table(counts_fp,header=T,sep='\t',quote="")
counts_raw[,1] <- gsub('_kaiju.*','',gsub('.*/','',counts_raw[,1]))

counts <- reshape(counts_raw, idvar='taxon_id', timevar='file', v.names='reads', direction='wide', drop=c('percent','taxon_name'))
rownames(counts) <- paste0('taxon_id_',counts[,1])
counts <- counts[,-1]
counts[is.na(counts)] <- 0
colnames(counts) <- sub('reads.','',colnames(counts))

counts <- counts[order(apply(apply(counts, 2, function(x) x / sum(x)),1,sum),decreasing=T),order(apply(counts,2,sum),decreasing = T)]
relabund <- apply(counts, 2, function(x) x / sum(x))

countsbin <- t(as.matrix(counts))
countsbin[countsbin > 0] <- 1

countsmod <- t(as.matrix(counts))
countsmod[countsmod==0] <- 0.5
logcountsmod <- log(countsmod)
logcountsmod_a <- t(apply(logcountsmod,1,function(x) x - mean(x)))


meta <- read.table(sample_metadata_fp, header=T, sep='\t')
meta$Gender <- as.factor(meta$Gender)
# there are three factor levels for 'group', and it makes most biological sense to compare to the baseline of 'Term'
meta$Group <- relevel(as.factor(meta$Group), 'Term')
meta$Individual <- as.factor(meta$Individual)
meta$DeliveryMode <- as.factor(meta$DeliveryMode)
meta$EnteralFeeds_2mo <- as.factor(meta$EnteralFeeds_2mo)
# some values of 'Breastmilk' have extra spaces in them for some reason, so we'll clean those up
levels(meta$EnteralFeeds_2mo)[gsub(' ', '', levels(meta$EnteralFeeds_2mo)) == 'Breastmilk'] <-
  'Breastmilk'


NS    = ncol(counts)
NF    = nrow(counts)

X_s <- model.matrix(~ Gender + Group + Individual, meta) ## standardize continuous variables before placing in model
X_s <- cbind(model.matrix(~ Gender, meta), model.matrix(~ 0 + Group, meta), model.matrix(~ 0 + Individual, meta), diag(1, NS)) ## standardize continuous variables before placing in model
X_s[,-1] <- apply(X_s[,-1], 2, function(x) x - mean(x))
rownames(X_s) <- meta$Sample
X_s <- X_s[colnames(counts),]

idx_s <- c(1, 2, rep(3, nlevels(meta$Group)), rep(4, nlevels(meta$Individual)), rep(5,NS))
NSB <- max(idx_s)
NB_s <- length(idx_s)

###
X_s_full <- diag(1,NB_s)
X_s_full[(NB_s-NS+1):NB_s,] <- X_s
X_s_qr <- qr(t(X_s_full[,1:(NB_s-NS)]))
X_s_full[1:(NB_s-NS),] <- qr.R(X_s_qr)[,order(X_s_qr$pivot)]
X_s_full_inv <- solve(X_s_full)
###


taxvec <- sapply(unique(counts_raw$taxon_id), function(x) {
  counts_raw$taxon_name[counts_raw$taxon_id==x & !is.na(counts_raw$taxon_id)][[1]]
})
names(taxvec) <- paste0('taxon_id_', unique(counts_raw$taxon_id))
temp <- strsplit(taxvec, ';')
most <- max(sapply(temp,length))
parsedtax <- lapply(temp,function(x) {length(x) <- most; return(x)})
taxtab <- do.call('rbind',parsedtax)
hi <- unique(as.vector(taxtab[,1:9]))
hi <- hi[sapply(hi, function(x) sum(x == taxtab[,1:9],na.rm=T)) > 1]

X_f <- cbind(1, sapply(hi, function(x) apply(taxtab,1,function(y) x %in% y)), diag(1,NF))
X_f[,-1] <- apply(X_f[,-1], 2, function(x) x-mean(x))
X_f <- X_f[rownames(counts),]

idx_f = c(1, rep(2,ncol(X_f)-1-NF), rep(3,NF))
NFB   = max(idx_f)
NB_f  = length(idx_f)

prior_scale_p <- sqrt(exp(mean(log(apply(countsbin,2,var)[apply(countsbin,2,var) > 0]))))
prior_scale_a <- sqrt(exp(mean(log(apply(counts,2,function(x) var(log(x[x>0])))))))


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
                K_s           = K_s,
                K_f           = K_f)

abundance_init <- X_s_full %*%  MASS::ginv(X_s_full[(NB_s-NS+1):NB_s,]) %*% logcountsmod %*% t(MASS::ginv(X_f_full[(NB_f-NF+1):NB_f,-1])) %*% t(X_f_full)[-1,-1]
prevalence_init <- X_s_full %*%  MASS::ginv(X_s_full[(NB_s-NS+1):NB_s,]) %*% countsbin %*% t(MASS::ginv(X_f_full[(NB_f-NF+1):NB_f,])) %*% t(X_f_full)

inits <- list(prevalence = prevalence_init,
              abundance = abundance_init,
              multinomial_nuisance = apply(logcountsmod,1,mean))

save.image(file.path(output_prefix,'zip_glm_setup.RData'))

cmdstanr::write_stan_json(standat[1:13], file.path(output_prefix,'zip_test_data.json'))
#cmdstanr::write_stan_json(standat[14:7], file.path(output_prefix,'zip_test_data_2.json'))
#cmdstanr::write_stan_json(standat[18:length(standat)], file.path(output_prefix,'zip_test_data_3.json'))

#system(paste0('sed \'$d\' ', file.path(output_prefix,'zip_test_data_1.json'), ' | sed \'$s/$/,/\' > ', file.path(output_prefix,'zip_test_data.json')))
#system(paste0('sed \'1d\' ', file.path(output_prefix,'zip_test_data_2.json'), ' | sed \'$d\' | sed \'$s/$/,/\' >> ', file.path(output_prefix,'zip_test_data.json')))
#system(paste0('sed \'1d\' ', file.path(output_prefix,'zip_test_data_3.json'), ' >> ', file.path(output_prefix,'zip_test_data.json')))


cmdstanr::write_stan_json(inits, file.path(output_prefix,'zip_test_inits.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], file.path(model_dir,'zip_glm')))

sampling_commands <- list(hmc = paste('./zip_glm',
                                      paste0('data file=',file.path(output_prefix,'zip_test_data.json')),
                                      paste0('init=',file.path(output_prefix,'zip_test_inits.json')),
                                      'output',
                                      paste0('file=',file.path(output_prefix,'zip_test_data_samples.csv')),
                                      paste0('refresh=', 1),
                                      'method=sample algorithm=hmc',
                                      'stepsize=0.01',
                                      'engine=nuts',
                                      'max_depth=6',
                                      'adapt t0=10',
                                      'delta=0.8',
                                      'kappa=0.75',
                                      'num_warmup=200',
                                      'num_samples=200',
                                      ('opencl platform=0 device=1')[opencl],
                                      sep=' '),
                          advi = paste('./zip_glm',
                                       paste0('data file=',file.path(output_prefix,'zip_test_data.json')),
                                       paste0('init=',file.path(output_prefix,'zip_test_inits.json')),
                                       'output',
                                       paste0('file=',file.path(output_prefix,'zip_test_data_samples.csv')),
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


