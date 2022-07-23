## in doMPI batches:
## zinbwave to get latent vars and sample intercepts ('size factors'). (maybe only do this on a subset of batches, and summarize the results so they're fixed for all batches of the below tests. mean size factor; procrustes average of latent vars? take first sample as reference, do procrustes rotation of all other samples to match it, take the elementwise averages across all samples, and finally perform PCA to rotate the result again to have orthogonal columns?)
## pscl::zeroinfl with model "~ 1 + offset(size_factor) + study + host_species + latent_vars + is_diseased + has_sctld" to get prevalence significance. still not sure how best to use 'disease exposed' or 'apparently healthy' type samples

args <- commandArgs(TRUE)

prefix   <- args[[1]]
sampledat <- args[[2]]
threshold <- as.numeric(args[[3]])
formula <- as.formula(args[[4]])
keycolumn <- as.numeric(args[[5]])
outdir <- args[[6]]

binary_only <- TRUE ## stick with this for now; leaving other code in case I want to refer to it in future

#maybe worth running with multiple chunks for estimates of W and offsets
files <- Sys.glob(paste0(prefix,'*'))
ref_subset <- sample(files,1)

counts <- as.matrix(read.table(gzfile(ref_subset), header=TRUE, row.names=1))
conditions <- read.table(sampledat, header=TRUE, row.names=1)

counts <- counts[,colnames(counts) %in% rownames(conditions)]
counts <- counts[apply(counts,1, \(x) sum(x>0)>0), ]

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

save.image(file.path(outdir,'02_id_diffs','zinbwave.RData'))

