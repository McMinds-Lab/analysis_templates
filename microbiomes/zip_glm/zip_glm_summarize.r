## get user-edited environmental variables output_prefix, counts_fp, etc
model_dir <- getwd()
source(file.path(model_dir, 'zip_glm.env'))
##

load(file.path(output_prefix,'zip_glm_setup.RData'))


stan.fit.var <- cmdstanr::read_cmdstan_csv(file.path(output_prefix,'zip_test_data_samples.csv'),
                                           format = 'draws_array',
                                           variables = c('sd_prevalence', 'sd_abundance'))

boxplot(stan.fit.var$draws[,1,-1,drop=T], las=2, cex.axis=0.5)



stan.fit.ls <- cmdstanr::read_cmdstan_csv(file.path(output_prefix,'zip_test_data_samples.csv'),
                                           format = 'draws_array',
                                           variables = 'L_s')
L_s <- matrix(apply(stan.fit.ls$draws,3,mean),ncol=K_s)
now <- princomp(L_s)
plot(now$scores, col = c('red', 'blue')[meta$DeliveryMode])
plot(now$scores, col = c('red', 'blue', 'green')[meta$EnteralFeeds_2mo])

stan.fit.bpi <- cmdstanr::read_cmdstan_csv(file.path(output_prefix,'zip_test_data_samples.csv'),
                                           format = 'draws_array',
                                           variables = 'beta_prevalence_i')

stan.fit.bpi.expanded <- array(stan.fit.bpi$draws, dim=c(dim(stan.fit.bpi$draws)[[1]], NB_s+K_s, NB_f+K_f))
dimnames(stan.fit.bpi.expanded) <- list(NULL, c(colnames(X_s), paste0('latent_s',1:K_s)), c(colnames(X_f), paste0('latent_f',1:K_f)))

stan.fit.bai <- cmdstanr::read_cmdstan_csv(file.path(output_prefix,'zip_test_data_samples.csv'),
                                           format = 'draws_array',
                                           variables = 'beta_abundance_i')
stan.fit.bai.expanded <- array(stan.fit.bai$draws, dim=c(dim(stan.fit.bai$draws)[[1]], NB_s+K_s, NB_f+K_f))
dimnames(stan.fit.bai.expanded) <- list(NULL, c(colnames(X_s), paste0('latent_s',1:K_s)), c(colnames(X_f), paste0('latent_f',1:K_f)))


save.image(file.path(output_prefix, 'results.RData'))

monteCarloP <- function(x, pn='p') {
  if(pn == 'n') {
    res <- (1 + sum(x >= 0, na.rm=TRUE)) / (1 + length(x))
  } else if(pn == 'p') {
    res <- (1 + sum(x <= 0, na.rm=TRUE)) / (1 + length(x))
  }
} # https://arxiv.org/pdf/1603.05766.pdf

bpi.pp <- apply(stan.fit.bpi.expanded, c(2,3), monteCarloP)
bpi.pn <- apply(stan.fit.bpi.expanded, c(2,3), monteCarloP, 'n')

bai.pp <- apply(stan.fit.bai.expanded, c(2,3), monteCarloP)
bai.pn <- apply(stan.fit.bai.expanded, c(2,3), monteCarloP, 'n')


bpi.pp[1:5,bpi.pp['GroupEarly and subsequent abx',] < 0.05]

bai.pp[1:5,bai.pp['GroupEarly and subsequent abx',] < 0.05]
bai.pn[1:5,bai.pn['GroupEarly and subsequent abx',] < 0.05]

relabund_0na <- relabund
relabund_0na[relabund_0na == 0] <- NA

meta$abundance_Escherichia <- apply(relabund_0na[X_f[,'Escherichia'] > 0,], 2, mean, na.rm=T)

beanplot::beanplot(
  meta$abundance_Escherichia ~ meta$Group,
  log = 'y',
  las = 2,
  cex.axis = 0.5
)


meta$richness_cellular <- apply(as.matrix(counts)[X_f[,'cellular organisms'] > 0,], 2, function(x) sum(x > 0))
meta$richness_cellular_norm <- meta$richness_cellular / apply(as.matrix(counts), 2, sum)

beanplot::beanplot(
  meta$richness_cellular_norm ~ meta$Group,
  log = 'y',
  las = 2,
  cex.axis = 0.5
)


meta$richness_bacteria <- apply(as.matrix(counts)[X_f[,'Bacteria'] > 0,], 2, function(x) sum(x > 0))
meta$richness_bacteria_norm <- meta$richness_bacteria / apply(as.matrix(counts),2,sum)

beanplot::beanplot(
  meta$richness_bacteria_norm ~ meta$Gender,
  log = 'y',
  las = 2,
  cex.axis = 0.5
)

meta$richness_Sfi11virus <- apply(as.matrix(counts)[X_f[,'Sfi11virus'] > 0,], 2, function(x) sum(x > 0))
meta$richness_Sfi11virus_norm <- meta$richness_Sfi11virus / apply(as.matrix(counts),2,sum)

beanplot::beanplot(
  meta$richness_Sfi11virus ~ meta$Group,
  las = 2,
  cex.axis = 0.5
)


stan.fit.draws <- stan.fit$post_warmup_draws

time <- stan.fit.draws[,,grep('^time\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
time <- apply(time,3,median)

finalMicrobeTree_time <- finalMicrobeTree
finalMicrobeTree_time$edge.length <- rev(time[sa[,2]])

var_prop_prevalence <- stan.fit.draws[,,grep('^var_prop_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
var_prop_prevalence <- apply(var_prop_prevalence,3,median)

var_prop_abundance <- stan.fit.draws[,,grep('^var_prop_abundance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
var_prop_abundance <- apply(var_prop_abundance,3,median)

sigma_prevalence <- stan.fit.draws[,,grep('^sigma_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_prevalence <- apply(sigma_prevalence,3,median)

sigma_abundance <- stan.fit.draws[,,grep('^sigma_abundance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_abundance <- apply(sigma_abundance,3,median)

sd_prevalence <- stan.fit.draws[,,grep('^sd_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sd_prevalence <- apply(sd_prevalence,3,median)

sd_abundance <- stan.fit.draws[,,grep('^sd_abundance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sd_abundance <- apply(sd_abundance,3,median)

theta_prevalence <- stan.fit.draws[,,grep('^theta_prevalence',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
theta_prevalence <- apply(theta_prevalence,3,median)

theta_abundance <- stan.fit.draws[,,grep('^theta_abundance',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
theta_abundance <- apply(theta_abundance,3,median)

inv_log_less_contamination <- stan.fit.draws[,,grep('^inv_log_less_contamination',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
inv_log_less_contamination <- apply(inv_log_less_contamination,3,median)

contaminant_overdisp <- stan.fit.draws[,,grep('^contaminant_overdisp',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
contaminant_overdisp <- apply(contaminant_overdisp,3,median)


delta_prevalence <- stan.fit.draws[,,grep('^delta_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
