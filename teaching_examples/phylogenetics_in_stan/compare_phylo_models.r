# https://rdrr.io/cran/PIGShift/src/R/brownian_motion_package.r#sym-OU.vcv
# https://rdrr.io/cran/geiger/src/R/utilities-phylo.R#sym-.ou.phylo
# https://discourse.mc-stan.org/t/rescaling-phylogenetic-tree-for-non-brownian-phylogenetic-models-in-brms/29637
# https://bookdown.org/content/4857/adventures-in-covariance.html#example-phylogenetic-distance

## code some simple examples of phylogenetic linear models in stan. progress from brownian model to OU model with fixed alpha to OU model with estimated alpha, then explore different parameterizations (simple vcv formula a la above code vs multilevel model with transformed branch lengths a la above code vs explicit evolutionary model with ancestral states.)

  
## define tree size, number of samples per tip, and mapping between samples and tree
N_tips <- 200
N_samples_per_tip <- 5
N_samples <- N_tips * N_samples_per_tip
idx_tips <- rep(1:N_tips, each = N_samples_per_tip)

## generate random tree
mytree <- ape::chronos(ape::rtree(N_tips))

## generate sample traits for non-phylogenetic (unstructured) independent variables (assume standardized)
N_iv_us <- 2
ivs_us <- sapply(1:N_iv_us, \(x) rnorm(N_samples))

## generate sample traits for phylogenetic independent variables (assume standardized)
N_iv_phylo <- 2
theta_iv <- rexp(2)
ivs_phylo <- sapply(theta_iv, \(x) castor::simulate_ou_model(mytree, 0, 1 / sqrt(2*x), x, include_nodes=FALSE)$tip_states) ## stationary sd = sigma / sqrt(2*theta), so this assumes sigma = 1 (i.e. standardized by sigma)

## combine model matrix
model_matrix <- cbind(ivs_us, ivs_phylo[idx_tips,])
colnames(model_matrix) <- paste0('iv_',1:(N_iv_us+N_iv_phylo))

## define sample trait betas (one of each category predefined to have no effect)
betas <- c(0,rnorm(N_iv_us-1),0,rnorm(N_iv_phylo-1))

## define strength of attraction for dependent variable, sigma of confounding phylogenetic effects, residual variance, and intercept
theta_dv <- 3
sigma_phy <- 2
sigma_resid <- 0.5
intercept <- rnorm(1,0,10)

## generate random data according to OU evolution with random intercept, and given tree, independent variables, and theta
y <- intercept + model_matrix %*% betas + castor::simulate_ou_model(mytree, 0, sigma_phy / sqrt(2*theta_dv), theta_dv, include_nodes=FALSE)$tip_states[idx_tips] + rnorm(N_samples,0,sigma_resid)

## combine data and model matrix
dat <- data.frame(y=y, tip=as.factor(mytree$tip.label[idx_tips]), model_matrix)

## analyze using phyr with geiger-generated cov
cov_geiger <- ape::vcv.phylo(geiger::rescale(mytree, 'OU', alpha = theta_dv, sigsq = 1))
res_phyr_geiger <- phyr::pglmm(y ~ iv_1 + iv_2 + iv_3 + iv_4, random.effects = list(re1 = list(cov_geiger[idx_tips,idx_tips])), data=dat)
## expect Std.Dev of re1 to be sigma_phy
summary(res_phyr_geiger)

## analyze using phyr with PIGShift-generated cov
cov_pigshift <- PIGShift::OU.vcv(mytree,theta_dv) 
res_phyr_pigshift <- phyr::pglmm(y ~ iv_1 + iv_2 + iv_3 + iv_4, random.effects = list(re1 = list(cov_pigshift[idx_tips,idx_tips])), data=dat)
summary(res_phyr_pigshift)

# cov_geiger and cov_pigshift are the same thing
all(round(cov_geiger,digits=5) == round(cov_pigshift,digits=5))

## analyze using geiger cov and stan model 1
## assumes working directory is 'phylogenetics_in_stan'
model_stan1_geiger <- cmdstanr::cmdstan_model('01_gaussian_ou_fixed_precalc.stan')
res_stan1_geiger <- model_stan1_geiger$sample(data = list(N_samples=N_samples, N_effects=ncol(model_matrix)+1, model_matrix=cbind(1,model_matrix), N_tips=N_tips, phy_cov=cov_geiger, idx_tips=idx_tips, y=as.vector(y)), parallel_chains=4)
res_stan1_geiger$summary(variables=c('beta','sigma','sigma_phy'))

## analyze using geiger cov and stan model 2 (should be same model, just calc the covariance within Stan)
model_stan2_geiger <- cmdstanr::cmdstan_model('02_gaussian_ou_fixed_stancalc.stan')
res_stan2_geiger <- model_stan2_geiger$sample(data = list(N_samples=N_samples, N_effects=ncol(model_matrix)+1, model_matrix=cbind(1,model_matrix), N_tips=N_tips, phy_dist=ape::cophenetic.phylo(mytree), phy_shared=phytools::findMRCA(mytree,type='height'), theta=theta_dv, idx_tips=idx_tips, y=as.vector(y)), parallel_chains=4)
res_stan2_geiger$summary(variables=c('beta','sigma','sigma_phy'))

## analyze using stan model 3 (using the covariance calculations that we demonstrated in the previous model, but now we don't have to guess theta beforehand!)
model_stan3 <- cmdstanr::cmdstan_model('03_gaussian_ou_est_theta_naive.stan')
res_stan3 <- model_stan3$sample(data = list(N_samples=N_samples, N_effects=ncol(model_matrix)+1, model_matrix=cbind(1,model_matrix), N_tips=N_tips, phy_dist=ape::cophenetic.phylo(mytree), phy_shared=phytools::findMRCA(mytree,type='height'), idx_tips=idx_tips, y=as.vector(y)), parallel_chains=4)
res_stan3$summary(variables=c('beta','sigma','sigma_phy','theta'))

## analyze using stan model 4 (put a prior on theta to avoid wacky estimates and improve efficiency of sampling)
model_stan4 <- cmdstanr::cmdstan_model('04_gaussian_ou_est_theta_naive_prior.stan')
res_stan4 <- model_stan4$sample(data = list(N_samples=N_samples, N_effects=ncol(model_matrix)+1, model_matrix=cbind(1,model_matrix), N_tips=N_tips, phy_dist=ape::cophenetic.phylo(mytree), phy_shared=phytools::findMRCA(mytree,type='height'), idx_tips=idx_tips, y=as.vector(y)), parallel_chains=4)
res_stan4$summary(variables=c('beta','sigma','sigma_phy','theta'))



