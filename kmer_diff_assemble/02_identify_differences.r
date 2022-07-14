## in doMPI batches:
## zinbwave to get latent vars and sample intercepts ('size factors'). (maybe only do this on a subset of batches, and somehow summarize the results so they're fixed for all batches of the below tests. mean size factor; procrustes average of latent vars?)
## pscl::zeroinfl with model "~ 1 + offset(size_factor) + study + host_species + latent_vars + is_diseased + has_sctld" to get prevalence significance. still not sure how best to use 'disease exposed' or 'apparently healthy' type samples
