functions {
    real zip_lpmf(int c, real p, real a) {
        if(c > 0) {
            return(log_inv_logit(p) + poisson_log_lpmf(c | a));
        } else {
            return(log_sum_exp(log1m_inv_logit(p), log_inv_logit(p) - exp(a))); // poisson lpmf with count 0 = log(exp(-exp(log_lambda)))
        }
    }
}
data {
    int NS;                                 // number of samples
    int NB_s;                               // number of sample factor levels
    int NSB;                                // number of sample factors
    array[NB_s] int idx_s;                  // mapping of sigmas to sample factor levels
    int rank_X_s;
    matrix[NS,rank_X_s] Q_s;                // Q of model matrix for samples (e.g. tissue compartments, duplicates, sequencing depth, etc.). must include intercept but not residuals
    matrix[NB_s,rank_X_s] R_inv_s;
    matrix[NB_s,NB_s-rank_X_s] Q_I_s;
    int NF;                                 // number of features (genes, taxa, etc)
    int NB_f;                               // number of feature factor levels
    int NFB;                                // number of feature factors
    array[NB_f] int idx_f;                  // mapping of sigmas to feature factor levels including residuals
    int rank_X_f;
    matrix[rank_X_f,NF] Q_f;                // Q of model matrix for features (e.g. taxonomy, functional annotations, gene length, etc.). must include intercept but not residuals
    matrix[rank_X_f,NB_f] R_inv_f;
    matrix[NB_f-rank_X_f,NB_f] Q_I_f;
    matrix[NB_f-rank_X_f,NB_f-1] Q_I_f2;
    array[NF,NS] int count;                 // observations
    real prior_scale_p;
    real prior_scale_a;
    int K_s;
    int K_f;
}
transformed data {
    array[NB_f-1] int idx_f2;
    array[NB_s+K_s] int idxk_s;
    array[NB_f+K_f] int idxk_f;
    array[NB_f+K_f-1] int idxk_f2;
    for(f in 2:NB_f) idx_f2[f-1] = idx_f[f] - 1;
    idxk_s[1:NB_s] = idx_s;
    for(s in 1:K_s) idxk_s[NB_s+s] = NSB + s;
    idxk_f[1:NB_f] = idx_f;
    for(f in 1:K_f) idxk_f[NB_f+f] = NFB + f;
    idxk_f2[1:(NB_f-1)] = idx_f2;
    for(f in 1:K_f) idxk_f2[NB_f-1+f] = NFB - 1 + f;
}
parameters {
    real<lower=0> global_scale_prevalence;
    real<lower=0> global_scale_abundance;
    vector<lower=0>[(NSB+K_s)*(NFB+K_f)-2] sd_prevalence_norm;
    vector<lower=0>[(NSB+K_s)*(NFB+K_f-1)-1] sd_abundance_norm;
    matrix[rank_X_s+K_s,rank_X_f+K_f] beta_prevalence_tilde;
    matrix[NB_s-rank_X_s,NB_f-rank_X_f] omega_prevalence_i;
    matrix[NB_s-rank_X_s,K_f] omega_prevalence_s;
    matrix[K_s,NB_f-rank_X_f] omega_prevalence_f;
    matrix[rank_X_s+K_s,rank_X_f-1+K_f] beta_abundance_tilde;
    matrix[NB_s-rank_X_s,NB_f-rank_X_f] omega_abundance_i;
    matrix[NB_s-rank_X_s,K_f] omega_abundance_s;
    matrix[K_s,NB_f-rank_X_f] omega_abundance_f;
    cholesky_factor_cov[NS,K_s] L_s;
    cholesky_factor_cov[NF,K_f] L_f;
    vector[NS] multinomial_nuisance;
}
transformed parameters {
    matrix<lower=0>[NSB+K_s,NFB+K_f] sd_prevalence = to_matrix(append_row(100 * global_scale_prevalence, append_row(sd_prevalence_norm, global_scale_prevalence)), NSB+K_s, NFB+K_f) * prior_scale_p;
    matrix<lower=0>[NSB+K_s,NFB+K_f-1] sd_abundance = to_matrix(append_row(sd_abundance_norm, global_scale_abundance), NSB+K_s, NFB+K_f-1) * prior_scale_a;
    matrix[NB_s+K_s,NB_f+K_f] beta_prevalence;
    matrix[NB_s+K_s,NB_f+K_f-1] beta_abundance;
    beta_prevalence[1:NB_s,1:NB_f]
        = R_inv_s * beta_prevalence_tilde[1:rank_X_s,1:rank_X_f] * R_inv_f
          + Q_I_s * omega_prevalence_i * Q_I_f;
    beta_prevalence[1:NB_s,(NB_f+1):]
        = R_inv_s * beta_prevalence_tilde[1:rank_X_s,(rank_X_f+1):]
          + Q_I_s * omega_prevalence_s;
    beta_prevalence[(NB_s+1):,1:NB_f]
        = beta_prevalence_tilde[(rank_X_s+1):,1:rank_X_f] * R_inv_f
          + omega_prevalence_f * Q_I_f;
    beta_prevalence[(NB_s+1):,(NB_f+1):]
        = beta_prevalence_tilde[(rank_X_s+1):,(rank_X_f+1):];

    beta_abundance[1:NB_s,1:(NB_f-1)]
        = R_inv_s * beta_abundance_tilde[1:rank_X_s,1:(rank_X_f-1)] * R_inv_f[2:,2:]
          + Q_I_s * omega_abundance_i * Q_I_f2;
    beta_abundance[1:NB_s,NB_f:]
        = R_inv_s * beta_abundance_tilde[1:rank_X_s,rank_X_f:]
          + Q_I_s * omega_abundance_s;
    beta_abundance[(NB_s+1):,1:(NB_f-1)]
        = beta_abundance_tilde[(rank_X_s+1):,1:(rank_X_f-1)] * R_inv_f[2:,2:]
          + omega_abundance_f * Q_I_f2;
    beta_abundance[(NB_s+1):,NB_f:]
        = beta_abundance_tilde[(rank_X_s+1):,rank_X_f:];
}
model {
    matrix[NS,rank_X_s+K_s] QL_s = append_col(Q_s,L_s);
    matrix[rank_X_f+K_f,NF] QL_f = append_row(Q_f,L_f');
    matrix[NS,NF] prevalence
        = QL_s * beta_prevalence_tilde * QL_f;
    matrix[NS,NF] abundance
        = QL_s * beta_abundance_tilde * QL_f[2:,]
          + rep_matrix(multinomial_nuisance, NF);
    // priors
    target += std_normal_lpdf(global_scale_prevalence);
    target += std_normal_lpdf(global_scale_abundance);
    target += normal_lpdf(sd_prevalence_norm | 0, global_scale_prevalence);
    target += normal_lpdf(sd_abundance_norm | 0, global_scale_abundance);
    target += normal_lpdf(to_vector(beta_prevalence) | 0, to_vector(sd_prevalence[idxk_s,idxk_f]));
    target += normal_lpdf(to_vector(beta_abundance) | 0, to_vector(sd_abundance[idxk_s,idxk_f2]));
    for(k in 1:K_s) target += std_normal_lpdf(L_s[k:,k]);
    for(k in 1:K_f) target += std_normal_lpdf(L_f[k:,k]);
    // likelihood
    for(f in 1:NF) {
        for(s in 1:NS) {
            target += zip_lpmf(count[f,s] | prevalence[s,f], abundance[s,f]);
        }
    }
}
