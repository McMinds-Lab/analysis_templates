functions {
    int num_zeros(int[] y) {
      int nz = 0;
      for (n in 1:size(y))
          nz += (y[n] == 0);
      return nz;
    }
    real zip_lpmf(int[] c, vector p, vector a, int[] i0, int[] in0) {
        vector[size(p)] lilp = log_inv_logit(p);
        real lp = sum(lilp[in0]) + poisson_log_lpmf(c[in0] | a[in0]);
        for(i in i0) {
            lp += log_sum_exp(log1m_inv_logit(p[i]), lilp[i] - exp(a[i])); // poisson lpmf with count 0 = log(exp(-exp(log_lambda)))
        }
        return(lp);
    }
}
data {
    int NS;                                 // number of samples
    int NB_s;                               // number of sample factor levels
    int NSB;                                // number of sample factors
    array[NB_s] int idx_s;                  // mapping of sigmas to sample factor levels
    matrix[NS,NB_s] X_s;                    // model matrix for samples (e.g. tissue compartments, duplicates, sequencing depth, etc.). must include intercept but not residuals
    int NF;                                 // number of features (genes, taxa, etc)
    int NB_f;                               // number of feature factor levels
    int NFB;                                // number of feature factors
    array[NB_f] int idx_f;                  // mapping of sigmas to feature factor levels including residuals
    matrix[NB_f,NF] X_f;                    // model matrix for features (e.g. taxonomy, functional annotations, gene length, etc.). must include intercept but not residuals
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
    array[NF*NS] int count_1d = to_array_1d(count);
    int N_zero = num_zeros(count_1d);
    array[N_zero] int i0;
    array[NF*NS-N_zero] int in0;
    int i0i = 1;
    int in0i = 1;
    for(f in 2:NB_f) idx_f2[f-1] = idx_f[f] - 1;
    idxk_s[1:NB_s] = idx_s;
    for(s in 1:K_s) idxk_s[NB_s+s] = NSB + 1;
    idxk_f[1:NB_f] = idx_f;
    for(f in 1:K_f) idxk_f[NB_f+f] = NFB + 1;
    idxk_f2[1:(NB_f-1)] = idx_f2;
    for(f in 1:K_f) idxk_f2[NB_f+f-1] = NFB;
    for(c in 1:size(count_1d)) {
        if(count_1d[c] == 0) {
            i0[i0i] = c;
            i0i += 1;
        } else {
            in0[in0i] = c;
            in0i += 1;
        }
    }
}
parameters {
    real<lower=0>                      global_scale_prevalence;
    real<lower=0>                      global_scale_abundance;
    vector<lower=0>[(NSB+2)*(NFB+2)-2] sd_prevalence_norm;
    vector<lower=0>[(NSB+2)*(NFB+1)-1] sd_abundance_norm;
    matrix[NB_s+K_s,NB_f+K_f]          beta_prevalence_i;
    matrix[NB_s+K_s,NF]                beta_prevalence_s;
    matrix[NS,NB_f+K_f]                beta_prevalence_f;
    matrix[NB_s+K_s,NB_f+K_f-1]        beta_abundance_i;
    matrix[NB_s+K_s,NF]                beta_abundance_s;
    matrix[NS,NB_f+K_f-1]              beta_abundance_f;
    matrix[NS,NF]                      residuals;
    vector[NS]                         multinomial_nuisance;
    cholesky_factor_cov[NS,K_s]        L_s;
    cholesky_factor_cov[NF,K_f]        L_f;
}
transformed parameters {
    matrix<lower=0>[NSB+2,NFB+2] sd_prevalence = to_matrix(append_row(100, append_row(sd_prevalence_norm, 1.0)), NSB+2, NFB+2) * (prior_scale_p * global_scale_prevalence);
    matrix<lower=0>[NSB+2,NFB+1] sd_abundance  = to_matrix(append_row(sd_abundance_norm, 1.0),                   NSB+2, NFB+1) * (prior_scale_a * global_scale_abundance);
}
model {
    matrix[NS,NB_s+K_s] XL_s = append_col(X_s,L_s);
    matrix[NB_f+K_f,NF] XL_f = append_row(X_f,L_f');
    matrix[NS,NF] prevalence
        =    XL_s * (beta_prevalence_s .* sd_prevalence[idxk_s, rep_array(NFB+2,NF)])
          + (XL_s * (beta_prevalence_i .* sd_prevalence[idxk_s, idxk_f])
             +      (beta_prevalence_f .* sd_prevalence[rep_array(NSB+2,NS), idxk_f])) * XL_f;
    matrix[NS,NF] abundance
        =    XL_s * (beta_abundance_s  .* sd_abundance[idxk_s, rep_array(NFB+1,NF)])
          + (XL_s * (beta_abundance_i  .* sd_abundance[idxk_s, idxk_f2])
             +      (beta_abundance_f  .* sd_abundance[rep_array(NSB+2,NS), idxk_f2])) * XL_f[2:,]
          + residuals                   * sd_abundance[NSB+2, NFB+1]
          + rep_matrix(multinomial_nuisance, NF);
    // priors
    target += std_normal_lpdf(global_scale_prevalence);
    target += std_normal_lpdf(global_scale_abundance);
    target += std_normal_lpdf(sd_prevalence_norm);
    target += std_normal_lpdf(sd_abundance_norm);
    target += std_normal_lpdf(to_vector(beta_prevalence_i));
    target += std_normal_lpdf(to_vector(beta_prevalence_s));
    target += std_normal_lpdf(to_vector(beta_prevalence_f));
    target += std_normal_lpdf(to_vector(beta_abundance_i));
    target += std_normal_lpdf(to_vector(beta_abundance_s));
    target += std_normal_lpdf(to_vector(beta_abundance_f));
    target += std_normal_lpdf(to_vector(residuals));
    for(k in 1:K_s) target += std_normal_lpdf(L_s[k:,k]);
    for(k in 1:K_f) target += std_normal_lpdf(L_f[k:,k]);
    // likelihood
    target += zip_lpmf(count_1d | to_vector(prevalence), to_vector(abundance), i0, in0);
}
