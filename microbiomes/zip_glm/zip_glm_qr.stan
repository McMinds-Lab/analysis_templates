functions {
    int num_zeros(array[] int y) {
      int nz = 0;
      for (n in 1:size(y))
          nz += (y[n] == 0);
      return nz;
    }
    real zip_lpmf(array[] int c, vector p, vector a, array[] int i0, array[] int in0) {
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
    matrix[NB_s-1,NB_s-1] X_sR_inv;         // inverse(qr_thin_R(X_s[,2:])), precomputed in R to use pivoting and be more stable
    int NF;                                 // number of features (genes, taxa, etc)
    int NB_f;                               // number of feature factor levels
    int NFB;                                // number of feature factors
    array[NB_f] int idx_f;                  // mapping of sigmas to feature factor levels including residuals
    matrix[NB_f,NF] X_f;                    // model matrix for features (e.g. taxonomy, functional annotations, gene length, etc.). must include intercept but not residuals
    matrix[NB_f-1,NB_f-1] X_fR_inv;         // inverse(qr_thin_R(X_f[2:,'))', precomputed in R to use pivoting and be more stable
    array[NF,NS] int count;                 // observations
    real prior_scale_p;
    real prior_scale_a;
    int K_s;
}
transformed data {
    array[NB_f-1] int idx_f2;
    array[NB_s+K_s] int idxk_s;
    array[NF*NS] int count_1d = to_array_1d(count);
    int N_zero = num_zeros(count_1d);
    array[N_zero] int i0;
    array[NF*NS-N_zero] int in0;
    int i0i = 1;
    int in0i = 1;
    for(f in 2:NB_f) idx_f2[f-1] = idx_f[f] - 1;
    idxk_s[1:NB_s] = idx_s;
    for(s in 1:K_s) idxk_s[NB_s+s] = NSB + 1;
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
    vector<lower=0>[(NSB+2)*(NFB+1)-2] sd_prevalence_norm;
    matrix<lower=0>[(NSB+2),NFB]       sd_abundance;
    vector<lower=0>[NS]                sd_resid_s;
    row_vector<lower=0>[NF]            sd_resid_f;
    matrix[NB_s+K_s,NF]                beta_prevalence_s;
    matrix[NB_s+K_s,NB_f]              beta_prevalence_i;
    matrix[NS,NB_f]                    beta_prevalence_f;
    matrix[NB_s+K_s,NF]                beta_abundance_s_tilde;
    matrix[NB_s+K_s,NB_f-1]            beta_abundance_i_tilde;
    matrix[NS,NB_f-1]                  beta_abundance_f_tilde;
    matrix[NS,NF]                      abundance;
    vector[NS]                         multinomial_nuisance;
    cholesky_factor_cov[NS,K_s]        L_s;
}
transformed parameters {
    matrix<lower=0>[NSB+2,NFB+1] sd_prevalence = to_matrix(append_row(100, append_row(sd_prevalence_norm, 1.0)), NSB+2, NFB+1) * (prior_scale_p * global_scale_prevalence);
    matrix[NB_s+K_s,NF]          beta_abundance_s;
    matrix[NB_s+K_s,NB_f-1]      beta_abundance_i;
    matrix[NS,NB_f-1]            beta_abundance_f;
    beta_abundance_s[1,]          =            beta_abundance_s_tilde[1,];
    beta_abundance_s[(NB_s+1):,]  =            beta_abundance_s_tilde[(NB_s+1):,];
    beta_abundance_s[2:NB_s,]     = X_sR_inv * beta_abundance_s_tilde[2:NB_s,];
    beta_abundance_i[1,]          =            beta_abundance_i_tilde[1,];
    beta_abundance_i[(NB_s+1):,]  =            beta_abundance_i_tilde[(NB_s+1):,];
    beta_abundance_i[2:NB_s,]     = X_sR_inv * beta_abundance_i_tilde[2:NB_s,] * X_fR_inv;
    beta_abundance_f              =            beta_abundance_f_tilde          * X_fR_inv;
}
model {
    matrix[NS,NB_s+K_s] XL_s = append_col(X_s,L_s);
    matrix[NS,NF] prevalence
        =    XL_s * (beta_prevalence_s .* sd_prevalence[idxk_s, rep_array(NFB+1,NF)])
          + (XL_s * (beta_prevalence_i .* sd_prevalence[idxk_s, idx_f])
             +      (beta_prevalence_f .* sd_prevalence[rep_array(NSB+2,NS), idx_f])) * X_f;
    matrix[NS,NF] abundance_predicted
        =    XL_s * beta_abundance_s
          + (XL_s * beta_abundance_i
             +      beta_abundance_f) * X_f[2:,]
          + rep_matrix(multinomial_nuisance, NF);
    abundance_predicted = diag_pre_multiply(sd_resid_s, diag_post_multiply(abundance_predicted, sd_resid_f));
    // priors
    target += std_normal_lpdf(global_scale_prevalence);
    target += std_normal_lpdf(sd_prevalence_norm);
    target += std_normal_lpdf(sd_resid_s);
    target += std_normal_lpdf(sd_resid_f);
    target += std_normal_lpdf(to_vector(beta_prevalence_s));
    target += std_normal_lpdf(to_vector(beta_prevalence_i));
    target += std_normal_lpdf(to_vector(beta_prevalence_f));
    target += normal_lpdf(sd_abundance[NSB+2,NFB] | 0, prior_scale_a);
    target += normal_lpdf(to_vector(sd_abundance)[:((NSB+2)*NFB-1)] | 0, sd_abundance[NSB+2,NFB]);
    target += normal_lpdf(to_vector(beta_abundance_s) | 0, to_vector(sd_abundance[idxk_s, rep_array(NFB,NF)]));
    target += normal_lpdf(to_vector(beta_abundance_i) | 0, to_vector(sd_abundance[idxk_s, idx_f2]));
    target += normal_lpdf(to_vector(beta_abundance_f) | 0, to_vector(sd_abundance[rep_array(NSB+2,NS), idx_f2]));
    target += normal_lpdf(to_vector(abundance) | to_vector(abundance_predicted), sd_abundance[NSB+2, NFB] * to_vector(sd_resid_s * sd_resid_f));
    for(k in 1:K_s) target += std_normal_lpdf(L_s[k:,k]);
    // likelihood
    target += zip_lpmf(count_1d | to_vector(prevalence), to_vector(abundance), i0, in0);
}
