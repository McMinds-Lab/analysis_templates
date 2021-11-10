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
    matrix[NB_s,NB_s] X_s_full_inv;              // completed model matrix for samples (e.g. tissue compartments, duplicates, sequencing depth, etc.). must include intercept and residuals
    int NF;                                 // number of features (genes, taxa, etc)
    int NB_f;                               // number of feature factor levels
    int NFB;                                // number of feature factors
    array[NB_f] int idx_f;                  // mapping of sigmas to feature factor levels including residuals
    matrix[NB_f,NB_f] X_f_full_inv;              // completed model matrix for features (e.g. taxonomy, functional annotations, gene length, etc.). must include intercept and residuals
    array[NF,NS] int count;                 // observations
    real prior_scale_p;
    real prior_scale_a;
}
transformed data {
    array[NB_f-1] int idx_f2;
    for(f in 2:NB_f) idx_f2[f-1] = idx_f[f] - 1;
}
parameters {
    real<lower=0> global_scale_prevalence;
    real<lower=0> global_scale_abundance;
    vector<lower=0>[NSB*NFB-2] sd_prevalence_norm;
    vector<lower=0>[NSB*(NFB-1)-1] sd_abundance_norm;
    matrix[NB_s,NB_f] prevalence;
    matrix[NB_s,NB_f-1] abundance;
    vector[NS] multinomial_nuisance;
}
transformed parameters {
    matrix<lower=0>[NSB,NFB] sd_prevalence = to_matrix(append_row(100 * global_scale_prevalence, append_row(sd_prevalence_norm, global_scale_prevalence)), NSB, NFB) * prior_scale_p;
    matrix<lower=0>[NSB,NFB-1] sd_abundance = to_matrix(append_row(sd_abundance_norm, global_scale_abundance), NSB, NFB-1) * prior_scale_a;
    matrix[NB_s,NB_f] beta_prevalence = X_s_full_inv * prevalence * X_f_full_inv;
    matrix[NB_s,NB_f-1] beta_abundance = abundance;
    beta_abundance[(NB_s-NS+1):,(NB_f-NF):] = beta_abundance[(NB_s-NS+1):,(NB_f-NF):] - rep_matrix(multinomial_nuisance, NF);
    beta_abundance = X_s_full_inv * beta_abundance * X_f_full_inv[2:,2:];
}
model {
    // priors
    target += std_normal_lpdf(global_scale_prevalence);
    target += std_normal_lpdf(global_scale_abundance);
    target += student_t_lpdf(sd_prevalence_norm | 5, 0, global_scale_prevalence);
    target += student_t_lpdf(sd_abundance_norm | 5, 0, global_scale_abundance);
    target += student_t_lpdf(to_vector(beta_prevalence) | 5, 0, to_vector(sd_prevalence[idx_s,idx_f]));
    target += student_t_lpdf(to_vector(beta_abundance) | 5, 0, to_vector(sd_abundance[idx_s,idx_f2]));
    // likelihood
    for(f in 1:NF) {
        int f2 = NB_f - NF + f;
        for(s in 1:NS) {
            int s2 = NB_s - NS + s;
            target += zip_lpmf(count[f,s] | prevalence[s2,f2], abundance[s2,f2-1]);
        }
    }
}
generated quantities {
    matrix[NB_s,NB_f] beta_prevalence_norm = beta_prevalence ./ sd_prevalence[idx_s,idx_f];
    matrix[NB_s,NB_f-1] beta_abundance_norm = beta_abundance ./ sd_abundance[idx_s,idx_f2];
}
