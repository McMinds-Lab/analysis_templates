#include functions.stan
data {
    int NS;                                 // number of samples
    int NB_s;                               // number of sample factor levels
    int NSB;                                // number of sample factors
    array[NB_s] int idx_s;                  // mapping of sigmas to sample factor levels
    matrix[NS,NB_s] X_s;                    // model matrix for samples (e.g. tissue compartments, duplicates, sequencing depth, etc.). must include intercept
    int NF;                                 // number of features (genes, taxa, etc)
    int NB_f;                               // number of feature factor levels
    int NFB;                                // number of feature factors
    array[NB_f] int idx_f;                  // mapping of sigmas to feature factor levels
    matrix[NB_f,NF] X_f;                    // model matrix for features (e.g. taxonomy, functional annotations, gene length, etc.). must include intercept
    array[NF,NS] int count;                 // observations
    real inv_log_max_contam;                // prior expectation of contamination rate
    real<lower=0> shape_gnorm;              // strength of prior pulling contamination toward zero
}
parameters {
    real<lower=0> global_scale;
    real<lower=0,upper=1> global_scale_prop;
    simplex[NSB*NFB] vp_prevalence;    // proportion of variance of sample effects
    simplex[NSB*NFB+1] vp_abundance;   // proportion of variance of sample effects
    matrix[NB_s,NB_f] beta_prevalence_raw;
    matrix[NB_s,NB_f] beta_abundance_raw;
    matrix[NS,NF] abundance_observed;
    vector[NS] multinomial_nuisance;
    real<upper=0> inv_log_less_contamination; // smaller = less average contamination
    real<lower=0> contaminant_overdisp;            // dispersion parameter for amount of contamination in true negative count observations
}
transformed parameters {
    real global_scale_prevalence = global_scale * sqrt(1-global_scale_prop);
    real global_scale_abundance = global_scale * sqrt(global_scale_prop);
    matrix[NSB,NFB] sd_prevalence = to_matrix(sqrt(vp_prevalence) * global_scale_prevalence, NSB, NFB);
    matrix[NSB,NFB] sd_abundance = to_matrix(sqrt(vp_abundance[1:(NSB*NFB)]) * global_scale_abundance, NSB, NFB);
    real sd_resid = sqrt(vp_abundance[NSB*NFB+1]) * global_scale_abundance;
    matrix[NB_s,NB_f] beta_prevalence = sd_prevalence[idx_s,idx_f] .* beta_prevalence_raw;
    matrix[NB_s,NB_f] beta_abundance = sd_abundance[idx_s,idx_f] .* beta_abundance_raw;
    real log_less_contamination = inv(inv_log_less_contamination);
}
model {
    // data wrangling
    matrix[NS,NF] prevalence = X_s * beta_prevalence * X_f;
    matrix[NS,NF] abundance_predicted = X_s * beta_abundance * X_f;
    row_vector[NF] abundance_contam
        = beta_abundance[1,] * X_f
          + log_inv_logit(beta_prevalence[1,] * X_f)
          + log_less_contamination;
    // priors
    target += student_t_lpdf(global_scale | 5, 0, 2.5);
    target += student_t_lpdf(to_vector(beta_prevalence_raw) | 5,0,1);
    target += student_t_lpdf(to_vector(beta_abundance_raw) | 5,0,1);
    target += generalized_std_normal_1_lpdf(inv_log_less_contamination / inv_log_max_contam | shape_gnorm);   // shrink amount of contamination in 'true zeros' toward zero
    target += lognormal_lpdf(contaminant_overdisp | 0, 0.1);                                               // shrink overdispersion of contaminant counts in 'true zeros' toward zero
    // likelihood
    for(f in 1:NF) {
        for(s in 1:NS) {
            target += log_sum_exp(log1m_inv_logit(prevalence[s,f])
                                  + student_t_lpdf(abundance_observed[s,f] |
                                                   5,
                                                   abundance_contam[f],
                                                   contaminant_overdisp * sd_resid), //estimated abundance if true negative
                                  log_inv_logit(prevalence[s,f])
                                  + student_t_lpdf(abundance_observed[s,f] |
                                                   5,
                                                   log_sum_exp(abundance_contam[f], abundance_predicted[s,f]),
                                                   sd_resid)); //estimated abundance if true positive
        }
    }
    target += poisson_log_lpmf(to_array_1d(count) |
                               to_vector(abundance_observed + rep_matrix(multinomial_nuisance, NF)));
}
