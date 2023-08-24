functions {
    int num_zeros(array[] int y) {
      int nz = 0;
      for (n in 1:size(y))
          nz += (y[n] == 0);
      return(nz);
    }
    real zip_lpmf(array[] int c, vector p, vector a, array[] int i0, array[] int in0) {
      vector[size(p)] lilp = log_inv_logit(p);
      real lp = sum(lilp[in0]) + poisson_log_lpmf(c[in0] | a[in0]);
      for(i in i0) {
          lp += log_sum_exp(log1m_inv_logit(p[i]), lilp[i] - exp(a[i])); // poisson lpmf with count 0 = log(exp(-exp(log_lambda)))
      }
      return(lp);
    }
    real calc_antagonism(real a, real b, real c, real x) {
      // a adjusts output scale
      // b adjusts scale of input distances
      // c influences shape of trough
      // x is phylogenetic distance
      return(antag = a * log(b * c * x) * exp(-b / c * x));
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
    real prior_scale_p;                     // expectation for scale of prevalence effects
    real prior_scale_a;                     // expectation for scale of abundance effects
    int K_s;                                // number of latent factors for feature covariance
    matrix[NF,NF] dist;                     // phylogenetic distance between features
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
    vector<lower=0>[(NSB+2)*NFB-1]     sd_abundance_norm;
    matrix[NB_s+K_s,NB_f]              beta_prevalence_i;
    matrix[NB_s+K_s,NF]                beta_prevalence_s;
    matrix[NS,NB_f]                    beta_prevalence_f;
    matrix[NB_s+K_s,NB_f-1]            beta_abundance_i;
    matrix[NB_s+K_s,NF]                beta_abundance_s;
    matrix[NS,NB_f-1]                  beta_abundance_f;
    matrix[NS,NF]                      residuals;               // 'overdispersion' of error for abundance
    vector[NS]                         multinomial_nuisance;    // converts poisson error into multinomial error
    cholesky_factor_cov[NS,K_s]        L_s;                     // sample low rank residual covariance
    row_vector[NB_f+NF]                beta_a_antagonism;       // feature phylogenetic antagonism parameter a
    row_vector[NB_f+NF]                beta_b_antagonism;       // feature phylogenetic antagonism parameter b
    row_vector[NB_f+NF]                beta_c_antagonism;       // feature phylogenetic antagonism parameter c
}
transformed parameters {
    matrix<lower=0>[NSB+2,NFB+1] sd_prevalence = to_matrix(append_row(100, append_row(sd_prevalence_norm, 1.0)), NSB+2, NFB+1) * (prior_scale_p * global_scale_prevalence);
    matrix<lower=0>[NSB+2,NFB]   sd_abundance  = to_matrix(append_row(sd_abundance_norm, 1.0),                   NSB+2, NFB)   * (prior_scale_a * global_scale_abundance);
}
model {
    matrix[NS,NB_s+K_s] XL_s = append_col(X_s,L_s);
    matrix[NS,NF] prevalence
      =    XL_s * (beta_prevalence_s .* sd_prevalence[idxk_s, rep_array(NFB+1,NF)])
        + (XL_s * (beta_prevalence_i .* sd_prevalence[idxk_s, idx_f])
           +      (beta_prevalence_f .* sd_prevalence[rep_array(NSB+2,NS), idx_f])) * X_f;
    matrix[NS,NF] abundance
      =    XL_s * (beta_abundance_s  .* sd_abundance[idxk_s, rep_array(NFB,NF)])
        + (XL_s * (beta_abundance_i  .* sd_abundance[idxk_s, idx_f2])
           +      (beta_abundance_f  .* sd_abundance[rep_array(NSB+2,NS), idx_f2])) * X_f[2:,]
        + residuals                  .* sd_abundance[NSB+2, NFB]
        + rep_matrix(multinomial_nuisance, NF);
    row_vector[NF] a_antagonism = beta_a_antagonism[1:NB_f] * X_f + beta_a_antagonism[(NB_f+1):];
    row_vector[NF] b_antagonism = beta_b_antagonism[1:NB_f] * X_f + beta_b_antagonism[(NB_f+1):];
    row_vector[NF] c_antagonism = beta_c_antagonism[1:NB_f] * X_f + beta_c_antagonism[(NB_f+1):];
    matrix[NS,NF] antagonism = rep_matrix(0,NS,NF);
    for(f in 1:NF) {
      for(s in 1:NS) {
        for(f2 in 1:(NF-1)) {
          // arrange indices to skip comparison to self
          if(f2 < f) {
            i = f2;
          } else {
            i = f2 + 1;
          }
          // if a comparison taxon is present, add its antagonism to the current taxon's probability (suddenly i realize that this won't work because a neighbor's presence will be completely sufficient to explain its presence and vice versa, so rest of model is essentially ignored. i really do need to work with covariance, not direct probability. or else use this format with latent value determined by mixture of probabilities based on all possibilities of the order of arrival of all microbes)
          // spectral_kernal <- function(nu,mu,dist) {exp(-2*pi^2*dist^2*nu)*cos(2*dist*pi*mu)}
          // 0.05,0.5,cophenetic(phy) is reasonable
          // 1,3,cophenetic(phy) is reasonable
          // 20,10,cophenetic(phy) is reasonable
          // should reparameterize so nu and mu share a scale, but unfortunately not linear. this way a prior could be set to keep a single wiggle at beginning, but just vary where that single wiggle happens
          if(count[i,s] > 0) {
            antagonism[s,f] += calc_antagonism(a_antagonism[f], b_antagonism[f], c_antagonism[f], dist[,f]);
          }
        }
      }
    }
    array[NS] cholesky_factor_cov[NF,NF] chol;
    for(s in NS) {
      for(f1 in NF) {
        for(f2 in 1:(f1-1)) {
          if(both present) {
            // antagonism of one added to other based on mixture probability of which arrived first
            // (presence of one might be entirely explained by presence of other, but one has to be explained by the rest of the model)
          } else if(one present) {
            // antagonism of present added to that of absent
          }
        }
      }
      
    }
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
    // likelihood
    target += zip_lpmf(count_1d | to_vector(prevalence), to_vector(abundance), i0, in0);
}
