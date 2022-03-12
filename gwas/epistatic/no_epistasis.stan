data {
  int NS; // number of samples
  int NU; // number of unique genetic lineages
  int NA; // number of alleles (could be same as number of snps, or could be coded to distinguish heterozygosity or two columns for if locus has three alleles...)
  int NV; // number of non-genetic predictors
  matrix[NS,NV] predictor_matrix; // probably normalized by sd per col
  matrix[NU,NA] snp_matrix;
  array[NS] int idx; // mapping of samples to lineages
  real sigma_prior;
  vector[NS] phenotype; // assuming continuous and unbounded for now
}
transformed data {
  matrix[NS,NV] predictor_matrix_norm = diag_post_multiply(predictor_matrix, inv_sqrt(columns_dot_self(predictor_matrix)));
  matrix[NU,NU] cor = tcrossprod(snp_matrix);
  cor = quad_form_diag(cor,inv_sqrt(diagonal(cor)));
}
parameters {
  real<lower=0> sigma;
  simplex[NV+2] variance_components;
  real intercept;
}
model {
  sigma ~ normal(0,sigma_prior);
  phenotype
    ~ multi_normal(intercept,
                   square(sigma)
                   * add_diag(tcrossprod(diag_post_multiply(predictor_matrix_norm,variance_components[1:NV]))
                              + variance_components[NV+1] * cor[idx,idx],
                              variance_components[NV+2] + 1e-9));
}
