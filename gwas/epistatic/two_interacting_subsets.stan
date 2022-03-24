// like ATOMM
data {
  int NS; // number of samples
  int NU1; // number of unique genetic categories (e.g. species, lineages, combinations of snps) in subset (genome) 1
  int NU2; // number of unique genetic categories (e.g. species, lineages, combinations of snps) in subset (genome) 2
  int NA1; // number of alleles in subset (genome) 1 (could be same as number of snps, or could be coded to distinguish heterozygosity or two columns for if locus has three alleles...)
  int NA2; // number of alleles in subset (genome) 2 (could be same as number of snps, or could be coded to distinguish heterozygosity or two columns for if locus has three alleles...)
  int NV; // number of non-genetic predictors
  matrix[NS,NV] predictor_matrix;
  matrix[NU1,NA1] snp_matrix_1;
  matrix[NU2,NA2] snp_matrix_2;
  array[NS] int idx_1;
  array[NS] int idx_2;
  real sigma_prior;
  vector[NS] phenotype; // assuming continuous and unbounded for now
}
transformed data {
  matrix[NS,NV] predictor_matrix_norm = diag_post_multiply(predictor_matrix, inv_sqrt(columns_dot_self(predictor_matrix))); // this normalization is probably notright if multi-categorical factors exist (because they should only have a single variance component estimated, and the normalization should be on the sum of all columns belonging to the factor)
  matrix[NU1,NU1] cor_1 = tcrossprod(snp_matrix_1);
  matrix[NU2,NU2] cor_2 = tcrossprod(snp_matrix_2);
  matrix[NS,NS] cor_1_2;
  cor_1 = quad_form_diag(cor_1,inv_sqrt(diagonal(cor_1)));
  cor_2 = quad_form_diag(cor_2,inv_sqrt(diagonal(cor_2)));
  cor_1_2 = cor_1[idx_1,idx_1] .* cor_2[idx_2,idx_2];
  cor_1_2 = quad_form_diag(cor_1_2,inv_sqrt(diagonal(cor_1_2)));
}
parameters {
  real<lower=0> sigma;
  simplex[NV+4] variance_components;
  real intercept;
}
model {
  sigma ~ normal(0,sigma_prior);
  phenotype
    ~ multi_normal(intercept,
                   square(sigma)
                   * add_diag(tcrossprod(diag_post_multiply(predictor_matrix_norm,variance_components[1:NV]))
                              + variance_components[NV+1] * cor_1
                              + variance_components[NV+2] * cor_2
                              + variance_components[NV+3] * cor_1_2
                              , variance_components[NV+4] + 1e-9));
}
