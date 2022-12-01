data {
  int N_samples;
  int N_effects;
  matrix[N_samples,N_effects] model_matrix;
  int N_tips;
  matrix[N_tips,N_tips] phy_cov;
  array[N_samples] int idx_tips;
  vector[N_samples] y;
}
transformed data {
  matrix[N_tips,N_tips] L = cholesky_decompose(phy_cov);
}
parameters {
  vector[N_effects] beta;
  real<lower=0> sigma_phy;
  vector[N_tips] beta_phy;
  real<lower=0> sigma;
}
model {
  beta_phy ~ multi_normal_cholesky(rep_vector(0,N_tips), sigma_phy * L);
  y ~ normal(model_matrix * beta + beta_phy[idx_tips], sigma);
}
