data {
  int N_samples;
  int N_effects;
  matrix[N_samples,N_effects] model_matrix;
  int N_tips;
  int N_nodes;
  array[N_nodes-1,2] int edge;
  vector[N_nodes-1] edge_lengths;
  array[N_samples] int idx_tips;
  vector[N_samples] y;
}
parameters {
  vector[N_effects] beta;
  real<lower=0> sigma_phy;
  vector[N_nodes] beta_phy;
  real<lower=0> sigma;
  real<lower=0> theta;
}
model {
  vector[N_nodes-1] ou_props = exp(-theta * edge_lengths); // precalculate values related to OU process 
  
  theta ~ exponential(1);
  
  beta_phy[N_tips+1] ~ normal(0, sigma_phy / sqrt(2*theta)); // root node is always N_tips + 1, and we will sample it directly from the stationary distribution
  beta_phy[edge[,2]] ~ normal(ou_props .* beta_phy[edge[,1]], sqrt(1 - square(ou_props)) * (sigma_phy / sqrt(2*theta))); // all other nodes and tips have simple univariate normal priors defined by the OU process
  
  y ~ normal(model_matrix * beta + beta_phy[idx_tips], sigma);
}
