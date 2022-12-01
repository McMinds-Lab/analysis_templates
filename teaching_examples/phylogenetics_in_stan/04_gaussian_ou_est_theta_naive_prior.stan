functions {
  matrix ou_cov(matrix dist, matrix shared, real theta) {
    // dist is the distance between tips i and j (in time)
    // shared is time between the root and the last common ancestor of i and j
    // theta is the OU strength of attraction
    real part1 = inv(2*theta); // precompute factor constant across entire matrix
    matrix[rows(dist),cols(dist)] V;
    for(j in 1:cols(dist)) {
      for(i in j:rows(dist)) {
        V[i,j] = part1 * exp(-theta * dist[i,j]) * (1-exp(-2 * theta * shared[i,j])); 
        V[j,i] = V[i,j];
      }
    }
    return(V);
  } 
  // uses formula (1) from https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/2041-210X.12201
  // rescale() in geiger relies on formula (2) from above (implicitly), as does OU.vcv() in PIGShift, thus neither are appropriate for non-ultrametric trees but both should be equivalent to the above for ultrametric trees
  // PIGShift's formulation seems to result in impossible covariance matrix when provided a non-ultrametric tree, whereas the rescale() method produces a proper covariance that is inappropriate for the model
}
data {
  int N_samples;
  int N_effects;
  matrix[N_samples,N_effects] model_matrix;
  int N_tips;
  matrix[N_tips,N_tips] phy_dist;
  matrix[N_tips,N_tips] phy_shared;
  array[N_samples] int idx_tips;
  vector[N_samples] y;
}
parameters {
  vector[N_effects] beta;
  real<lower=0> sigma_phy;
  vector[N_tips] beta_phy;
  real<lower=0> sigma;
  real<lower=0> theta;
}
transformed parameters {
  matrix[N_tips,N_tips] phy_cov = ou_cov(phy_dist, phy_shared, theta);
}
model {
  theta ~ exponential(1);
  beta_phy ~ multi_normal(rep_vector(0,N_tips), sigma_phy^2 * phy_cov);
  y ~ normal(model_matrix * beta + beta_phy[idx_tips], sigma);
}
