functions {
  real ancestral_state_ou_lp(vector s2, vector mu, vector ou_props, real diff_raw) {
    vector[2] s2m = s2 ./ square(ou_props);
    vector[2] mum = mu ./ ou_props;
    target += normal_lpdf(mu[2:] | mu[1], sqrt(sum(s2))); // fairly sure this is incorrect for nonbinary nodes
    return(sqrt(prod(s2m) / sum(s2m)) * diff_raw + sum(mum ./ s2m) / sum(inv(s2m)));
  }
}
data {
  int N_samples;
  int N_effects;
  matrix[N_samples,N_effects] model_matrix;
  int N_tips;
  int N_nodes; // not equal to phylo$Nnode, but rather phylo$Nnode + length(phylo$tip_label)
  array[N_nodes-1,2] int edge; // important that this be in 'postorder' using ape::reorder.phylo(phylo, 'postorder')
  vector[N_nodes-1] edge_lengths;
  array[N_samples] int idx_tips;
  vector[N_samples] y;
}
transformed data {
  int N_anc_nodes = N_nodes - N_tips;
  array[N_anc_nodes,2] int idx_idx;
  int current = edge[1,1];
  int anc = 1;
  idx_idx[1,1] = 1;
  for(n in 1:(N_nodes-1)) {
    if(edge[n,1] != current) {
      idx_idx[anc,2] = n-1;
      anc += 1;
      idx_idx[anc,1] = n;
      current = edge[n,1];
    }
  }
  idx_idx[anc,2] = (N_nodes-1);
}
parameters {
  vector[N_effects] beta;
  real<lower=0> sigma_phy;
  vector[N_nodes] beta_phy_raw;
  real<lower=0> sigma;
  real<lower=0> theta;
}
transformed parameters {
  vector[N_nodes-1] ou_props = exp(-theta * edge_lengths); // precalculate values related to OU process
  vector[N_nodes-1] sigma2 = 1 - square(ou_props);
  vector[N_nodes] beta_phy = beta_phy_raw;
  for(n in 1:N_anc_nodes) {
    beta_phy[edge[idx_idx[n,1],1]] 
      = ancestral_state_ou_lp(sigma2[idx_idx[n,1]:idx_idx[n,2]], 
                              beta_phy[edge[idx_idx[n,1]:idx_idx[n,2],2]],
                              ou_props[idx_idx[n,1]:idx_idx[n,2]],
                              beta_phy[edge[idx_idx[n,1],1]]);
  }
  beta_phy = sigma_phy / sqrt(2*theta) * beta_phy;
}
model {
  
  theta ~ exponential(1);
  
  beta_phy_raw ~ std_normal();
  beta_phy[N_tips+1] ~ normal(0, sigma_phy / sqrt(2*theta));
  target += log(sigma_phy / sqrt(2*theta));

  y ~ normal(model_matrix * beta + beta_phy[idx_tips], sigma);
  
}
