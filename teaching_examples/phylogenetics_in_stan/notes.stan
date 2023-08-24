functions {
  real two_brownian_obs_lpdf(real anc, vector desc, vector time) {
    real constant = log(inv(2*pi()));
    real term1 = 0.5 * log(time[1]) - square(desc[1]-anc) / (2*time[1]);
    real term2 = 0.5 * log(time[2]) - square(desc[2]-anc) / (2*time[2]);
    return(constant + term1 + term2);
  }
  real two_brownian_obs_2_lpdf(real anc, vector desc, vector time) {
    real constant = log(inv(2*pi()));
    real mu = inv(time) * desc / sum(inv(time));
    real s2 = prod(time) / sum(time);
    real term1 = 0.5 * log(sum(time)) - square(desc[2]-desc[1]) / (2*sum(time));
    real term2 = 0.5 * log(s2) - square(mu-anc) / (2*s2);
    return(constant + term1 + term2);
  }
  real two_ou_obs_lpdf(real anc, vector desc, vector ou_prop) {
    real constant = log(inv(2*pi()));
    vector[2] mu = ou_prop * anc;
    vector[2] sigma2 = (1 - square(ou_prop));
    real term1 = 0.5 * log(sigma2[1]) - square(desc[1]-mu[1]) / (2*sigma2[1]);
    real term2 = 0.5 * log(sigma2[2]) - square(desc[2]-mu[2]) / (2*sigma2[2]);
    return(constant + term1 + term2);
  }
  real two_ou_obs_2_1_lpdf(real anc, vector desc, vector ou_prop) {
    real constant = log(inv(2*pi()));
    vector[2] sigma2 = (1 - square(ou_prop));
    real term1 = 0.5 * log(sigma2[1]) - square(desc[1] - ou_prop[1] * anc) / (2*sigma2[1]);
    real term2 = 0.5 * log(sigma2[2]) - square(desc[2] - ou_prop[2] * anc) / (2*sigma2[2]);
    return(constant + term1 + term2);
  }
  real two_ou_obs_2_lpdf(real anc, vector desc, vector ou_prop) {
    real constant = log(inv(2*pi()));
    vector[2] descM = desc ./ ou_prop;
    vector[2] sigma2 = (1 - square(ou_prop)) ./ square(ou_prop);

    real mu = inv(sigma2) * descM / sum(inv(sigma2));
    real s2 = prod(sigma2) / sum(sigma2);
    real term1 = 0.5 * log(sum(sigma2)) - square(descM[2]-descM[1]) / (2*sum(sigma2));
    real term2 = 0.5 * log(s2) - square(mu-anc) / (2*s2);
    return(constant + term1 + term2);
  }
}
