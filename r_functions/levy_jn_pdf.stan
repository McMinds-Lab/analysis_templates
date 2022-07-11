functions {
  real inverse_gaussian_cf(real time, real delta, real alpha, real xi) {
    return time * delta * (square(alpha) - sqrt(square(alpha) + square(xi)));
  }
  real levy_branch_sd() {
    return integrate_1d();
  }
}
//doi: 10.1073/pnas.1710920114
//https://discourse.mc-stan.org/t/implementing-levy-alpha-stable-distribution/3565
//https://mc-stan.org/docs/2_29/functions-reference/functions-1d-integrator.html
//https://mc-stan.org/docs/2_29/stan-users-guide/integrate-1d.html
//https://academic.oup.com/sysbio/article/62/2/193/1668359
