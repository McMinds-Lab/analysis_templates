rougamma <- function(n,shape,scale,time,sigma2r_mean,start) {

  ## generate random variates from an Ornstein-Uhlenbeck process where the rate of evolution varies over time but the stationary distribution is fixed
  ## n = number of variates to sample
  ## shape = shape of gamma distribution determining overall rate of evolution over a unit of time
  ## scale = scale of the stationary distribution
  ## time = length of time between the starting state and the sampling state
  ## sigma2r_mean = the mean squared rate of evolution /relative to the stationary variance/
  ## start = the starting value from which the distribution evolves

  ## integration of sigma^2 over a time period can be modeled as the sum of many small deviations, of a number proportional to time, each drawn from a gamma distribution
  ## the sum of gamma distributions is another gamma with shape equal to the sum of the shapes of the input
  ## a sum of N deviations each with shape S and expectation E has distribution gamma(N*S,S/E), which can be generalized to a continuous integral with real-valued N (such as time)

    sigma2r_t <- rgamma(n, time*shape, shape/sigma2r_mean)

  ## note that in this model, the expectation for total evolution (sigma2r_t) increases proportionally with time, but so does the relative precision of our estimate
  ## small time periods have more variable and more skewed possibilities, with more 'outliers'
  ## as the base shape goes to infinity, variability of evolutionary rate goes to 0, and the model converges to standard OU evolution with sigma2r_t fixed to sigma2r_mean * t

  ## stationary variance = sigma^2/(2*theta)
  ## if stationary variance is fixed to x, then theta = sigma^2/(2x)
  ## if we have a measure of sigma2 that is already relative to x (sigma2r = sigma^2 / x), then theta = 0.5 * sigma2r
  ## because theta is a simple scaling of sigma2r, the integral of theta over time is a simple scaling of the integral of sigma2r over time:

    theta_t <- 0.5 * sigma2r_t

  ## the deterministic portion of OU evolution (when the attractor is 0, such that we are modeling deviances from the stationary mean) can be calculated as a simple fraction of the starting value, with the fraction defined as:

    oup <- exp(-theta_t)

  ## thus the expectation of the process is:

    mu <- oup * start

  ## the expected deviation from this value is:

    dev <- scale * sqrt(1-oup^2)

  ## leaving a simple draw from a normal distribution to generate our samples:

    return(rnorm(n,mu,dev))

  ## this should be trivially extendable to multivariate OU if it is assumed that all variables' rates of evolution (relative to their stationary variance, and thus thetas) are the same
  ## more work is needed to allow different variables to evolve at different relative rates

}
