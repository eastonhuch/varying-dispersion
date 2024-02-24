#include<RcppArmadillo.h>
#include <stdexcept>
#include <iostream>

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double integrate_dln_scalar(double mu, double sigma, double k, bool calc_deriv=false, bool wrt_mu=false, double sd_multiplier=20.0) {
  // Numeric integration to approximate moments and their derivatives
  // This will not be accurate when there is high skew/kurtosis
  double sigma2 = sigma * sigma;
  double ln_mean = exp(mu + sigma2 / 2.0);
  double ln_var = (exp(sigma2) - 1.0) * exp(2.0 * mu + sigma2);
  double ln_sd = sqrt(ln_var);
  int y_min = floor(ln_mean - sd_multiplier * ln_sd);
  if (y_min < 0) y_min = 0;
  int y_max = ceil(ln_mean + sd_multiplier * ln_sd);
  double z_lower, z_upper, dist_lower, dist_upper, factor_lower, factor_upper, y_k;
  double result = 0, inc = 0, inc_last = 0;
  if (y_max < 0) {
    throw std::overflow_error("moments of discrete log-normal are too large");
  }
  int y_inc = ceil((y_max - y_min) / 10000.0); // skip values if needed
  // std::cout << "y_min: " << y_min << std::endl;
  // std::cout << "y_max: " << y_max << std::endl;
  // std::cout << "y_inc: " << y_inc << std::endl;
  if (calc_deriv && wrt_mu) {
    factor_lower = 1.0 / sigma;
    factor_upper = 1.0 / sigma;
  }
  if (!calc_deriv) {
    factor_lower = 1.0;
    factor_upper = 1.0;
  }
  for (int y = y_min; y < y_max; y+=y_inc) {
    // dist_lower
    if (y <= 0) {
      dist_lower = 0;
    } else {
      z_lower = (log(y) - mu) / sigma;
      if (calc_deriv) {
        dist_lower = -R::dnorm(z_lower, 0.0, 1.0, false);
      } else {
        dist_lower = R::pnorm(z_lower, 0.0, 1.0, true, false);
      }
    }

    // dist_upper
    z_upper = (log(y + 1.0) - mu) / sigma;
    if (calc_deriv) {
      dist_upper = -R::dnorm(z_upper, 0.0, 1.0, false);
    } else {
      dist_upper = R::pnorm(z_upper, 0.0, 1.0, true, false);
    }

    // factors
    if (calc_deriv && (!wrt_mu)) {
      factor_lower = z_lower / sigma;
      factor_upper = z_upper / sigma;
    }

    // increment
    y_k = pow(y, k);
    inc_last = inc;
    inc = y_k * (factor_upper * dist_upper - factor_lower * dist_lower);
    result += (inc_last + inc)/2.0;
  }
  
  return y_inc * result;
}

//[[Rcpp::export]]
Rcpp::NumericVector integrate_dln(Rcpp::NumericVector mu, Rcpp::NumericVector sigma, double k, bool calc_deriv=false, bool wrt_mu=false, double sd_multiplier=20.0) {
  if (mu.size() != sigma.size()) {
    throw std::invalid_argument("mu and sigma must have equal length");
  }
  Rcpp::NumericVector result(mu.size());
  for (int i = 0; i < mu.size(); i++) {
    result[i] = integrate_dln_scalar(mu[i], sigma[i], k, calc_deriv, wrt_mu, sd_multiplier);
  }
  return result;
}
  

/*** R
# approximate_moment <- function(mu, sigma, k=1, plot_dist=TRUE) {
#   sigma2 <- sigma^2
#   ln_mean <- exp(mu + sigma2/2)
#   ln_var <- (exp(sigma2) - 1) * exp(2*mu + sigma2)
#   ln_sd <- sqrt(ln_var)
#   y_min <- max(0, floor(ln_mean - 20*ln_sd))
#   y_max <- ceiling(ln_mean + 20*ln_sd)
#   ys <- seq(y_min, y_max)
#   ps <- pnorm((log(ys+1)-mu)/sigma) - pnorm((log(ys)-mu)/sigma)
#   if (plot_dist) plot(ys, ps, type="l")
#   sum(ys^k * ps)
# }
# mu <- 6
# sigma <- 0.2
# k <- 1
# delta <- 1e-8
# r_est <- approximate_moment(mu, sigma, k)
# r_est_delta <- approximate_moment(mu, sigma+delta, k)
# (r_est_delta - r_est) / delta
# integrate_dln_scalar(mu, sigma, k, calc_deriv=TRUE, wrt_mu=FALSE)
# integrate_dln(rep(mu, 2), rep(sigma, 2), k, calc_deriv=TRUE, wrt_mu=FALSE)
*/
