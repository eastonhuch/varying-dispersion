# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
require(mpcmp)
require(cmp)
source("cmp-helpers.R")

# Generate some data
set.seed(1)
n <- 100L
x1 <- rnorm(n)  # We'll keep the design fixed across simulations
x2 <- rnorm(n)
x3 <- rnorm(n)

# NOTE: The mpcmp package is a little temperamental
# You might run into errors if you make the coefficients much larger
beta_true <- c(3, 0.05, -0.1)  # Regression coefficients for the mean
log_mu <- beta_true[1] + beta_true[2]*x1 - beta_true[3]*x2
mu <- exp(log_mu)
log_nu <- 0.1 - 0.2*x3  # Determines dispersion; depends only on x3
nu <- exp(log_nu)

lambda <- comp_lambdas(mu, nu)$lambda  # Helper from mpcmp
# This calls a C function from the cmp package that efficiently calculates the moments
var_true <- .C(
  "cmp_cum_all", llambda=as.double(log(lambda)), lnu=as.double(log_nu),
  flag=as.double(1L), n=n, mean=rep(0, n), var=rep(0, n))$var
sd_true <- sqrt(var_true)  # True standard deviation that we're trying to estimate

# Prepare to loop
reps <- 50L
mu_estimates <- matrix(FALSE, nrow=n, ncol=reps) # Setting types for the matrices
mu_estimates[] <- NA # Setting default to NA in case of model-fitting failures
sd_estimates <- mu_covered

mu_covered <- mu_estimates
mu_interval_widths <- matrix(0, nrow=n, ncol=reps)
mu_interval_widths[] <- NA
sd_covered <- mu_covered
sd_interval_widths <- mu_interval_widths

for (i in seq(reps)) {
  print(i)
  y <- mpcmp::rcomp(n, mu=mu, nu=nu)
  dat <- data.frame(x1, x2, x3, y)
  
  # Fit MPCMP model
  try({
    mod_cmp <- fit_cmp(y ~ x1 + x2, y ~ x3, data=dat)
    mu_estimates[,i] <- mod_cmp$fitted_values
    sd_estimates[,i] <- mod_cmp$sd_estimates
    mu_covered[,i] <- (mod_cmp$fitted_lower_bounds <= mu) & (mu <= mod_cmp$fitted_upper_bounds)
    mu_interval_widths[,i] <- mod_cmp$fitted_interval_widths
    sd_covered[,i] <- (mod_cmp$sd_lower_bounds <= sd_true) & (sd_true <= mod_cmp$sd_upper_bounds)
    sd_interval_widths[,i] <- mod_cmp$sd_interval_widths
  })
}

nominal_coverage <- 0.95

# Performance of point estimates
mu_rmse <- rowMeans((mu_estimates - mu)^2)
mu_avg_rmse <- mean(mu_rmse)
mu_avg_rmse

sd_rmse <- rowMeans((sd_estimates - sd_true)^2)
sd_avg_rmse <- mean(sd_rmse)
sd_avg_rmse

# Performance of confidence intervals for fitted values
mu_observation_coverage <- rowMeans(mu_covered)
mu_marginal_coverage <- mean(mu_observation_coverage)
mu_coverage_rMSE <- sqrt(mean((mu_observation_coverage - nominal_coverage)^2))
mean(mu_interval_widths)
mu_marginal_coverage
mu_coverage_rMSE

# Performance of confidence intervals for standard deviation
sd_observation_coverage <- rowMeans(sd_covered)
sd_marginal_coverage <- mean(sd_observation_coverage)
sd_coverage_rMSE <- sqrt(mean((sd_observation_coverage - nominal_coverage)^2))
mean(sd_interval_widths)
sd_marginal_coverage
sd_coverage_rMSE