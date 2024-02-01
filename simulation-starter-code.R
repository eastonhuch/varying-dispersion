# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
require(mpcmp)
require(cmp)
source("cmp-helpers.R")
source("gpp.R")

# Generate some data
set.seed(1)
n <- 100L
x1 <- rnorm(n)  # We'll keep the design fixed across simulations
x2 <- rnorm(n)
x3 <- rnorm(n)

# NOTE: The mpcmp package is a little temperamental
# You might run into errors if you make the coefficients much larger
beta_true <- c(3, 0.05, -0.1)  # Regression coefficients for the mean
X <- model.matrix(~ x1 + x2)
log_mu <- c(X %*% beta_true)
mu <- exp(log_mu)
log_nu <- 0.1 - 0.2*x3  # Determines dispersion; depends only on x3
nu <- exp(log_nu)

lambda <- comp_lambdas(mu, nu)$lambda  # Helper from mpcmp
# This calls a C function from the cmp package that efficiently calculates the moments
var_true <- .C(
  "cmp_cum_all", llambda=as.double(log(lambda)), lnu=as.double(log_nu),
  flag=as.double(1L), n=n, mean=rep(0, n), var=rep(0, n))$var
sd_true <- sqrt(var_true)  # True standard deviation that we're trying to estimate

# GP-P control parameters
tol <- 1e-12
max_iter <- 200
phi_method <- "joint"
stephalving_max <- 200
phi_start <- 0

# Prepare to loop
reps <- 400L
method_names <- c("MPCMP", "GP-1")
num_methods <- length(method_names)
mu_estimates <- array(0, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
mu_estimates[] <- NA # Setting default to NA in case of model-fitting failures
sd_estimates <- mu_estimates
mu_interval_widths <- mu_estimates
sd_interval_widths <- mu_interval_widths

mu_covered <- array(FALSE, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
mu_covered[] <- NA
sd_covered <- mu_covered

for (i in seq(reps)) {
  print(i)
  y <- mpcmp::rcomp(n, mu=mu, nu=nu)
  dat <- data.frame(x1, x2, x3, y)
  
  # Fit MPCMP model
  try({
    mod_cmp <- fit_cmp(y ~ x1 + x2, y ~ x3, data=dat)
    mu_estimates[,i,"MPCMP"] <- mod_cmp$fitted_values
    sd_estimates[,i,"MPCMP"] <- mod_cmp$sd_estimates
    mu_covered[,i,"MPCMP"] <- (mod_cmp$fitted_lower_bounds <= mu) & (mu <= mod_cmp$fitted_upper_bounds)
    mu_interval_widths[,i,"MPCMP"] <- mod_cmp$fitted_interval_widths
    sd_covered[,i,"MPCMP"] <- (mod_cmp$sd_lower_bounds <= sd_true) & (sd_true <= mod_cmp$sd_upper_bounds)
    sd_interval_widths[,i,"MPCMP"] <- mod_cmp$sd_interval_widths
  })
  
  # Fit GP-P model
  try({
    mod_pois <- glm(y ~ x1 + x2, family=poisson, data=dat)
    beta_start <- mod_pois$coefficients
    gp1 <- gpp(
      y, X, betastart=beta_start, phistart=phi_start, P=1, 
      tol=tol, max_iter=max_iter, phi_method=phi_method, verbose=F,
      stephalving_max=stephalving_max)
    mu_estimates[,i,"GP-1"] <- gp1$fitted_values
    sd_estimates[,i,"GP-1"] <- gp1$sd_estimates
    mu_covered[,i,"GP-1"] <- (gp1$fitted_lower_bounds <= mu) & (mu <= gp1$fitted_upper_bounds)
    mu_interval_widths[,i,"GP-1"] <- gp1$fitted_interval_widths
    sd_covered[,i,"GP-1"] <- (gp1$sd_lower_bounds <= sd_true) & (sd_true <= gp1$sd_upper_bounds)
    sd_interval_widths[,i,"GP-1"] <- gp1$sd_interval_widths
  })
}

nominal_coverage <- 0.95

# Performance of point estimates
mu_rmse <- sqrt(apply((mu_estimates - mu)^2, MARGIN=c(1,3), mean))
mu_avg_rmse <- colMeans(mu_rmse)
mu_avg_rmse

sd_rmse <- sqrt(apply((sd_estimates - sd_true)^2, MARGIN=c(1,3), mean))
sd_avg_rmse <- colMeans(sd_rmse)
sd_avg_rmse

# Performance of confidence intervals for fitted values
mu_observation_coverage <- apply(mu_covered, MARGIN=c(1,3), mean)
mu_marginal_coverage <- colMeans(mu_observation_coverage)
mu_coverage_rMSE <- sqrt(colMeans((mu_observation_coverage - nominal_coverage)^2))
mu_avg_interval_width <- apply(mu_interval_widths, MARGIN=c(3), mean)

mu_marginal_coverage
mu_coverage_rMSE
mu_avg_interval_width

# Performance of confidence intervals for standard deviation
sd_observation_coverage <- apply(sd_covered, MARGIN=c(1,3), mean)
sd_marginal_coverage <- colMeans(sd_observation_coverage)
sd_coverage_rMSE <- sqrt(colMeans((sd_observation_coverage - nominal_coverage)^2))
sd_avg_interval_width <- apply(sd_interval_widths, MARGIN=c(3), mean)

sd_marginal_coverage
sd_coverage_rMSE
sd_avg_interval_width