# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
require(mpcmp)
require(cmp)
require(scales)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")

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
alpha_true <- c(0.1, -0.2)  # Not directly comparable across all models
# alpha_true <- c(0, 0)  # Equidispersion for fast simulation (e.g., to check for errors)
Z <- model.matrix(~ x3)  # Determines dispersion; depends only on x3
log_nu <- c(Z %*% alpha_true)
nu <- exp(log_nu)

lambda <- comp_lambdas(mu, nu)$lambda  # Helper from mpcmp
# This calls a C function from the cmp package that efficiently calculates the moments
var_true <- .C(
  "cmp_cum_all", llambda=as.double(log(lambda)), lnu=as.double(log_nu),
  flag=as.double(1L), n=n, mean=rep(0, n), var=rep(0, n))$var
sd_true <- sqrt(var_true)  # True standard deviation that we're trying to estimate

# Control parameters
tol <- 1e-8
max_iter <- 100
phi_method <- "joint"
stephalving_max <- 10

# Prepare to loop
reps <- 100L
method_names <- c("MPCMP", "GP-1", "EPL")
# method_names <- c("GP-1", "EPL")  # Removing a method from the list removes it from the simulation
num_methods <- length(method_names)
mu_estimates <- array(0, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
mu_estimates[] <- NA # Setting default to NA in case of model-fitting failures
sd_estimates <- mu_estimates
mu_interval_widths <- mu_estimates
sd_interval_widths <- mu_estimates
fitting_times <- matrix(0, nrow=reps, ncol=num_methods, dimnames=list(NULL, method_names))
fitting_times[] <- NA

mu_covered <- array(FALSE, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
mu_covered[] <- NA
sd_covered <- mu_covered

for (i in seq(reps)) {
  print(i)
  y <- mpcmp::rcomp(n, mu=mu, nu=nu)
  # y <- rpois(n, mu)
  dat <- data.frame(x1, x2, x3, y)
  
  # Fit MPCMP model
  if ("MPCMP" %in% method_names ) try({
    start_time <- Sys.time()
    mod_cmp <- fit_cmp(y ~ x1 + x2, y ~ x3, data=dat)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    fitting_times[i,"MPCMP"] <- elapsed_time
    mu_estimates[,i,"MPCMP"] <- mod_cmp$fitted_values
    sd_estimates[,i,"MPCMP"] <- mod_cmp$sd_estimates
    mu_covered[,i,"MPCMP"] <- (mod_cmp$fitted_lower_bounds <= mu) & (mu <= mod_cmp$fitted_upper_bounds)
    mu_interval_widths[,i,"MPCMP"] <- mod_cmp$fitted_interval_widths
    sd_covered[,i,"MPCMP"] <- (mod_cmp$sd_lower_bounds <= sd_true) & (sd_true <= mod_cmp$sd_upper_bounds)
    sd_interval_widths[,i,"MPCMP"] <- mod_cmp$sd_interval_widths
  })
  
  # Fit GP-P model
  if ("GP-1" %in% method_names ) try({
    start_time <- Sys.time()
    mod_quasipois <- glm(y ~ x1 + x2, family=quasipoisson(), data=dat)
    beta_start <- mod_quasipois$coefficients
    phi_start <- summary(mod_quasipois)$dispersion
    mod_gp1 <- gpp(
      y, X, betastart=beta_start, phistart=phi_start, P=1, 
      tol=tol, max_iter=max_iter, phi_method=phi_method, verbose=F,
      stephalving_max=stephalving_max)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    fitting_times[i,"GP-1"] <- elapsed_time
    mu_estimates[,i,"GP-1"] <- mod_gp1$fitted_values
    sd_estimates[,i,"GP-1"] <- mod_gp1$sd_estimates
    mu_covered[,i,"GP-1"] <- (mod_gp1$fitted_lower_bounds <= mu) & (mu <= mod_gp1$fitted_upper_bounds)
    mu_interval_widths[,i,"GP-1"] <- mod_gp1$fitted_interval_widths
    sd_covered[,i,"GP-1"] <- (mod_gp1$sd_lower_bounds <= sd_true) & (sd_true <= mod_gp1$sd_upper_bounds)
    sd_interval_widths[,i,"GP-1"] <- mod_gp1$sd_interval_widths
  })
  
  # Fit EPL model
  if ("EPL" %in% method_names ) try({
    start_time <- Sys.time()
    mod_quasipois <- glm(y ~ x1 + x2, family=quasipoisson(), data=dat)
    beta_start <- mod_quasipois$coefficients
    alpha_start <- rep(0, ncol(Z))
    names(alpha_start) <- colnames(Z)
    alpha_start[1] <- log(summary(mod_quasipois)$dispersion)
    mod_epl <- epl(y, X, Z, betastart = beta_start, alphastart = alpha_start,
                   max_iter=max_iter, stephalving_maxiter=stephalving_max, verbose=FALSE)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    fitting_times[i,"EPL"] <- elapsed_time
    mu_estimates[,i,"EPL"] <- mod_epl$fitted_values
    sd_estimates[,i,"EPL"] <- mod_epl$sd_estimates
    mu_covered[,i,"EPL"] <- (mod_epl$fitted_lower_bounds <= mu) & (mu <= mod_epl$fitted_upper_bounds)
    mu_interval_widths[,i,"EPL"] <- mod_epl$fitted_interval_widths
    sd_covered[,i,"EPL"] <- (mod_epl$sd_lower_bounds <= sd_true) & (sd_true <= mod_epl$sd_upper_bounds)
    sd_interval_widths[,i,"EPL"] <- mod_epl$sd_interval_widths
  })
}

nominal_coverage <- 0.95
format_float <- label_number(0.001)
format_percent <- label_percent(0.1)

# Fitting times
avg_fitting_times <- colMeans(fitting_times)
consolidated_results <- data.frame(`Avg Elapsed Time (Seconds)`=format_float(avg_fitting_times), check.names=FALSE)
boxplot(fitting_times, log="y", main="Model-fitting Time by Method", ylab="Avg Elapsed Time (Seconds)")

# Mean
mu_rmse <- sqrt(apply((mu_estimates - mu)^2, MARGIN=c(1,3), mean))
mu_avg_rmse <- colMeans(mu_rmse)
mu_observation_coverage <- apply(mu_covered, MARGIN=c(1,3), mean)
mu_marginal_coverage <- colMeans(mu_observation_coverage)
mu_coverage_rMSE <- sqrt(colMeans((mu_observation_coverage - nominal_coverage)^2))
mu_avg_interval_width <- apply(mu_interval_widths, MARGIN=c(3), mean)

consolidated_results$`Mean: Avg rMSE` <- format_float(mu_avg_rmse)
consolidated_results$`Mean: Avg Interval Width` <- format_float(mu_avg_interval_width)
consolidated_results$`Mean: Marginal Coverage` <- format_percent(mu_marginal_coverage)
consolidated_results$`Mean: Coverage rMSE` <- format_float(mu_coverage_rMSE)

# SD
sd_rmse <- sqrt(apply((sd_estimates - sd_true)^2, MARGIN=c(1,3), mean))
sd_avg_rmse <- colMeans(sd_rmse)
sd_observation_coverage <- apply(sd_covered, MARGIN=c(1,3), mean)
sd_marginal_coverage <- colMeans(sd_observation_coverage)
sd_coverage_rMSE <- sqrt(colMeans((sd_observation_coverage - nominal_coverage)^2))
sd_avg_interval_width <- apply(sd_interval_widths, MARGIN=c(3), mean)

consolidated_results$`SD: Avg rMSE` <- format_float(sd_avg_rmse)
consolidated_results$`SD: Avg Interval Width` <- format_float(sd_avg_interval_width)
consolidated_results$`SD: Marginal Coverage` <- format_percent(sd_marginal_coverage)
consolidated_results$`SD: Coverage rMSE` <- format_float(sd_coverage_rMSE)

# Make table pretty
final_table <- t(consolidated_results)
View(final_table)
