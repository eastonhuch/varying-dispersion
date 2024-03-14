# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
require(tidyverse)
require(xtable)
require(mpcmp)
require(cmp)
require(scales)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

# Generate some data
set.seed(1)
n <- 400L
x1 <- rnorm(n)  # We'll keep the design fixed across simulations
x2 <- rnorm(n)

# NOTE: The mpcmp package is a little temperamental
# You might run into errors if you make the coefficients much larger
beta_true <- c(3, 0.05, -0.1, 0.02)  # Regression coefficients for the mean
y <- rep(0, n)  # We'll create this later
master_formula <- y ~ x1*x2
X <- model.matrix(master_formula)
log_mu <- c(X %*% beta_true)
mu <- exp(log_mu)
alpha_true <- c(0.1, 0, -0.2, 0.05)  # Not directly comparable across all models
# alpha_true <- c(0, 0)  # Equidispersion for fast simulation (e.g., to check for errors)
Z <- model.matrix(master_formula)
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
reps <- 400L
method_names <- c("MPCMP", "GP-1", "EPL", "DLN-Newton", "DLN-EM")
# method_names <- c("GP-1", "EPL")  # Removing a method from the list removes it from the simulation
num_methods <- length(method_names)
mu_estimates <- array(0, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
mu_estimates[] <- NA # Setting default to NA in case of model-fitting failures
sd_estimates <- mu_estimates
mu_interval_widths <- mu_estimates
sd_interval_widths <- mu_estimates
fitting_times <- matrix(0, nrow=reps, ncol=num_methods, dimnames=list(NULL, method_names))
fitting_times[] <- NA
used_score_for_cov <- matrix(FALSE, nrow=reps, ncol=num_methods, dimnames=list(NULL, method_names))
used_score_for_cov[] <- NA

mu_covered <- array(FALSE, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
mu_covered[] <- NA
sd_covered <- mu_covered

for (i in seq(reps)) {
  print(i)
  y <- mpcmp::rcomp(n, mu=mu, nu=nu)
  # y <- rpois(n, mu)
  dat <- data.frame(x1, x2, y)
  
  # Fit MPCMP model
  if ("MPCMP" %in% method_names ) try({
    start_time <- Sys.time()
    mod_cmp <- fit_cmp(master_formula, master_formula, data=dat)
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
    mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
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
    mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
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
  
  # Fit DLN-Newton model
  if ("DLN-Newton" %in% method_names ) try({
    start_time <- Sys.time()
    mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
    beta_start <- mod_quasipois$coefficients
    alpha_start <- rep(0, ncol(Z))
    names(alpha_start) <- colnames(Z)
    alpha_start[1] <- log(sqrt(summary(mod_quasipois)$dispersion))
    mod_dln_newton <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton",
                   max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    fitting_times[i,"DLN-Newton"] <- elapsed_time
    mu_estimates[,i,"DLN-Newton"] <- mod_dln_newton$fitted_values
    sd_estimates[,i,"DLN-Newton"] <- mod_dln_newton$sd_estimates
    mu_covered[,i,"DLN-Newton"] <- (mod_dln_newton$fitted_lower_bounds <= mu) & (mu <= mod_dln_newton$fitted_upper_bounds)
    mu_interval_widths[,i,"DLN-Newton"] <- mod_dln_newton$fitted_interval_widths
    sd_covered[,i,"DLN-Newton"] <- (mod_dln_newton$sd_lower_bounds <= sd_true) & (sd_true <= mod_dln_newton$sd_upper_bounds)
    sd_interval_widths[,i,"DLN-Newton"] <- mod_dln_newton$sd_interval_widths
    used_score_for_cov[i,"DLN-Newton"] <- mod_dln_newton$covariance_via_score
  })
  
  # Fit DLN-EM model
  if ("DLN-EM" %in% method_names ) try({
    start_time <- Sys.time()
    mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
    beta_start <- mod_quasipois$coefficients
    alpha_start <- rep(0, ncol(Z))
    names(alpha_start) <- colnames(Z)
    alpha_start[1] <- log(sqrt(summary(mod_quasipois)$dispersion))
    mod_dln_em <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="EM",
                          max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    fitting_times[i,"DLN-EM"] <- elapsed_time
    mu_estimates[,i,"DLN-EM"] <- mod_dln_em$fitted_values
    sd_estimates[,i,"DLN-EM"] <- mod_dln_em$sd_estimates
    mu_covered[,i,"DLN-EM"] <- (mod_dln_em$fitted_lower_bounds <= mu) & (mu <= mod_dln_em$fitted_upper_bounds)
    mu_interval_widths[,i,"DLN-EM"] <- mod_dln_em$fitted_interval_widths
    sd_covered[,i,"DLN-EM"] <- (mod_dln_em$sd_lower_bounds <= sd_true) & (sd_true <= mod_dln_em$sd_upper_bounds)
    sd_interval_widths[,i,"DLN-EM"] <- mod_dln_em$sd_interval_widths
    used_score_for_cov[i,"DLN-EM"] <- mod_dln_newton$covariance_via_score
  })
}

nominal_coverage <- 0.95
format_float <- label_number(0.001)
format_percent <- label_percent(0.1)
n_bootstraps <- 400L
print_val_with_se <- function(val, se, percent=FALSE) {
  digits <- abs(min(c(floor(log10(min(se))), 0)))
  if (percent) digits <- max(c(digits - 2, 0))
  multiplier <- ifelse(percent, 100, 1)
  fmt <- paste0("%.", digits, "f")
  val_digits <- sprintf(fmt, val*multiplier)
  se_digits <- sprintf(fmt, se*multiplier)
  if (percent) {
    res <- paste0(val_digits, "% (", se_digits, "%)")
  } else {
    res <- paste0(val_digits, " (", se_digits, ")")
  }
  res
}

# Fitting times
avg_fitting_times <- colMeans(fitting_times)
avg_fitting_times_ses <- apply(fitting_times, 2, sd) / sqrt(reps)
avg_fitting_times_chr <- print_val_with_se(avg_fitting_times, avg_fitting_times_ses)
consolidated_results <- data.frame(`Avg Elapsed Time (Seconds)`=avg_fitting_times_chr, check.names=FALSE)
row.names(consolidated_results) <- method_names
boxplot(fitting_times, log="y", main="Model-fitting Time by Method", ylab="Elapsed Time (Seconds)")

# Mean
mu_rmse <- sqrt(apply((mu_estimates - mu)^2, MARGIN=c(2,3), mean))
mu_avg_rmse <- colMeans(mu_rmse)
mu_avg_rmse_ses <- apply(mu_rmse, 2, sd) / sqrt(reps)
mu_rep_coverage <- apply(mu_covered, MARGIN=c(2,3), mean)
mu_marginal_coverage <- colMeans(mu_observation_coverage)
mu_marginal_coverage_ses <- apply(mu_observation_coverage, 2, sd) / sqrt(reps)
mu_rep_interval_widths <- apply(mu_interval_widths, MARGIN=c(2, 3), mean)
mu_avg_interval_widths <- colMeans(mu_rep_interval_widths)
mu_avg_interval_widths_ses <- apply(mu_rep_interval_widths, 2, sd) / sqrt(reps)

mu_observation_coverage <- apply(mu_covered, MARGIN=c(1,3), mean)
mu_coverage_rMSE <- sqrt(colMeans((mu_observation_coverage - nominal_coverage)^2))
mu_coverage_rMSE_ses <- apply(mu_covered, MARGIN=3, function(M) {
  bootstrap_coverage_rmses <- rep(0, n_bootstraps)
  for (j in seq(n_bootstraps)) {
    M_boot <- M[,sample(ncol(M), replace=TRUE)]
    bootstrap_observation_coverages <- rowMeans(M_boot)
    bootstrap_coverage_rmses[j] <- sqrt(mean((bootstrap_observation_coverages - nominal_coverage)^2))
  }
  sd(bootstrap_coverage_rmses)
})

consolidated_results$`Mean: Avg rMSE` <- print_val_with_se(mu_avg_rmse, mu_avg_rmse_ses)
consolidated_results$`Mean: Avg Interval Width` <- print_val_with_se(mu_avg_interval_widths, mu_avg_interval_widths_ses)
consolidated_results$`Mean: Marginal Coverage` <- print_val_with_se(mu_marginal_coverage, mu_marginal_coverage_ses, percent=TRUE)
consolidated_results$`Mean: Coverage rMSE` <- print_val_with_se(mu_coverage_rMSE, mu_coverage_rMSE_ses, percent=TRUE)

# SD
sd_rmse <- sqrt(apply((sd_estimates - sd_true)^2, MARGIN=c(2,3), mean))
sd_avg_rmse <- colMeans(sd_rmse)
sd_avg_rmse_ses <- apply(sd_rmse, 2, sd) / sqrt(reps)
sd_rep_coverage <- apply(sd_covered, MARGIN=c(2,3), mean)
sd_marginal_coverage <- colMeans(sd_observation_coverage)
sd_marginal_coverage_ses <- apply(sd_observation_coverage, 2, sd) / sqrt(reps)
sd_rep_interval_widths <- apply(sd_interval_widths, MARGIN=c(2, 3), mean)
sd_avg_interval_widths <- colMeans(sd_rep_interval_widths)
sd_avg_interval_widths_ses <- apply(sd_rep_interval_widths, 2, sd) / sqrt(reps)

sd_observation_coverage <- apply(sd_covered, MARGIN=c(1,3), mean)
sd_coverage_rMSE <- sqrt(colMeans((sd_observation_coverage - nominal_coverage)^2))
sd_coverage_rMSE_ses <- apply(sd_covered, MARGIN=3, function(M) {
  bootstrap_coverage_rmses <- rep(0, n_bootstraps)
  for (j in seq(n_bootstraps)) {
    M_boot <- M[,sample(ncol(M), replace=TRUE)]
    bootstrap_observation_coverages <- rowMeans(M_boot)
    bootstrap_coverage_rmses[j] <- sqrt(mean((bootstrap_observation_coverages - nominal_coverage)^2))
  }
  sd(bootstrap_coverage_rmses)
})

consolidated_results$`SD: Avg rMSE` <- print_val_with_se(sd_avg_rmse, sd_avg_rmse_ses)
consolidated_results$`SD: Avg Interval Width` <- print_val_with_se(sd_avg_interval_widths, sd_avg_interval_widths_ses)
consolidated_results$`SD: Marginal Coverage` <- print_val_with_se(sd_marginal_coverage, sd_marginal_coverage_ses, percent=TRUE)
consolidated_results$`SD: Coverage rMSE` <- print_val_with_se(sd_coverage_rMSE, sd_coverage_rMSE_ses, percent=TRUE)

# Plot avg SD vs. x2
sd_order <- order(sd_true)
plot(x2[sd_order], sd_true[sd_order], xlab=expression(x[3]), ylab="SD", ylim=c(2.5, 7))
points(x2[sd_order], sd_avg[sd_order, "DLN-Newton"], col=2)
points(x2[sd_order], sd_avg[sd_order, "MPCMP"], col=3)

# How often did we use score-based covariance?
colSums(used_score_for_cov)

# Make table pretty
final_table <- t(consolidated_results)
View(final_table)

xtable(
  final_table,
  label="tab:main_sim",
  align=paste0("l", paste0(rep("r", num_methods), collapse="")),
  caption=paste0(
    "Method comparison for data simulated from MPCMP (n=", n, ").",
    " Avg rMSE averages over the covariate values.",
    " We fix covariates across repetitions to assess conditional coverage.",
    " Coverage rMSE calculates an rMSE across observation-level coverages (averaging over repetitions).")
) %>%
  print(
    include.colnames=TRUE,
    include.rownames=TRUE,
    hline.after=c(-1, 0, 1, 5, 9)
  )
