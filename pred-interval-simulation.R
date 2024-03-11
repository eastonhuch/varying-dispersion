# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
require(MASS)
require(tidyverse)
require(mpcmp)
require(cmp)
require(scales)
require(xtable)
source("discrete-log-normal.R")

# Set parameters
set.seed(1)
master_formula <- y ~ x1
y <- x1 <- rep(0, 10)  # Just for calculating p
Z <- X <- model.matrix(master_formula)
p <- ncol(X)
beta_idx <- seq(p)
alpha_idx <- p + beta_idx
nominal_coverage <- 0.95  # This is hard-coded in the model-fitting functions
beta_intercept_prior_mean <- 4
alpha_intercept_prior_mean <- -1
theta_prior_mean <- rep(0, 2*p)
theta_prior_mean[1] <- beta_intercept_prior_mean
theta_prior_mean[p+1] <- alpha_intercept_prior_mean

intercept_prior_sd <- 0.2
coef_prior_sd <- 0.05
theta_prior_precision <- diag(2*p) * coef_prior_sd^(-2)
theta_prior_precision[1, 1] <- intercept_prior_sd^(-2)
theta_prior_precision[p+1, p+1] <- intercept_prior_sd^(-2)
theta_prior_cov <- chol2inv(chol(theta_prior_precision))

# Control parameters
tol <- 1e-8
max_iter <- 100
phi_method <- "joint"
stephalving_max <- 10

# Arrays for results
n_vals <- c(10, 30L, 100L, 400L)
num_n_vals <- length(n_vals)
reps <- 40L
n_bootstraps <- 400L
method_names <- c("Plug-in", "Asymp. Bayes", "Full Bayes")
num_methods <- length(method_names)

marginal_coverages <- matrix(0, nrow=num_n_vals, ncol=num_methods, dimnames = list(n_vals, method_names))
marginal_coverages[] <- NA
marginal_coverage_ses <- marginal_coverages
coverage_rmses <- marginal_coverages
coverage_rmses_ses <- marginal_coverages
interval_widths <- marginal_coverages
interval_width_ses <- marginal_coverages
interval_width_medians <- marginal_coverages
interval_width_medians_ses <- marginal_coverages

# Starting values
beta_start <- theta_prior_mean[beta_idx]
names(beta_start) <- colnames(X)
alpha_start <- theta_prior_mean[alpha_idx]
names(alpha_start) <- colnames(Z)

# Loop
for (n in n_vals) {
  cat("n: ", n, "\n")
  
  # Create data
  x1 <- rnorm(n)  # We'll keep the design fixed across simulations
  x2 <- x1 + rnorm(n, sd=0.05)
  y <- rep(0, n)  # We'll create this later
  X <- model.matrix(master_formula)
  Z <- model.matrix(master_formula)
  
  y_star_covered <- array(FALSE, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
  y_star_covered[] <- NA
  y_star_interval_widths <- array(0, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
  y_star_interval_widths[] <- NA
  
  for (i in seq(reps)) {
    print(i)
    theta_true <- mvrnorm(1, theta_prior_mean, theta_prior_cov)
    beta_true <- theta_true[beta_idx]
    alpha_true <- theta_true[alpha_idx]
    z_mu <- c(X %*% beta_true)
    z_sigma <- c(exp(Z %*% alpha_true))
    y <- floor(exp(rnorm(n, mean=z_mu, sd=z_sigma)))
    y_star <- floor(exp(rnorm(n, mean=z_mu, sd=z_sigma)))
    dat <- data.frame(x1, y)
    
    # Plug-in Method
    if ("Plug-in" %in% method_names ) try({
      # Fit model
      mod_plug_in <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton", pred_interval_method="Plug-in",
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
      
      # Store results
      y_star_covered[,i,"Plug-in"] <- (mod_plug_in$pred_lower_bounds <= y_star) & (y_star <= mod_plug_in$pred_upper_bounds)
      y_star_interval_widths[,i,"Plug-in"] <- mod_plug_in$pred_interval_widths
    })
    
    # Asymp. Bayes method
    if ("Asymp. Bayes" %in% method_names ) try({
      # Fit model
      mod_approx_bayes <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton", pred_interval_method="Asymp. Bayes",
                         max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
      
      # Store results
      y_star_covered[,i,"Asymp. Bayes"] <- (mod_approx_bayes$pred_lower_bounds <= y_star) & (y_star <= mod_approx_bayes$pred_upper_bounds)
      y_star_interval_widths[,i,"Asymp. Bayes"] <- mod_approx_bayes$pred_interval_widths
    })
    
    # Full Bayes method
    if ("Full Bayes" %in% method_names ) try({
      # Fit model
      mod_full_bayes <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton", pred_interval_method="Full Bayes",
                            prior_mean=theta_prior_mean, prior_precision=theta_prior_precision,
                            max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
      
      # Store results
      y_star_covered[,i,"Full Bayes"] <- (mod_full_bayes$pred_lower_bounds <= y_star) & (y_star <= mod_full_bayes$pred_upper_bounds)
      y_star_interval_widths[,i,"Full Bayes"] <- mod_full_bayes$pred_interval_widths
    })
  }
  
  # Marginal coverage
  rep_coverage <- apply(y_star_covered, MARGIN=c(2,3), function(x) mean(x, na.rm=TRUE))
  marginal_coverages[as.character(n),] <- colMeans(rep_coverage, na.rm=TRUE)
  marginal_coverage_ses[as.character(n),] <- apply(rep_coverage, 2, function(x) sd(x, na.rm=TRUE)) / sqrt(reps)
  
  # Coverage rMSE
  observation_coverage <- apply(y_star_covered, MARGIN=c(1,3), function(x) mean(x, na.rm=TRUE))
  coverage_rmses[as.character(n),] <- sqrt(colMeans((observation_coverage - nominal_coverage)^2, na.rm=TRUE))
  coverage_rmses_ses[as.character(n),] <- apply(y_star_covered, MARGIN=3, function(M) {
    bootstrap_coverage_rmses <- rep(0, n_bootstraps)
    for (j in seq(n_bootstraps)) {
      M_boot <- M[,sample(ncol(M), replace=TRUE)]
      bootstrap_observation_coverages <- rowMeans(M_boot)
      bootstrap_coverage_rmses[j] <- sqrt(mean((bootstrap_observation_coverages - nominal_coverage)^2))
    }
    sd(bootstrap_coverage_rmses)
  })
  
  # Interval widths
  # rep_interval_widths <- apply(y_star_interval_widths, MARGIN=c(2,3), function(x) mean(x, na.rm=TRUE))
  # interval_widths[as.character(n),] <- colMeans(rep_interval_widths, na.rm=TRUE)
  # interval_width_ses[as.character(n),] <- apply(rep_interval_widths, 2, function(x) sd(x, na.rm=TRUE)) / sqrt(reps)
  rep_interval_width_medians <- apply(y_star_interval_widths, MARGIN=c(2,3), function(x) median(x, na.rm=TRUE))
  interval_width_medians[as.character(n),] <- apply(rep_interval_width_medians, 2, function(x) median(x, na.rm=TRUE))
  interval_width_medians_ses[as.character(n),] <- apply(rep_interval_width_medians, 2, function(x) {
    medians <- rep(0, n_bootstraps)
    for (j in seq(n_bootstraps)) {
      x_boot <- x[sample(length(x), replace=TRUE)]
      medians[j] <- median(x_boot)
    }
    sd(medians)
  })
}

# Display results
format_float <- label_number(0.001)
format_float_1d <- label_number(0.1)
format_percent <- label_percent(0.1)
result_df <- data.frame(n=c("Method", paste0("n=", n_vals)))

# Marginal Coverage
for (method_name in method_names) {
  result_df <- result_df %>% mutate(
    !!(paste0(method_name, ": Marginal Coverage")) := c(method_name, paste0(format_float(marginal_coverages[,method_name]), " (", format_float(marginal_coverage_ses[,method_name]), ")"))
  )
}

# Coverage rMSE
for (method_name in method_names) {
  result_df <- result_df %>% mutate(
    !!(paste0(method_name, ": Coverage rMSE")) := c(method_name, paste0(format_float(coverage_rmses[,method_name]), " (", format_float(coverage_rmses_ses[,method_name]), ")"))
  )
}

# Interval Length
for (method_name in method_names) {
  result_df <- result_df %>% mutate(
    !!(paste0(method_name, ": Median Interval Length")) := c(method_name, paste0(format_float_1d(interval_width_medians[,method_name]), " (", format_float_1d(interval_width_medians_ses[,method_name]), ")"))
  )
}

row.names(result_df) <- result_df$n
result_df$n <- NULL
# View(result_df)

result_df_t <- t(result_df)
metric_names <- c(
  "\\multirow{3}{*}{\\parbox{5em}{Marginal Coverage}}", "", "",
  "\\multirow{3}{*}{\\parbox{5em}{Coverage rMSE}}", "", "",
  "\\multirow{3}{*}{\\parbox{5em}{Median Length}}", "", ""
)
result_df_t <- cbind(Summary=metric_names, result_df_t)

# metric_names <- c("Marginal Coverage", "Coverage rMSE", "Median Length")
xtable(
  result_df_t,
  label="tab:pred_sim",
  caption="Coverage and median interval lengths for three prediction interval methods. The Plug-in and Asymp. Bayes methods correspond with those explained in Section \\ref{sec:log_normal_inference}. The prior for the fully Bayesian method matches the data generating distribution."
  ) %>%
  print(
    include.colnames=TRUE,
    include.rownames=FALSE,
    sanitize.text.function=identity,
    hline.after=c(-1, 0, 3, 6, 9)
    # add.to.row=list(pos=list(0), command=paste0(" ", paste0('& \\multicolumn{', num_methods, '}{c}{', metric_names, '} ', collapse=''), '\\\\')),
  )
