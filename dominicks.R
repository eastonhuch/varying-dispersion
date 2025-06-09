# Dominick's Analysis
require(tidyverse)
require(xtable)
require(scales)
require(splines)
require(lubridate)
require(MASS)
require(glmmTMB)
require(mvtnorm)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

# Read in data
omitted_upcs <- c(1780013217)  # Sudden drops causing model fitting issues
dat <- read.csv("dominicks-data/all-categories.csv") %>%
  filter(!(upc %in% omitted_upcs))
dat$date <- as.Date(dat$date)
dat$yday <- yday(dat$date)
dat$yday_radians <- dat$yday / 366 * 2 * pi
upcs <- sort(unique(dat$upc))
num_upcs <- length(upcs)

# Model specification
get_formula <- function(spline_dof, n_trig_terms) {
  # Return formula given terms
  if (spline_dof >= 1) {
    X_spline_formula <- formula <- as.formula(paste0("sales ~ ns(days_std, ", spline_dof, ")"))
  } else {
    X_spline_formula <- sales ~ 1
  }
  if (n_trig_terms >= 1) {
    for (k in seq(n_trig_terms)) {
      new_terms <- bquote(. ~ . + cos(.(k) * yday_radians) + sin(.(k) * yday_radians))
      X_spline_formula <- update(X_spline_formula, new_terms)
    }
  }
  X_spline_formula
}

get_dat_upc <- function(dat_subset, u) {
  # Get filtered data
  dat_subset_upc <- filter(dat_subset, upc == u)
  dat_subset_upc$days_std <- scale(as.numeric(dat_subset_upc$date))
  dat_subset_upc
}

get_Xy <- function(dat_subset, u, spline_dof, n_trig_terms) {
  # Get outcome vector (y) and design matrix (X)
  my_formula <- get_formula(spline_dof, n_trig_terms)
  dat_subset_upc <- get_dat_upc(dat_subset, u)
  X <- model.matrix(my_formula, data=dat_subset_upc)
  y <- model.response(model.frame(my_formula, data=dat_subset_upc))
  result_list <- list(X=X, y=y)
  result_list
}

get_upc_bic <- function(dat_subset, u, X_spline_dof, X_n_trig_terms) {
  # Get BIC for NB model with UPC=u
  Xy_list <- get_Xy(dat_subset, u, X_spline_dof, X_n_trig_terms)
  X <- Xy_list$X
  y <- Xy_list$y
  mod_nb <- glmmTMB(y ~ 0 + X, family = nbinom2)  # OR: glm.nb(y ~ 0 + X)
  BIC(mod_nb)
}

get_avg_bic <- function(dat_subset, dofs) {
  # Get average BIC across UPCs
  try_get_upc_bic <- function(u) {
    tryCatch(
      {get_upc_bic(dat_subset, u, dofs[1], dofs[2])},
      error = function(e) NA
    )
  }
  bics <- sapply(upcs, function(u) try_get_upc_bic(u))
  mean(bics, na.rm=TRUE)
}

get_next_model <- function(dat_subset, dofs) {
  # Return next model based on BIC
  p <- length(dofs)
  dofs_mat <- matrix(rep(dofs, p), nrow=p) + diag(p)
  bics <- apply(dofs_mat, 2, function(d) get_avg_bic(dat_subset, d))
  best_idx <- which.min(bics)
  best_bic <- bics[best_idx]
  best_dofs <- dofs_mat[,best_idx]
  list(
    bics=bics,
    dofs_mat=dofs_mat,
    best_bic=best_bic,
    best_dofs=best_dofs
  )
}
dat_non_holdout <- dat %>%
  filter(!in_holdout) %>%
  mutate(days_std = scale(date))
last_date <- max(dat_non_holdout$date)
first_prediction_date <- last_date - (51 * 7)  # Use 52 weeks for prediction evaluation
dat_pre_prediction <- filter(dat, date < first_prediction_date)

# Use these functions to choose the best model for X
dofs_curr <- c(0, 0)
bic_curr <- Inf
bic_delta <- -1
while (bic_delta < 0) {
  dofs_prev <- dofs_curr
  bic_prev <- bic_curr
  r <- get_next_model(dat_pre_prediction, dofs_curr)
  dofs_curr <- r$best_dofs
  bic_curr <- r$best_bic
  bic_delta <- bic_curr - bic_prev
  print(bic_curr)
}

X_spline_dof <- dofs_prev[1]
X_n_trig_terms <- dofs_prev[2]
Z_spline_dof <- 0  # We'll keep the formula for Z simple to avoid overfitting
Z_n_trig_terms <- 1
X_formula <- get_formula(X_spline_dof, X_n_trig_terms)
Z_formula <- get_formula(Z_spline_dof, Z_n_trig_terms)

# We'll use these to make some plots
date_xaxt_locations <- as.Date(c("1989-10-01", "1990-01-01", "1990-04-01", "1990-07-01", "1990-10-01", "1991-01-01", "1991-04-01", "1991-07-01", "1991-10-01", "1992-01-01", "1992-04-01", "1992-07-01", "1992-10-01", "1993-01-01", "1993-04-01", "1993-07-01", "1993-10-01"))
date_xaxt_label_locations <- as.Date(c("1990-01-01", "1990-07-01", "1991-01-01", "1991-07-01", "1992-01-01", "1992-07-01", "1993-01-01", "1993-07-01"))
get_month_year <- function(d) {
  d_year <- as.character(year(d))
  d_year_2 <- substr(d_year, 3, 4)
  d_month <- as.character(month(d))
  paste0(d_month , "/", d_year_2)
}
date_xaxt_label_names <- sapply(date_xaxt_label_locations, get_month_year)

# Is there evidence of underdispersion?
min_phis <- c()
avg_phis <- c()
full_reduced_pvals <- c()
epl_pvals <- c()
for (u in upcs) {
  print(u)
  dat_non_holdout_upc <- filter(dat_non_holdout, upc==u)
  X <- model.matrix(X_formula, dat_non_holdout_upc)
  y <- model.response(model.frame(X_formula, data=dat_non_holdout_upc))
  Z <- model.matrix(Z_formula, dat_non_holdout_upc)
  mod_qp <- glm(y ~ 0 + X, family=quasipoisson())
  beta_start <- mod_qp$coefficients
  alpha_start <- rep(0, ncol(Z))
  names(alpha_start) <- colnames(Z)
  alpha_start[1] <- log(summary(mod_qp)$dispersion)
  mod_epl <- epl(y, X, Z, betastart = beta_start, alphastart = alpha_start,
                 tol=1e-4, max_iter=100, stephalving_maxiter=100, verbose=TRUE)
  alpha_trig <- mod_epl$alpha[-1]
  cov_alpha_trig <- mod_epl$cov_alpha[-1, -1]
  chi2_epl <- c(alpha_trig %*% solve(cov_alpha_trig, alpha_trig))
  epl_pval <- pchisq(chi2_epl, df=2, lower.tail=FALSE)
  epl_pvals <- c(epl_pvals, epl_pval)
  min_phis <- c(min_phis, min(mod_epl$phi))
  avg_phis <- c(avg_phis, mean(mod_epl$phi))
  
  beta_start_dln <- rep(0, ncol(X))
  beta_start_dln[1] <- log(mean(y))
  alpha_start_dln <- rep(0, ncol(Z))
  alpha_start_dln[1] <- log(sd(y) / mean(y))
  mod_dln_em <- dln(
    y, X, Z, betastart=beta_start_dln, alphastart=alpha_start_dln, method="EM",
    # pred_interval_method="None", skip_inference=TRUE,
    pred_interval_method="Asymp. Bayes", skip_inference=FALSE, Xnew=X, Znew=Z,
    max_iter=100, stephalving_maxiter=100, 
    tol=1e-4, verbose=TRUE)
  dev_full <- mod_dln_em$dev
  
  Z_reduced <- matrix(1, nrow=nrow(X), ncol=1)
  mod_dln_em_reduced <- dln(
    y, X, Z_reduced, betastart=beta_start_dln, alphastart=alpha_start_dln[1], method="EM",
    pred_interval_method="None", max_iter=100, stephalving_maxiter=100, 
    tol=1e-4, skip_inference=TRUE, verbose=TRUE)
  dev_reduced <- mod_dln_em_reduced$dev
  dev_diff <- dev_reduced - dev_full
  dev_diff_dof <- ncol(Z) - 1
  dev_diff_pval <- pchisq(dev_diff, dev_diff_dof, lower.tail=FALSE)
  full_reduced_pvals <- c(full_reduced_pvals, dev_diff_pval)
  
  # Plot time series
  pdf(paste0("./figures/upc/", u, "-sales.pdf"), width=5, height=3)
  par(mai=c(0.9, 0.9, 0.1, 0), mgp=c(2.3, 1, 0), cex=1.35)
  plot(NULL, type="n",
       xlim=c(min(dat_non_holdout$date), max(dat_non_holdout$date)),
       ylim=c(0, max(mod_dln_em$pred_upper_bounds)*1.5),
       # main=paste("UPC:", u),
       xlab="Date", ylab="Sales", xaxt="n")
  axis(1, at=date_xaxt_locations, labels=rep("", length(date_xaxt_locations)))
  axis(1, at=date_xaxt_label_locations, labels=date_xaxt_label_names)
  polygon(
    c(dat_non_holdout_upc$date, rev(dat_non_holdout_upc$date)),
    c(mod_dln_em$pred_lower_bounds, rev(mod_dln_em$pred_upper_bounds)),
    density=200,
    col="lightgray"
  )
  lines(dat_non_holdout_upc$date, dat_non_holdout_upc$sales)
  legend("topleft", bg="white", col=c("lightgray", "black"), lty=c(NA, 1), pch=c(15, NA), legend=c("95% PI", "Observed"))
  dev.off()
  
  pdf(paste0("./figures/upc/", u, "-phi.pdf"), width=5, height=3)
  par(mai=c(0.9, 0.9, 0.1, 0), mgp=c(2.3, 1, 0), cex=1.35)
  plot(NULL, type="n",
       xlim=c(min(dat_non_holdout$date), max(dat_non_holdout$date)),
       ylim=c(min(mod_epl$log_phi_lower_bounds), max(mod_epl$log_phi_upper_bounds)),
       # main=paste("UPC:", u),
       xlab="Date", ylab=expression(paste("log(", phi[t], ")")), xaxt="n")
  polygon(
    c(dat_non_holdout_upc$date, rev(dat_non_holdout_upc$date)),
    c(mod_epl$log_phi_lower_bounds, rev(mod_epl$log_phi_upper_bounds)),
    density=200,
    col="lightgray"
  )
  lines(dat_non_holdout_upc$date, mod_epl$log_phi)
  abline(h=0)
  axis(1, at=date_xaxt_locations, labels=rep("", length(date_xaxt_locations)))
  axis(1, at=date_xaxt_label_locations, labels=date_xaxt_label_names)
  legend("topleft", bg="white", col=c("lightgray", "black"), lty=c(NA, 1), pch=c(15, NA), legend=c("95% CI", "Estimate"))
  dev.off()
}

# How many exhibit underdispersion during at least one time point?
sum(min_phis < 1)  # 2

# We reject constant dispersion for most products
sum(full_reduced_pvals < 0.05)  # 243
hist(full_reduced_pvals)

# Similar results for EPL model
sum(epl_pvals < 0.05)
hist(epl_pvals)

# Compute week-ahead predictions for 1993
# The varying-dispersion NB model is commented out because it doesn't always fit
unique_dates <- unique(dat_non_holdout$date)
prediction_dates <- unique_dates[unique_dates >= first_prediction_date]
num_prediction_dates <- length(prediction_dates)
analyze_upc_prediction_date <- function(u, prediction_date, alpha_error=0.05, verbose=FALSE) {
  selected_upc <- u
  dat_upc <- dat_non_holdout %>% filter(
    upc == selected_upc,
    date <= prediction_date)
  is_train <- dat_upc$date < prediction_date
  is_test <- !is_train
  if (sum(is_test) != 1) stop("Test set should have one observation")
  dat_upc_train <- dat_upc[is_train,]
  dat_upc_test <- dat_upc[is_test,]
  
  lm_X <- lm(X_formula, data=dat_upc_train)
  y_train <- model.response(model.frame(X_formula, data=dat_upc_train))
  y_test <- model.response(model.frame(X_formula, data=dat_upc_test))
  X_train <- model.matrix(lm_X, data=dat_upc_train)
  X_test <- model.matrix(lm_X, data=dat_upc_test)
  
  lm_Z <- lm(Z_formula, data=dat_upc_train)
  Z_train <- model.matrix(lm_Z, data=dat_upc_train)
  Z_test <- model.matrix(lm_Z, data=dat_upc_test)
  
  qp_start_time <- Sys.time()
  mod_qp <- glm(X_formula, data=dat_upc_train, family=quasipoisson())
  mod_qp_pred <- predict(mod_qp, newdata=dat_upc_test, type="response")
  qp_end_time <- Sys.time()
  qp_elapsed_time <- qp_end_time - qp_start_time
  
  nb_start_time <- Sys.time()
  mod_nb <- suppressWarnings(glmmTMB(y_train ~ 0 + X_train, family = nbinom2))
  params_nb <- mod_nb$fit$par
  dim_params <- length(params_nb)
  mod_nb_pred <- exp(c(X_test %*% params_nb[-dim_params]))
  cov_nb <- chol2inv(chol(mod_nb$obj$he()))
  n_samples <- 1000
  sampled_params_nb <- mvrnorm(n_samples, params_nb, cov_nb)
  sampled_means_nb <- exp(c(sampled_params_nb[,-dim_params] %*% c(X_test)))
  sampled_phis_nb <- exp(sampled_params_nb[, dim_params])
  sampled_vars_nb <- sampled_means_nb * (1 + sampled_means_nb/sampled_phis_nb)
  sampled_sizes_nb <- sampled_means_nb^2 / (sampled_vars_nb - sampled_means_nb) 
  sampled_outcomes_nb <- rnbinom(n_samples, size=sampled_sizes_nb, mu=sampled_means_nb)
  nb_lower_bound <- quantile(sampled_outcomes_nb, probs=alpha_error/2, na.rm=TRUE)
  nb_upper_bound <- quantile(sampled_outcomes_nb, probs=1-alpha_error/2, na.rm=TRUE)
  nb_end_time <- Sys.time()
  nb_elapsed_time <- nb_end_time - nb_start_time
  
  # nb_vardisp_start_time <- Sys.time()
  # mod_nb_vardisp <- suppressWarnings(glmmTMB(y_train ~ 0 + X_train, dispformula=~0+Z_train, family=nbinom2))
  # summary_nb_vardisp <- summary(mod_nb_vardisp)
  # coef_nb_vardisp <- summary_nb_vardisp$coefficients$cond[,"Estimate"]
  # dim_mean_model <- ncol(X)
  # dim_var_model <- ncol(Z)
  # mod_nb_vardisp_pred <- exp(c(X_test %*% coef_nb_vardisp))
  # zero_mat <- matrix(0, nrow=dim_mean_model, ncol=dim_var_model)
  # vcov_nb_vardisp <- vcov(mod_nb_vardisp)
  # cov_nb_vardisp <- rbind(
  #   cbind(vcov_nb_vardisp$cond, zero_mat),
  #   cbind(t(zero_mat), vcov_nb_vardisp$disp))
  # # cov_nb_vardisp <- chol2inv(chol(mod_nb_vardisp$obj$he()))  # Estimated cov via Hessian
  # sampled_params_nb_vardisp <- mvrnorm(n_samples, mod_nb_vardisp$fit$parfull, cov_nb_vardisp)
  # sampled_means_nb_vardisp <- exp(c(sampled_params_nb_vardisp[,seq(dim_mean_model)] %*% c(X_test)))
  # sampled_phis_nb_vardisp <- exp(sampled_params_nb_vardisp[,dim_mean_model + seq(dim_var_model)] %*% c(Z_test))
  # sampled_vars_nb_vardisp <- sampled_means_nb_vardisp * (1 + sampled_means_nb_vardisp/sampled_phis_nb_vardisp)
  # sampled_sizes_nb_vardisp <- sampled_means_nb_vardisp^2 / (sampled_vars_nb_vardisp - sampled_means_nb_vardisp) 
  # sampled_outcomes_nb_vardisp <- rnbinom(n_samples, size=sampled_sizes_nb_vardisp, mu=sampled_means_nb_vardisp)
  # nb_vardisp_lower_bound <- quantile(sampled_outcomes_nb_vardisp, probs=alpha_error/2, na.rm=TRUE)
  # nb_vardisp_upper_bound <- quantile(sampled_outcomes_nb_vardisp, probs=1-alpha_error/2, na.rm=TRUE)
  # nb_vardisp_end_time <- Sys.time()
  # nb_vardisp_elapsed_time <- nb_vardisp_end_time - nb_vardisp_start_time
  
  beta_start <- mod_qp$coefficients
  alpha_start <- rep(0, ncol(Z_train))
  names(alpha_start) <- colnames(Z_train)
  alpha_start[1] <- log(summary(mod_qp)$dispersion)
  max_iter <- 100
  tol <- 1e-4
  stephalving_max <- 100
  mod_epl <- epl(y_train, X_train, Z_train, betastart = beta_start, alphastart = alpha_start, Xnew=X_test, Znew=Z_test,
                 tol=tol, max_iter=max_iter, stephalving_maxiter=stephalving_max, verbose=verbose)
  
  beta_start_dln <- rep(0, ncol(X_train))
  beta_start_dln[1] <- log(mean(y_train))
  alpha_start_dln <- rep(0, ncol(Z_train))
  alpha_start_dln[1] <- log(sd(y_train) / mean(y_train))
  mod_dln_em <- dln(y_train, X_train, Z_train, betastart = beta_start_dln, alphastart = alpha_start_dln, method="EM",
                    pred_interval_method="Asymp. Bayes", Xnew=X_test, Znew=Z_test,
                    tol=tol, max_iter=max_iter, stephalving_maxiter=stephalving_max, verbose=verbose)
  
  results <- c(
    qp_error=unname(mod_qp_pred - y_test),
    nb_error=unname(mod_nb_pred - y_test),
    # nb_vardisp_error=unname(mod_nb_vardisp_pred - y_test),
    epl_error=unname(mod_epl$fitted_values - y_test),
    dln_error=unname(mod_dln_em$fitted_values - y_test),
    nb_covered=unname((nb_lower_bound <= y_test) & (y_test <= nb_upper_bound)),
    # nb_vardisp_covered=unname((nb_vardisp_lower_bound <= y_test) & (y_test <= nb_vardisp_upper_bound)),
    dln_covered=unname((mod_dln_em$pred_lower_bounds <= y_test) & (y_test <= mod_dln_em$pred_upper_bounds)),
    nb_length=unname(nb_upper_bound - nb_lower_bound),
    # nb_vardisp_length=unname(nb_vardisp_upper_bound - nb_vardisp_lower_bound),
    dln_length=unname(mod_dln_em$pred_upper_bounds - mod_dln_em$pred_lower_bounds),
    qp_elapsed_time=unname(qp_elapsed_time),
    nb_elapsed_time=unname(nb_elapsed_time),
    # nb_vardisp_elapsed_time=unname(nb_vardisp_elapsed_time),
    epl_elapsed_time=unname(mod_epl$elapsed_time),
    dln_elapsed_time=unname(mod_dln_em$elapsed_time)
  )
  results
}
method_names <- c("QP", "NB", "EPL", "DLN")  # c("QP", "NB", "NB-Var", "EPL", "DLN")
num_methods <- length(method_names)
results_dimnames <- list(
  method=method_names,
  upc=as.character(upcs),
  prediction_date=as.character(prediction_dates))
all_errors <- array(0, dim=c(num_methods, num_upcs, num_prediction_dates), dimnames=results_dimnames)
all_elapsed_time <- all_errors
pred_method_names <- c("NB", "DLN")  # c("NB", "NB-Var", "DLN")
num_pred_method_names <- length(pred_method_names)
pred_results_dimnames <- list(
  method=pred_method_names,
  upc=as.character(upcs),
  prediction_date=as.character(prediction_dates))
all_covered <- array(FALSE, dim=c(num_pred_method_names, num_upcs, num_prediction_dates), dimnames=pred_results_dimnames)
all_lengths <- array(0, dim=c(num_pred_method_names, num_upcs, num_prediction_dates), dimnames=pred_results_dimnames) 
alpha_error <- 0.05
for (u_idx in seq(num_upcs)) {
  u <- upcs[u_idx]
  print(as.character(u))
  gc()
  for (d_idx in seq(num_prediction_dates)) {
    d <- prediction_dates[d_idx]
    results_u_d <- analyze_upc_prediction_date(u, d, alpha_error=alpha_error)
    all_errors["QP", u_idx, d_idx] <- results_u_d["qp_error"]
    all_errors["NB", u_idx, d_idx] <- results_u_d["nb_error"]
    # all_errors["NB-Var", u_idx, d_idx] <- results_u_d["nb_vardisp_error"]
    all_errors["EPL", u_idx, d_idx] <- results_u_d["epl_error"]
    all_errors["DLN", u_idx, d_idx] <- results_u_d["dln_error"]
    all_elapsed_time["QP", u_idx, d_idx] <- results_u_d["qp_elapsed_time"]
    all_elapsed_time["NB", u_idx, d_idx] <- results_u_d["nb_elapsed_time"]
    # all_elapsed_time["NB-Var", u_idx, d_idx] <- results_u_d["nb_vardisp_elapsed_time"]
    all_elapsed_time["EPL", u_idx, d_idx] <- results_u_d["epl_elapsed_time"]
    all_elapsed_time["DLN", u_idx, d_idx] <- results_u_d["dln_elapsed_time"]
    all_covered["NB", u_idx, d_idx] <- results_u_d["nb_covered"]
    # all_covered["NB-Var", u_idx, d_idx] <- results_u_d["nb_vardisp_covered"]
    all_covered["DLN", u_idx, d_idx] <- results_u_d["dln_covered"]
    all_lengths["NB", u_idx, d_idx] <- results_u_d["nb_length"]
    # all_lengths["NB-Var", u_idx, d_idx] <- results_u_d["nb_vardisp_length"]
    all_lengths["DLN", u_idx, d_idx] <- results_u_d["dln_length"]
  }
}

# Analyze the results
paste_ci <- function(x, l, u, num_digits, name="Metric") {
  formatter <- function(y) format(round(y, digits=num_digits), big.mark=",", trim=TRUE)
  x_formatted <- formatter(x)
  ci_formatted <- paste0(" (", formatter(l), ", ", formatter(u), ")")
  result <- cbind(x_formatted, ci_formatted)
  colnames(result) <- c(name, "(95% CI)")
  result
}

summarize_metric <- function(x, num_digits=2, alpha_error=0.05, f=identity, name="Metric") {
  ndims <- length(dim(x))
  if (ndims == 2) {
    margin1 <- 1
    margin2 <- 1
  } else if (ndims == 3) {
    margin1 <- c(1, 2)
    margin2 <- 1
  } else {
    stop("x must have 2 or 3 dimensions")
  }
  x_means <- f(apply(x, MARGIN=margin1, FUN=mean))
  m <- apply(x_means, MARGIN=margin2, FUN=mean)
  sample_sizes <- apply(x_means, MARGIN=margin2, FUN=length)
  se <- apply(x_means, MARGIN=margin2, FUN=sd) / sqrt(sample_sizes)
  z_star <- qnorm(1 - alpha_error/2)
  m_lower <- m - z_star * se
  m_upper <- m + z_star * se
  result <- paste_ci(m, m_lower, m_upper, num_digits, name)
  result
}

# Error table
bias_chr <- summarize_metric(all_errors, num_digits=0, alpha=alpha_error, name="Bias")
rmse_chr <- summarize_metric(all_errors^2, num_digits=0, alpha=alpha_error, f=sqrt, name="rMSE")
elapsed_time_chr <- summarize_metric(all_elapsed_time, num_digits=4, alpha=alpha_error, f=sqrt, name="Time [sec]")
bias_rmse_df <- cbind(
  Method=method_names,
  bias_chr,
  rmse_chr,
  elapsed_time_chr)
print(bias_rmse_df)
xtable(
  bias_rmse_df,
  label="tab:errors",
  align="llrlrlrl",
  caption="Bias and rMSE for the week-ahead predictions on the DFF data set. The values in parentheses indicate 95\\% Monte Carlo confidence intervals.") %>%
  print.xtable(include.rownames = FALSE)

# Coverage
coverage_chr <- summarize_metric(all_covered, num_digits=3, alpha=alpha_error, name="Coverage")
length_chr <- summarize_metric(all_lengths, num_digits=0, alpha=alpha_error, name="Avg Length")
ci_df <- cbind(
  Method=pred_method_names,
  coverage_chr,
  length_chr)
print(ci_df)
xtable(
  ci_df,
  label="tab:pred-coverage",
  align="llrlrl",
  caption="Coverage and average lengths for 95\\% prediction intervals (95\\% CI)") %>%
  print.xtable(include.rownames = FALSE)

# Difference between models
summarize_metric_diffs <- function(x, num_digits=0, alpha=0.05, f=identity, metric_name="Metric") {
  x_means <- f(apply(x, MARGIN=c(1, 2), FUN=mean))
  dim1_names <- dimnames(x_means)[[1]]
  dim1 <- dim(x_means)[1]
  dim2 <- dim(x_means)[2]
  x_means1 <- array(0, dim=c(dim1, dim1, dim2))
  x_means2 <- x_means1
  for (i in seq(dim1)) {
    x_means1[i,,] <- x_means
    x_means2[,i,] <- x_means
  }
  diffs <- x_means2 - x_means1
  diffs_mean <- apply(diffs, MARGIN=c(1,2), FUN=mean)
  diffs_mean_se <- apply(diffs, MARGIN=c(1,2), FUN=sd) / sqrt(dim2)
  z_star <- qnorm(1 - alpha_error/2)
  diffs_mean_lower <- diffs_mean - z_star * diffs_mean_se
  diffs_mean_upper <- diffs_mean + z_star * diffs_mean_se
  
  formatter <- function(y) format(round(y, digits=num_digits), big.mark=",", trim=TRUE)
  diffs_mean_formatted <- formatter(diffs_mean)
  colnames(diffs_mean_formatted) <- dim1_names
  diffs_mean_formatted[as.logical(diag(dim1))] <- "-"
  diffs_mean_lower_formatted <- formatter(diffs_mean_lower)
  diffs_mean_upper_formatted <- formatter(diffs_mean_upper)
  cis_formatted <- matrix(
    paste0("(", diffs_mean_lower_formatted, ", ", diffs_mean_upper_formatted, ")"),
    nrow=nrow(diffs_mean))
  cis_formatted[as.logical(diag(dim1))] <- ""
  colnames(cis_formatted) <- rep("(95% CI)", dim1)
  result <- cbind(diffs_mean_formatted, cis_formatted)
  result_idx <- c(0, dim1) + rep(seq(dim1), each=2)
  result <- result[, result_idx]
  result <- cbind(
    Metric=metric_name,
    Method=dim1_names,
    result)
  result
}

# Error diffs
diff_bias_chr <- summarize_metric_diffs(all_errors, num_digits=0, alpha=alpha, metric_name="Bias")
diff_rmse_chr <- summarize_metric_diffs(all_errors^2, num_digits=0, alpha=alpha, f=sqrt, metric_name="rMSE")
diff_df <- rbind(diff_bias_chr, diff_rmse_chr)
diff_df
xtable(
  diff_df,
  label="tab:error-diff",
  align="lllrlrlrlrl",
  caption="Pairwise differences between methods in terms of bias and rMSE (95\\% CIs).") %>%
  print.xtable(include.rownames = FALSE)

# Coverage diffs
diff_coverage_chr <- summarize_metric_diffs(all_covered, num_digits=3, alpha=alpha, metric_name="Coverage")
diff_length_chr <- summarize_metric_diffs(all_lengths, num_digits=3, alpha=alpha, metric_name="Length")
diff_ci_df <- rbind(diff_coverage_chr, diff_length_chr)
diff_ci_df
xtable(diff_ci_df) %>%
  print.xtable(include.rownames = FALSE)

log2_avg_phis <- log2(avg_phis)
pdf("./figures/phis.pdf", width=5, height=3)
par(mai=c(0.9, 0.9, 0.1, 0), mgp=c(2.3, 1, 0), cex=1.35)
hist(
  log2_avg_phis,
  # main=expression(paste("Histogram of ", hat(phi))),
  main="",
  xlab=expression(hat(phi)),
  ylab="# UPCs",
  ylim=c(0, 40),
  xaxt="n",
  yaxt="n")
phis_xaxt_locations <- seq(0, ceiling(max(log2_avg_phis)), 3)
phis_xaxt_labels <- format(2^phis_xaxt_locations, big.mark=",", trim=TRUE)
axis(1, at=phis_xaxt_locations, labels=phis_xaxt_labels)
phis_yaxt_locations <- c(0, 10, 20, 30, 40)
axis(2, at=phis_yaxt_locations)
dev.off()

# Look at plots for these UPCs
min_phi_idx <- which.min(min_phis)
u_min_phi <- upcs[min_phi_idx]
print(u_min_phi)
max_phi_idx <- which.max(avg_phis)
u_max_phi <- upcs[max_phi_idx]
print(u_max_phi)

# Plot example problematic UPC/date for varying dispersion NB model
selected_upc <- "4400000693"
prediction_date <- as.Date("1992-10-22")
dat_upc <- dat_non_holdout %>% filter(
  upc == selected_upc,
  date <= prediction_date)
is_train <- dat_upc$date < prediction_date
is_test <- !is_train
if (sum(is_test) != 1) stop("Test set should have one observation")
dat_upc_train <- dat_upc[is_train,]
dat_upc_test <- dat_upc[is_test,]

lm_X <- lm(X_formula, data=dat_upc_train)
y_train <- model.response(model.frame(X_formula, data=dat_upc_train))
y_test <- model.response(model.frame(X_formula, data=dat_upc_test))
X_train <- model.matrix(lm_X, data=dat_upc_train)
X_test <- model.matrix(lm_X, data=dat_upc_test)

lm_Z <- lm(Z_formula, data=dat_upc_train)
Z_train <- model.matrix(lm_Z, data=dat_upc_train)
Z_test <- model.matrix(lm_Z, data=dat_upc_test)

mod_nb_vardisp <- glmmTMB(y_train ~ 0 + X_train, dispformula=~0+Z_train, family=nbinom2)
summary_nb_vardisp <- summary(mod_nb_vardisp)
beta_nb_vardisp <- summary_nb_vardisp$coefficients$cond[,"Estimate"]
alpha_nb_vardisp <- summary_nb_vardisp$coefficients$disp[,"Estimate"]
dim_mean_model <- ncol(X)
dim_var_model <- ncol(Z)
mus_nb_vardisp <- c(exp(X_train %*% beta_nb_vardisp))
phis_nb_vardisp <- c(exp(Z_train %*% alpha_nb_vardisp))
vars_nb_vardisp <- mus_nb_vardisp * (1 + mus_nb_vardisp/phis_nb_vardisp)
sizes_nb_vardisp <- mus_nb_vardisp^2 / (vars_nb_vardisp - mus_nb_vardisp) 
nb_vardisp_lower_bound <- qnbinom(0.025, size=sizes_nb_vardisp, mu=mus_nb_vardisp)
nb_vardisp_upper_bound <- qnbinom(0.975, size=sizes_nb_vardisp, mu=mus_nb_vardisp)

pdf(paste0("./figures/upc-error/", selected_upc, "-fit.pdf"), width=5, height=3)
par(mai=c(0.9, 0.9, 0.1, 0), mgp=c(2.3, 1, 0), cex=1.35)
plot(NULL, type="n",
     xlim=c(min(dat_upc_train$date), max(dat_upc_train$date)),
     ylim=c(0, 350),
     # main="Sales vs. Model Fit",
     xlab="Date", ylab="Sales",
     xaxt="n"
)
polygon(
  c(dat_upc_train$date, rev(dat_upc_train$date)),
  c(nb_vardisp_lower_bound, nb_vardisp_upper_bound),
  density=200,
  col="lightgray"
)
lines(dat_upc$date, dat_upc$sales, type="l", col="black")
legend("bottomleft", bg="white", col=c("lightgray", "black"), lty=c(NA, 1), pch=c(15, NA), legend=c("95% Confidence", "Observed"))
axis(1, at=date_xaxt_locations, labels=rep("", length(date_xaxt_locations)))
axis(1, at=date_xaxt_label_locations, labels=date_xaxt_label_names)
dev.off()

pdf(paste0("./figures/upc-error/", selected_upc, "-dispersion.pdf"), width=5, height=3)
par(mai=c(0.9, 0.9, 0.1, 0), mgp=c(2.3, 1, 0), cex=1.35)
disp_nb_vardisp <- vars_nb_vardisp / mus_nb_vardisp
plot(NULL, type="n", xaxt="n",
     xlim=c(min(dat_upc_train$date), max(dat_upc_train$date)),
     ylim=c(0.8, 22), xlab="Date",
     # main="Estimated Dispersion over Time",
     ylab="Estimated Dispersion", log="y")
abline(h=1, col="lightgray", lty=2)
lines(dat_upc_train$date, disp_nb_vardisp)
axis(1, at=date_xaxt_locations, labels=rep("", length(date_xaxt_locations)))
axis(1, at=date_xaxt_label_locations, labels=date_xaxt_label_names)
dev.off()

pdf(paste0("./figures/upc-error/", selected_upc, "-log-size.pdf"), width=5, height=3)
par(mai=c(0.9, 0.9, 0.1, 0), mgp=c(2.3, 1, 0), cex=1.35)
plot(NULL, type="n", xaxt="n",
     xlim=c(min(dat_upc_train$date), max(dat_upc_train$date)),
     ylim=c(min(sizes_nb_vardisp), max(sizes_nb_vardisp)), xlab="Date",
     # main="Estimated Size Parameter over Time",
     ylab="Estimated Size", log="y")
abline(h=1, col="lightgray", lty=2)
lines(dat_upc_train$date, sizes_nb_vardisp)
axis(1, at=date_xaxt_locations, labels=rep("", length(date_xaxt_locations)))
axis(1, at=date_xaxt_label_locations, labels=date_xaxt_label_names)
dev.off()
