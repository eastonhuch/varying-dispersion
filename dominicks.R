# Dominick's Analysis
require(tidyverse)
require(xtable)
require(scales)
require(splines)
require(lubridate)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

# Read in data
dat <- read.csv("./dominicks.csv")
unique(dat$upc)
dat$date <- as.Date(dat$date)
dat$days_std <- scale(as.numeric(dat$date))
dat$yday <- yday(dat$date)
dat$yday_radians <- dat$yday / 366 * 2 * pi

# Look at some example time series
upcs <- sort(unique(dat$upc))
for (u in upcs) {
  sales_example <- filter(dat, upc == u)
  plot(sales_example$date, sales_example$sales, type="l",
       main=u, log="y", xlab="Date", ylab="Sales")
  Sys.sleep(0.2)
}

# Compute week-ahead predictions for 1993
first_prediction_date <- as.Date("1993-01-01")
unique_dates <- unique(dat$date)
prediction_dates <- unique_dates[unique_dates >= first_prediction_date]
num_prediction_dates <- length(prediction_dates)
num_upcs <- length(upcs)
analyze_upc_prediction_date <- function(upc, prediction_date) {
  selected_upc <- upc
  dat_upc <- dat %>% filter(
    upc == selected_upc,
    date <= prediction_date
  )
  is_train <- dat_upc$date < prediction_date
  is_test <- !is_train
  if (sum(is_test) != 1) stop("Test set should have one observation")
  dat_upc_train <- dat_upc[is_train,]
  dat_upc_test <- dat_upc[is_test,]
  
  X_formula <- sales ~
    ns(days_std, 10) +
    cos(yday_radians) +
    sin(yday_radians) +
    cos(2 * yday_radians) +
    sin(2 * yday_radians)
  lm_X <- lm(X_formula, data=dat_upc_train)
  y_train <- pull(dat_upc_train, "sales")
  y_test <- pull(dat_upc_test, "sales")
  X_train <- model.matrix(lm_X, data=dat_upc_train)
  X_test <- model.matrix(lm_X, data=dat_upc_test)
  
  Z_formula <- sales ~
    ns(days_std, 3) +
    cos(yday_radians) +
    sin(yday_radians)
  lm_Z <- lm(Z_formula, data=dat_upc_train)
  Z_train <- model.matrix(lm_Z, data=dat_upc_train)
  Z_test <- model.matrix(lm_Z, data=dat_upc_test)
  
  quasi_start_time <- Sys.time()
  mod_quasipois <- glm(y_train ~ 0 + X_train, family=quasipoisson())
  quasi_end_time <- Sys.time()
  quasi_elapsed_time <- quasi_end_time - quasi_start_time
  
  beta_start <- mod_quasipois$coefficients
  alpha_start <- rep(0, ncol(Z_train))
  names(alpha_start) <- colnames(Z_train)
  alpha_start[1] <- log(summary(mod_quasipois)$dispersion)
  max_iter <- 100
  tol <- 1e-4
  stephalving_max <- 50
  mod_epl <- epl(y_train, X_train, Z_train, betastart = beta_start, alphastart = alpha_start, Xnew=X_test, Znew=Z_test,
                 tol=tol, max_iter=max_iter, stephalving_maxiter=stephalving_max)
  
  beta_start_dln <- rep(0, ncol(X_train))
  beta_start_dln[1] <- log(mean(y_train))
  alpha_start_dln <- rep(0, ncol(Z_train))
  alpha_start_dln[1] <- log(sd(y_train) / mean(y_train))
  mod_dln_em <- dln(y_train, X_train, Z_train, betastart = beta_start_dln, alphastart = alpha_start_dln, method="EM",
                    pred_interval_method="Asymp. Bayes", Xnew=X_test, Znew=Z_test,
                    max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=FALSE)
  
  results <- c(
    epl_error = y_test - mod_epl$fitted_values,
    dln_error = y_test - mod_dln_em$fitted_values,
    dln_covered = (mod_dln_em$pred_lower_bounds <= y_test) & (y_test <= mod_dln_em$pred_upper_bounds)
  )
  
  results
}
results_dimnames <- list(
  upc=as.character(upcs),
  prediction_date=as.character(prediction_dates))
epl_errors <- matrix(0, nrow=num_upcs, ncol=num_prediction_dates, dimnames=results_dimnames)
dln_errors <- matrix(0, nrow=num_upcs, ncol=num_prediction_dates, dimnames=results_dimnames)
dln_covered <- matrix(FALSE, nrow=num_upcs, ncol=num_prediction_dates, dimnames=results_dimnames)

for (u_idx in seq(num_upcs)) {
  u <- upcs[u_idx]
  print(as.character(u))
  for (d_idx in seq(num_prediction_dates)) {
    d <- prediction_dates[d_idx]
    results_u_d <- analyze_upc_prediction_date(u, d)
    epl_errors[u_idx, d_idx] <- results_u_d["epl_error"]
    dln_errors[u_idx, d_idx] <- results_u_d["dln_error"]
    dln_covered[u_idx, d_idx] <- results_u_d["dln_covered"]
  }
}

# Analyze the results
# Error rate and normal quantile
alpha <- 0.05

paste_ci <- function(x, l, u, num_digits) {
  formatter <- function(y) sprintf(paste0("%.", num_digits, "f"), y)
  paste0(
    formatter(x),
    " (",
    formatter(l),
    ", ",
    formatter(u),
    ")"
  ) 
}

summarize_metric <- function(x, y=NULL, num_digits=2, alpha=0.05, f=identity) {
  x_means <- f(colMeans(x))
  if (is.null(y)) {
    diff_means <- x_means
  } else {
    y_means <- f(colMeans(y))
    diff_means <- x_means - y_means
  }
  diff_mean <- mean(diff_means)
  diff_se <- sd(diff_means) / sqrt(length(diff_means))
  z_star <- qnorm(1 - alpha/2)
  diff_mean_lower <- diff_mean - z_star * diff_se
  diff_mean_upper <- diff_mean + z_star * diff_se
  result <- paste_ci(diff_mean, diff_mean_lower, diff_mean_upper, num_digits)
  result
}

# EPL Error
epl_bias_chr <- summarize_metric(epl_errors, num_digits=0, alpha=alpha)
epl_rmse_chr <- summarize_metric(epl_errors^2, num_digits=0, alpha=alpha, f=sqrt)

# DLN Error
dln_bias_chr <- summarize_metric(dln_errors, num_digits=0, alpha=alpha)
dln_rmse_chr <- summarize_metric(dln_errors^2, num_digits=0, alpha=alpha, f=sqrt)

# Difference between models
diff_bias_chr <- summarize_metric(dln_errors, epl_errors, num_digits=0, alpha=alpha)
diff_rmse_chr <- summarize_metric(dln_errors^2, epl_errors^2, num_digits=0, alpha=alpha, f=sqrt)

# Bias, rMSE table
bias_rmse_df <- data.frame(
  method=c("EPL", "DLN", "Difference"),
  bias=c(epl_bias_chr, dln_bias_chr, diff_bias_chr),
  rmse=c(epl_rmse_chr, dln_rmse_chr, diff_rmse_chr))
print(bias_rmse_df)
xtable(bias_rmse_df) %>%
  print.xtable(include.rownames = FALSE)

# Coverage
coverage_chr <- summarize_metric(dln_covered)
cat("Coverage for DLN: ", coverage_chr, "\n")
