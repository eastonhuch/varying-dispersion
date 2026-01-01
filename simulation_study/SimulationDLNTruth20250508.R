################################
# Simulation to compare mean and sd estimation
# when data is generated using DLN instead of COM-Pois
# Dec 23, 2025
# Increasing precision of COM-Pois estimation
################################

require(mpcmp) #install_github("thomas-fung/mpcmp")
require(cmp) #install_github("SuneelChatla/cmp")
require(scales)
require(MASS)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

# Generate some data

myns <- c(100L, 250L, 500L, 1000L)
for(n in myns){
  set.seed(20240730+n)
  #n <- 400L
  x1 <- rnorm(n)  # We'll keep the design fixed across simulations
  x2 <- rnorm(n)

  # NOTE: The mpcmp package can be a bit temperamental
  # You might run into errors if you make the coefficients much larger
  beta_true <- c(2.85, .025, -.025, .01) #c(3, 0.05, -0.1, 0.02) # Regression coefficients for the mean
  y <- rep(0, n)  # We'll create this later
  master_formula <- y ~ x1*x2
  X <- model.matrix(master_formula)
  muz <- c(X%*%beta_true)
  alpha_true <- c(-1.4, 0, -0.1, 0.05) #c(0.1, 0, -0.2, 0.05) # Not directly comparable across all models; these values of alpha (alpha0=1; combined with the values of beta) result in underdispersion; to create overdispersed data, use alpha0=-1 or alpha0=-0.5
  Z <- model.matrix(master_formula)
  log_nu <- c(Z %*% alpha_true)
  nu <- exp(log_nu)
  sd_truez <- nu
  var_truez <- sd_truez^2

  muy <- integrate_dln(muz, sd_truez, 1)
  ey2 <- integrate_dln(muz, sd_truez, 2)
  vary <- ey2 - (muy^2)
  sdy <- sqrt(vary)
  mu <- muy
  sd_true <- sdy

  # Control parameters
  tol <- 1e-10
  max_iter <- 100
  phi_method <- "joint"
  stephalving_max <- 50

  # Prepare to loop
  reps <- 100L
  method_names <- c("Pois", "QP", "NB", "MPCMP", "GP-1", "EPL", "DLN-N", "DLN-EM")
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
    zlatent <- rnorm(n, muz, sd=sd_truez)
    y <- floor(exp(zlatent))
    dat <- data.frame(x1, x2, y)

    # Fit Pois model
    if("Pois" %in% method_names){
      start_time <- Sys.time()
      mod_pois <- glm(master_formula, family="poisson")
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      fitting_times[i,"Pois"] <- elapsed_time
      mu_estimates[,i,"Pois"] <- mod_pois$fitted
      sd_estimates[,i,"Pois"] <- sqrt(mod_pois$fitted)
      #muhat = exp(etahat); etahat = x*betahat => cov(muhat) = xt*cov(betahat)*x
      semu <- X%*%summary(mod_pois)$cov.scaled%*%t(X)
      lowerbnd <- exp(predict(mod_pois) - qnorm(.975)*sqrt(diag(semu)))
      upperbnd <- exp(predict(mod_pois) + qnorm(.975)*sqrt(diag(semu)))
      mu_covered[,i,"Pois"] <- (lowerbnd <= mu) & (mu <= upperbnd)
      mu_interval_widths[,i,"Pois"] <- upperbnd - lowerbnd
      sd_covered[,i,"Pois"] <- (sqrt(lowerbnd) <= sd_true) & (sd_true <= sqrt(upperbnd))
      sd_interval_widths[,i,"Pois"] <- sqrt(upperbnd) - sqrt(lowerbnd)
    }

    # Fit quasi-Pois model
    if("QP" %in% method_names){
      start_time <- Sys.time()
      mod_qp <- glm(master_formula, family="quasipoisson")
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      fitting_times[i,"QP"] <- elapsed_time
      mu_estimates[,i,"QP"] <- mod_qp$fitted
      sd_estimates[,i,"QP"] <- sqrt(mod_qp$fitted)
      #muhat = exp(etahat); etahat = x*betahat => cov(muhat) = xt*cov(betahat)*x
      semu <- X%*%summary(mod_qp)$cov.scaled%*%t(X)
      lowerbnd <- exp(predict(mod_qp) - qnorm(.975)*sqrt(diag(semu)))
      upperbnd <- exp(predict(mod_qp) + qnorm(.975)*sqrt(diag(semu)))
      mu_covered[,i,"QP"] <- (lowerbnd <= mu) & (mu <= upperbnd)
      mu_interval_widths[,i,"QP"] <- upperbnd - lowerbnd
      sd_covered[,i,"QP"] <- (sqrt(lowerbnd) <= sd_true) & (sd_true <= sqrt(upperbnd))
      sd_interval_widths[,i,"QP"] <- sqrt(upperbnd) - sqrt(lowerbnd)
    }

    # Fit NB model
    if("NB" %in% method_names) try({
      start_time <- Sys.time()
      mod_nb <- glm.nb(master_formula)
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      fitting_times[i,"NB"] <- elapsed_time
      mu_estimates[,i,"NB"] <- mod_nb$fitted
      sd_estimates[,i,"NB"] <- sqrt(mod_nb$fitted)
      #muhat = exp(etahat); etahat = x*betahat => cov(muhat) = xt*cov(betahat)*x
      semu <- X%*%summary(mod_nb)$cov.scaled%*%t(X)
      lowerbnd <- exp(predict(mod_nb) - qnorm(.975)*sqrt(diag(semu)))
      upperbnd <- exp(predict(mod_nb) + qnorm(.975)*sqrt(diag(semu)))
      mu_covered[,i,"NB"] <- (lowerbnd <= mu) & (mu <= upperbnd)
      mu_interval_widths[,i,"NB"] <- upperbnd - lowerbnd
      sd_lbnd <- lowerbnd*(lowerbnd+mod_nb$theta)/mod_nb$theta
      sd_ubnd <- upperbnd*(upperbnd+mod_nb$theta)/mod_nb$theta
      sd_covered[,i,"NB"] <- (sqrt(sd_lbnd) <= sd_true) & (sd_true <= sqrt(sd_ubnd))
      sd_interval_widths[,i,"NB"] <- sqrt(sd_ubnd) - sqrt(sd_lbnd)
    })

    # Fit MPCMP model
    if ("MPCMP" %in% method_names ) try({
      start_time <- Sys.time()
      mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-6), error=function(e) e)
      if(!is.null(mod_cmp$message)){
        mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-6, beta_true=beta_true), error=function(e) e)
      }
      if(!is.null(mod_cmp$message)){
        mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-6, beta_true=beta_true), error=function(e) e)
      }
      if(!is.null(mod_cmp$message)){
        mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-6, beta_true=beta_true), error=function(e) e)
      }
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
    if ("DLN-N" %in% method_names ) try({
      start_time <- Sys.time()
      mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
      beta_start <- mod_quasipois$coefficients
      alpha_start <- rep(0, ncol(Z))
      names(alpha_start) <- colnames(Z)
      alpha_start[1] <- log(sqrt(summary(mod_quasipois)$dispersion))
      mod_dln_newton <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton",
                            max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=FALSE)
      alpha_start[1] <- mod_dln_newton$alpha[1]
      mod_dln_newton <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton",
                            max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=FALSE)
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      fitting_times[i,"DLN-N"] <- elapsed_time
      mu_estimates[,i,"DLN-N"] <- mod_dln_newton$fitted_values
      sd_estimates[,i,"DLN-N"] <- mod_dln_newton$sd_estimates
      mu_covered[,i,"DLN-N"] <- (mod_dln_newton$fitted_lower_bounds <= mu) & (mu <= mod_dln_newton$fitted_upper_bounds)
      mu_interval_widths[,i,"DLN-N"] <- mod_dln_newton$fitted_interval_widths
      sd_covered[,i,"DLN-N"] <- (mod_dln_newton$sd_lower_bounds <= sd_true) & (sd_true <= mod_dln_newton$sd_upper_bounds)
      sd_interval_widths[,i,"DLN-N"] <- mod_dln_newton$sd_interval_widths
      used_score_for_cov[i,"DLN-N"] <- mod_dln_newton$covariance_via_score
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
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=FALSE)
      alpha_start[1] <- mod_dln_em$alpha[1]
      mod_dln_em <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="EM",
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=FALSE)
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
  format_float <- label_number(0.01)
  format_float_with_se <- function(m, se) {
    paste0(
      format_float(m), " (",
      format_float(se), ")"
    )
  }
  format_percent <- label_percent(1)
  format_percent_with_se <- function(m, se) {
    paste0(
      format_percent(m), " (",
      format_percent(se), ")"
    )
  }

  # Fitting times
  avg_fitting_times <- colMeans(fitting_times)
  avg_fitting_times_se <- apply(fitting_times, 2, function(x) sd(x)/sqrt(length(x)))
  consolidated_results <- data.frame(`Avg Elapsed Time (Seconds)`=format_float_with_se(avg_fitting_times, avg_fitting_times_se), check.names=FALSE)
  boxplot(fitting_times, log="y", main="Model-fitting Time by Method", ylab="Elapsed Time (Seconds)")
  
  # NA values
  is_na <- is.na(mu_estimates) | is.na(mu_covered) | is.na(sd_estimates) | is.na(sd_covered)
  na_rep <- apply(is_na, c(2, 3), mean)
  na_avg <- apply(na_rep, 2, mean)
  na_avg_se <- apply(na_rep, 2, function(x) sd(x) / sqrt(length(x)))
  consolidated_results$`Percent Failed` <- format_percent_with_se(na_avg, na_avg_se)

  # Mean
  mu_bias_rep <- apply(mu_estimates - mu, MARGIN=c(2,3), function(x) mean(x, na.rm=TRUE))
  mu_avg_bias <- apply(mu_bias_rep, MARGIN=2, function(x) mean(x, na.rm=TRUE))
  mu_avg_bias_se <- apply(mu_bias_rep, MARGIN=2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
  mu_rmse_rep <- apply((mu_estimates - mu)^2, MARGIN=c(2,3), function(x) sqrt(mean(x, na.rm=TRUE)))
  mu_avg_rmse <- apply(mu_rmse_rep, MARGIN=2, function(x) mean(x, na.rm=TRUE))
  mu_avg_rmse_se <- apply(mu_rmse_rep, MARGIN=2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
  mu_rep_coverage <- apply(mu_covered, MARGIN=c(2,3), function(x) mean(x, na.rm=TRUE))
  mu_marginal_coverage <- apply(mu_rep_coverage, 2, function(x) mean(x, na.rm=TRUE))
  mu_marginal_coverage_se <- apply(mu_rep_coverage, 2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
  get_coverage_rMSE <- function(covered_array) {
    observation_coverage <- apply(covered_array, MARGIN=c(1,3), function(x) mean(x, na.rm=TRUE))
    coverage_rMSE <- apply((observation_coverage - nominal_coverage)^2, 2, function(x) sqrt(mean(x, na.rm=TRUE)))
    coverage_rMSE
  }
  mu_coverage_rMSE <- get_coverage_rMSE(mu_covered)
  n_bootstrap <- 1000
  mu_coverage_rMSE_bootstrap <- sapply(seq(n_bootstrap), function(i) {
    idx <- sample(seq(reps), reps, replace=TRUE)
    mu_covered_boot <- mu_covered[,idx,]
    coverage_rMSE <- get_coverage_rMSE(mu_covered_boot)
    coverage_rMSE
  })
  mu_coverage_rMSE_se <- apply(mu_coverage_rMSE_bootstrap, 1, sd)
  mu_interval_widths_rep <- apply(mu_interval_widths, c(2,3), mean)
  mu_avg_interval_width <- apply(mu_interval_widths_rep, 2, function(x) mean(x, na.rm=TRUE))
  mu_avg_interval_width_se <- apply(mu_interval_widths_rep, 2, function(x) sd(x, na.rm=TRUE))

  consolidated_results$`Mean: Bias` <- format_float_with_se(mu_avg_bias, mu_avg_bias_se)
  consolidated_results$`Mean: RMSE` <- format_float_with_se(mu_avg_rmse, mu_avg_rmse_se)
  consolidated_results$`Mean: CI Width` <- format_float_with_se(mu_avg_interval_width, mu_avg_interval_width_se)
  consolidated_results$`Mean: Marginal Coverage` <- format_percent_with_se(mu_marginal_coverage, mu_marginal_coverage_se)
  consolidated_results$`Mean: Coverage RMSE` <- format_float_with_se(mu_coverage_rMSE, mu_coverage_rMSE_se)

  # SD
  sd_bias_rep <- apply(sd_estimates - sd_true, MARGIN=c(2,3), function(x) mean(x, na.rm=TRUE))
  sd_avg_bias <- apply(sd_bias_rep, MARGIN=2, function(x) mean(x, na.rm=TRUE))
  sd_avg_bias_se <- apply(sd_bias_rep, MARGIN=2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
  sd_rmse_rep <- apply((sd_estimates - sd_true)^2, MARGIN=c(2,3), function(x) sqrt(mean(x, na.rm=TRUE)))
  sd_avg_rmse <- apply(sd_rmse_rep, MARGIN=2, function(x) mean(x, na.rm=TRUE))
  sd_avg_rmse_se <- apply(sd_rmse_rep, MARGIN=2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
  sd_rep_coverage <- apply(sd_covered, MARGIN=c(2,3), function(x) mean(x, na.rm=TRUE))
  sd_marginal_coverage <- apply(sd_rep_coverage, 2, function(x) mean(x, na.rm=TRUE))
  sd_marginal_coverage_se <- apply(sd_rep_coverage, 2, function(x) sd(x, na.rm=TRUE)/sqrt(length(x)))
  sd_coverage_rMSE <- get_coverage_rMSE(sd_covered)
  sd_coverage_rMSE_bootstrap <- sapply(seq(n_bootstrap), function(i) {
    idx <- sample(seq(reps), reps, replace=TRUE)
    sd_covered_boot <- sd_covered[,idx,]
    coverage_rMSE <- get_coverage_rMSE(sd_covered_boot)
    coverage_rMSE
  })
  sd_coverage_rMSE_se <- apply(sd_coverage_rMSE_bootstrap, 1, sd)
  sd_interval_widths_rep <- apply(sd_interval_widths, c(2,3), mean)
  sd_avg_interval_width <- apply(sd_interval_widths_rep, 2, function(x) mean(x, na.rm=TRUE))
  sd_avg_interval_width_se <- apply(sd_interval_widths_rep, 2, function(x) sd(x, na.rm=TRUE))
  
  consolidated_results$`SD: Bias` <- format_float_with_se(sd_avg_bias, sd_avg_bias_se)
  consolidated_results$`SD: RMSE` <- format_float_with_se(sd_avg_rmse, sd_avg_rmse_se)
  consolidated_results$`SD: CI Width` <- format_float_with_se(sd_avg_interval_width, sd_avg_interval_width_se)
  consolidated_results$`SD: Marginal Coverage` <- format_percent_with_se(sd_marginal_coverage, sd_marginal_coverage_se)
  consolidated_results$`SD: Coverage RMSE` <- format_float_with_se(sd_coverage_rMSE, sd_coverage_rMSE_se)

  # How often did we use score-based covariance?
  colSums(used_score_for_cov)

  # Make table pretty
  rownames(consolidated_results) <- method_names
  #View(consolidated_results)

  save.image(file=paste0("SimDLNTruth20251223", n, ".Rdata"))
}