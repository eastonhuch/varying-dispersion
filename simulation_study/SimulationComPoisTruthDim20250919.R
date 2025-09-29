################################
# Simulation to compare mean and sd estimation
# when data is generated using COM-Pois
# September 19, 2025
# In this version, I'm adding additional covariates
################################

require(mpcmp) #install_github("thomas-fung/mpcmp")
require(cmp) #install_github("SuneelChatla/cmp")
require(scales)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

# Generate some data

myns <- c(50L, 100L, 250L, 500L, 1000L)
for(n in myns){
  set.seed(20240501+n)
  #n <- 400L
  x1 <- rnorm(n)  # We'll keep the design fixed across simulations
  x2 <- rnorm(n)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rpois(n, 1)

  # NOTE: The mpcmp package can be a bit temperamental
  # You might run into errors if you make the coefficients much larger
  beta_true <- c(3, 0.05, -0.1, 0.02, 0.04, -0.01) # Regression coefficients for the mean
  y <- rep(0, n)  # We'll create this later
  master_formula <- y ~ x1 + x2 + I(x1*x2) + x3 + x4
  X <- model.matrix(master_formula)
  log_mu <- c(X %*% beta_true)
  mu <- exp(log_mu)
  alpha_true <- c(0.1, 0, -0.2, 0.05, -0.2, 0.1) # Not directly comparable across all models; these values of alpha (alpha0=1; combined with the values of beta) result in underdispersion; to create overdispersed data, use alpha0=-1 or alpha0=-0.5
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
  tol <- 1e-10
  max_iter <- 100
  phi_method <- "joint"
  stephalving_max <- 50

  # Prepare to loop
  reps <- 1000L
  method_names <- c("Pois", "QP", "NB","MPCMP", "GP-1", "EPL", "DLN-N", "DLN-EM")
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

    # Fit Pois model
    if("Pois" %in% method_names){
      start_time <- Sys.time()
      mod_pois <- glm(master_formula, family="poisson", data=dat)
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
      mod_nb <- glm.nb(master_formula, data=dat)
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
      sd_lbnd <- lowerbnd*(lowerbnd+mod_nb$theta)/mod_nb$theta #this is the variance even though it's labeled as sd
      sd_ubnd <- upperbnd*(upperbnd+mod_nb$theta)/mod_nb$theta
      sd_covered[,i,"NB"] <- (sqrt(sd_lbnd) <= sd_true) & (sd_true <= sqrt(sd_ubnd))
      sd_interval_widths[,i,"NB"] <- sqrt(sd_ubnd) - sqrt(sd_ubnd)
    })

    # Fit MPCMP model
    if ("MPCMP" %in% method_names ) {
      start_time <- Sys.time()
      mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-4), error=function(e) e)
      if(!is.null(mod_cmp$message)){
        mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-4, beta_true=beta_true), error=function(e) e)
      }
      if(!is.null(mod_cmp$message)){
        mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-4, beta_true=beta_true), error=function(e) e)
      }
      if(!is.null(mod_cmp$message)){
        mod_cmp <- tryCatch(fit_cmp(master_formula, master_formula, data=dat, delta=1e-4, beta_true=beta_true), error=function(e) e)
      }
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      fitting_times[i,"MPCMP"] <- elapsed_time
      if(is.null(mod_cmp$message)){
      mu_estimates[,i,"MPCMP"] <- mod_cmp$fitted_values
      sd_estimates[,i,"MPCMP"] <- mod_cmp$sd_estimates
      mu_covered[,i,"MPCMP"] <- (mod_cmp$fitted_lower_bounds <= mu) & (mu <= mod_cmp$fitted_upper_bounds)
      mu_interval_widths[,i,"MPCMP"] <- mod_cmp$fitted_interval_widths
      sd_covered[,i,"MPCMP"] <- (mod_cmp$sd_lower_bounds <= sd_true) & (sd_true <= mod_cmp$sd_upper_bounds)
      sd_interval_widths[,i,"MPCMP"] <- mod_cmp$sd_interval_widths
      }
    }

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
  format_float <- label_number(0.001)
  format_percent <- label_percent(0.1)

  # Fitting times
  avg_fitting_times <- colMeans(fitting_times)
  consolidated_results <- data.frame(`Avg Elapsed Time (Seconds)`=format_float(avg_fitting_times), check.names=FALSE)
  boxplot(fitting_times, log="y", main="Model-fitting Time by Method", ylab="Elapsed Time (Seconds)")

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
  sd_avg <- apply(sd_estimates, MARGIN=c(1,3), mean)
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

  # Plot avg SD vs. x2
  sd_order <- order(sd_true)
  plot(x2[sd_order], sd_true[sd_order], xlab=expression(x[2]), ylab="SD", ylim=c(2.5, 7))
  points(x2[sd_order], sd_avg[sd_order, "DLN-N"], col=2)
  points(x2[sd_order], sd_avg[sd_order, "MPCMP"], col=3)

  # How often did we use score-based covariance?
  colSums(used_score_for_cov)

  # Make table pretty
  final_table <- t(consolidated_results)
  #View(final_table)

  save.image(file=paste0("SimComPoisTruthDim20250919",n,".Rdata"))
}