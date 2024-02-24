library(Rcpp)
sourceCpp("discrete-log-normal-helpers.cpp")

ll_latentnorm <- function(y, X, Z, beta, alpha) {
  mu <- drop(X %*% beta)
  sigma <- exp(drop(Z %*% alpha))
  z_bar <- (log1p(y) - mu) / sigma 
  z_underbar <- (log(y) - mu) / sigma
  ll <- pnorm(z_underbar, log.p=TRUE) + log(exp(pnorm(z_bar, log.p=TRUE) - pnorm(z_underbar, log.p=TRUE)) - 1)
  mean(ll)
}

# This is the hardest part!
get_kappa <- function(z_bar, z_underbar, q) {
  if (any(z_underbar >= z_bar)) stop("z_bar must be greater than z_underbar")
  
  # -Inf for z_underbar
  kappa_q_neg_inf <- exp(dnorm(z_bar, log=TRUE) - pnorm(z_bar, log.p=TRUE))
  
  # Typical case
  # Denominator
  z_underbar_pos <- z_underbar > 0
  z_underbar_neg <- !z_underbar_pos  # Technically includes 0 too, but I deal with this separately
  
  log_cdf_bar <- pnorm(z_bar, log.p=TRUE)
  log_cdf_underbar <- pnorm(z_underbar, log.p=TRUE)
  log_den_neg <- log_cdf_underbar + log(exp(log_cdf_bar - log_cdf_underbar) - 1)

  log_cdf_neg_bar <- pnorm(-z_bar, log.p=TRUE)
  log_cdf_neg_underbar <- pnorm(-z_underbar, log.p=TRUE)
  log_dem_pos <- log_cdf_neg_bar + log(exp(log_cdf_neg_underbar - log_cdf_neg_bar) - 1)
  
  log_den <- rep(0, length(z_bar))
  log_den[z_underbar_pos] <- log_dem_pos[z_underbar_pos]
  log_den[z_underbar_neg] <- log_den_neg[z_underbar_neg]
  
  # Numerator
  z_bar_q_pos <- (z_bar^q) > 0
  z_underbar_q_pos <- (z_underbar^q) > 0
  both_pos <- z_bar_q_pos & z_underbar_q_pos
  both_neg <- (!z_bar_q_pos) & (!z_underbar_q_pos)
  sign_switch <- z_bar_q_pos & (!z_underbar_q_pos)
  log_mod_bar <- q * log(abs(z_bar)) + dnorm(z_bar, log=TRUE)
  log_mod_underbar <- q * log(abs(z_underbar)) + dnorm(z_underbar, log=TRUE)
  bar_mod_larger <- log_mod_bar > log_mod_underbar
  keep_order <- sign_switch | (both_pos & bar_mod_larger) | (both_neg & (!bar_mod_larger))
  sign <- (-1)^(!keep_order)
  log_mod_b <- keep_order * log_mod_bar + (!keep_order) * log_mod_underbar
  log_mod_a <- (!keep_order) * log_mod_bar + keep_order * log_mod_underbar
  mod_b_greater <- log_mod_b > log_mod_a
  log_mod_c <- (!mod_b_greater) * log_mod_b + mod_b_greater * log_mod_a
  log_mod_d <- mod_b_greater * log_mod_b + (!mod_b_greater) * log_mod_a
  const <- (-1)^(!sign_switch)
  log_num <- log_mod_c + log(exp(log_mod_d - log_mod_c) + const)
  kappa_q <- sign * exp(log_num - log_den)
  
  # 0 for z_bar
  kappa_q_0_z_bar <- -exp(log_mod_underbar - log_den)
  
  # 0 for z_underbar
  kappa_q_0_z_underbar <- exp(log_mod_bar - log_den)
  
  # Put these together
  z_underbar_neg_inf <- z_underbar == -Inf
  kappa_q[z_underbar_neg_inf] <- kappa_q_neg_inf[z_underbar_neg_inf]
  
  z_bar_0 <- z_bar == 0
  kappa_q[z_bar_0] <- kappa_q_0_z_bar[z_bar_0]
  
  z_underbar_0 <- z_underbar == 0
  kappa_q[z_underbar_0] <- kappa_q_0_z_underbar[z_underbar_0]
  
  # Return result
  kappa_q
}

# Tests for the above function
# z_underbar <- c(-Inf, -100, -7, -1, 0, 6, 20, 99)
# z_bar <-      c(-4. ,  -99, -6,  0, 1, 7, 21, 100)
# q <- 3
# get_kappa(z_bar, z_underbar, q)
# (z_bar^q * dnorm(z_bar) - z_underbar^q * dnorm(z_underbar)) /
#   (pnorm(z_bar) - pnorm(z_underbar))

gradutils <- function(y, X, Z, beta, alpha) {
  mu <- drop(X %*% beta)
  sigma <- exp(drop(Z %*% alpha))
  
  z_bar <- (log1p(y) - mu) / sigma 
  z_underbar <- (log(y) - mu) / sigma
  
  kappa_0 <- get_kappa(z_bar, z_underbar, 0)
  kappa_1 <- get_kappa(z_bar, z_underbar, 1)
  kappa_2 <- get_kappa(z_bar, z_underbar, 2)
  kappa_3 <- get_kappa(z_bar, z_underbar, 3)
  
  list(
    kappa_0=kappa_0,
    kappa_1=kappa_1,
    kappa_2=kappa_2,
    kappa_3=kappa_3,
    mu=mu,
    sigma=sigma
  )
}

em_gradutils <- function(Z, sigma, v, alpha) {
  grad <- t(Z) %*% (v/sigma^2 - 1)
  hess <- -2 * crossprod(sqrt(v)/sigma * Z)
  list(
    grad=grad,
    hess=hess
  )
}

dln <- function(
  y, X, Z, betastart, alphastart, method="Newton",
  max_iter = 100, stephalving_maxiter=10, tol=1e-8, verbose=TRUE) {
  
  # Initialize some values
  result_list <- list()
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  beta_idx <- seq(p)
  alpha_idx <- seq(p+1, p+q)
  thetastart <- c(betastart, alphastart)
  
  # A bunch of helper functions
  get_log_like <- function(beta, alpha) ll_latentnorm(y, X, Z, beta, alpha)
  get_dev <- function(beta, alpha) -2 * get_log_like(beta, alpha)  # Doesn't include additive factor for saturated model

  get_grad <- function(beta, alpha) {
    gradutils_result <- gradutils(y, X, Z, beta, alpha)
    grad_beta <- -colMeans(gradutils_result$kappa_0 / gradutils_result$sigma * X)
    grad_alpha <- -colMeans(gradutils_result$kappa_1 * Z)
    c(grad_beta, grad_alpha)
  }
  
  get_hess <- function(beta, alpha) {
    # Calculate kappas, mu, and sigma
    gradutils_result <- gradutils(y, X, Z, beta, alpha)
    kappa_0 <- gradutils_result$kappa_0
    kappa_1 <- gradutils_result$kappa_1
    kappa_2 <- gradutils_result$kappa_2
    kappa_3 <- gradutils_result$kappa_3
    mu <- gradutils_result$mu
    sigma <- gradutils_result$sigma
    
    k_beta <- (kappa_0^2 + kappa_1) / sigma^2
    k_alpha <- kappa_1 * (kappa_1 - 1) + kappa_3
    k_beta_alpha <- (kappa_2 + kappa_0 * (kappa_1 - 1)) / sigma
    
    hess11 <- -crossprod(sqrt(k_beta) * X) / n
    hess12 <- -(t(X) %*% (k_beta_alpha * Z)) / n
    hess22 <- -crossprod(sqrt(k_alpha) * Z) / n
    hess <- rbind(
      cbind(hess11, hess12),
      cbind(t(hess12), hess22)
    )
    hess
  }
  
  # Prepare to loop
  dev_last <- Inf
  dev <- get_dev(betastart, alphastart)
  beta <- betastart
  alpha <- alphastart
  theta <- c(betastart, alphastart)
  
  # Model-fitting loop
  iter <- 0
  if (verbose) cat("Iteration: ", iter, ", Deviance: ", dev, "\n", sep="")
  while (((dev_last - dev) > tol) && (iter <= max_iter)) {
    # Looping mechanics
    iter <- iter + 1
    dev_last <- dev
    
    if (method == "Newton") {
      grad <- get_grad(beta, alpha)
      hess <- get_hess(beta, alpha)
      inc <- solve(hess, grad)
      dev <- Inf
      step <- 1
      stephalving_iter <- 0
      while (((dev >= dev_last) || (!is.finite(dev))) && (stephalving_iter <= stephalving_maxiter)) {
        theta_new <- theta - step * inc
        dev <- get_dev(theta_new[beta_idx], theta_new[alpha_idx])
        if (verbose && (stephalving_iter > 0)) {
          cat(
            "Step-halving Iterations: ", stephalving_iter,
            ", Deviance: ", dev,
            "\n", sep="")
        }
        stephalving_iter <- stephalving_iter + 1
        step <- step / 2
      }
      theta <- theta_new
      beta <- theta[beta_idx]
      alpha <- theta[alpha_idx]
    } else if (method == "EM") {
      # Update beta
      gradutils_result <- gradutils(y, X, Z, beta, alpha)
      e1 <- gradutils_result$mu - gradutils_result$sigma * gradutils_result$kappa_0
      beta <- solve(
        crossprod(X / gradutils_result$sigma),
        t(X) %*% (e1 / gradutils_result$sigma^2))
      
      # Update alpha
      gradutils_result <- gradutils(y, X, Z, beta, alpha)  # Need to update these again
      mu <- gradutils_result$mu
      sigma <- gradutils_result$sigma
      kappa_0 <- gradutils_result$kappa_0
      kappa_1 <- gradutils_result$kappa_1
      e1 <- mu - sigma * kappa_0
      e2 <- (sigma^2 - sigma^2*kappa_1+mu^2 - 2*mu*sigma*kappa_0) # Need to double-check this
      v <- e2 - 2*e1*mu + mu^2
      em_gradutils_result <- em_gradutils(Z, sigma, v, alpha)
      inc <- solve(em_gradutils_result$hess, em_gradutils_result$grad)
      dev <- Inf
      step <- 1
      stephalving_iter <- 0
      while (((dev >= dev_last) || (!is.finite(dev))) && (stephalving_iter <= stephalving_maxiter)) {
        alpha_new <- alpha - step * inc
        dev <- get_dev(beta, alpha_new)
        if (verbose && (stephalving_iter > 0)) {
          cat(
            "Step-halving Iterations: ", stephalving_iter,
            ", Deviance: ", dev,
            "\n", sep="")
        }
        stephalving_iter <- stephalving_iter + 1
        step <- step / 2
      }
      alpha <- alpha_new
      theta <- c(beta, alpha)
    } else {
      stop("method must be 'Newton' or 'EM'")
    }
    if (verbose) cat("Iteration: ", iter, ", Deviance: ", dev, "\n", sep="")
  }
  
  # Add results to list
  theta <- c(beta, alpha)
  names(theta) <- names(thetastart)
  result_list$beta <- beta
  result_list$alpha <- alpha
  result_list$theta <- theta
  hess <- get_hess(beta, alpha)
  result_list$cov_theta <- chol2inv(chol(-hess))
  result_list$cov_beta <- result_list$cov_theta[beta_idx, beta_idx]
  result_list$cov_alpha <- result_list$cov_theta[alpha_idx, alpha_idx]
  
  # Fitted values
  z_mu <- drop(exp(X %*% beta))
  z_sigma <- drop(exp(Z %*% alpha))
  result_list$fitted_values <- integrate_dln(z_mu, z_sigma, 1)
  mean_mu_derivs <- integrate_dln(z_mu, z_sigma, 1, calc_deriv=TRUE, wrt_mu=TRUE)
  mean_sigma_derivs <- integrate_dln(z_mu, z_sigma, 1, calc_deriv=TRUE, wrt_mu=FALSE)
  mean_beta_grads <- mean_mu_derivs * z_mu * X
  mean_alpha_grads <- mean_sigma_derivs * z_sigma * Z
  mean_theta_grads <- cbind(mean_beta_grads, mean_alpha_grads)
  fitted_ses <- sqrt(rowSums((mean_theta_grads %*% result_list$cov_theta) * mean_theta_grads))
  result_list$fitted_lower_bounds <- result_list$fitted_values - 1.96 * fitted_ses
  result_list$fitted_upper_bounds <- result_list$fitted_values + 1.96 * fitted_ses
  result_list$fitted_interval_widths <- result_list$fitted_upper_bounds - result_list$fitted_lower_bounds
  
  # Standard deviations
  fitted_second_moments <- integrate_dln(z_mu, z_sigma, 2)
  result_list$sd_estimates <- fitted_second_moments - result_list$fitted_values^2
  moment2_mu_derivs <- integrate_dln(z_mu, z_sigma, 2, calc_deriv=TRUE, wrt_mu=TRUE)
  moment2_sigma_derivs <- integrate_dln(z_mu, z_sigma, 2, calc_deriv=TRUE, wrt_mu=FALSE)
  sd_beta_grads <- 0.5 * result_list$sd_estimates * (moment2_mu_derivs - 2 * result_list$fitted_values * mean_mu_derivs) * z_mu * X
  sd_alpha_grads <- 0.5 * result_list$sd_estimates * (moment2_sigma_derivs - 2 * result_list$fitted_values * mean_sigma_derivs) * z_sigma * Z
  sd_theta_grads <- cbind(sd_beta_grads, sd_alpha_grads)
  sd_ses <- sqrt(rowSums((sd_theta_grads %*% result_list$cov_theta) * sd_theta_grads))
  result_list$sd_lower_bounds <- result_list$sd_estimates - 1.96 * sd_ses
  result_list$sd_upper_bounds <- result_list$sd_estimates + 1.96 * sd_ses
  result_list$sd_interval_widths <- result_list$sd_upper_bounds - result_list$sd_lower_bounds
  
  result_list
}

mod_dln_newton <- dln(y, X, Z, beta_start, alpha_start, method="Newton", max_iter=100, stephalving_maxiter=10, tol=1e-16, verbose=TRUE)
mod_dln_em <- dln(y, X, Z, beta_start, alpha_start, method="EM", max_iter=100, stephalving_maxiter=10, tol=1e-16, verbose=TRUE)

# Same estimates!!!
mod_dln_newton$theta
mod_dln_em$theta

# Use this for testing
# betastart <- beta_start
# alphastart <- alpha_start
# max_iter <- 100
# stephalving_maxiter <- 10
# tol <- 1e-8
# verbose <- TRUE
