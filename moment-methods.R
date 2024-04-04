# Code for implementing the moment-based methods

# epl: pseudolikelihood model
epl <- function(y, X, Z, betastart, alphastart, Xnew=NULL, Znew=NULL, tol=1e-8, max_iter=100, stephalving_maxiter=10, verbose=FALSE) {
  start_time <- Sys.time()
  result_list <- list()
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  beta_idx <- seq(p)
  alpha_idx <- seq(p+1, p+q)
  
  get_mu <- function(beta) c(exp(X %*% beta))
  get_phi <- function(alpha) c(exp(Z %*% alpha))
  
  get_row_grads <- function(beta, alpha) {
    # Helpers
    mu <- get_mu(beta)
    phi <- get_phi(alpha)
    error <- y - mu
    r2 <- error^2 / mu
    
    # Individual row_grads
    row_grads1 <- (error / phi) * X
    row_grads2 <- (r2 / phi - 1) * Z
    
    # Assemble them
    row_grads <- cbind(row_grads1, row_grads2)
    row_grads
  }
  
  get_grad <- function(beta, alpha) {
    row_grads <- get_row_grads(beta, alpha)
    grad <- colMeans(row_grads)
    grad
  }
  
  get_hess <- function(beta, alpha) {
    mu <- get_mu(beta)
    phi <- get_phi(alpha)
    error <- y - mu
    r2 <- error^2 / mu
    
    # 11
    weight11 <- sqrt(mu / phi)
    hess11 <- -crossprod(weight11 * X)
      
    # 12
    hess12 <- -(t((error / phi) * X) %*% Z)
      
    # 21
    hess21 <- -t(((y^2 - mu^2) / (phi * mu)) * Z) %*% X
      
    # 22
    weight22 <- sqrt(r2 / phi)
    hess22 <- -crossprod(weight22 * Z)
      
    # Arrange in block matrix
    hess <- rbind(
      cbind(hess11, hess12),
      cbind(hess21, hess22)) / n
    hess
  }
  
  get_dev <- function(beta, alpha) {
    mu <- get_mu(beta)
    phi <- get_phi(alpha)
    error <- y - mu
    r2 <- error^2 / mu
    dev <- mean(r2/phi) + mean(log(2 * pi * phi * mu))
    dev
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
    
    # Update params
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
    if (verbose) cat("Iteration: ", iter, ", Deviance: ", dev, "\n", sep="")
  }
  
  if (iter > max_iter) warning("Hit max_iter and failed to converge")
  end_fit_time <- Sys.time()
  if (verbose) print("Finished model fitting; beginning inference")
  
  # Sandwich estimator
  row_grads <- get_row_grads(beta, alpha) / n
  hess <- get_hess(beta, alpha)
  meat <- crossprod(row_grads)
  half_sandwich <- solve(hess, t(chol(meat)))
  sandwich <- tcrossprod(half_sandwich)
  
  # Save basic results
  result_list$beta <- beta
  result_list$alpha <- alpha
  result_list$theta <- theta
  result_list$cov_beta <- sandwich[beta_idx, beta_idx]
  result_list$cov_alpha <- sandwich[alpha_idx, alpha_idx]
  result_list$cov_theta <- sandwich
  
  # Fitted values
  if (is.null(Xnew)) Xnew <- X
  if (is.null(Znew)) Znew <- Z
  mu <- c(exp(Xnew %*% beta))
  result_list$mu <- mu
  result_list$fitted_values <- mu
  
  # Confidence intervals for fitted values
  fitted_ses <- sqrt(rowSums((Xnew %*% result_list$cov_beta) * Xnew)) * mu
  result_list$fitted_lower_bounds <- mu - 1.96 * fitted_ses
  result_list$fitted_upper_bounds <- mu + 1.96 * fitted_ses
  result_list$fitted_interval_widths <- result_list$fitted_upper_bounds - result_list$fitted_lower_bounds
  
  # Estimated standard deviations
  phi <- c(exp(Znew %*% alpha))
  var_estimates <- phi * mu
  sd_estimates <- sqrt(var_estimates)
  result_list$phi <- phi
  result_list$var_estimates <- var_estimates
  result_list$sd_estimates <- sd_estimates
  
  # Confidence intervals for standard deviations
  W <- cbind(Xnew, Znew)
  sd_grads <- result_list$sd_estimates * W
  sd_ses <- 0.5 * sqrt(rowSums((sd_grads %*% result_list$cov_theta) * sd_grads))
  result_list$sd_lower_bounds <- result_list$sd_estimates - 1.96 * sd_ses
  result_list$sd_upper_bounds <- result_list$sd_estimates + 1.96 * sd_ses
  result_list$sd_interval_widths <- result_list$sd_upper_bounds - result_list$sd_lower_bounds
  if (verbose) print("Finished inference")

  # Timing
  end_time <- Sys.time()
  result_list$fit_time <- end_fit_time - start_time
  result_list$inference_time <- end_time - end_fit_time
  result_list$elapsed_time <- end_time - start_time

  result_list
}

# Z <- model.matrix(~x3)
# mod_quasipois <- glm(y ~ x1 + x2, data=dat, family=quasipoisson())
# phi0 <- summary(mod_quasipois)$dispersion
# alpha_start <- rep(0, ncol(Z))
# alpha_start[1] <- log(phi0)
# names(alpha_start) <- c("(Intercept)", "x3")
# mod_epl <- epl(y, X, Z, betastart = beta_start, alphastart = alpha_start, verbose=TRUE)