library(Matrix)
library(MASS)

ll_latentnorm <- function(y, X, Z, beta, alpha) {
  mu <- drop(X %*% beta)
  sigma <- exp(drop(Z %*% alpha))
  z_bar <- (log1p(y) - mu) / sigma 
  z_underbar <- (log(y) - mu) / sigma
  ll <- pnorm(z_underbar, log.p=TRUE) + log(exp(pnorm(z_bar, log.p=TRUE) - pnorm(z_underbar, log.p=TRUE)) - 1)
  return(ll)
}

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

z_underbar <- c(-100, -7, -1, 0, 6, 20, 99)
z_bar <-      c( -99, -6,  0, 1, 7, 21, 100)
q <- 2
get_kappa(z_bar, z_underbar, q)
(z_bar^q * dnorm(z_bar) - z_underbar^q * dnorm(z_underbar)) /
  (pnorm(z_bar) - pnorm(z_underbar))

gradutils <- function(y, X, beta, alpha) {
  mu <- drop(X %*% beta)
  sigma <- exp(drop(X %*% alpha))
  
  z_bar <- (log1p(y) - mu) / sigma 
  z_underbar <- (log(y) - mu) / sigma
  
  kappa_0 <- get_kappa(z_bar, z_underbar, 0)
  kappa_1 <- get_kappa(z_bar, z_underbar, 1)
  kappa_2 <- get_kappa(z_bar, z_underbar, 2)
  kappa_3 <- get_kappa(z_bar, z_underbar, 3)
  
  return(list(kappa_0, kappa_1, kappa_2, kappa_3, mu, sigma))
}

vec_matrix_multiply <- function(a, B) {
  return(t(apply(B, 1, function(x) x * a)))
}

em_gradutils <- function(W, sigma, c, alpha, return_hessian = FALSE) {
  sigma_neg_2 <- sigma^(-2)
  grad <- t(W) %*% (sigma_neg_2 * c - 1) - alpha
  hessian <- NULL
  if (return_hessian) {
    W_sqrt_k <- vec_matrix_multiply(sqrt(c)/sigma, W)
    hessian <- -2 * t(W_sqrt_k) %*% W_sqrt_k
  }
  return(list(grad, hessian))
}

discrete_lognormal <- function(endog, exog, start_params = NULL, method = "EM", maxiter = 100, use_hessian = FALSE, step_size = 1e-4, tol = 1e-6, penalty = 0) {
  nparams <- 22
  nloglikeobs <- function(params) {
    beta <- params[1:11]
    alpha <- params[12:22]
    ll <- ll_latentnorm(endog, exog, beta, alpha)
    params_alt <- params
    params_alt[1] <- 0
    return(-ll - penalty * sum(params_alt^2) / length(endog))
  }
  
  score <- function(params) {
    beta <- params[1:11]
    alpha <- params[12:22]
    gradutils_result <- gradutils(endog, exog, beta, alpha)
    kappa_0 <- gradutils_result[[1]]
    kappa_1 <- gradutils_result[[2]]
    kappa_2 <- gradutils_result[[3]]
    kappa_3 <- gradutils_result[[4]]
    mu <- gradutils_result[[5]]
    sigma <- gradutils_result[[6]]
    grad_beta <- -matrix(kappa_0 / sigma, nrow = length(endog), ncol = 1) %*% exog - penalty * 2 * beta
    grad_alpha <- -kappa_1 %*% exog - penalty * 2 * alpha
    return(c(grad_beta, grad_alpha))
  }
  
  hessian <- function(params) {
    beta <- params[1:11]
    alpha <- params[12:22]
    gradutils_result <- gradutils(endog, exog, beta, alpha)
    kappa_0 <- gradutils_result[[1]]
    kappa_1 <- gradutils_result[[2]]
    kappa_2 <- gradutils_result[[3]]
    kappa_3 <- gradutils_result[[4]]
    mu <- gradutils_result[[5]]
    sigma <- gradutils_result[[6]]
    k_beta <- (kappa_0^2 + kappa_1) / sigma^2
    k_alpha <- kappa_1 * (kappa_1 - 1) + kappa_3
    k_beta_alpha <- (kappa_2 + kappa_0 * (kappa_1 - 1)) / sigma
    H_beta <- Matrix(0, nrow = 11, ncol = 11)
    H_alpha <- Matrix(0, nrow = 11, ncol = 11)
    H_beta_alpha <- Matrix(0, nrow = 11, ncol = 11)
    for (i in 1:nrow(exog)) {
      x <- matrix(exog[i,], nrow = 1)
      xxT <- tcrossprod(x)
      H_beta <- H_beta - k_beta[i] * xxT
      H_alpha <- H_alpha - k_alpha[i] * xxT
      H_beta_alpha <- H_beta_alpha - k_beta_alpha[i] * xxT
    }
    H_all <- rbind(cbind(H_beta, H_beta_alpha), cbind(t(H_beta_alpha), H_alpha))
    penalty_matrix <- penalty * 2 * diag(nparams)
    penalty_matrix[1, 1] <- 0
    return(H_all - penalty_matrix)
  }
  
  if (is.null(start_params)) {
    start_params <- rep(0, nparams)
    start_params[1] <- log(mean(endog))
  }
  
  if (method == "EM") {
    loss <- mean(nloglikeobs(start_params))
    beta <- start_params[1:11]
    alpha <- start_params[12:22]
    W <- exog
    penalty_alpha <- penalty * diag(11)
    WtW_plus_penalty <- t(W) %*% W + penalty_alpha
    penalty_beta <- penalty_alpha
    penalty_beta[1, 1] <- 0
    converged <- FALSE
    for (i in 1:maxiter) {
      loss_last <- loss
      gradutils_result <- gradutils(endog, exog, beta, alpha)
      e1 <- mu <- gradutils_result[[5]] - sigma <- gradutils_result[[6]] * kappa_0 <- gradutils_result[[1]]
      e2 <- (sigma^2 - sigma^2 * kappa_1 + mu^2 - 2 * mu * sigma * kappa_0)
      X_sqrt_w <- vec_matrix_multiply(1 / sigma, exog)
      XtSiX <- t(X_sqrt_w) %*% X_sqrt_w
      XtSiX <- XtSiX + penalty_beta
      XtSie1 <- t(exog) %*% (sigma^(-2) * e1)
      beta <- solve(XtSiX, XtSie1)
      c <- e2 - 2 * e1 * mu + mu^2
      if (use_hessian) {
        em_gradutils_result <- em_gradutils(W, sigma, c, alpha, return_hessian = TRUE)
        grad <- em_gradutils_result[[1]]
        hessian <- em_gradutils_result[[2]]
        grad <- grad - alpha
        hessian <- hessian - penalty_alpha
        alpha <- alpha - solve(hessian, grad)
      } else {
        grad <- em_gradutils(W, sigma, c, alpha)
        d <- -grad
        prop_increase <- 0.5
        step_multiplier <- 0.5
        curr_step_size <- step_size
        f_start <- sum(nloglikeobs(c(beta, alpha)))
        while (TRUE) {
          alpha_prop <- alpha - curr_step_size * d
          f_stop <- sum(nloglikeobs(c(beta, alpha_prop)))
          required_change <- prop_increase * curr_step_size * sum(grad * d)
          if (f_stop - f_start <= required_change) {
            break
          }
          curr_step_size <- curr_step_size * step_multiplier
        }
      }
      alpha_prop <- alpha - curr_step_size * d
      alpha <- alpha_prop
      loss <- mean(nloglikeobs(c(beta, alpha)))
      obj <- loss_last - loss
      if (abs(obj) < tol) {
        converged <- TRUE
        break
      }
    }
    if (!converged) {
      stop("Hit maxiter and failed to converge")
    }
    params <- c(beta, alpha)
    return(list(params = params, loss = loss, converged = converged))
  } else {
    stop("Unsupported method")
  }
}