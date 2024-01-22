# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
# require(mpcmp)
# require(cmp)

# This function fits the mean-parameterized CMP model via mpcmp
# delta is used to numerically approximate the gradient for the Delta method
fit_cmp <- function(formula_mu, formula_nu, data, delta=0.0001) {
  # Fit model
  mod_cmp <- mpcmp::glm.cmp(
    formula=formula_mu,
    formula_nu=formula_nu,
    data=data
  )
  
  # Extract model components
  X <- model.matrix(mod_cmp)$x
  Z <- model.matrix(mod_cmp)$s
  n <- nrow(X)
  p <- length(mod_cmp$coefficients)
  p_beta <- length(mod_cmp$coefficients_beta)
  p_gamma <- length(mod_cmp$coefficients_gamma)
  beta_idx <- seq(p_beta)
  gamma_idx <- seq(p_beta+1, p)
  
  # Calculate standard errors for fitted values
  fitted_ses <- sapply(seq(n), function(j) {
    x <- X[j,]
    mod_cmp$fitted_values[j] * sqrt(c(t(x) %*% mod_cmp$variance_beta %*% x))
    # NOTE: Need fitted value because we're applying delta method
  })
  
  # Create confidence intervals for fitted values
  mod_cmp$fitted_lower_bounds <- mod_cmp$fitted_values - 1.96 * fitted_ses
  mod_cmp$fitted_upper_bounds <- mod_cmp$fitted_values + 1.96 * fitted_ses
  mod_cmp$fitted_interval_widths <- mod_cmp$fitted_upper_bounds - mod_cmp$fitted_lower_bounds

  # Coverage for standard deviations
  sd_grads <- matrix(0, nrow=n, ncol=p)
  lambda_hats <- comp_lambdas(mod_cmp$fitted_values, mod_cmp$nu)$lambda
  var_hats <- .C(
    "cmp_cum_all", llambda=as.double(log(lambda_hats)), lnu=as.double(log(mod_cmp$nu)),
    flag=as.double(1L), n=n, mean=rep(0, n), var=rep(0, n))$var
  sd_hats <- sqrt(var_hats)
  get_sd_grad_j <- function(j) {
    coefs <- mod_cmp$coefficients
    coefs[j] <- coefs[j] + delta
    mu_hats_j <- exp(X %*% coefs[beta_idx])
    nu_hats_j <- exp(Z %*% coefs[gamma_idx])
    lambda_hats_j <- comp_lambdas(mu_hats_j, nu_hats_j)$lambda
    var_hats_j <- .C(
      "cmp_cum_all", llambda=as.double(log(lambda_hats_j)), lnu=as.double(log(nu_hats_j)),
      flag=as.double(1L), n=n, mean=rep(0, n), var=rep(0, n))$var
    sd_hats_j <- sqrt(var_hats_j)
    grad_j <- (sd_hats_j - sd_hats) / delta
    grad_j
  }
  for (j in seq(p)) {
    sd_grads[,j] <- get_sd_grad_j(j)
  }
  
  vcov_full <- matrix(0, nrow=p, ncol=p)
  vcov_full[beta_idx, beta_idx] <- mod_cmp$variance_beta
  vcov_full[gamma_idx, gamma_idx] <- mod_cmp$variance_gamma
  sd_ses <- sapply(seq(n), function(k) {
    grad_k <- sd_grads[k,]
    sqrt(t(grad_k) %*% vcov_full %*% grad_k)
  })
  
  mod_cmp$sd_lower_bounds <- sd_hats - 1.96 * sd_ses
  mod_cmp$sd_upper_bounds <- sd_hats + 1.96 * sd_ses
  mod_cmp$sd_interval_widths <- mod_cmp$sd_upper_bounds - mod_cmp$sd_lower_bounds
  
  mod_cmp
}