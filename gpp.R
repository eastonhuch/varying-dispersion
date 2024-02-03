# Requirement(s)
require(RMKdiscrete) # Used for random sampling from GP-P (not necessary to fit model)
require(cmna) # Used for bisection function (only necessary for moment estimator)

################################################################################
# Generalized Poisson functions for (theta, delta) parameterization

repeat_if_necessary <- function(val, n) {
  # If val is length 1, then repeat it n times
  if (length(val) == 1) rep(val, n) else val
}

check_gp_parameter_space <- function(theta, delta) {
  # Returns TRUE if the parameter space is value
  # Otherwise, issues a warning and returns FALSE

  # Check on individual parameters
  theta_is_valid <- theta > 0
  delta_is_valid <- abs(delta) <= 1

  # Check combination of parameter values
  delta_is_nonnegative <- delta >= 0
  lower_bound <- rep(-1, length(delta))
  neg_theta_over_4 <- -theta/4
  low_theta <- neg_theta_over_4 > -1
  lower_bound[low_theta] <- neg_theta_over_4[low_theta]
  delta_above_lower_bound <- lower_bound < delta
  combination_is_valid <- delta_is_nonnegative | delta_above_lower_bound

  # Bring all conditions together
  is_valid <- theta_is_valid & delta_is_valid & combination_is_valid
  if (! all(is_valid)) warning("Some parameter values are invalid")

  # Return
  is_valid
}

get_max_y <- function(theta, delta) {
  # Calculates vector of highest permissible y values
  # NOTE: This function assumes that the parameter values are valid
  # Use check_gp_parameter_space to validate whether that's the case
  m_raw <- -theta/delta
  m_raw_floor <- floor(m_raw)
  m_for_negative_delta <- m_raw_floor - as.integer(m_raw_floor == m_raw)
  n <- length(theta)
  delta_repped <- repeat_if_necessary(delta, n)
  m <- ifelse(delta_repped >= 0, Inf, m_for_negative_delta)
  m
}

check_gp_y_values <- function(y, theta, delta) {
  # Returns TRUE if all of the y values are valid
  # Otherwise, it issues a warning and returns FALSE
  m <- get_max_y(theta, delta)
  ys_are_valid <- y <= m
  any_problems <- (!all(ys_are_valid)) || any(is.na(ys_are_valid))
  if (any_problems) warning("Some y values are above the support")
  ys_are_valid
}

dgpois <- function(y, theta, delta, log=FALSE, omit_constant=FALSE) {
  # Calculates the probability mass function for the generalized Poisson distribution

  # Check parameter space and support
  parameters_are_valid <- check_gp_parameter_space(theta, delta)
  ys_are_valid <- check_gp_y_values(y=y, theta=theta, delta=delta)
  everything_is_valid <- parameters_are_valid & ys_are_valid

  # Suppress additional warning for this step
  old_warning_value <- getOption("warn")
  options(warn=-1)
  val1 <- log(theta) + (y-1)*log(theta + delta*y) - theta - delta*y
  val1[is.na(val1)] <- -Inf
  options(warn=old_warning_value)

  # Now transform
  val2 <- if (omit_constant) val1 else val1 - lfactorial(y)
  val3 <- ifelse(everything_is_valid, val2, -Inf)
  val4 <- if (log) val3 else exp(val3)
  val4
}

pgpois <- function(y, theta, delta){
  # Calculates cumulative distribution function for generalized Poisson distribution
  n <- length(y)
  theta_repped <- repeat_if_necessary(theta, n)
  delta_repped <- repeat_if_necessary(delta, n)
  sum_probs <- function(i) {
    probs <- dgpois(seq(0, y[i]), theta=theta_repped[i], delta=delta_repped[i], log=FALSE, propto=FALSE)
    sum(probs)
  }
  sapply(seq_along(y), sum_probs)
}

rgpois <- function(n, theta, delta) {
  # Draws random samples from the generalized Poisson distribution
  rLGP(n=n, theta=theta, lambda=delta)
}

################################################################################
# GP-P functions, parameterized by (mu, phi, P)
# NOTE: These functions transform parameters back to the original (theta, delta) parameter space

gpp_to_standard_gp <- function(mu, phi, P) {
  # Transforms (mu, phi, P) to (theta, delta)
  mu_to_P_minus_one <- mu^(P-1)
  denominator <- 1 + phi * mu_to_P_minus_one
  theta <- mu / denominator
  delta <- phi * mu_to_P_minus_one / denominator
  cbind(theta, delta)
}

dgpoisP <- function(y, mu, phi, P, log=FALSE, omit_constant=FALSE) {
  # Calculates probability mass function for GP-P distribution
  theta_delta <- gpp_to_standard_gp(mu=mu, phi=phi, P=P)
  dgpois(y=y, theta=theta_delta[,1], delta=theta_delta[,2], log=log, omit_constant=omit_constant)
}

pgpoisP <- function(y, mu, phi, P){
  # Calculates cumulative distribution function for GP-P distribution
  theta_delta <- gpp_to_standard_gp(mu=mu, phi=phi, P=P)
  pgpois(y=y, theta=theta_delta[,1], delta=theta_delta[,2])
}

rgpoisP <- function(n, mu, phi, P){
  # Draws random samples from the GP-P distribution
  theta_delta <- gpp_to_standard_gp(mu=mu, phi=phi, P=P)
  rgpois(n=n, theta=theta_delta[,1], delta=theta_delta[,2])
}

################################################################################
# Helper functions for fitting GP-P regression model
get_objective <- function(l_new, l_old) (l_new - l_old) / abs(l_old)

update_working_estimates <- function(working_list) {
  # Updates eta, mu, theta, delta, loss, and objective
  working_list[["eta_no_offset"]] <- working_list[["X"]] %*% working_list[["beta"]]
  working_list[["eta"]] <- working_list[["eta_no_offset"]] + working_list[["offset"]]
  working_list[["mu"]] <- exp(working_list[["eta"]])
  theta_delta <- gpp_to_standard_gp(
    working_list[["mu"]],
    working_list[["phi"]],
    working_list[["P"]])
  working_list[["theta"]] <- theta_delta[,1]
  working_list[["delta"]] <- theta_delta[,2]
  log_likelihood <- sum(dgpoisP(
    y=working_list[["y"]],
    mu=working_list[["mu"]],
    phi=working_list[["phi"]],
    P=working_list[["P"]],
    log=TRUE))
  penalty_term <- 0.5 * as.numeric(
    t(working_list[["regularized_parameters"]]) %*%
      working_list[["penalty"]] %*%
      working_list[["regularized_parameters"]])
  working_list[["loss"]] <- penalty_term - log_likelihood
  working_list[["objective"]] <- ifelse(
    is.null(working_list[["loss_old"]]),
    Inf,
    get_objective(working_list[["loss"]], working_list[["loss_old"]]))
  working_list
}

update_old_params <- function(working_list) {
  # Stores old versions of parameters for future reference
  params <- c(
    "beta",
    "phi",
    "scored_parameters",
    "regularized_parameters",
    "loss")
  for (param in params) working_list[[paste0(param, "_old")]] <- working_list[[param]]
  working_list
}

update_score_information <- function(working_list) {
  # Updates the score and Fisher information

  # Pull values out of list for simplicity
  y <- working_list[["y"]]
  X <- working_list[["X"]]
  mu <- working_list[["mu"]]
  phi <- working_list[["phi"]]
  P <- working_list[["P"]]

  # Helpers
  piece2 <- 1 + phi*mu^(P-1)
  wii <- c((mu + 2 * {(P-2)*phi*(mu^{P-1})}^2 / {1 + 2*phi*(mu^{P-2})}) /
      (piece2^2))
  dldmu <- 1/mu +
    (y-1) *
      (1 + {P-1}*phi*y*mu^{P-2}) /
      (mu + phi*y*mu^{P-1}) -
    (1 + 2*{P-1}*phi*y*mu^{P-2}) /
      piece2 +
    ({(P-1)*phi*mu^(P-1)} * {1 + phi*y*mu^(P-2)}) /
      (piece2^2)
  denom <- (piece2^2) * (1 + 2*phi*mu^{P-2})

  # Score (U)
  working_list[["U_beta"]] <- c(t(dldmu*mu) %*% X)
  if (working_list[["phi_method"]] == "joint") {
    working_list[["U_phi"]] <- sum(
      mu^{(P-2)}*y*(y-1) / {1 + phi*(mu^{P-2})*y} -
        mu^{(P-1)}*{y-mu} / {piece2^2} -
        mu^{(P-1)}*y / piece2)
    working_list[["U"]] <- c(working_list[["U_beta"]], working_list[["U_phi"]])
  } else {
    working_list[["U"]] <- working_list[["U_beta"]]
  }
  working_list[["U_star"]] <- working_list[["U"]] - as.numeric(
    working_list[["penalty"]] %*%
    working_list[["regularized_parameters"]]
  )[seq(working_list[["num_scored_parameters"]])]

  # Fisher information (J)
  X_sqrt_wii <- sqrt(wii) * X
  working_list[["J_beta"]] <- crossprod(X_sqrt_wii)
  if (working_list[["phi_method"]] == "joint") {
    working_list[["J_phi"]] <- sum(2*mu^{2*(P-1)} / denom)
    working_list[["J_beta_phi"]] <-t(X) %*% ({P-2}*phi*mu^{2*(P-1)} / denom)
    working_list[["J"]] <- rbind(
      cbind(working_list[["J_beta"]], working_list[["J_beta_phi"]]),
      c(working_list[["J_beta_phi"]], working_list[["J_phi"]])
    )
  } else {
    working_list[["J"]] <- working_list[["J_beta"]]
  }
  working_list[["J_star"]] <- working_list[["J"]] + working_list[["penalty"]][
    seq(working_list[["num_scored_parameters"]]),
    seq(working_list[["num_scored_parameters"]])
  ]

  # Increment to coefficients
  working_list[["scored_parameters_inc"]] <- solve(working_list[["J_star"]], working_list[["U_star"]])

  # Return the updated list
  working_list
}

update_coef <- function(working_list, step_size=1) {
  # Updates coefficients using step_size
  working_list[["scored_parameters"]] <- working_list[["scored_parameters_old"]] + step_size * working_list[["scored_parameters_inc"]]
  k <- working_list[["k"]]
  working_list[["beta"]] <- working_list[["scored_parameters"]][1:k]
  if (working_list[["phi_method"]] == "joint") {
    working_list[["phi"]] <- working_list[["scored_parameters"]][k+1]
  }
  regularized_parameters <- working_list[["beta"]]
  if (working_list[["use_regularization_for_phi"]]) {
    regularized_parameters <- c(regularized_parameters, working_list[["phi"]])
  }
  working_list[["regularized_parameters"]] <- regularized_parameters
  working_list <- update_working_estimates(working_list)
  working_list
}

check_objective_for_step_halving <- function(objective, tol) {
  objective > tol || is.null(objective) || is.na(objective)
}

gr_log_search <- function(f, lower, upper=1e4, tol=1e-8) {
  # Performs Golden-section search on the log scale
  # Adapted from Wikipedia article on the Golden-section search:
  # https://en.wikipedia.org/wiki/Golden-section_search

  # Starting values
  gr <- (sqrt(5) + 1) / 2
  a <- 0
  b <- log(upper + 1 - lower)
  g <- function(x) f(exp(x) + lower - 1)
  c <- b - (b - a) / gr
  d <- a + (b - a) / gr

  # Search until tolerance is met
  while (abs(b - a) > tol) {
    if (is.infinite(g(d)) || is.na(g(d)) || g(c) < g(d)) {
      b <- d
    } else{
      a <- c
    }
    c <- b - (b - a) / gr
    d <- a + (b - a) / gr
  }

  # Transform back to original scale and return result
  result <- exp(b) + lower - 1
  result
}

################################################################################
# Fit GP-P via Fisher scoring
gpp <- function(
  y, X, betastart, phistart, P, tol=1e-8, max_iter=100, phi_method="joint",
  stephalving_max=10, penalty=NULL, regularize_intercept=FALSE, offset=NULL,
  phi_max=1e4, verbose=FALSE) {
  # Fits the GP-P regression model using one of a variety of methods
  # Possible values for phi_method include: "fixed", "joint", "moment", "separate"

  k <- length(betastart)
  num_params <- ifelse(phi_method == "fixed", k, k+1)
  scored_parameters <- betastart
  if (phi_method  == "joint") scored_parameters <- c(scored_parameters, phistart)
  num_scored_parameters <- length(scored_parameters)
  use_regularization_for_phi <- phi_method %in% c("joint", "separate")
  consider_step_halving <- phi_method %in% c("fixed", "joint")
  regularized_parameters <- betastart
  if (use_regularization_for_phi) regularized_parameters <- c(regularized_parameters, phistart)
  num_regularized_parameters <- length(regularized_parameters)
  n <- length(y)
  if (is.null(offset)) offset <- rep(0, n)

  working_list <- list(
    y=y,
    X=X,
    n=n,
    beta=betastart,
    phi=phistart,
    scored_parameters=scored_parameters,
    num_scored_parameters=num_scored_parameters,
    regularized_parameters=regularized_parameters,
    consider_step_halving=consider_step_halving,
    offset=offset,
    P=P,
    k=k,
    num_params=num_params,
    objective=Inf,
    phi_method=phi_method,
    use_regularization_for_phi=use_regularization_for_phi
  )

  if (is.null(penalty)) {
    working_list$penalty <- matrix(0, nrow=num_regularized_parameters, ncol=num_regularized_parameters)
  } else if (!is.numeric(penalty)) {
    stop("penalty is not NULL and not numeric; what type is it?")
  } else if (length(penalty) == 1) {
    penalty_vec <- c(
      ifelse(regularize_intercept, penalty, 0),
      rep(penalty, num_regularized_parameters-1)
    )
    working_list$penalty <- diag(penalty_vec)
  } else if (is.vector(penalty)) {
    working_list$penalty <- diag(penalty)
  } else if (is.matrix(penalty)) {
    working_list$penalty <- penalty
  } else {
    stop("penalty must be NULL, numeric vector, or numeric matrix")
  }

  iters <- 0

	while ( abs(working_list[["objective"]]) > tol && iters < max_iter) {
	  iters <- iters + 1
	  if (verbose) cat("Iteration", iters, "\n")

	  # Update working_list in steps
    working_list <- update_working_estimates(working_list)
    working_list <- update_old_params(working_list)
    working_list <- update_score_information(working_list)
    if (verbose) cat("Loss", working_list[["loss"]], "\n")
    working_list <- update_coef(working_list)

    # Step halving
    do_step_halving <- FALSE
    if (working_list[["consider_step_halving"]]) {
      do_step_halving <- check_objective_for_step_halving(
        working_list$objective, tol)
    }

    if (do_step_halving) {
      step_counter <- 1
      continue_step_halving <- TRUE
      while (continue_step_halving) { # Step-halving paper uses -tol, not 0
        # Break if we're at stephalving_max
        if (step_counter > stephalving_max) {
          warning("Hit stephalving_max; exiting loop")
          working_list <- update_coef(working_list, 0)
          break
        }

        # Update parameters
        step_size <- 2^(-step_counter)
        if (verbose) {
          cat("Step halving iteration", step_counter, "\n")
          cat("Objective: ", working_list$objective, "\n")
        }
        working_list <- update_coef(working_list, step_size)
        step_counter <- step_counter + 1
        continue_step_halving <- check_objective_for_step_halving(
          working_list$objective, tol)
      }
    } # End step halving

		# Univariate optimization of phi
    moment_equation <- function(phi) {
      num <- (working_list[["y"]] - working_list[["mu"]])^2
      den <- (1 + phi * working_list[["mu"]]^(P-1))^2 * working_list[["mu"]]
      sum(num/den) + working_list[["k"]] - working_list[["n"]]
    }
    if (phi_method == "moment") {
      working_list[["phi"]] <- bisection(moment_equation, -2^(-P), 1e2, tol=1e-8)
      working_list <- update_working_estimates(working_list)
    } else if (phi_method == "separate") {
      get_loss <- function(phi, working_list) {
        working_list[["phi"]] <- phi
        working_list <- update_working_estimates(working_list)
        working_list[["loss"]]
      }
      working_list[["phi"]] <- gr_log_search(
        function(phi) get_loss(phi, working_list),
        lower=-2^(-working_list[["P"]]),
        upper=phi_max,
        tol=tol)
      working_list <- update_working_estimates(working_list)
    }

	  if (verbose) {
	    cat("phi", working_list$phi, "\n")
	    cat("Objective: ", working_list$objective, "\n")
	  }
	} # End optimization

  # Update Fisher info again
  working_list <- update_score_information(working_list)
  result_list <- working_list[c(
    "beta",
    "phi",
    "J_beta",
    "J_phi",
    "J_beta_phi",
    "J_star",
    "loss",
    "mu",
    "P"
  )]
  result_list[["iters"]] <- iters
  result_list[["fitted_values"]] <- result_list[["mu"]]
  full_cov <- chol2inv(chol(result_list[["J_star"]]))
  if (nrow(full_cov) == k) {
    result_list[["cov_beta"]] <- full_cov
    result_list[["cov_theta"]] <- full_cov
  } else if (nrow(full_cov) == (k+1)) {
    result_list[["cov_beta"]] <- full_cov[1:k, 1:k]
    result_list[["cov_theta"]] <- full_cov
  } else {
    stop("nrow(full_cov) should be k or k+1")
  }
  
  # Calculate standard errors for fitted values
  fitted_ses <- sapply(seq(n), function(j) {
    x <- X[j,]
    result_list[["mu"]][j] * sqrt(c(t(x) %*% result_list[["cov_beta"]] %*% x))
    # NOTE: Need fitted value because we're applying delta method
  })
  
  # Create confidence intervals for fitted values
  result_list[["fitted_lower_bounds"]] <- result_list[["mu"]] - 1.96 * fitted_ses
  result_list[["fitted_upper_bounds"]] <- result_list[["mu"]] + 1.96 * fitted_ses
  result_list[["fitted_interval_widths"]] <- result_list[["fitted_upper_bounds"]] - result_list[["fitted_lower_bounds"]]
  
  # Calculated fitted standard deviations
  result_list[["var_estimates"]] <- (1 + result_list[["phi"]] * result_list[["mu"]]^(P-1))^2 * result_list[["mu"]]
  result_list[["sd_estimates"]] <- sqrt(result_list[["var_estimates"]])
  
  # Add confidence intervals for standard deviation if joint method is used
  if (phi_method == "joint") {
    # Calculate standard errors for standard deviations
    sd_ses <- sapply(seq(n), function(j) {
      x_j <- X[j,]
      mu_j <- result_list[["mu"]][j]
      phi <- result_list[["phi"]]
      grad <- rep(0, k+1)
      # Calculate the below using Wolfram Alpha: d/d\theta \sqrt{[1 + \phi * exp(x*\theta)^(p-1)]^2 * exp(x*\theta)}
      grad[1:k] <- x_j * (
        (2*(P-1)*phi*mu_j^P*(phi*mu_j^(P-1) + 1) + mu_j*(phi*mu_j^(P-1) + 1)^2) /
          (2 * sqrt(mu_j * (phi*mu_j^(P-1) + 1)^2))
      )
      # Calculate the below using Wolfram Alpha: d/d\phi \sqrt{[1 + \phi * exp(x*\theta)^(p-1)]^2 * exp(x*\theta)}
      grad[k+1] <- mu_j^P * sqrt((mu_j + phi*mu_j^P)^2 / mu_j) / (mu_j + phi*mu_j^P)
      
      # Delta method
      sqrt(c(t(grad) %*% result_list[["cov_theta"]] %*% grad))
    })

    # Create confidence intervals for standard deviations
    result_list[["sd_lower_bounds"]] <- result_list[["sd_estimates"]] - 1.96 * sd_ses
    result_list[["sd_upper_bounds"]] <- result_list[["sd_estimates"]] + 1.96 * sd_ses
    result_list[["sd_interval_widths"]] <- result_list[["sd_upper_bounds"]] - result_list[["sd_lower_bounds"]] 
  }
  
  # Return result list
	result_list
}