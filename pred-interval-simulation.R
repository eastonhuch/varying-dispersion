# devtools::install_github("SuneelChatla/cmp")
# devtools::install_github("thomas-fung/mpcmp") 
require(mpcmp)
require(cmp)
require(scales)
source("discrete-log-normal.R")

# Set parameters
set.seed(1)
beta_true <- c(4, 0.05)
alpha_true <- c(0, 0) # Setting this to zero so that we can simulate faster
master_formula <- y ~ x1
nominal_coverage <- 0.95  # This is hard-coded in the model-fitting functions

# Control parameters
tol <- 1e-8
max_iter <- 100
phi_method <- "joint"
stephalving_max <- 10

# Arrays for results
n_vals <- c(10L, 25L, 100L, 400L)
num_n_vals <- length(n_vals)
reps <- 100L
method_names <- c("Plug-in", "Approx. Bayes", "Full Bayes")
num_methods <- length(method_names)

marginal_coverages <- matrix(0, nrow=num_n_vals, ncol=num_methods, dimnames = list(n_vals, method_names))
marginal_coverages[] <- NA
marginal_coverage_ses <- marginal_coverages
coverage_rmses <- marginal_coverages
interval_widths <- marginal_coverages
interval_width_ses <- marginal_coverages
interval_width_medians <- marginal_coverages

# Loop
for (n in n_vals) {
  cat("n: ", n, "\n")
  
  # Create data
  x1 <- rnorm(n)  # We'll keep the design fixed across simulations
  y <- rep(0, n)  # We'll create this later
  X <- model.matrix(master_formula)
  log_mu <- c(X %*% beta_true)
  mu <- exp(log_mu)
  Z <- model.matrix(master_formula)
  log_nu <- c(Z %*% alpha_true)
  nu <- exp(log_nu)
  lambda <- comp_lambdas(mu, nu)$lambda  # Helper from mpcmp
  # This calls a C function from the cmp package that efficiently calculates the moments
  var_true <- .C(
    "cmp_cum_all", llambda=as.double(log(lambda)), lnu=as.double(log_nu),
    flag=as.double(1L), n=as.integer(n), mean=rep(0, n), var=rep(0, n))$var
  sd_true <- sqrt(var_true)  # True standard deviation that we're trying to estimate
  y_star_covered <- array(FALSE, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
  y_star_covered[] <- NA
  y_star_interval_widths <- array(0, dim=c(n, reps, num_methods), dimnames = list(NULL, NULL, method_names))
  y_star_interval_widths[] <- NA
  
  for (i in seq(reps)) {
    print(i)
    # y <- mpcmp::rcomp(n, mu=mu, nu=nu)
    # y_star <- mpcmp::rcomp(n, mu=mu, nu=nu)
    y <- rpois(n, mu)
    y_star <- rpois(n, mu)
    dat <- data.frame(x1, y)
    
    # Plug-in Method
    if ("Plug-in" %in% method_names ) try({
      # Fit model
      mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
      beta_start <- mod_quasipois$coefficients
      alpha_start <- rep(0, ncol(Z))
      names(alpha_start) <- colnames(Z)
      alpha_start[1] <- log(sqrt(summary(mod_quasipois)$dispersion))
      mod_plug_in <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="EM", pred_interval_method="Plug-in",
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
      
      # Store results
      y_star_covered[,i,"Plug-in"] <- (mod_plug_in$pred_lower_bounds <= y_star) & (y_star <= mod_plug_in$pred_upper_bounds)
      y_star_interval_widths[,i,"Plug-in"] <- mod_plug_in$pred_interval_widths
    })
    
    # Approx. Bayes method
    if ("Approx. Bayes" %in% method_names ) try({
      # Fit model
      mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
      beta_start <- mod_quasipois$coefficients
      alpha_start <- rep(0, ncol(Z))
      names(alpha_start) <- colnames(Z)
      alpha_start[1] <- log(sqrt(summary(mod_quasipois)$dispersion))
      mod_bayes <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="EM", pred_interval_method="Approx. Bayes",
                         max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
      
      # Store results
      y_star_covered[,i,"Approx. Bayes"] <- (mod_bayes$pred_lower_bounds <= y_star) & (y_star <= mod_bayes$pred_upper_bounds)
      y_star_interval_widths[,i,"Approx. Bayes"] <- mod_bayes$pred_interval_widths
    })
    
    # Full Bayes method
    if ("Full Bayes" %in% method_names ) try({
      # Fit model
      mod_quasipois <- glm(master_formula, family=quasipoisson(), data=dat)
      beta_start <- mod_quasipois$coefficients
      alpha_start <- rep(0, ncol(Z))
      names(alpha_start) <- colnames(Z)
      alpha_start[1] <- log(sqrt(summary(mod_quasipois)$dispersion))
      mod_bayes <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="EM", pred_interval_method="Full Bayes",
                       max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)
      
      # Store results
      y_star_covered[,i,"Full Bayes"] <- (mod_bayes$pred_lower_bounds <= y_star) & (y_star <= mod_bayes$pred_upper_bounds)
      y_star_interval_widths[,i,"Full Bayes"] <- mod_bayes$pred_interval_widths
    })
  }
  
  # Marginal coverage
  rep_coverage <- apply(y_star_covered, MARGIN=c(2,3), mean)
  marginal_coverages[as.character(n),] <- colMeans(rep_coverage)
  marginal_coverage_ses[as.character(n),] <- apply(rep_coverage, 2, sd) / sqrt(n)
  
  # Coverage rMSE
  observation_coverage <- apply(y_star_covered, MARGIN=c(1,3), mean)
  coverage_rmses[as.character(n),] <- sqrt(colMeans((observation_coverage - nominal_coverage)^2))
  
  # Interval widths
  rep_interval_widths <- apply(y_star_interval_widths, MARGIN=c(2,3), mean)
  interval_widths[as.character(n),] <- colMeans(rep_interval_widths)
  interval_width_ses[as.character(n),] <- apply(rep_interval_widths, 2, sd) / sqrt(n)
  interval_width_medians[as.character(n),] <- apply(y_star_interval_widths, MARGIN=c(3), median)
}

# Check results quickly
marginal_coverages
coverage_rmses
interval_widths
interval_width_medians

# Display results
format_float <- label_number(0.001)
format_percent <- label_percent(0.1)

marginal_coverages_text <- paste0(format_float(marginal_coverages), "(", format_float(marginal_coverage_ses), ")")
interval_widths_text <- paste0(format_float(interval_widths), "(", format_float(interval_width_ses), ")")

