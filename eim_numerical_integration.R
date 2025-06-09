# Libraries and directory
require(lattice)
require(latticeExtra)
require(scales)
setwd("/path/to/varying-dispersion")
source("./gpp.R")

# Set simulation parameters
P <- 1 # Can adjust this but might need to increase max_y
mu_max <- 10
mu_inc <- 0.1
phi_inc <- 2^-P/20
mus <- seq(mu_inc, mu_max, mu_inc)
phis <- seq(-2^-P, 2^-P, phi_inc)
param_grid <- expand.grid(mus, phis)
colnames(param_grid) <- c("mu", "phi")
den <- 1 + param_grid$phi * param_grid$mu^(P-1)
thetas <- param_grid$mu / den
deltas <- param_grid$phi * param_grid$mu^(P-1) / den
param_grid$is_viable <- check_gp_parameter_space(thetas, deltas)
# We get a warning because some of these combinations are invalid

# Viable parameter values are shown in blue
levelplot(is_viable ~ mu * phi, data=param_grid,
          main="Viable Parameter Choices",
          xlab=expression(mu), ylab=expression(varphi))

# Function for calculating analytical EIM
get_analytical_eim <- function(mu, phi) {
  multiplier <- 1 / (1 + phi * mu^(P-1))^2
  den <- 1 + 2*phi*mu^(P-2)
  J_11 <- 1/mu + 2*((P-2) * phi * mu^(P-2))^2 / den
  J_12 <- 2*(P-2)*phi*mu^(2*P-3) / den
  J_22 <- 2*mu^(2*(P-1)) / den
  J <- multiplier * matrix(c(J_11, J_12, J_12, J_22), nrow=2, ncol=2)
  J
}

# Function for numerically approximating EIM
get_numerical_eim <- function(mu, phi, max_y) {
  y <- seq(0, max_y)
  
  # Score for mu
  temp <- 1 + phi*mu^(P-1)
  score_mu <- 1/mu +
    (y-1) *
    (1 + {P-1}*phi*y*mu^{P-2}) /
    (mu + phi*y*mu^{P-1}) -
    (1 + 2*{P-1}*phi*y*mu^{P-2}) /
    temp +
    ({(P-1)*phi*mu^(P-1)} * {1 + phi*y*mu^(P-2)}) /
    (temp^2)
  
  # Score for phi
  score_phi <- mu^{(P-2)}*y*(y-1) / {1 + phi*(mu^{P-2})*y} -
    mu^{(P-1)}*{y-mu} / {temp^2} -
    mu^{(P-1)}*y / temp
  
  # Probabilities
  probs <- dgpoisP(y, mu, phi, P, log=FALSE)
  use_entry <- probs > 0 & (!is.na(probs)) & is.finite(probs) &
    (!is.na(score_mu)) & is.finite(score_mu) &
    (!is.na(score_phi)) & is.finite(score_phi)
  probs <- probs[use_entry]
  
  # Join scores and estimate EIM
  scores <- cbind(score_mu, score_phi)[use_entry,]
  mean_score <- as.numeric(colSums(probs * scores))
  weighted_scores <- sqrt(probs) * scores
  numerical_eim <- crossprod(weighted_scores) - tcrossprod(mean_score)
  # NOTE: The assumptions required for the last term to be zero are not met
  # We include an estimate of it for completeness, but is has limited influence
  # on the results
  
  numerical_eim
}

# Function for calculating EIM error
# NOTE: You may need to increase max_y for some values of P
get_eim_error <- function(params, max_y=1000) {
  # params is a 2-d vector whose first element is mu and second, phi
  mu <- params[1]
  phi <- params[2]
  analytical_eim <- get_analytical_eim(mu, phi)
  numerical_eim <- get_numerical_eim(mu, phi, max_y)
  l1_error <- sum(abs(analytical_eim - numerical_eim))
  l1_numerical_eim <- sum(abs(numerical_eim))
  l1_error / l1_numerical_eim
}

# Calculate errors
param_grid$eim_error <- apply(param_grid[, c("mu", "phi")], 1, get_eim_error)
param_grid$eim_error[!param_grid$is_viable] <- NA
param_grid$percent_eim_error <- param_grid$eim_error * 100
param_grid$log_eim_error <- log(param_grid$eim_error)

# Final plot of error
pdf(paste0("./figures/eim_error.pdf"), width=6, height=4)
finite_elements <- function(x) {
  include <- (!is.na(x)) & is.finite(x)
  x[include]
}
min_log_eim_error <- min(finite_elements(param_grid$log_eim_error))
max_log_eim_error <- max(finite_elements(param_grid$log_eim_error))
labels <- seq(min_log_eim_error, max_log_eim_error, length.out=15)
if (P == 1) labels <- seq(log(1e-17), 0, length.out=18)
labels_transformed <- format(exp(labels), digits=2)
levelplot(
  log_eim_error ~ mu * phi,
  data=param_grid,
  at=labels,
  colorkey=list(
    at=labels,
    labels=labels_transformed),
  main=paste0("EIM Error for the GP-", P, " Model"),
  xlab=expression(mu), ylab=expression(varphi))
dev.off()
