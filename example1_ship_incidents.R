# Libraries and directory
library(MASS)
library(COUNT)
setwd("/path/to/gpp")
source("./gpp.R")

# Read in and process data
dat <- read.csv("./ship_incidents.csv")
# Data is available from "Generalized Linear Models" by McCullagh & Nelder
remove_row <- dat$Necessarily.Empty.Cell | dat$Accidentally.Empty.Cell
dat <- dat[!remove_row,] # Filter down to ships in service
my_formula <- formula(
  Number.of.damage.incidents ~
  Ship.Type + 
  Year.of.construction +
  Period.of.operation)
X <- model.matrix(my_formula, data=dat)
y <- dat$Number.of.damage.incidents
offset <- log(dat$Aggregate.months.service)

# Starting values from Poisson regression
Poisfit <- glm(y~0+X, family=poisson, offset=offset)
betastart <- Poisfit$coefficients
logLik(Poisfit) # -68.28077
deviance(Poisfit)
phistart <- 0

# Control parameters for model-fitting functions
tol <- 1e-12
max_iter <- 200
phi_method <- "joint"
stephalving_max <- 200

# Fit GP-1
gp1 <- gpp(
  y, X, betastart=betastart, phistart=phistart, offset=offset, P=1, 
  tol=tol, max_iter=max_iter, phi_method=phi_method, verbose=T,
  stephalving_max=stephalving_max)
gp1_loglikelihood <- -gp1$loss # -67.65261
saturated_gp1_loglikelihood <- sum(dgpoisP(y[y>0], y[y>0], gp1$phi, 1, log=TRUE, omit_constant=FALSE)) # Saturated Loglik: -52.24568
gp1_dev <- 2 * (saturated_gp1_loglikelihood - gp1_loglikelihood)
gp1_dev # Deviance for the GP-1 model
# Full-reduced model test: Fail to reject
2 * pchisq(2 * (gp1_loglikelihood - logLik(Poisfit)), 1, lower.tail=FALSE)

# Fit model with different values of P
Ps <- seq(-2, 5, by=0.01)
loglikelihoods <- c()
avg_dispersion_parameter <- c()
phis <- c()

for (P in Ps) {
  mod <- gpp(y, X, betastart=betastart, phistart=phistart, offset=offset, P=P, tol=tol, max_iter=max_iter, phi_method=phi_method, verbose=T, stephalving_max=stephalving_max)
  loglikelihoods <- c(loglikelihoods, -mod$loss)
  avg_dispersion_parameter <- c(avg_dispersion_parameter, mean((1 + mod$phi * mod$mu^(P-1))^2))
  phis <- c(phis, mod$phi)
}
pdf("ships.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(Ps, loglikelihoods, type="l",
     main=expression(paste("(a) Log likelihood values for varying P")),
     xlab="P", ylab="Log likelihood")
plot(Ps, phis, type="l",
     main=expression(paste("(b) Estimated values of ", varphi, " for varying P")),
     xlab="P", ylab=expression(paste("Estimated value of ", varphi)))
abline(h=0, col=2)
dev.off()

# Test best model compared to GP-1
2 * pchisq(2 * (max(loglikelihoods) - gp1_loglikelihood), 1, lower.tail=FALSE) # Fail to reject

# Test best model compared to Poisson
2 * pchisq(2 * (max(loglikelihoods) - logLik(Poisfit)), 1, lower.tail=FALSE) # Fail to reject

# Best model
Ps[which.max(loglikelihoods)] # Best P is 0.26
gp_best <- gpp(y, X, betastart=betastart, phistart=phistart, offset=offset, P=0.26, tol=tol, max_iter=max_iter, phi_method=phi_method, verbose=T, stephalving_max=stephalving_max)
gp_best$phi # Estimated value of phi
1/sqrt(gp_best$J_phi) # SE of phi
J_beta_inv <- solve(gp_best$J_beta) # For SEs
data.frame(beta=gp_best$beta, se=sqrt(diag(J_beta_inv))) # Coefficient table