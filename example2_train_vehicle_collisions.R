setwd("path/to/gpp")
source('./gpp.R')

# Read in data
dat <- read.csv("./train_vehicle_collisions.csv")
# Data is available in .xlsx format here:
# https://ceprofs.civil.tamu.edu/dlord/Highway_Safety_Analytics_and_Modeling.htm

# Data transformations
my_formula <- formula(Crashes ~ LNAADT + ADRT + PCA + TDD + PSH)
X <- model.matrix(my_formula, data=dat)
num_predictors <- ncol(X)
X[,2:num_predictors] <- scale(X[,2:num_predictors])
y <- dat$Crashes

# Starting values from Poisson regression
Poisfit <- glm(y~0+X, family=poisson)
betastart <- Poisfit$coefficients
phistart <- 0

# Control parameters for model-fitting functions
tol <- 1e-8
max_iter <- 100
phi_method <- "joint"
stephalving_max <- 100

# Fit GP-1
gp1 <- gpp(
  y, X, betastart=betastart, phistart=phistart,
  P=1, tol=tol, max_iter=max_iter, phi_method=phi_method,
  verbose=T, stephalving_max=stephalving_max)

# Regularized fit
penalty <- c(0, rep(1e2, num_predictors-1), 0)
gp1_penalized <- gpp(
  y, X, betastart=betastart_regularized, phistart=phistart,
  P=1, tol=tol, penalty=penalty, max_iter=max_iter, phi_method=phi_method, stephalving_max=stephalving_max,
  verbose=T)

# Show how it fits with different penalties
betastart_regularized <- c(log(mean(dat$Crashes)), rep(0, num_predictors-1))
penalties <- 10^seq(-2, 6, by=0.1)
loglikelihoods <- c()
avg_dispersion_parameter <- c()
phis <- c()

for (penalty in penalties) {
  penalty <- c(0, rep(penalty, num_predictors-1), 0)
  mod <- gpp(
    y, X, betastart=betastart_regularized, phistart=phistart, P=1, tol=tol,
    penalty=penalty, max_iter=max_iter, phi_method=phi_method, verbose=T, 
    stephalving_max=stephalving_max)
  loglikelihoods <- c(loglikelihoods, -mod$loss)
  avg_dispersion_parameter <- c(avg_dispersion_parameter, mean((1 + mod$phi)^2))
  phis <- c(phis, mod$phi)
}

# Plot results
pdf("collisions.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(penalties, loglikelihoods, type="l", log="x",
     main=expression("(a) Log likelihood values for varying penalty"),
     xlab="penalty", ylab="Log likelihood")
plot(penalties, phis, type="l", log="x",
     main=expression(paste("(b) Estimated values of ", varphi, " for varying penalty")),
     xlab="penalty", ylab=expression(paste("Estimated value of ", varphi)))
abline(h=0, col=2)
dev.off()

