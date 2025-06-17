#######################
# Plot of data
#######################
require(mpcmp) #install_github("thomas-fung/mpcmp")
require(cmp) #install_github("SuneelChatla/cmp")
require(scales)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

myfilepath <- "/path/to/simulation_study/Figs/"

library(ggplot2)

n <- 100L
set.seed(2)
x1 <- rnorm(n)  # We'll keep the design fixed across simulations
x2 <- rnorm(n)

# NOTE: The mpcmp package can be a bit temperamental
# You might run into errors if you make the coefficients much larger
beta_true <- c(3, 0.05, -0.1, 0.02) # Regression coefficients for the mean
y <- rep(0, n)  # We'll create this later
master_formula <- y ~ x1*x2
X <- model.matrix(master_formula)
log_mu <- c(X %*% beta_true)
mu <- exp(log_mu)
alpha_true <- c(0.1, 0, -0.2, 0.05) # Not directly comparable across all models
Z <- model.matrix(master_formula)
log_nu <- c(Z %*% alpha_true)
nu <- exp(log_nu)

lambda <- comp_lambdas(mu, nu)$lambda  # Helper from mpcmp
# This calls a C function from the cmp package that efficiently calculates the moments
var_true <- .C(
  "cmp_cum_all", llambda=as.double(log(lambda)), lnu=as.double(log_nu),
  flag=as.double(1L), n=n, mean=rep(0, n), var=rep(0, n))$var
sd_true <- sqrt(var_true)  # True standard deviation that we're trying to estimate

y <- mpcmp::rcomp(n, mu=mu, nu=nu)
dat <- data.frame(x1, x2, y, mu, var_true, sd_true, overind = (var_true/mu > 1))

mycols <- gray(c(.1, .5)) #c("#d8b365", "#5ab4ac")

ggplot(data=dat, aes(x1,y)) +geom_point(aes(color=overind)) + labs(color=("Overdispersed")) + scale_color_manual(values = mycols)
ggsave(filename=paste0(myfilepath, "x1y.pdf"), width=4, height=3, units="in")
ggplot(data=dat, aes(x2,y))+geom_point(aes(color=overind)) + labs(color=("Overdispersed")) + scale_color_manual(values = mycols)
ggsave(filename=paste0(myfilepath, "x2y.pdf"), width=4, height=3, units="in")
ggplot(data=dat, aes(mu,var_true))+geom_point(aes(color=overind)) + labs(color=("Overdispersed"), y="variance", x="mean") + scale_color_manual(values = mycols)
ggsave(filename=paste0(myfilepath, "disp.pdf"), width=4, height=3, units="in")