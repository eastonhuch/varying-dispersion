# Test gpp on Twitter data
setwd("~/Documents/research/GPP/gpp/")
require(tidyverse)
require(splines)
require(mvtnorm)
require(microbenchmark)
require(COMPoissonReg)
source("./gpp.R")

# Transform data
twitter_dat <- read.csv("../twitter-data/twcs/twcs.csv")
twitter_dat$is_apple <- grepl("@AppleSupport", twitter_dat$text, ignore.case=TRUE)
apple_timestamps <- twitter_dat %>%
  filter(is_apple) %>%
  pull(created_at) %>%
  strptime(format="%a %b %e %H:%M:%S %z %Y", tz="UTC")
apple_rounded_timestamps <- format(round(apple_timestamps, units="hours"), format="%Y-%m-%d %H:%M", tz="UTC")
apple_rounded_timestamp_counts <- c(table(apple_rounded_timestamps))
timestamp_df <- data.frame(
  time=strptime(names(apple_rounded_timestamp_counts), format="%Y-%m-%d %H:%M", tz="UTC"),
  count=apple_rounded_timestamp_counts
)
row.names(timestamp_df) <- 1:nrow(timestamp_df)
timestamp_df <- timestamp_df %>% 
  filter(time >= strptime("2017-10-05 00:00", format="%Y-%m-%d %H:%M", tz="UTC")) %>%
  arrange(time)
  
# Visualize data
plot(timestamp_df$time, timestamp_df$count, type="l",
     xlab="Hour", ylab="Num Tweets")

# Apply gpp to full data set
timestamp_df$day <- as.Date(timestamp_df$time)
timestamp_df$hour <- as.numeric(format(timestamp_df$time, format="%H", tz="UTC"))
timestamp_df$hour_0_2pi <- timestamp_df$hour / 23 * 2 * pi
timestamp_df$weekday <- format(timestamp_df$time, format="%A", tz="UTC")
timestamp_df$time_numeric <- c(scale(as.numeric(timestamp_df$time)))
x_formula <- count ~
  weekday +
  ns(time_numeric, df=8) +
  I(sin(hour_0_2pi)) +
  I(cos(hour_0_2pi)) +
  I(sin(hour_0_2pi/2)) +
  I(cos(hour_0_2pi/2)) +
  I(sin(hour_0_2pi/3)) +
  I(cos(hour_0_2pi/3))
X <- model.matrix(x_formula, data=timestamp_df)
mod_pois <- glm(timestamp_df$count ~ 0 + X, family=poisson)
beta_init <- mod_pois$coefficients
phi_init <- 0
tol <- 1e-8
max_iter <- 100
phi_method <- "joint"
stephalving_max <- 100

gp1 <- gpp(
  timestamp_df$count, X, betastart=beta_init, phistart=phi_init,
  P=1, tol=tol, max_iter=max_iter, phi_method=phi_method,
  verbose=T, stephalving_max=stephalving_max)
gp1_log_like <- -gp1$loss

plot(timestamp_df$time, timestamp_df$count, type="l",
     xlab="Hour", ylab="Num Tweets")
fitted_values <- exp(X %*% gp1$beta)
lines(timestamp_df$time, fitted_values, col=2)
legend("topleft", legend=c("Data", "Fitted Values"), lty=1, col=1:2)

# Compare to COM-Poisson
microbenchmark({gpp(
  timestamp_df$count, X, betastart=beta_init, phistart=phi_init,
  P=1, tol=tol, max_iter=max_iter, phi_method=phi_method,
  verbose=F, stephalving_max=stephalving_max)},
  times=10)

cmp_start_time <- Sys.time()
cmp <- glm.cmp(x_formula, data=timestamp_df, init=get.init(beta_init, 0))
cmp_end_time <- Sys.time()
cmp_elapsed_time <- cmp_end_time - cmp_start_time
cmp_elapsed_time
cmp_log_like <- logLik(cmp)

print(gp1_log_like)
print(cmp_log_like)

# Apply one day at a time to identify outliers
first_forecast_date <- as.Date(strptime("2017-10-16", format="%Y-%m-%d", tz="UTC"))
last_date <- as.Date(max(timestamp_df$time))
num_forecast_days <- as.integer(last_date - first_forecast_date)

max_x <- 500
kde_quantile <- function(x, probs) {
  length_x <- length(x)
  density_x <- density(x, n=max_x+1, from=0, to=max_x)
  masses_raw <- density_x$y
  masses <- masses_raw / sum(masses_raw)
  cum_masses <- cumsum(masses)
  quantiles <- sapply(probs, function(p) {
    below_p <- cum_masses < p
    if (p < 0.5) {
      if (sum(below_p) > 0) {
        rev(density_x$x[below_p])[1]
      } else {
        0
      }
    } else {
      if (sum(!below_p) > 0) {
        density_x$x[!below_p][1] 
      } else {
        max_x
      }
    }
  })
  quantiles
}

n_posterior_samples <- 1000
for (i in seq(0, num_forecast_days)) {
  start_date <- first_forecast_date + i
  filtered_timestamp_df <- filter(timestamp_df, time < start_date)
  mod_lm <- lm(x_formula, data=filtered_timestamp_df)
  filtered_X <- model.matrix(mod_lm, data=filtered_timestamp_df)
  p <- ncol(filtered_X)
  filtered_mod_pois <- glm(filtered_timestamp_df$count ~ 0 + filtered_X, family=poisson)
  filtered_beta_init <- filtered_mod_pois$coefficients
  filtered_gp1 <- gpp(
    filtered_timestamp_df$count, filtered_X, betastart=beta_init, phistart=phi_init,
    P=1, tol=tol, max_iter=max_iter, phi_method=phi_method,
    verbose=T, stephalving_max=stephalving_max)
  
  new_day_timestamp_df <- filter(timestamp_df, day == start_date)
  new_day_counts <- new_day_timestamp_df$count
  new_day_X <- model.matrix(mod_lm, data=new_day_timestamp_df)
  new_day_timestamp_df$pred <- exp(new_day_X %*% filtered_gp1$beta)
  sampled_params <- rmvnorm(
    n_posterior_samples,
    c(filtered_gp1$beta, filtered_gp1$phi),
    sigma=chol2inv(chol(filtered_gp1$J_star))
  )
  sampled_means <- exp(new_day_X %*% t(sampled_params[,1:p]))
  sampled_phis <- sampled_params[,p+1]
  sampled_outcomes <- matrix(
    rgpoisP(n=length(sampled_means), mu=c(sampled_means), phi=c(sampled_phis), P=1),
    nrow=nrow(sampled_means),
    ncol=ncol(sampled_means)
  )
  cis <- apply(sampled_outcomes, MARGIN=1, FUN=function(x) kde_quantile(x, probs=c(0.00001, 0.99999)))
  new_day_timestamp_df$in_interval <- (cis[1,] <= new_day_counts) & (new_day_counts <= cis[2,])
  timestamp_df[timestamp_df$day == start_date, "in_interval"] <- new_day_timestamp_df$in_interval
}

# Plot outliers at the end
timestamp_df$is_anomaly <- !timestamp_df$in_interval
plot(timestamp_df$time, timestamp_df$count, xlab="Hour", ylab="Num Tweets", type="l")
points(timestamp_df$time[timestamp_df$is_anomaly], timestamp_df$count[timestamp_df$is_anomaly], col=2, cex=0.3)
legend("topleft", legend=c("Data", "Anomalies"), lty=1, col=1:2)
