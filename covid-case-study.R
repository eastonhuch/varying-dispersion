# COVID Case Study
require(tidyverse)
require(xtable)
require(mpcmp)
require(cmp)
require(scales)
require(splines)
source("cmp-helpers.R")
source("gpp.R")
source("moment-methods.R")
source("discrete-log-normal.R")

# Read in data
covid_raw <- read.csv("./covid-data/covid.csv")
covid <- covid_raw %>% filter(
  sub_region=="Western Europe",
  date >= "2020-07-01",
  date <  "2022-09-01"
)
unique(covid$country_name)
covid$date <- as.Date(covid$date)
covid$days_std <- scale(as.numeric(covid$date))

# Look at some example time series
countries <- sort(unique(covid$country_name))
for (c in countries) {
  covid_example <- filter(covid, country_name==c)
  plot(covid_example$date, covid_example$new_confirmed, type="l",
       main=c, xlab="Date", ylab="New Cases")
  Sys.sleep(0.2)
}

# Create design matrices
cutoff_date <- as.Date("2022-08-18")
is_train <- covid$date < cutoff_date
is_test <- !is_train
covid_train <- covid[is_train,]
X_formula <- new_confirmed ~
  as.factor(location_key) +
  as.factor(location_key)*ns(days_std, 40) +
  as.factor(location_key)*as.factor(weekdays(date))
lm_X <- lm(X_formula, data=covid_train)
X <- model.matrix(lm_X, data=covid_train)
Xnew <- model.matrix(lm_X, data=covid)
Z_formula <- new_confirmed ~
  as.factor(location_key) +
  as.factor(location_key)*ns(days_std, 10) +
  as.factor(location_key)*as.factor(weekdays(date))
lm_Z <- lm(Z_formula, data=covid_train)
Z <- model.matrix(lm_Z, data=covid_train)
Znew <- model.matrix(lm_Z, data=covid)
ynew <- covid$new_confirmed
y <- ynew[is_train]

mod_quasipois <- glm(y ~ 0 + X, family=quasipoisson())
beta_start <- mod_quasipois$coefficients
alpha_start <- rep(0, ncol(Z))
names(alpha_start) <- colnames(Z)
alpha_start[1] <- log(summary(mod_quasipois)$dispersion)
max_iter <- 100
stephalving_max <- 10
mod_epl <- epl(y, X, Z, betastart = beta_start, alphastart = alpha_start, Xnew=Xnew, Znew=Znew,
               max_iter=max_iter, stephalving_maxiter=stephalving_max, verbose=FALSE)

mod_approx_bayes <- dln(y, X, Z, betastart = beta_start, alphastart = alpha_start, method="Newton",
                        pred_interval_method="Asymp. Bayes", Xnew=Xnew, Znew=Znew,
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=FALSE)


# Let's look at model fit
for (c in countries) {
  c_idx <- (covid$country_name == c) & (covid$date >= "2022-05-01")
  c_dates <- covid$date[c_idx]
  lower <- mod_epl$fitted_values[c_idx] - 2*mod_epl$sd_estimates[c_idx]
  upper <- mod_epl$fitted_values[c_idx] + 2*mod_epl$sd_estimates[c_idx]
  lower_rev_upper <- c(lower, rev(upper))
  plot(c_dates, covid$new_confirmed[c_idx], type="n", ylim=range(lower_rev_upper),
       main=c, xlab="Date", ylab="New Cases", cex=0.2)
  abline(v=cutoff_date, col=1, lty=2)
  text(cutoff_date-10, max(upper), "Train")
  text(cutoff_date+9, max(upper), "Test")
  polygon(
    c(c_dates, rev(c_dates)),
    c(lower, rev(upper)),
    density=400,
    col="lightgray"
  )
  lines(c_dates, covid$new_confirmed[c_idx], type="l")
  Sys.sleep(1)
}
