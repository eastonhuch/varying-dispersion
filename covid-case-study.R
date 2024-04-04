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
countries_with_incomplete_data <- c(
  "Belarus",
  "Denmark",
  "Finland",
  "Greece",
  "Hungary",
  "Moldova",
  "Romania",
  "Spain",
  "Macedonia",
  "Ukraine",
  "San Marino"
)
covid <- read.csv("./covid-data/covid.csv") %>% filter(
  region=="Europe",
  !(country_name %in% countries_with_incomplete_data), # Data issues
  # sub_region=="Western Europe",
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
       main=c, log="y",xlab="Date", ylab="New Cases")
  Sys.sleep(0.2)
}

# Create design matrices
first_visual_date <- as.Date("2022-05-01")
in_visuals <- covid$date >= first_visual_date
covid_visuals <- covid[in_visuals,]

first_test_date <- as.Date("2022-08-25")
is_train <- covid$date < first_test_date
is_test <- !is_train
covid_train <- covid[is_train,]

X_formula <- new_confirmed ~
  as.factor(location_key) +
  as.factor(location_key)*ns(days_std, 60) +
  as.factor(location_key)*as.factor(weekdays(date))
lm_X <- lm(X_formula, data=covid_train)
X <- model.matrix(lm_X, data=covid_train)
X_visuals <- model.matrix(lm_X, data=covid_visuals)

Z_formula <- new_confirmed ~
  as.factor(location_key) +
  as.factor(location_key)*ns(days_std, 20) +
  as.factor(location_key)*as.factor(weekdays(date))
lm_Z <- lm(Z_formula, data=covid_train)
Z <- model.matrix(lm_Z, data=covid_train)
Z_visuals <- model.matrix(lm_Z, data=covid_visuals)

y <- covid_train$new_confirmed
y_visuals <- covid_visuals$new_confirmed

quasi_start_time <- Sys.time()
mod_quasipois <- glm(y ~ 0 + X, family=quasipoisson())
quasi_end_time <- Sys.time()
quasi_elapsed_time <- quasi_end_time - quasi_start_time
print(quasi_elapsed_time)

beta_start <- mod_quasipois$coefficients
alpha_start <- rep(0, ncol(Z))
names(alpha_start) <- colnames(Z)
alpha_start[1] <- log(summary(mod_quasipois)$dispersion)
rm(mod_quasipois)
gc()
max_iter <- 100
tol <- 1e-4  # The fit changes very little after this point
stephalving_max <- 50
mod_epl <- epl(y, X, Z, betastart = beta_start, alphastart = alpha_start, Xnew=X_visuals, Znew=Z_visuals,
               tol=tol, max_iter=max_iter, stephalving_maxiter=stephalving_max, verbose=TRUE)
cat("Fit time: ", mod_epl$fit_time, "\n")
cat("Inference time: ", mod_epl$inference_time, "\n")
cat("Elapsed time: ", mod_epl$elapsed_time, "\n")
gc()

# The EM algorithm works great!
beta_start_dln <- rep(0, ncol(X))
beta_start_dln[1] <- log(mean(y))
alpha_start_dln <- rep(0, ncol(Z))
alpha_start_dln[1] <- log(sd(y) / mean(y))

mod_dln_em <- dln(y, X, Z, betastart = beta_start_dln, alphastart = alpha_start_dln, method="EM",
                        pred_interval_method="Asymp. Bayes", Xnew=X_visuals, Znew=Z_visuals,
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=TRUE)
cat("Fit time: ", mod_dln_em$fit_time, "\n")
cat("Inference time: ", mod_dln_em$inference_time, "\n")
cat("Elapsed time: ", mod_dln_em$elapsed_time, "\n")
gc()

# This one is prone to getting stuck
# mod_dln_newton <- dln(y, X, Z,
#                       betastart = mod_dln_em$beta + rnorm(length(mod_dln_em$beta), sd=1e-2),
#                       alphastart = mod_dln_em$alpha + rnorm(length(mod_dln_em$alpha), sd=1e-2),
#                       method="Newton", pred_interval_method="Asymp. Bayes", Xnew=X_visuals, Znew=Z_visuals,
#                       max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=1e-8, verbose=TRUE)
# cat("Fit time: ", mod_dln_em$fit_time, "\n")
# cat("Inference time: ", mod_dln_em$inference_time, "\n")
# cat("Elapsed time: ", mod_dln_em$elapsed_time, "\n")

# Let's look at model fit
for (c in countries) {
  par(mfrow=c(1,2))
  c_idx <- covid_visuals$country_name == c
  c_dates <- covid_visuals$date[c_idx]
  
  # EPL
  lower <- mod_epl$fitted_values[c_idx] - 2*mod_epl$sd_estimates[c_idx]
  upper <- mod_epl$fitted_values[c_idx] + 2*mod_epl$sd_estimates[c_idx]
  lower_rev_upper <- c(lower, rev(upper))
  plot(c_dates, covid_visuals$new_confirmed[c_idx], type="n", ylim=range(lower_rev_upper),
       main=paste("EPL:", c), xlab="Date", ylab="New Cases", cex=0.2)
  abline(v=first_test_date, col=1, lty=2)
  text(first_test_date-10, max(upper), "Train")
  text(first_test_date+9, max(upper), "Test")
  polygon(
    c(c_dates, rev(c_dates)),
    c(lower, rev(upper)),
    density=400,
    col="lightgray"
  )
  lines(c_dates, covid_visuals$new_confirmed[c_idx], type="l")
  
  # DLN
  lower <- mod_dln_em$pred_lower_bounds[c_idx]
  upper <- mod_dln_em$pred_upper_bounds[c_idx]
  lower_rev_upper <- c(lower, rev(upper))
  plot(c_dates, covid_visuals$new_confirmed[c_idx], type="n", ylim=range(lower_rev_upper),
       main=paste("DLN:", c), xlab="Date", ylab="New Cases", cex=0.2)
  abline(v=first_test_date, col=1, lty=2)
  text(first_test_date-10, max(upper), "Train")
  text(first_test_date+9, max(upper), "Test")
  polygon(
    c(c_dates, rev(c_dates)),
    c(lower, rev(upper)),
    density=400,
    col="lightgray"
  )
  lines(c_dates, covid_visuals$new_confirmed[c_idx], type="l")
  
  # Wait
  Sys.sleep(0.2)
}

# How does the out-of-sample coverage look?
covid_visuals$y_pred_lower <- mod_dln_em$pred_lower_bounds
covid_visuals$y_pred_upper <- mod_dln_em$pred_upper_bounds
covid_visuals$is_test <- covid_visuals$date >= first_test_date
covid_test <- covid_visuals[covid_visuals$is_test,]
covid_test$y_covered <- (covid_test$y_pred_lower <= covid_test$new_confirmed) &
  (covid_test$new_confirmed <= covid_test$y_pred_upper)
coverage_by_country <- covid_test %>%
  group_by(country_name) %>%
  summarise(coverage=mean(y_covered))
marginal_coverage <- mean(coverage_by_country$coverage)
marginal_coverage_se <- sd(coverage_by_country$coverage) / sqrt(nrow(coverage_by_country))
cat(
  "Marginal Coverage: ",
  percent(marginal_coverage, accuracy = 0.1),
  " (", percent(marginal_coverage_se, accuracy = 0.1), ")",
  sep="")