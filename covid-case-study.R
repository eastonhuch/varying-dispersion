# COVID Case Study
require(tidyverse)
require(xtable)
# require(mpcmp)
# require(cmp)
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
  date >= "2020-07-01",
  date <  "2022-09-01"
)
unique(covid$country_name)
covid$date <- as.Date(covid$date)
covid$days_std <- scale(as.numeric(covid$date))

# Look at some example time series
countries <- sort(unique(covid$country_name))
location_keys <- sort(unique(covid$location_key))
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
  as.factor(paste0(location_key, ":", weekdays(date))) +
  as.factor(location_key):ns(days_std, 50)
X_mod <- lm(X_formula, data=covid_train)
X <- model.matrix(X_mod, data=covid_train)
X_visuals <- model.matrix(X_mod, data=covid_visuals)

Z_formula <- new_confirmed ~ 
  as.factor(paste0(location_key, ":", weekdays(date))) +
  as.factor(location_key):ns(days_std, 12)
Z_mod <- lm(Z_formula, data=covid_train)
Z <- model.matrix(Z_mod, data=covid_train)
Z_visuals <- model.matrix(Z_mod, data=covid_visuals)

y <- covid_train$new_confirmed
y_visuals <- covid_visuals$new_confirmed

# Fit Quasipoisson model
quasi_start_time <- Sys.time()
mod_quasipois <- glm(y ~ 0 + X, family=quasipoisson())
mean(is.na(mod_quasipois$coefficients))
mod_quasipois$coefficients[is.na(mod_quasipois$coefficients)]
quasi_end_time <- Sys.time()
quasi_elapsed_time <- quasi_end_time - quasi_start_time
cat("Elapsed time: "); print(quasi_elapsed_time)

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
cat("Fit time: "); print(mod_epl$fit_time)
cat("Inference time: "); print(mod_epl$inference_time)
cat("Elapsed time: "); print(mod_epl$elapsed_time)
gc()

# The EM algorithm works great!
beta_start_dln <- rep(0, ncol(X))
beta_start_dln[1] <- log(mean(y))
alpha_start_dln <- rep(0, ncol(Z))
alpha_start_dln[1] <- log(sd(y) / mean(y))

mod_dln_em <- dln(y, X, Z, betastart = beta_start_dln, alphastart = alpha_start_dln, method="EM",
                        pred_interval_method="Asymp. Bayes", Xnew=X_visuals, Znew=Z_visuals,
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=TRUE)
cat("Fit time: "); (mod_dln_em$fit_time)
cat("Inference time: "); print(mod_dln_em$inference_time)
cat("Inference time: "); print(mod_dln_em$pred_interval_time)
cat("Elapsed time: "); print(mod_dln_em$elapsed_time)
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
plot_country_fits <- function(country_vec) {
  for (c in country_vec) {
    pdf(paste0("./figures/country/", gsub(" ", "-", c), ".pdf") , width=4, height=4)
    par(mai=c(0.5, 0.9, 0.6, 0.1), cex.lab=1.4, cex.axis=1.1)
    c_idx <- covid_visuals$country_name == c
    c_dates <- covid_visuals$date[c_idx]
    
    # EPL
    # lower <- mod_epl$fitted_values[c_idx] - 2*mod_epl$sd_estimates[c_idx]
    # upper <- mod_epl$fitted_values[c_idx] + 2*mod_epl$sd_estimates[c_idx]
    # lower_rev_upper <- c(lower, rev(upper))
    # plot(c_dates, covid_visuals$new_confirmed[c_idx], type="n", ylim=range(lower_rev_upper),
    #      main=paste("EPL:", c), xlab="Date", ylab="New Cases", cex=0.2)
    # abline(v=first_test_date, col=1, lty=2)
    # text(first_test_date-5, max(upper), adj=c(1, 1), "Train-Test\nCutoff")
    # polygon(
    #   c(c_dates, rev(c_dates)),
    #   c(lower, rev(upper)),
    #   density=400,
    #   col="lightgray"
    # )
    # lines(c_dates, covid_visuals$new_confirmed[c_idx], type="l")
    
    # DLN
    lower <- mod_dln_em$pred_lower_bounds[c_idx]
    upper <- mod_dln_em$pred_upper_bounds[c_idx]
    lower_rev_upper <- c(lower, rev(upper))
    plot(c_dates, covid_visuals$new_confirmed[c_idx], type="n", ylim=range(lower_rev_upper),
         main="", xlab="", ylab="New Cases", cex=0.2)
    # text(first_test_date-5, max(upper), adj=c(1, 1), "Train-Test\nCutoff")
    polygon(
      c(c_dates, rev(c_dates)),
      c(lower, rev(upper)),
      density=400,
      col="lightgray"
    )
    lines(c_dates, covid_visuals$new_confirmed[c_idx], type="l")
    abline(v=first_test_date, col=1, lty=2)
    dev.off()
  }
}
plot_country_fits(countries)

# Save fitted values/preds to dataframe
covid_visuals$epl_fitted_values <- mod_epl$fitted_values
covid_visuals$epl_var_estimates <- mod_epl$var_estimates
covid_visuals$dln_fitted_values <- mod_dln_em$fitted_values
covid_visuals$dln_var_estimates <- mod_dln_em$sd_estimates^2
covid_visuals$y_pred_lower <- mod_dln_em$pred_lower_bounds
covid_visuals$y_pred_upper <- mod_dln_em$pred_upper_bounds
covid_visuals$is_test <- covid_visuals$date >= first_test_date
covid_test <- covid_visuals[covid_visuals$is_test,]
covid_test$y_true <- covid_test$new_confirmed

# How does the out-of-sample coverage look?
covid_test$y_covered <- (covid_test$y_pred_lower <= covid_test$y_true) &
  (covid_test$y_true <= covid_test$y_pred_upper)
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

# Comparison of EPL vs. DLN across some predictive metrics
country_bias_mae <- covid_test %>%
  mutate(
    epl_error = epl_fitted_values - y_true,
    dln_error = dln_fitted_values - y_true
  ) %>%
  group_by(country_name) %>%
  summarise(
    epl_bias=mean(epl_error),
    epl_mae=mean(abs(epl_error)),
    dln_bias=mean(dln_error),
    dln_mae=mean(abs(dln_error))
  )

avg_by_country <- covid_visuals %>%
  filter(!is_test) %>%
  group_by(country_name) %>%
  summarise(mean_new_confirmed=mean(new_confirmed))
country_bias_mae_percents <- inner_join(country_bias_mae, avg_by_country, by="country_name") %>%
  mutate(
    epl_mae_percent = epl_mae / mean_new_confirmed,
    epl_bias_percent = epl_bias / mean_new_confirmed,
    dln_mae_percent = dln_mae / mean_new_confirmed,
    dln_bias_percent = dln_bias / mean_new_confirmed,
    diff_mae_percent = dln_mae_percent - epl_mae_percent,
    diff_bias_percent = dln_bias_percent - epl_bias_percent
  )

pdf("./figures/case-mae.pdf", width=4, height=4)
par(mai=c(0.8, 0.9, 0.1, 0.1), cex.lab=1.4, cex.axis=1.1)
boxplot(country_bias_mae_percents$epl_mae_percent, country_bias_mae_percents$dln_mae_percent,
        xlab="Method", ylab="% MAE", names=c("EPL", "DLN"), yaxt="n"
)
yaxt_mae <- seq(0, 1.25, 0.25)
axis(2, at=yaxt_mae, lab=percent(yaxt_mae))
dev.off()

pdf("./figures/case-bias.pdf", width=4, height=4)
par(mai=c(0.8, 0.9, 0.1, 0.1), cex.lab=1.4, cex.axis=1.1)
boxplot(country_bias_mae_percents$epl_bias_percent, country_bias_mae_percents$dln_bias_percent,
        xlab="Method", ylab="% Bias", names=c("EPL", "DLN"), yaxt="n"
)
yaxt_bias <- seq(-0.25, 0.75, 0.25)
axis(2, at=yaxt_bias, lab=percent(yaxt_bias))
dev.off()

pdf("./figures/case-diff.pdf", width=4, height=4)
par(mai=c(0.8, 0.9, 0.1, 0.1), cex.lab=1.4, cex.axis=1.1)
boxplot(country_bias_mae_percents$diff_mae_percent, country_bias_mae_percents$diff_bias_percent,
        xlab="Performance Metric", ylab="Difference (DLN - EPL)", names=c("% MAE", "% Bias"), yaxt="n"
)
yaxt_diff <- seq(-0.15, 0.3, 0.15)
axis(2, at=yaxt_diff, lab=percent(yaxt_diff))
dev.off()

# Overall metrics
epl_bias_percent <- mean(country_bias_mae_percents$epl_bias_percent)
epl_bias_percent_se <- sd(country_bias_mae_percents$epl_bias_percent) / nrow(country_bias_mae_percents)
dln_bias_percent <- mean(country_bias_mae_percents$dln_bias_percent)
dln_bias_percent_se <- sd(country_bias_mae_percents$dln_bias_percent) / nrow(country_bias_mae_percents)
diff_bias_percent <- mean(country_bias_mae_percents$diff_bias_percent)
diff_bias_percent_se <- sd(country_bias_mae_percents$diff_bias_percent) / nrow(country_bias_mae_percents)

epl_mae_percent <- mean(country_bias_mae_percents$epl_mae_percent)
epl_mae_percent_se <- sd(country_bias_mae_percents$epl_mae_percent) / nrow(country_bias_mae_percents)
dln_mae_percent <- mean(country_bias_mae_percents$dln_mae_percent)
dln_mae_percent_se <- sd(country_bias_mae_percents$dln_mae_percent) / nrow(country_bias_mae_percents)
diff_mae_percent <- mean(country_bias_mae_percents$diff_mae_percent)
diff_mae_percent_se <- sd(country_bias_mae_percents$diff_mae_percent) / nrow(country_bias_mae_percents)

percent_accuracy <- 0.1
make_result_vec <- function(b, b_se, m, m_se, fit_time="N/A", acc=0.1) {
  if (!is.character(fit_time)) {
    fit_time <- sprintf("%.2f", as.numeric(fit_time))
  }
  c(
    paste0(percent(b, accuracy=acc), " (", percent(b_se, accuracy=acc), ")"),
    paste0(percent(m, accuracy=acc), " (", percent(m_se, accuracy=acc), ")"),
    fit_time
  )
}
epl_result_vec <- make_result_vec(
  epl_bias_percent, epl_bias_percent_se,
  epl_mae_percent, epl_mae_percent_se,
  mod_epl$elapsed_time)
dln_result_vec <- make_result_vec(
  dln_bias_percent, dln_bias_percent_se,
  dln_mae_percent, dln_mae_percent_se,
  mod_dln_em$elapsed_time)
diff_result_vec <- make_result_vec(
  diff_bias_percent, diff_bias_percent_se,
  diff_mae_percent, diff_mae_percent_se)

result_table <- data.frame(
  EPL=epl_result_vec,
  DLN=dln_result_vec,
  Difference=diff_result_vec
)
row.names(result_table) <- c("Bias %", "MAE %", "Elapsed Time (Minutes)")

xtable(
  result_table,
  align="lrrr",
  label="tab:pred_case",
  caption="Performance comparison of the extended pseudo-likelihood (EPL) and discrete log-normal (DLN) models in predicting COVID case counts for 27 European countries. The DLN method exhibits slightly more bias, lower MAE, and a shorter model-fitting time. Figure \ref{fig:case_performance} plots the country-level data."
) %>%
  print()