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
last_date <- as.Date("2022-08-25")
covid <- read.csv("./covid-data/covid.csv") %>% filter(
  region=="Europe",
  !(country_name %in% countries_with_incomplete_data), # Data issues
  # country_name %in% c("Croatia", "Lithuania", "France"), # For testing
  # sub_region=="Western Europe",
  date >= "2020-07-01",
  date <  last_date
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

# Create dates
first_visual_date <- as.Date("2022-05-01")
in_visuals <- covid$date >= first_visual_date
covid_visuals <- covid[in_visuals,]
first_test_date <- as.Date("2022-08-11") # Try making this earlier to see what effect it has
in_test <- covid$date >= first_test_date
avg_by_loc <- covid[!in_test,] %>%
  group_by(location_key) %>%
  summarise(mean_new_confirmed = mean(new_confirmed))
covid_with_avg <- inner_join(covid, avg_by_loc, by="location_key")
covid_test <- covid_with_avg[in_test,]

# Perform cross validation
x_dfs <- seq(30, 80, 10)
z_dfs <- seq(4, 20, 4)
# x_dfs <- seq(40, 60, 10) # For fast testing
# z_dfs <- seq(8, 10, 2)
dfs <- expand.grid(x_dfs, z_dfs)
colnames(dfs) <- c("x", "z")
df_row <- dfs[1,]
cv <- function(df_row) {
  print(df_row)
  x_df <- df_row["x"]
  z_df <- df_row["z"]
  first_validation_date <- first_test_date
  last_validation_date <- first_validation_date + 6
  while (last_validation_date < last_date) {
    # Subset data
    in_train <- covid$date < first_validation_date
    in_val <- (covid_test$date >= first_validation_date) & (covid_test$date <= last_validation_date)
    covid_train <- covid[in_train,]
    # covid_val <- covid_test[in_val,]
    
    for (c in countries) {
      # Subset dataframes
      covid_train_c <- filter(covid_train, country_name == c)
      y <- covid_train_c$new_confirmed
      in_val_c <- in_val & (covid_test$country_name == c)
      covid_val_c <- covid_test[in_val_c,]
      
      # Create design matrices
      X_formula <- new_confirmed ~ as.factor(weekdays(date)) + ns(days_std, x_df)
      X_mod <- lm(X_formula, data=covid_train_c)
      X <- model.matrix(X_mod, data=covid_train_c)
      X_val <- model.matrix(X_mod, data=covid_val_c)
      
      Z_formula <- new_confirmed ~ as.factor(weekdays(date)) + ns(days_std, z_df)
      Z_mod <- lm(Z_formula, data=covid_train_c)
      Z <- model.matrix(Z_mod, data=covid_train_c)
      Z_val <- model.matrix(Z_mod, data=covid_val_c)
      
      # Fit models
      mod_quasipois <- glm(y ~ 0 + X, family=quasipoisson())
      beta_start <- mod_quasipois$coefficients
      alpha_start <- rep(0, ncol(Z))
      names(alpha_start) <- colnames(Z)
      alpha_start[1] <- log(summary(mod_quasipois)$dispersion)
      
      max_iter <- 100
      tol <- 1e-8
      stephalving_max <- 50
      mod_epl <- epl(y, X, Z, betastart = beta_start, alphastart = alpha_start, Xnew=X_val, Znew=Z_val,
                     tol=tol, max_iter=max_iter, stephalving_maxiter=stephalving_max, verbose=FALSE)

      beta_start_dln <- rep(0, ncol(X))
      beta_start_dln[1] <- log(mean(y))
      alpha_start_dln <- rep(0, ncol(Z))
      alpha_start_dln[1] <- log(sd(y) / mean(y))
      mod_dln_em <- dln(y, X, Z, betastart=beta_start_dln, alphastart=alpha_start_dln, method="EM",
                        pred_interval_method="Asymp. Bayes", Xnew=X_val, Znew=Z_val,
                        max_iter=max_iter, stephalving_maxiter=stephalving_max, tol=tol, verbose=FALSE)      
      
      # Add predictions to dataframe
      covid_test[in_val_c, "epl_pred"] <- mod_epl$fitted_values
      covid_test[in_val_c, "dln_pred"] <- mod_dln_em$fitted_values
    }

    # Increment dates
    first_validation_date <- first_validation_date + 7
    last_validation_date <- first_validation_date + 6
  }
  
  # Check for invalid predictions
  # These can occur due to overflow issues when computing the fitted values for the DLN
  epl_pred_is_valid <- is.finite(covid_test$epl_pred)
  num_invalid_epl_preds <- sum(!epl_pred_is_valid)
  if (num_invalid_epl_preds > 0) {
    print(paste(num_invalid_epl_preds, "invalid EPL predictions"))
  }
  dln_pred_is_valid <- is.finite(covid_test$dln_pred)
  num_invalid_dln_preds <- sum(!dln_pred_is_valid)
  if (num_invalid_dln_preds > 0) {
    print(paste(num_invalid_dln_preds, "invalid DLN predictions"))
  }
  
  # Compute standardized mean absolute errors (SMAEs) and return them
  epl_smae <- mean((abs(covid_test$new_confirmed - covid_test$epl_pred) / covid_test$mean_new_confirmed)[epl_pred_is_valid])
  dln_smae <- mean((abs(covid_test$new_confirmed - covid_test$dln_pred) / covid_test$mean_new_confirmed)[dln_pred_is_valid])
  c(epl=epl_smae, dln=dln_smae)
}

# Run cross-validation
cv_results <- apply(dfs, 1, cv)

# Visualize CV results
heatmap_data <- as.data.frame(dfs)
heatmap_data$epl_smae <- cv_results[1,]
heatmap_data$dln_smae <- cv_results[2,]
ggplot(heatmap_data, aes(x, z, fill=log(epl_smae))) + geom_tile()
ggplot(heatmap_data, aes(x, z, fill=log(dln_smae))) + geom_tile()
plot(heatmap_data$x, heatmap_data$epl_smae, log="y") # df=50 for X
plot(heatmap_data$x, heatmap_data$dln_smae, log="y") # Same for DLN
heatmap_data_x_50 <- filter(heatmap_data, x==50)
heatmap_data_x_50_vals <- unlist(heatmap_data_x_50[, c("epl_smae", "dln_smae")])
plot(heatmap_data_x_50$z, heatmap_data_x_50$epl_smae, log="y", type="l", ylim=range(heatmap_data_x_50_vals),
     xlab="DF for Z", ylab="SMAE")
lines(heatmap_data_x_50$z, heatmap_data_x_50$dln_smae, col=2)

# Final model
# Degrees of freedom for X: 50
# Degrees of freedom for Z: 12