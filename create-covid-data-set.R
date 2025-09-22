# Create COVID case study data
require(tidyverse)

# Read in and process data
# This lookup can be downloaded from lukes/ISO-3166-Countries-with-Regional-Codes on GitHub
country_lookup <- read.csv("./covid-data/country-lookup.csv") %>% select(
  name,
  alpha.2,
  region,
  sub.region,
  intermediate.region
)
country_lookup$country_code <- country_lookup$alpha.2
colnames(country_lookup) <- gsub("[.]", "_", colnames(country_lookup))
country_lookup$alpha_2 <- NULL
index <- read.csv("./covid-data/index.csv")
demographics <- read.csv("./covid-data/demographics.csv")
economy <- read.csv("./covid-data/economy.csv")
health <- read.csv("./covid-data/health.csv")
epidemiology <- read.csv("./covid-data/epidemiology.csv")
epidemiology_select <- epidemiology %>% select(
  location_key,
  date,
  new_confirmed
)

# Filter down to countries of interest
countries <- index %>% filter(
    aggregation_level == 0,
    !is.na(location_key)
  ) %>% left_join(
    demographics,
    by="location_key"
  ) %>% mutate(
    log_population = log(population),
    pct_pop_urban = population_urban / population,
    pct_pop_60_plus = (population_age_60_69 + population_age_70_79 + population_age_80_and_older) / population,
  ) %>% filter(
    !is.na(pct_pop_urban),
    !is.na(pct_pop_60_plus),
    # !is.na(human_development_index)
  ) %>% left_join(
    economy,
    by="location_key"
  ) %>% mutate(
    log_gdp_per_capita_usd = log(gdp_per_capita_usd)
  ) %>% filter(
    !is.na(log_gdp_per_capita_usd)
  ) %>% left_join(
    health,
    by="location_key"
  ) %>% mutate(
    pct_health_spending = health_expenditure_usd / gdp_per_capita_usd
  ) %>% filter(
    !is.na(life_expectancy),
    !is.na(pct_health_spending)
  )

# Explore these columns and make sure there aren't any crazy values
hist(countries$log_population)
hist(countries$pct_pop_urban)
hist(countries$pct_pop_60_plus)
hist(countries$log_gdp_per_capita_usd)
hist(countries$life_expectancy)
hist(countries$pct_health_spending)
countries <- filter(countries, pct_health_spending < 0.4)  # Filter out Sudan
hist(countries$pct_health_spending)

# Join in case count data
covid <- countries %>% left_join(
  epidemiology_select,
  by="location_key"
) %>% filter(
  !is.na(new_confirmed),
  new_confirmed >= 0
) %>% left_join(
  country_lookup,
  by="country_code"
) %>% arrange(
  location_key,
  date
) %>% select(
  location_key,
  country_name,
  region,
  sub_region,
  intermediate_region,
  date,
  new_confirmed,
  log_population,
  pct_pop_urban,
  pct_pop_60_plus,
  log_gdp_per_capita_usd,
  life_expectancy,
  pct_health_spending
)
write.csv(covid, file="./covid-data/covid.csv", row.names=FALSE)
