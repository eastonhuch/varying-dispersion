# Imports
require(IndexNumR)
require(tidyverse)

# Read in raw data
categories <- c(
  cer="Cereals",
  coo="Cookies",
  cra="Crackers"
)
last_date <- as.Date("1993-09-16")  # No missing values before this date
process_category <- function(category) {
  # Read in raw data
  dat_raw <- dominicksData(
    categories[category],
    paste0("dominicks-data/", category, "/", "w", category, ".csv"),
    paste0("dominicks-data/", category, "/", "upc", category, ".csv"))
  
  # Sum to date, upc level
  dat_date_upc <- dat_raw %>%
    mutate(date=start) %>%
    group_by(date, upc) %>%
    summarise(sales = sum(quantity)) %>%
    mutate(in_holdout = date > last_date)
  
  # Inspect UPCs and filter down to subset
  dat_upc <- dat_date_upc %>%
    filter(!in_holdout) %>%
    group_by(upc) %>%
    summarize(
      sales = sum(sales),
      first_date = min(date),
      last_date = max(date),
      days = length(date)) %>%
    filter(days == max(days))
  included_upcs <- dat_upc$upc
  
  # Filter data
  dat_date_upc_filtered <- dat_date_upc %>%
    filter(upc %in% included_upcs) %>%
    arrange(date, upc)
  
  # Write data
  write.csv(
    dat_date_upc_filtered,
    paste0("dominicks-data/", category, "-processed.csv"),
    row.names=FALSE)
}

for (category in names(categories)) {
  print(category)
  process_category(category)
}

for (category in names(categories)) {
  dat_c <- read.csv(paste0("dominicks-data/", category, "-processed.csv"))
  if (category == names(categories)[1]) {
    dat <- dat_c
  } else {
    dat <- rbind(dat, dat_c)
  }
}
write.csv(dat, "dominicks-data/all-categories.csv", row.names=FALSE)