# Imports
require("IndexNumR")
require(tidyverse)

# Read in raw data
dat_raw <- dominicksData("Bottled Juices", "wbjc.csv", "upcbjc.csv")

# Sum to date, upc level
dat_date_upc <- dat_raw %>%
  mutate(date=start) %>%
  group_by(date, upc) %>%
  summarise(sales = sum(quantity))

# Explore data by date
dat_date <- dat_date_upc %>%
  group_by(date) %>%
  summarise(sales=sum(sales)) %>%
  arrange(date)
plot(
  dat_date$date, dat_date$sales,
  type="l", ylim=c(0, max(dat_date$sales)))
last_date <- as.Date("1993-09-16")  # No missing values before this date
dat_date_to_last_date <- filter(dat_date, date <= last_date)
dat_date_upc_to_last_date <- filter(dat_date_upc, date <= last_date)

# Inspect UPCs and filter down to subset
dat_upc <- dat_date_upc_to_last_date %>%
  group_by(upc) %>%
  summarize(
    sales = sum(sales),
    first_date = min(date),
    last_date = max(date)) %>%
  mutate(days = as.integer(last_date - first_date)) %>%
  filter(days == max(days)) %>%
  filter(rank(-sales) <= 40)
hist(dat_upc$sales)
min(dat_upc$sales)
included_upcs <- dat_upc$upc

# Filter data
dat_date_upc_filtered <- dat_date_upc_to_last_date %>%
  filter(upc %in% included_upcs) %>%
  arrange(date, upc)
min(dat_date_upc_filtered$sales)
# No sales of zero

# Check for missing values
date_scaffold <- dat_date_to_last_date[, c("date")]
upc_scaffold <- data.frame(upc=included_upcs)
full_scaffold <- cross_join(date_scaffold, upc_scaffold)
if (nrow(full_scaffold) != nrow(dat_date_upc_filtered)) {
  stop("Some dates are missing in dat_date_upc_filtered")
}

# Save processed data
write.csv(dat_date_upc_filtered, "./dominicks.csv", row.names=FALSE)

# Plot all UPCs at once
num_weeks <- nrow(dat_date_to_last_date)
sales_y_lim <- c(
  min(dat_date_upc_filtered$sales),
  max(dat_date_upc_filtered$sales))
plot(
  dat_date_to_last_date$date, rep(0.5, num_weeks),
  ylim = sales_y_lim,
  log="y",
  main="Sales Over Time by UPC",
  xlab="Date", ylab="Sales"
)
for (i in seq_along(included_upcs)) {
  u <- included_upcs[i]
  dat_single_upc <- filter(dat_date_upc_filtered, upc == u)
  lines(dat_single_upc$date, dat_single_upc$sales, col=i)
}

