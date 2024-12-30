############################## FRONT MATTER ##############################
#   SCRIPT:     gen_oil_vars.R
#   AUTHOR:     Jacob Bradt (jbradt@g.harvard.edu)
#   NOTES:      Script to construct oil supply shock time series

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, fixest, sf, units, zoo, vars, jsonlite, fredr, ggplot2)

# Set project root:
i_am("R/gen_oil_vars.R")

# API keys:
fredr::fredr_set_key("588eb023c04a4e8fd31b6b828b6ddbd6")

# STEP 1: Download and process global oil production data ----

# Import EIA global oil production data series:
production <- jsonlite::fromJSON(paste(readLines(here("data/eia_intl_crude_prod.txt")), collapse=""))
production <- data.table(production$data)
names(production) <- c("year_mon", "crude_prod")

# Generate monthly percent change in global production:
production[, year_mon := as.yearmon(year_mon, "%Y%m")]
production <- production[order(year_mon)]
production[, crude_prod := as.numeric(crude_prod)]
production[, d_prod := ((crude_prod - lag(crude_prod)) / lag(crude_prod))*100 ]
production[, crude_prod := NULL]

# STEP 2: Download and process oil price data ----

# Download acquisition cost of imported crude oil for US refiners from EIA:
price <- jsonlite::fromJSON(paste(readLines(here("data/eia_refiner_costs.txt")), collapse=""))
price <- data.table(price$data)
names(price) <- c("year_mon", "refiner_price")
price[, year_mon := as.yearmon(year_mon, "%Y%m")]

# Deflate refiner import acquisition cost using CPI:
cpi <-
  data.table(fredr::fredr(
    "CPIAUCSL",
    observation_start = as.Date(paste0(min(year(price$year)), "-01-01")),
    observation_end = as.Date(paste0(max(year(price$year)), "-12-31"))
  ))
cpi[, year_mon := as.yearmon(date)]
cpi <- cpi[, .(cpi = mean(value)), by = year_mon]
price <- merge(price, cpi)
base_cpi <- cpi[year_mon == "Jan 2011"]$cpi
price[, refiner_price := as.numeric(refiner_price) * (base_cpi / cpi)]
price[, cpi := NULL]
rm(cpi, base_cpi)

# STEP 3: Download index of real economic activity ----

# Download index of global real economic activity:
rea <- data.table(fredr::fredr("IGREA"))
rea <- rea[,c(1,3)]
setnames(rea, old = "value", new = "rea")

# Format REA date data:
rea[, year_mon := as.yearmon(date)]
rea[,date := NULL]

# STEP 4: Merge full data for SVAR ----

# Merge index of global real economic activity:
full_data <- merge(production, rea)

# Merge production and price data:
full_data <- merge(full_data, price)

# Create ts object:
full_data <- full_data[year(year_mon) >= 1985]
full_data_ts <- ts(full_data[, 2:4], frequency = 12, start = c(1985, 1), end = c(2023,10))

# STEP 5: Estimate SVAR and save estimated structural shocks----

# Estimate reduced form VAR
var_est <- VAR(full_data_ts, p = 24, type = "none")

# Estimate structural coefficients
A <- diag(1, 3)
A[lower.tri(A)] <- NA
svar_est <- SVAR(var_est, Amat = A, estmethod ="direct")

# Calculate supply shocks:
supply_shock <- data.table(t(solve(svar_est$A) %*% t(resid(var_est))))

# Add temporal data:
supply_shock[, year_mon := full_data[25:nrow(full_data)]$year_mon]

# Save estimated structural shocks:
fwrite(supply_shock, here("data/oil_struc_shocks.csv"))
rm(A, full_data, full_data_ts, price, production, rea, svar_est, var_est)

# STEP 6: Plot estimated structural shocks----

# Aggregate to the quarter level
supply_shock[, (c("qtr", "year")) := list(as.yearqtr(year_mon), year(year_mon))]
supply_shock <- supply_shock %>%
  group_by(qtr) %>%
  summarise(d_prod = mean(d_prod),
            rea = mean(rea),
            refiner_price = mean(refiner_price))

# Plot structural shocks:
fig.supply <- ggplot(supply_shock) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_line(aes(qtr, d_prod)) + 
  theme_minimal() + 
  scale_x_yearqtr(format = '%Y Q%q') +
  theme(plot.subtitle = element_text(hjust = -0.05), legend.title = element_blank(), 
        axis.line.x.bottom = element_line(), axis.line.y.left = element_line()) +
  xlab(element_blank()) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Estimated Oil Supply Shock")
ggsave(here( "output/fig_shock_oilsupply_wide.jpg"), plot = fig.supply, width = 6, height = 3, units = "in", dpi = 500)
ggsave(here( "output/fig_shock_oilsupply.jpg"), plot = fig.supply, width = 6, height = 4, units = "in", dpi = 500)
fig.rea <- ggplot(supply_shock) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_line(aes(qtr, rea)) + 
  theme_minimal() + 
  scale_x_yearqtr(format = '%Y Q%q') +
  theme(plot.subtitle = element_text(hjust = -0.05), legend.title = element_blank(), 
        axis.line.x.bottom = element_line(), axis.line.y.left = element_line()) +
  xlab(element_blank()) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Estimated Aggregate Demand Shock")
ggsave(here( "output/fig_shock_aggdem_wide.jpg"), plot = fig.rea, width = 6, height = 3, units = "in", dpi = 500)
fig.oil.dem <- ggplot(supply_shock) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_line(aes(qtr, refiner_price)) + 
  theme_minimal() + 
  scale_x_yearqtr(format = '%Y Q%q') +
  theme(plot.subtitle = element_text(hjust = -0.05), legend.title = element_blank(), 
        axis.line.x.bottom = element_line(), axis.line.y.left = element_line()) +
  xlab(element_blank()) +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Estimated Oil Demand Shock")
ggsave(here( "output/fig_shock_oildem_wide.jpg"), plot = fig.oil.dem, width = 6, height = 3, units = "in", dpi = 500)
rm(fig.oil.dem, fig.rea, fig.supply, supply_shock)
