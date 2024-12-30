############################## FRONT MATTER ##############################
#   SCRIPT:     gen_commuting_time.R
#   AUTHOR:     Jacob Bradt (jbradt@g.harvard.edu)
#   NOTES:      Script to construct origin-specific commuting time for
#   replicating English et al (2018) Deepwater Horizon study

# Load dependencies:
pacman::p_load(here, data.table, tidyverse, fixest, sf, units, zoo, fredr, tidycensus)

# Set project root:
i_am("R/gen_commuting_time.R")

# STEP 1: Import and format respondent data ----

# Import and format local respondents data:
respondents_local <- fread(here("data/respondents_local.csv"))
respondents_local <- distinct(respondents_local[, c("caseid", "cutoffym", "zipcode")])

# Import and format national respondents data:
respondents_national <- fread(here("data/respondents_national.csv"))
respondents_national <- distinct(respondents_national[, c("caseid", "cutoffym", "zipcode")])

# Merge national and local origins data:
respondents <- rbindlist(list(respondents_local, respondents_national), fill = T)
rm(respondents_local, respondents_national)

# STEP 2: Import commuting data ----

# Import commuting data:
commute <- fread(here("data/nhgis0028_ds201_20135_zcta.csv"))

# Get commuting table:
commute_tab <- commute[, grepl("UFHE", names(commute)), with = F]
commute_tab[, (colnames(commute_tab)) := lapply(.SD, function(x){x / UFHE001}), .SDcols = colnames(commute_tab)]

# Calculate average commuting time:
commute_vals <- c(seq(2.5,45, 5.0), 52, 75, 90)
commute[, mean_commute := rowSums(as.matrix(commute_tab[, 2:ncol(commute_tab)]) %*% commute_vals) ] 
rm(commute_tab, commute_vals)

# STEP 3: Merge commute times to respondent data ----

# Merge commute times to respondent data:
respondents <- merge(respondents, commute[, c("ZCTA5A", "mean_commute")], by.x = "zipcode", by.y = "ZCTA5A")

# Save respondent commute times:
fwrite(respondents, here("data/respondents_commutes.csv"))
rm(respondents)

# STEP 4: Plot distribution of commuting times ----

# Plot commute times:
fig_commute_time <- ggplot(commute, aes(mean_commute)) +
  geom_histogram(fill = "#009E73", color = "black", alpha = 0.75) +
  theme_minimal() + 
  #scale_x_yearqtr(format = '%Y Q%q') +
  theme(plot.subtitle = element_text(hjust = -0.05), legend.title = element_blank(), 
        axis.line.x.bottom = element_line(), axis.line.y.left = element_line()) +
  xlab("One-way Average Commuting Time (min.)") +
  ylab(element_blank()) + 
  labs(title = element_blank(), subtitle = "Count")
ggsave(here( "output/fig_commute_time.jpg"), plot = fig_commute_time, width = 6, height = 5, units = "in", dpi = 500)
rm(commute, fig_commute_time)
