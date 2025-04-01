# Load Required Libraries
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(sf)
library(stringr)
library(plm)
library(spdep)
library(fastmatrix)

# Load Data Files
setwd("/Users/kexinxie/Downloads/GitHub/SpatialSpilloverEffect_Measles") 
mean_cost_inci_exp_county_all <- fread("data/mean_costs_inci_all_expt_county_all.csv")
median_cost_inci_exp_county_all <- fread("data/median_costs_inci_all_expt_county_all.csv")
Q3_cost_inci_exp_county_all <- fread("data/Q3_costs_inci_all_expt_county_all.csv")
nighty_cost_inci_exp_county_all <- fread("data/nighty_costs_inci_all_expt_county_all.csv")

county_chara <- read.csv("data/simulation_county_characteristic.csv")
hh_chara <- read.csv("data/simulation_household_characteristic.csv")
va.shp <- st_read("data/Virginia_City_County_Census_Boundaries/SDE_USDC_CENSUS_VA_COUNTY.shp")

commuting <- fread("data/VA-commuting-flows-2020.csv") %>%
  mutate(
    Re_County_FIPS_Code = sprintf("%03d", Re_County_FIPS_Code),
    Re_State_FIPS_Code = sprintf("%02d", Re_State_FIPS_Code),
    Wo_County_FIPS_Code = sprintf("%03d", Wo_County_FIPS_Code),
    Wo_State_FIPS_Code = sprintf("%02d", Wo_State_FIPS_Code),
    county_fips_Re = str_c(Re_State_FIPS_Code, Re_County_FIPS_Code),
    county_fips_Wo = str_c(Wo_State_FIPS_Code, Wo_County_FIPS_Code)
  )

commuting_virginia <- commuting %>%
  select(county_fips_Re, county_fips_Wo, Workers_in_Commuting_Flow) %>%
  spread(key = county_fips_Wo, value = Workers_in_Commuting_Flow) %>%
  column_to_rownames(var = "county_fips_Re") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate_all(as.numeric)

# Function to prepare data based on selected model
prepare_model_data <- function(cost_data) {
  
  model_data <- cost_data %>%
    left_join(county_chara, by = "county_fips") %>%
    left_join(hh_chara, by = "county_fips") %>%
    mutate(
      n_mmr = case_when(
        alpha == 0 ~ n_alpha_0,
        alpha == 5 ~ n_alpha_5,
        alpha == 15 ~ n_alpha_15,
        alpha == 25 ~ n_alpha_25
      ),
      incidence_rate = incidence / pop,
      p_mmr = n_mmr / pop,
      totalCost.rate.log = log(totalCost / pop)
    ) %>%
    left_join(va.shp %>% mutate(GEOID = as.integer(GEOID)), by = c("county_fips" = "GEOID")) %>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
    mutate(expt=rep(1:30, 133)) %>%
    arrange(expt, county_fips)

  return(model_data)
}

# Choose model and prepare data
if (model == "mean") {
  model_cost_inci_exp_county_all <- prepare_model_data(mean_cost_inci_exp_county_all)
} else if (model == "median") {
  model_cost_inci_exp_county_all <- prepare_model_data(median_cost_inci_exp_county_all)
} else if (model == "Q3") {
  model_cost_inci_exp_county_all <- prepare_model_data(Q3_cost_inci_exp_county_all)
} else if (model == "nighty") {
  model_cost_inci_exp_county_all <- prepare_model_data(nighty_cost_inci_exp_county_all)
}

# Format and panelize data
model_cost_inci_exp_county_all <- model_cost_inci_exp_county_all %>% arrange(expt, county_fips)

unique.county <- model_cost_inci_exp_county_all %>%
  group_by(county_fips) %>%
  distinct(county_fips, .keep_all = TRUE)

panel_cost_inci_exp_county_all <- plm::pdata.frame(
  as.data.frame(model_cost_inci_exp_county_all[, c(1:27, 48)]),
  index = c("county_fips", "expt"),
  drop.index = FALSE,
  drop.const.series = TRUE
)

# Construct commuting-based adjacency matrix
pop <- unique.county$pop
adjmat.commuting.only <- commuting_virginia / sqrt(pop %*% t(pop))
adjmat.commuting.only <- adjmat.commuting.only %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
diag(adjmat.commuting.only) <- 0
adjmat.commuting.only <- as.matrix(abs(adjmat.commuting.only))

# Convert adjacency matrix to listw and Kronecker product for panel
va.shp.listw.commuting.only <- mat2listw(adjmat.commuting.only, style = "W")
va.shp.adjmat.commuting.only <- listw2mat(va.shp.listw.commuting.only)

pool.va.shp.adjmat.commuting.only <- kronecker.prod(diag(30), va.shp.adjmat.commuting.only)
pool.va.shp.listw.commuting.only <- mat2listw(pool.va.shp.adjmat.commuting.only, style = "W")

# Final data objects
data <- panel_cost_inci_exp_county_all
pdata <- data %>% arrange(expt, county_fips)

w <- va.shp.adjmat.commuting.only
W <- pool.va.shp.adjmat.commuting.only
w_listw <- mat2listw(w, style = "W")
W_listw <- mat2listw(W, style = "W")

# Factor conversions (optional)
if (tau_factor) {
  levels <- c("0.5")
  data$tau <- factor(data$tau, ordered = FALSE)
  pdata$tau <- factor(pdata$tau, ordered = FALSE)
  unique.county$tau <- factor(unique.county$tau, ordered = FALSE)
  data$tau <- relevel(data$tau, ref = levels)
  pdata$tau <- relevel(pdata$tau, ref = levels)
  unique.county$tau <- relevel(unique.county$tau, ref = levels)
}

if (vhi_factor) {
  levels <- c("75")
  data$vhi <- factor(data$vhi, ordered = FALSE)
  pdata$vhi <- factor(pdata$vhi, ordered = FALSE)
  unique.county$vhi <- factor(unique.county$vhi, ordered = FALSE)
  data$vhi <- relevel(data$vhi, ref = levels)
  pdata$vhi <- relevel(pdata$vhi, ref = levels)
  unique.county$vhi <- relevel(unique.county$vhi, ref = levels)
}

# Cleanup
rm(county_chara, hh_chara, commuting, model_cost_inci_exp_county_all,
   panel_cost_inci_exp_county_all, pop, adjmat.commuting.only,
   va.shp.listw.commuting.only, va.shp.adjmat.commuting.only,
   pool.va.shp.adjmat.commuting.only, pool.va.shp.listw.commuting.only)
