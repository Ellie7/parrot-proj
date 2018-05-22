#### Tamora's script for extracting the data from the excel sheet
# (mean_hatch n_hatch, se_hatch, mean_nestling_surv, n_nestling_surv, se_nestling_surv)

library(tidyverse)
library(readxl)


## Set working directory to data directory first
parrot_lh_data <- read_excel("Life History 2006-2014.xlsx")
glimpse(parrot_lh_data)

## Add columns to calculate:
## - proportion hatched = brood size at hatching / total # eggs
## - nestling survival = # fledged / brood size at hatching
parrot_lh_data <- parrot_lh_data %>% 
  mutate(prop_hatched = `brood size at hatching`/`total eggs`,
         nestling_survival = ifelse(`brood size at hatching` > 0, 
                                    fledged/`brood size at hatching`, NA))

## Andrew's SE function (modified to remove NAs)
se_fnc <- function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}

## Summarise to get yearly means and sample sizes
lh_summary <- parrot_lh_data %>% 
  group_by(year) %>% 
  summarise(mean_hatch = mean(prop_hatched, na.rm = TRUE), 
            n_hatch = sum(!is.na(prop_hatched)),
            se_hatch = se_fnc(prop_hatched),
            mean_nestling_surv = mean(nestling_survival, na.rm = TRUE), 
            n_nestling_surv = sum(!is.na(nestling_survival)),
            se_nestling_surv = se_fnc(nestling_survival))

## Summarise to get total means and sample sizes
total_summary <- parrot_lh_data %>% 
  summarise(mean_hatch = mean(prop_hatched, na.rm = TRUE), 
            n_hatch = sum(!is.na(prop_hatched)),
            se_hatch = se_fnc(prop_hatched),
            mean_nestling_surv = mean(nestling_survival, na.rm = TRUE), 
            n_nestling_surv = sum(!is.na(nestling_survival)),
            se_nestling_surv = se_fnc(nestling_survival))