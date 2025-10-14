# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2025-05-27 12:34:23
# | LastEditTime: 2025-08-19 21:44:55
# | FilePath: \R_scripts\1_data_cleaning.R
# | Description: 
# | ---------------------------------------

library(tidyverse)
library(readxl)
library(naniar)

fread <- data.table::fread

set.seed(2025)

# ---------------------------------------- load data

data_raw <- read_excel("D:/OneDrive/After_work/2025/IFX_efficiency/data/IFX-20250806 - 副本.xlsx")

# ---------------------------------------- missing over 40%

variables_missing_over_40 <- data_raw %>%
    miss_var_summary() %>%
    filter(pct_miss >= 40) %>%
    pull(variable)

data_tidy <- data_raw %>%
    select(-all_of(variables_missing_over_40), -c(data_source, monitoring_date, Montreal_type)) %>%
    mutate(
        age = as.numeric(age),
        ADA = case_when(
            ADA == "+" ~ 1,
            ADA == "-" ~ 0,
            TRUE ~ NA
        ),
        CDAI_before = ifelse(CDAI_before == "/", NA, CDAI_before),
        monitoring_c = ifelse(monitoring_c == ">45", "45", monitoring_c)
    ) %>%
    mutate(
        across(c(CDAI_before, monitoring_c, RBC, `D-dimer`, APTT, TT, PT, Fg), as.numeric),
        across(c(gender, `fecal calprotectin`:MTX, L:behavior_P, dose, Age), as.factor)
    ) %>%
    rename(
        "Montreal_age" = "Age",
        "Montreal_L" = "L",
        "Montreal_B" = "behavior_B",
        "Montreal_P" = "behavior_P",
    ) %>%
    mutate(
        CDAI_change = CDAI_after - CDAI_before,
        group = ifelse(CDAI_change < -70, 1, 0)
    ) %>%
    bind_cols(arrange(miss_case_summary(.), case)) %>%
    filter(pct_miss < 40) %>%
    filter(!is.na(group)) %>%
    rename(
        "fecal_calprotectin" = "fecal calprotectin",
        "D_dimer" = "D-dimer",
    )

# ----------------------------------------

case_missing_over_40 <- data_tidy %>%
    filter(pct_miss >= 40)

# ---------------------------------------- handle missing data

var_simple_impute <- data_tidy %>%
    miss_var_summary() %>%
    filter(pct_miss > 0 & pct_miss < 5) %>%
    pull(variable)

data_simply_impute <- data_tidy %>%
    mutate(
        CRP_before = impute_median(CRP_before),
        D_dimer = impute_median(D_dimer),
        APTT = impute_median(APTT),
        TT = impute_median(TT),
        PT = impute_median(PT),
        Fg = impute_median(Fg),
        WBC = impute_median(WBC),
        height = impute_median(height),
        weight = impute_median(weight),
        RBC = impute_median(RBC),
        monitoring_c = impute_median(monitoring_c),
        Montreal_age = impute_median(Montreal_age),
        Montreal_L = impute_median(Montreal_L),
        Montreal_B = impute_median(Montreal_B),
        Montreal_P = impute_median(Montreal_P)
    ) %>%
    mutate(
        dose = impute_mode(dose),
        fecal_calprotectin = impute_mode(fecal_calprotectin)
    ) %>%
    select(-c(name, case, n_miss, pct_miss))

# ---------------------------------------- 

var_multiple_impute <- data_tidy %>%
    miss_var_summary() %>%
    filter(pct_miss >= 5 & pct_miss < 30) %>%
    pull(variable)

library(mice)

impute_strategy <- mice(data_simply_impute, variables = var_multiple_impute,m = 5, maxit = 4, method = "rf", seed = 2025)

data_complete_1 <- complete(impute_strategy, action = 1) %>%
    mutate(imputation = 1)

data_complete_2 <- complete(impute_strategy, action = 2) %>%
    mutate(imputation = 2)

data_complete_3 <- complete(impute_strategy, action = 3) %>%
    mutate(imputation = 3)

data_complete_4 <- complete(impute_strategy, action = 4) %>%
    mutate(imputation = 4)

data_complete_5 <- complete(impute_strategy, action = 5) %>%
    mutate(imputation = 5)

data_imputed <- bind_rows(data_complete_1, data_complete_2, data_complete_3, data_complete_4, data_complete_5) %>%
    group_by(order) %>%
    summarise(
        across(where(is.numeric), mean, na.rm = TRUE),
        across(where(is.factor), ~ factor(levels(.)[which.max(table(.))]))
    )

# ---------------------------------------- exclude

patient_CDAI_before_less_than_70 <- data_imputed %>%
    filter(CDAI_before < 70)

data_tidy <- data_imputed %>%
    filter(!order %in% patient_CDAI_before_less_than_70$order) %>%
    select(-c(order, imputation)) %>%
    mutate(group = as.factor(group))

save(data_tidy, file = "output/1_data_tidy.RData")

# ! end 