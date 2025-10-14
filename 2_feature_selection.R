# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2025-08-10 10:38:00
# | LastEditTime: 2025-09-12 13:06:28
# | FilePath: \R_scripts\2_feature_selection.R
# | Description:
# | ---------------------------------------

library(tidyverse)
library(gtsummary)
library(rstatix)
library(nortest)
library(glmnet)
library(ggfortify)

# ----------------------------------------

load("output/1_data_tidy.RData")

data_tidy <- data_tidy %>%
    select(-c(CDAI_change, CDAI_after, monitoring_c))

# ---------------------------------------- descriptive statistics
## ---------------------------------------- numeric variables

numeric_var <- data_tidy %>%
    select_if(is.numeric) %>%
    names()

### ---------------------------------------- normality test

ad.test.multi <- function(x) {
    ad.test(x) %>%
        broom::tidy()
}

normality <- map_df(
    data_tidy[numeric_var],
    ad.test.multi
) %>%
    bind_cols(variable = numeric_var) %>%
    select(-method) %>%
    rename("p_norm" = "p.value")

### ---------------------------------------- variance test

variance <- map_df(
    data_tidy[numeric_var],
    function(x) {
        levene_test(data_tidy, x ~ data_tidy$group)
    }
) %>%
    select(p_vari = p) %>%
    cbind(variable = numeric_var)

### ---------------------------------------- t test

t_test_res <- data_tidy %>%
    select(all_of(numeric_var), group) %>%
    pivot_longer(
        cols = all_of(numeric_var),
        names_to = "variable"
    ) %>%
    group_by(variable) %>%
    t_test(value ~ group, var.equal = TRUE) %>%
    select(variable, p_t_test = p)

### ---------------------------------------- wilcox test

wilcox_test_res <- data_tidy %>%
    select(all_of(numeric_var), group) %>%
    pivot_longer(
        cols = all_of(numeric_var),
        names_to = "variable"
    ) %>%
    group_by(variable) %>%
    wilcox_test(value ~ group) %>%
    select(variable, p_wilcox_test = p)

p_numb <- normality %>%
    left_join(variance, by = "variable") %>%
    left_join(t_test_res, by = "variable") %>%
    left_join(wilcox_test_res, by = "variable") %>%
    group_by(variable) %>%
    mutate(p_final = if_else(p_norm >= 0.05 & p_vari >= 0.05, p_t_test, p_wilcox_test)) %>%
    add_significance(p.col = "p_final")

## ---------------------------------------- nominal variables

nominal_var <- data_tidy %>%
    select_if(is.factor) %>%
    select(-group) %>%
    names()

p_chisq <- data_tidy %>%
    select(all_of(nominal_var), group) %>%
    pivot_longer(
        cols = -group,
        names_to = "variable"
    ) %>%
    group_by(variable) %>%
    do(chisq_test(.$group, .$value)) %>%
    select(variable, p_chisq = p)

## ---------------------------------------- fisher

freq <- freq_table(data_tidy, group, all_of(nominal_var))

fisher_map <- function(x) {
    temp <- freq_table(data_tidy, group, x) %>%
        select(-prop) %>%
        pivot_wider(
            names_from = "group",
            values_from = "n"
        ) %>%
        replace(is.na(.), 0)
    p_fisher <- temp %>%
        select(all_of(unique(data_tidy$group))) %>%
        fisher_test(simulate.p.value = TRUE)
    min_n <- temp %>%
        pivot_longer(
            cols = all_of(unique(data_tidy$group))
        ) %>%
        arrange(value) %>%
        slice(1) %>%
        select(value) %>%
        cbind(variable = x)
    cbind(p_fisher, min_n)
}

p_fisher <- map_df(
    all_of(nominal_var),
    fisher_map
) %>%
    select(p_fisher = p, min_n = value, variable)

p_norm <- tibble(variable = all_of(nominal_var)) %>%
    left_join(p_chisq, by = "variable") %>%
    left_join(p_fisher, by = "variable") %>%
    group_by(variable) %>%
    mutate(p_final = ifelse(nrow(data_tidy) > 40 & min_n >= 5, p_chisq, p_fisher))

p_numb
p_norm

p_value <- bind_rows(
    p_numb %>%
        select(variable, p_final),
    p_norm %>%
        select(variable, p_final)
)

# ---------------------------------------- get summary

theme_gtsummary_language("en", big.mark = "")

total_summary <- data_tidy %>%
    relocate(age, gender, height, weight, CDAI_before, CRP_before, dose) %>%
    mutate(
        fecal_calprotectin = factor(fecal_calprotectin, levels = c(0, 1, 2)),
        Montreal_age = factor(Montreal_age, levels = c(1, 2, 3)),
        Montreal_L = factor(Montreal_L, levels = c(1, 2, 3, 4, 5, 6)),
        Montreal_B = factor(Montreal_B, levels = c(1, 2, 3, 4)),
        Montreal_P = factor(Montreal_P, levels = c(0, 1))
    ) %>%
    tbl_summary(
        by = group,
        statistic = list(
            all_continuous2() ~ c("{mean} \u00B1 {sd}", "{median} ({p25}, {p75})", "{p_miss}"),
            all_categorical() ~ "{n} ({p}%)"
        ),
        digits = list(
            all_continuous() ~ 1,
            all_categorical() ~ c(0, 1)
        ),
        value = list(
            ADA ~ "1",
            Azathioprine ~ "1",
            Mesalazine ~ "1",
            Thalidomide ~ "1",
            MTX ~ "1"
        )
    ) %>%
    as_tibble() %>%
    rename("variable" = `**Characteristic**`) %>%
    left_join(p_value) %>%
    mutate(p_final = ifelse(p_final >= 0.001, round(p_final, 3), "< 0.001")) %>%
    mutate(across(where(is.numeric), as.character)) %>%
    mutate(variable = case_when(
        variable == "age" ~ "Age (year)",
        variable == "gender" ~ "Gender",
        variable == "height" ~ "Height(cm)",
        variable == "weight" ~ "Weight(kg)",
        variable == "CDAI_before" ~ "CDAI",
        variable == "CRP_before" ~ "CRP",
        variable == "dose" ~ "IFX dose(mg)",
        variable == "D_dimer" ~ "D-dimer(μg/L)",
        variable == "Fg" ~ "Fibrinogen(g/L)",
        variable == "fecal calprotectin" ~ "Fecal_calprotectin",
        variable == "Montreal_age" ~ "Montreal Classification(age)",
        variable == "Montreal_L" ~ "Montreal Classification(L)",
        variable == "Montreal_B" ~ "Montreal Classification(B)",
        variable == "Montreal_P" ~ "Montreal Classification(P)",
        .default = variable
    ))

total_summary %>%
    replace(is.na(.), "") %>%
    write_csv(file = "output/2_total_summary.csv")

# ---------------------------------------- lasso regression
## ---------------------------------------- select best lambda

set.seed(2025)

var_include_1 <- p_value %>%
    filter(p_final <= 0.1) %>%
    pull(variable)

data_lasso <- data_tidy %>%
    select(all_of(var_include_1), group)

lasso_var <- data_lasso %>%
    select(-group) %>%
    as.matrix()

lasso_pred <- data_lasso %>%
    select(group) %>%
    mutate(group = as.numeric(group)) %>%
    as.matrix()

lasso_cv <- cv.glmnet(lasso_var, lasso_pred, alpha = 1, nfold = 20, nlambda = 400, family = "binomial", type.measure = "auc")

lasso_cv_plot <- autoplot(
    lasso_cv,
    label = F,
    color = "steelblue",
    alpha = 0.5
) +
    theme_bw() +
    xlab("Log Lambda") +
    ylab("AUROC") +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.position = "none"
    )

tiff(filename = "plot/lasso_cv_plot.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

lasso_cv_plot

dev.off()

## ---------------------------------------- lasso regression with best lambda

lambda_min <- lasso_cv$lambda.min
lambda_1se <- lasso_cv$lambda.1se

best_lambda <- lasso_cv$lambda.1se

best_lambda

lasso_fit <- glmnet(x = lasso_var, y = lasso_pred, alpha = 1, lambda = best_lambda, fimaly = "binomial", type.measure = "auc")

coef(lasso_fit) %>%
    round(5)

## ---------------------------------------- trace plot

library(ggrepel)

trace_mod <- glmnet(x = lasso_var, y = lasso_pred, alpha = 1, fimaly = "binomial", type.measure = "auc")

trace_plot_labels <- coef(lasso_fit) %>%
    cbind(trace_mod$beta[, dim(trace_mod$beta)[2]]) %>%
    as.matrix() %>%
    as.data.frame() %>%
    set_names(c("s0", "y")) %>%
    slice(-1) %>%
    arrange(desc(y)) %>%
    rownames_to_column() %>%
    mutate(
        label = case_when(
            rowname == "Fg" ~ "Fibrinogen",
            rowname == "monitoring_c" ~ "IFX concentration",
            rowname == "CDAI_before" ~ "CDAI",
            rowname == "Montreal_age" ~ "Age at diagnosis",
            rowname == "CRP_before" ~ "CRP",
            rowname == "D_dimer" ~ "D-dimer",
            rowname == "age" ~ "Age",
            rowname == "weight" ~ "Weight",
            rowname == "height" ~ "Height",
            rowname == "ALB" ~ "Albumin",
            .default = rowname
        )
    ) %>%
    mutate(
        x = rep(-8.05, dim(trace_mod$beta)[1]),
        x_point = rep(-8.05, dim(trace_mod$beta)[1])
    )

trace_plot <- autoplot(trace_mod, label = F, size = 0.8, xvar = "lambda") +
    geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "gray") +
    geom_vline(aes(xintercept = log(lambda_1se)), linetype = 2, linewidth = 0.5, color = "black") +
    geom_vline(aes(xintercept = log(lambda_min)), linetype = 2, linewidth = 0.5, color = "black") +
    geom_point(
        data = trace_plot_labels,
        aes(x = x_point, y = y)
    ) +
    geom_text_repel(
        data = trace_plot_labels,
        aes(label = label, x = x, y = y),
        color = "black",
        max.overlaps = 50,
        nudge_x = -4,
        direction = "y",
        force = 5,
        hjust = 0,
        size = 4
    ) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.position = "none"
    )

tiff(filename = "plot/trace_plot.tiff", width = 15, height = 10, res = 300, units = "in", compression = "lzw")

trace_plot

dev.off()

## ----------------------------------------

library(cowplot)

lasso_grid_plot <- plot_grid(lasso_cv_plot, trace_plot, ncol = 1, align = "hv", labels = "AUTO", label_size = 20)

tiff(filename = "plot/lasso_grid_plot.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

lasso_grid_plot

dev.off()

# ---------------------------------------- output final data

final_var <- trace_plot_labels %>%
    filter(s0 != 0) %>%
    pull(rowname)

data_model <- data_tidy %>%
    select(all_of(final_var), group)

save(data_model, file = "output/data_model.RData")

# ! end
