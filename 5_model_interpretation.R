# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2024-10-15 21:13:11
# | LastEditTime: 2025-10-12 14:00:46
# | FilePath: \R_scripts\5_model_interpretation.R
# | Description:
# | ---------------------------------------

library(DALEX)
library(DALEXtra)
library(randomForest)
library(tidymodels)
library(probably)
library(stacks)

load(file = "output/final_mod.RData")

load(file = "output/data_model.RData")

load(file = "output/testset.RData")

data_model <- data_model |>
    as_tibble()

# ---------------------------------------- calibration curve

plot_calibration_curve <- IFX_test %>%
    mutate(predict(final_mod, new_data = ., type = "prob")) %>%
    cal_plot_breaks(truth = group, estimate = .pred_0, include_rug = FALSE, num_breaks = 8) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold")
    )

plot_calibration_curve

tiff(filename = "plot/calibration curve.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_calibration_curve

dev.off()

# ---------------------------------------- 

explainer_rf <- DALEXtra::explain_tidymodels(
    model = final_mod,
    data = data_model,
    y = data_model$group,
    label = "Random Forest"
)

# 2---------------------------------------- in dataset level

variable_label <- c(
    "ESR" = "ESR",
    "CDAI_before" = "CDAI",
    "Montreal_age" = "Montreal age",
    "CRP_before" = "CRP",
    "RBC" = "RBC"
)

plot_rf_profile <- model_profile(
    explainer_rf,
    variable = names(data_model[, 1:5])
) %>%
    plot() +
    facet_wrap(
        ~`_vname_`,
        labeller = as_labeller(variable_label),
        scales = "free_x",
        ncol = 2
    ) +
    labs(y = "Probability of Clinical Response") +
    theme_bw() +
    theme(
        title = element_blank()
    ) +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.spacing = unit(1.5, "lines")
    )

plot_rf_profile

tiff(filename = "plot/plot_rf_profile.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

plot_rf_profile

dev.off()

# 2---------------------------------------- in person level

bd_rf <- predict_parts(explainer_rf, new_observation = data_model[, -6][1, ])

plot_rf_bd <- plot(
    bd_rf,
    vnames = c("Intercept", "CDAI = 83.4", "RBC = 3.89", "CRP = 55.5", "ESR = 23", "Age (16 - 40)", "Prediction"),
    title = NULL,
    subtitle = NULL,
    digits = 3,
    add_contributions = FALSE
) +
    geom_text(aes(y = right_side), size = 5, color = "black", nudge_y = 0.12) +
    xlab("Probability of Clinical Response") +
    theme(
        title = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_blank(),
        panel.spacing = unit(1, "lines")
    )

plot_rf_bd

tiff(filename = "plot/plot_rf_bd.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

plot_rf_bd

dev.off()

# ----------------------------------------

library(cowplot)

plot_interpretation <- plot_grid(plot_calibration_curve, vip_plot, plot_rf_profile, plot_rf_bd, ncol = 2, align = "hv", labels = "AUTO", label_size = 20)

plot_interpretation

tiff(filename = "plot/plot interpretation.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_interpretation

dev.off()

# ! end

