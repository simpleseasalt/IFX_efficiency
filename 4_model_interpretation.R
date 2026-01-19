# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2024-10-15 21:13:11
# | LastEditTime: 2026-01-18 23:07:12
# | FilePath: \R_scripts\4_model_interpretation.R
# | Description:
# | ---------------------------------------

library(DALEX)
library(DALEXtra)
library(randomForest)
library(tidymodels)
library(probably)
library(stacks)
library(themis)
library(ggthemes)
library(shapviz)
library(kernelshap)

load(file = "output/final_mod.RData")

load(file = "output/data_model.RData")

load(file = "output/testset.RData")

load(file = "output/IFX_workspace.RData")

# ---------------------------------------- final roc

data_stack_roc <- bind_rows(
    stack_train_pred %>% mutate(Dataset = "Train"),
    stack_test_pred %>% mutate(Dataset = "Test")
) %>%
    mutate(Dataset = as.factor(Dataset, levels = c("Train", "Test")))

plot_stack_roc <- data_stack_roc %>%
    group_by(Dataset) %>%
    roc_curve(
        group,
        .pred_0
    ) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = Dataset)) +
    geom_text(
        aes(x = 0.6, y = 0.25, label = paste("AUROC\nTrain 0.935\nTest 0.897")),
        color = "black",
        size = 5,
        hjust = 0,
        check_overlap = T
        ) +
    geom_path(lwd = 1) +
    geom_abline(lty = 3) +
    theme_bw() +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"),
        legend.position = "none"
    )

plot_stack_roc

tiff(filename = "plot/stack roc.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_stack_roc

dev.off()

# ---------------------------------------- calibration curve

plot_calibration_curve <- IFX_test %>%
    mutate(predict(final_mod, new_data = ., type = "prob")) %>%
    cal_plot_breaks(truth = group, estimate = .pred_0, include_rug = FALSE, num_breaks = 8) +
    labs(
        x = "Predicted Probability",
        y = "Observed Probability"
    ) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"),
        legend.position = "none"
    )

plot_calibration_curve

tiff(filename = "plot/calibration curve.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_calibration_curve

dev.off()

# ---------------------------------------- 

explainer_final <- DALEXtra::explain_tidymodels(
    model = final_mod,
    data = data_model,
    y = data_model$group,
    label = "Random Forest"
)

##---------------------------------------- in dataset level

variable_label <- c(
    "ESR" = "ESR",
    "CDAI_before" = "CDAI",
    "Montreal_age" = "Montreal age",
    "CRP_before" = "CRP",
    "RBC" = "RBC"
)

stack_profile <- model_profile(
    explainer_final,
    variable = names(IFX_train[, c("CRP_before", "CDAI_before", "ESR", "RBC")])
)

plot_profile <- stack_profile %>%
    plot() +
    facet_wrap(
        ~`_vname_`,
        labeller = as_labeller(variable_label),
        scales = "free_x",
        ncol = 2
    ) +
    coord_cartesian(ylim = c(0, 1.0)) +
    scale_y_continuous(breaks = seq(0.0, 1, 0.2)) +
    labs(y = "Probability of Clinical Response") +
    theme_bw() +
    theme(
        title = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = NA)
    )

plot_profile

tiff(filename = "plot/plot_stack_profile.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

plot_profile

dev.off()

## --------------------------------------- Montreal age

plot_profile_age <- model_profile(
    explainer_final,
    variable = names(IFX_train[, c("CRP_before", "CDAI_before", "ESR", "RBC")]),
    type = "partial",
    groups = "Montreal_age"
) %>%
    plot(geom = "profiles") +
    facet_wrap(
        ~`_vname_`,
        labeller = as_labeller(variable_label),
        scales = "free_x",
        ncol = 2
    ) +
    coord_cartesian(ylim = c(0, 1.0)) +
    scale_y_continuous(breaks = seq(0.0, 1, 0.2)) +
    scale_color_discrete(
        labels = c(
            "Montreal age (< 16)",
            "Montreal age (16 -40)",
            "Montreal age (≥ 40)"
        ) # 自定义标签
    ) +
    labs(y = "Probability of Clinical Response") +
    theme_bw() +
    theme(
        title = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.spacing = unit(2, "lines"),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA)
    )

plot_profile_age

tiff(filename = "plot/plot_profile_age.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_profile_age

dev.off()

##---------------------------------------- in person level

example_patient <- predict_parts(explainer_final, new_observation = data_model[, -6][1, ])

plot_explain_patient <- plot(
    example_patient,
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
        panel.spacing = unit(0.5, "lines")
    )

plot_explain_patient

tiff(filename = "plot/plot_rf_bd.tiff", width = 10, height = 6, res = 300, units = "in", compression = "lzw")

plot_explain_patient

dev.off()

# ----------------------------------------

library(cowplot)

plot_top_col <- plot_grid(plot_calibration_curve, plot_explain_patient, nrow = 2, align = , label_size = 20, labels = "AUTO", rel_heights = 1, rel_widths = 1)

plot_top_col

plot_interpretation <- plot_grid(plot_top_col, plot_profile, nrow = 1,  labels = c("", "C"), label_size = 20, rel_widths = c(1, 1.5))

plot_interpretation

tiff(filename = "plot/plot interpretation.tiff", width = 15, height = 10, res = 300, units = "in", compression = "lzw")

plot_interpretation

dev.off()

# ---------------------------------------- shap
## --------------------------------------- importance
X_explain <- data_model %>%
    select(-group)

calculate_shap <- final_mod %>%
    permshap(X = X_explain, type = "prob") %>%
    shapviz()

shap_value <- calculate_shap$.pred_0

plot_shap_importance <- shap_value %>%
    sv_importance(kind = "beeswarm") +
    theme_bw() +
    theme(
        title = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA)
    )

tiff(filename = "plot/shap importance.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_shap_importance

dev.off()

## --------------------------------------- waterfall

shap_value %>%
    sv_waterfall(
        row_id = 10,
        annotation_size = 6
        ) +
    theme_bw() +
    theme(
        title = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA)
    )




# ! end


