# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2025-08-11 12:58:07
# | LastEditTime: 2025-10-12 15:12:32
# | FilePath: \R_scripts\3_model_construction.R
# | Description: 
# | ---------------------------------------

library(tidyverse)
library(tidymodels)
library(future)
library(themis)
library(stacks)

tidymodels_prefer()

set.seed(2025)

# ---------------------------------------- load data

load("output/1_data_tidy.RData")

load("output/data_model.RData")

# ---------------------------------------- data split

IFX_split <- initial_split(data_tidy, prop = 0.70, strata = group)

IFX_train <- training(IFX_split)

IFX_test <- testing(IFX_split)

variable_include <- names(data_model)

IFX_recipe <- recipe(group ~ CDAI_before + ESR + Montreal_age + CRP_before + RBC, data = IFX_train) %>%
    step_normalize(all_numeric()) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_smote(group, over_ratio = 1, seed = 2025)


save(IFX_train, file = "output/trainset.RData")
save(IFX_test, file = "output/testset.RData")

IFX_recipe_prep <- IFX_recipe %>%
    prep()

save(IFX_recipe_prep, file = "output/IFX_recipe.RData")

# ---------------------------------------- model settings

mod_plr <- logistic_reg(penalty = tune()) %>%
    set_engine("glmnet") %>%
    set_mode("classification")

mod_svm <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")

mod_rf <- rand_forest(min_n = tune(), trees = tune()) %>%
    set_engine("randomForest") %>%
    set_mode("classification")

mod_xgb <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    min_n = tune(),
    sample_size = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
) %>%
    set_engine("xgboost") %>%
    set_mode("classification")

## --------------------------------------- metrics

metrics_result <- metric_set(roc_auc, accuracy, precision, sensitivity, specificity, recall, f_meas)

## ---------------------------------------- cross validaion

IFX_cv <- vfold_cv(IFX_train, strata = "group", repeats = 1, v = 10)

## --------------------------------------- define multiprocess

plan(multisession, workers = parallel::detectCores() - 2)

## ---------------------------------------- workflow

IFX_wf <- workflow_set(
    preproc = list(
        recipe = IFX_recipe
    ),
    models = list(
        SVM = mod_svm,
        PLR = mod_plr,
        RF = mod_rf,
        XGB = mod_xgb
    ),
    cross = FALSE
) %>%
    mutate(
        wflow_id = case_when(
            wflow_id == "recipe_SVM" ~ "SVM",
            wflow_id == "recipe_PLR" ~ "PLR",
            wflow_id == "recipe_RF" ~ "RF",
            wflow_id == "recipe_XGB" ~ "XGBoost",
        )
    )

tune_control <- control_bayes(
    allow_par = TRUE,
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
)

grid_result <- workflow_map(
    IFX_wf,
    fn = "tune_bayes",
    resamples = IFX_cv,
    metrics = metrics_result,
    verbose = FALSE,
    iter = 25,
    initial = 10,
    seed = 2025,
    control = tune_control
)

autoplot(
    grid_result,
    rank_metric = "roc_auc"
) +
    theme_bw()

updated_mods_plot <- autoplot(
    grid_result,
    metric = "roc_auc",
    rank_metric = "roc_auc"
) +
    labs(x = "Model Rank", y = "AUROC") +
    scale_shape_discrete(guide = "none") +
    scale_color_discrete(name = "Models", label = c("XGBoost", "PLR", "RF", "SVM")) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"),
        legend.position = "right"
    )

updated_mods_plot

tiff(filename = "plot/updated_mods_plot.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

updated_mods_plot

dev.off()

## ---------------------------------------- export result
### ---------------------------------------- glm

PLR_tune <- grid_result %>%
    extract_workflow_set_result("PLR") %>%
    unnest(cols = .metrics) %>%
    filter(.metric == "roc_auc") %>%
    group_by(penalty) %>%
    summarise(
        mean = mean(.estimate),
        sd = sd(.estimate)
    )

write.csv(PLR_tune, "output/PLR_tune.csv")

### ---------------------------------------- svm

SVM_tune <- grid_result %>%
    extract_workflow_set_result("SVM") %>%
    unnest(cols = .metrics) %>%
    filter(.metric == "roc_auc") %>%
    group_by(cost) %>%
    summarise(
        mean = mean(.estimate),
        sd = sd(.estimate)
    )

write.csv(SVM_tune, "output/SVM_tune.csv")

### ---------------------------------------- rf

RF_tune <- grid_result %>%
    extract_workflow_set_result("RF") %>%
    unnest(cols = .metrics) %>%
    filter(.metric == "roc_auc") %>%
    group_by(min_n, trees) %>%
    summarise(
        mean = mean(.estimate),
        sd = sd(.estimate)
    )

write.csv(RF_tune, "output/RF_tune.csv")

### ---------------------------------------- xgb

XGB_tune <- grid_result %>%
    extract_workflow_set_result("XGBoost") %>%
    unnest(cols = .metrics) %>%
    filter(.metric == "roc_auc") %>%
    group_by(min_n, tree_depth, learn_rate, loss_reduction) %>%
    summarise(
        mean = mean(.estimate),
        sd = sd(.estimate)
    )

write.csv(XGB_tune, "output/XGB_tune.csv")

## ---------------------------------------- plot AUROC

mod_best_roc <- grid_result %>%
    collect_metrics() %>%
    filter(.metric == "roc_auc") %>%
    group_by(wflow_id) %>%
    arrange(-mean) %>%
    slice(1) %>%
    mutate(best_model = paste(wflow_id, .config, sep = "_"))

updated_result <- grid_result %>%
    collect_metrics() %>%
    mutate(best_model = paste(wflow_id, .config, sep = "_")) %>%
    filter(best_model %in% mod_best_roc$best_model)

updated_result %>%
    select(wflow_id, metrics = .metric, mean, std_err) %>%
    write_csv("output/updated_result.csv")

best_mods <- updated_result %>%
    mutate(model_name = paste(.config, model, sep = "_")) %>%
    pull(model_name)

prediction_in_best_mods <- grid_result %>%
    workflowsets::collect_predictions() %>%
    mutate(model_name = paste(.config, model, sep = "_")) %>%
    filter(
        model_name %in% best_mods
    )

updated_roc_plot <- prediction_in_best_mods %>%
    group_by(wflow_id) %>%
    roc_curve(
        group,
        .pred_0
    ) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = wflow_id)) +
    geom_path(lwd = 1) +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_bw() +
    scale_color_discrete(name = "Models") +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"),
        legend.position = "right"
    )

tiff(filename = "plot/updated_roc_plot.tiff", width = 10, height = 6.6, res = 300, units = "in", compression = "lzw")

updated_roc_plot

dev.off()

# ---------------------------------------- fit on test data

PLR_fit <- finalize_workflow(
    extract_workflow(grid_result, id = "PLR"),
    select_best(
        grid_result[grid_result$wflow_id == "PLR", "result"][[1]][[1]],
        metric = "roc_auc"
    )
) %>%
    last_fit(
        split = IFX_split,
        metrics = metrics_result
    )

XGB_fit <- finalize_workflow(
    extract_workflow(grid_result, id = "XGBoost"),
    select_best(
        grid_result[grid_result$wflow_id == "XGBoost", "result"][[1]][[1]],
        metric = "roc_auc"
    )
) %>%
    last_fit(
        split = IFX_split,
        metrics = metrics_result
    )

RF_fit <- finalize_workflow(
    extract_workflow(grid_result, id = "RF"),
    select_best(
        grid_result[grid_result$wflow_id == "RF", "result"][[1]][[1]],
        metric = "roc_auc"
    )
) %>%
    last_fit(
        split = IFX_split,
        metrics = metrics_result
    )

SVM_fit <- finalize_workflow(
    extract_workflow(grid_result, id = "SVM"),
    select_best(
        grid_result[grid_result$wflow_id == "SVM", "result"][[1]][[1]],
        metric = "roc_auc"
    )
) %>%
    last_fit(
        split = IFX_split,
        metrics = metrics_result
    )

last_fit_list <- list(
    PLR_fit,
    RF_fit,
    SVM_fit,
    XGB_fit
)

last_fit_res <- last_fit_list %>%
    map(collect_metrics) %>%
    rlist::list.stack() %>%
    mutate(model = c(
        rep("PLR", nrow(.) / 4),
        rep("RF", nrow(.) / 4),
        rep("SVM", nrow(.) / 4),
        rep("XGBoost", nrow(.) / 4)
    )) %>%
    select(-c(.estimator, .config)) %>%
    pivot_wider(
        names_from = .metric,
        values_from = .estimate
    )

write_csv(last_fit_res, "output/final_result_in_testset.csv")

## --------------------------------------- combine roc plot

auroc_train <- grid_result %>%
    collect_metrics() %>%
    filter(.metric == "roc_auc") %>%
    group_by(wflow_id) %>%
    arrange(desc(mean)) %>%
    dplyr::slice(1)

combine_roc_label <- data.frame(
    group = auroc_train$wflow_id,
    roc_train = round(auroc_train$mean, 3),
    roc_test = round(as.numeric(last_fit_res$roc_auc), 3)
) %>%
    as_tibble() %>%
    mutate(across(c(roc_train, roc_test), format, nsmall = 3)) %>%
    mutate(
        label = paste("AUROC\n", "Trainset", .$roc_train, "\n", "Testset", .$roc_test)
    ) %>%
    rename("wflow_id" = "group")

data_train_roc <- prediction_in_best_mods %>%
    group_by(wflow_id) %>%
    roc_curve(
        group,
        .pred_0
    ) %>%
    mutate(
        split = "train"
    )

data_test_roc <- last_fit_list %>%
    map(collect_predictions) %>%
    map(roc_curve, group, .pred_0) %>%
    data.table::rbindlist(idcol = "wflow_id") %>%
    as_tibble() %>%
    mutate(
        wflow_id = case_when(
            wflow_id == 1 ~ "PLR",
            wflow_id == 2 ~ "XGBoost",
            wflow_id == 3 ~ "RF",
            wflow_id == 4 ~ "SVM",
        ),
        split = "test",
    )

final_roc_plot <- bind_rows(data_train_roc, data_test_roc) %>%
    left_join(combine_roc_label) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = split)) +
    geom_path(lwd = 1) +
    geom_abline(lty = 3) +
    facet_wrap(~wflow_id) +
    geom_text(
        aes(x = 0.6, y = 0.25, label = label),
        color = "black",
        size = 5,
        hjust = 0,
        check_overlap = T
    ) +
    coord_equal() +
    theme_bw() +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    scale_color_discrete(name = "Dataset", label = c("Test", "Train")) +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        panel.spacing.x = unit(1, "lines")
    )

final_roc_plot

tiff(filename = "plot/final_roc_plot.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

final_roc_plot

dev.off()

## ---------------------------------------- confuse matrix
### ---------------------------------------- in trainset

train_fit <- prediction_in_best_mods %>%
    select(c(wflow_id, group, .pred_0)) %>%
    mutate(split = "train")

train_fit %>%
    mutate(.pred_res = factor(if_else(.pred_0 >= 0.5, "0", "1"))) %>%
    group_by(wflow_id) %>%
    conf_mat(group, .pred_res) %>%
    as.data.frame()

# ---------------------------------------- stack model

stack_wf <- workflow_set(
    preproc = list(
        recipe = IFX_recipe
    ),
    models = list(
        PLR = mod_plr,
        RF = mod_rf
    ),
    cross = FALSE
) %>%
    mutate(
        wflow_id = case_when(
            wflow_id == "recipe_PLR" ~ "PLR",
            wflow_id == "recipe_RF" ~ "RF"
        )
    )

stack_result <- workflow_map(
    stack_wf,
    fn = "tune_bayes",
    resamples = IFX_cv,
    metrics = metrics_result,
    verbose = FALSE,
    iter = 25,
    initial = 10,
    control = tune_control
)

IFX_stack <- stacks() %>%
    add_candidates(extract_workflow_set_result(grid_result, "PLR")) %>%
    add_candidates(extract_workflow_set_result(grid_result, "RF"))

IFX_stack_fit <- IFX_stack %>%
    blend_predictions(
        metric = metrics_result
    ) %>%
    fit_members()

stack_train_pred <- predict(IFX_stack_fit, IFX_train) %>%
    bind_cols(predict(IFX_stack_fit, IFX_train, type = "prob")) %>%
    bind_cols(IFX_train %>% select(group)) %>%
    mutate(wflow_id = "Stack model")

stack_test_pred <- predict(IFX_stack_fit, IFX_test) %>%
    bind_cols(predict(IFX_stack_fit, IFX_test, type = "prob")) %>%
    bind_cols(IFX_test %>% select(group)) %>%
    mutate(wflow_id = "Stack model")

train_metrics <- metrics_result(stack_train_pred, truth = group, estimate = .pred_class, .pred_0) %>%
    mutate(
        data_type = "train",
        wflow_id = "Stack model"
    )

test_metrics <- metrics_result(stack_test_pred, truth = group, estimate = .pred_class, .pred_0) %>%
    mutate(
        data_type = "test",
        wflow_id = "Stack model"
    )

all_metrics <- bind_rows(train_metrics, test_metrics) %>%
    pivot_wider(
        names_from = ".metric",
        values_from = ".estimate"
    ) %>%
    select(-.estimator) %>%
    mutate(across(accuracy:roc_auc, ~ format(., digits = 3)))

all_metrics

write.csv(all_metrics, file = "output/result_stack_model.csv")

# --------------------------------------- output all result
## --------------------------------------- metrics

metrics_train_all <- bind_rows(
    updated_result %>%
        select(wflow_id, .metric, mean) %>%
        group_by(wflow_id) %>%
        pivot_wider(
            names_from = ".metric",
            values_from = "mean"
        ),
    train_metrics %>%
        select(wflow_id, .metric, mean = .estimate) %>%
        pivot_wider(
            names_from = ".metric",
            values_from = "mean"
        )
) %>%
    mutate(Dataset = "Train")

metrics_test_all <- bind_rows(
    last_fit_res %>%
        rename("wflow_id" = "model"),
    test_metrics %>%
        select(wflow_id, .metric, mean = .estimate) %>%
        pivot_wider(
            names_from = ".metric",
            values_from = "mean"
        )
) %>%
    mutate(Dataset = "Test")

final_result <- bind_rows(
    metrics_train_all,
    metrics_test_all
) %>%
    relocate(
        wflow_id,
        Dataset
    ) %>%
    relocate(
        f_meas, roc_auc,
        .after = last_col()
    )

final_result %>%
    rename(
        "Model" = "wflow_id",
        "Accuracy" = "accuracy",
        "Precision" = "precision",
        "Sensitivity" = "sensitivity",
        "Specificity" = "specificity",
        "Recall" = "recall",
        "F score" = "f_meas",
        "AUROC" = "roc_auc"
    ) %>%
    mutate(across(c(Accuracy:AUROC), format, nsmall = 3, digits = 3)) %>%
    write_csv(final_result, "output/final_result.csv")

## --------------------------------------- roc plot
### -------------------------------------- trainset

roc_labels <- final_result %>%
    select(wflow_id, Dataset, roc_auc) %>%
    pivot_wider(
        names_from = Dataset,
        values_from = roc_auc
    ) %>%
    mutate(across(c(Train, Test), format, nsmall = 3, digits = 3)) %>%
    mutate(
        label_train = paste(
            wflow_id,
            Train
        ),
        label_test = paste(
            wflow_id,
            Test
        )
    )

label_train <- str_c(
    "AUROC\n",
    str_c(roc_labels$label_train, collapse = "\n"),
    collapse = ""
    ) %>%
    as_tibble()

plot_roc_train <- prediction_in_best_mods %>%
    select(wflow_id, names(stack_train_pred)) %>%
    bind_rows(stack_train_pred) %>%
    group_by(wflow_id) %>%
    roc_curve(
        group,
        .pred_0
    ) %>%
    left_join(roc_labels) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = wflow_id)) +
    geom_path(lwd = 1) +
    geom_abline(lty = 3) +
    geom_text(
        data = label_train,
        aes(x = 0.6, y = 0.25, label = value),
        color = "black",
        size = 6,
        hjust = 0,
        check_overlap = T
    ) +
    coord_equal() +
    theme_bw() +
    scale_color_discrete(name = "Models") +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"),
        legend.position = "right"
    )

plot_roc_train

tiff(filename = "plot/plot_roc_train.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_roc_train

dev.off()

## --------------------------------------- testset

label_test <- str_c(
    "AUROC\n",
    str_c(roc_labels$label_test, collapse = "\n"),
    collapse = ""
) %>%
    as_tibble()

plot_roc_test <- prediction_in_best_mods %>%
    select(wflow_id, names(stack_test_pred)) %>%
    bind_rows(stack_test_pred) %>%
    group_by(wflow_id) %>%
    roc_curve(
        group,
        .pred_0
    ) %>%
    left_join(roc_labels) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = wflow_id)) +
    geom_path(lwd = 1) +
    geom_abline(lty = 3) +
    geom_text(
        data = label_test,
        aes(x = 0.6, y = 0.25, label = value),
        color = "black",
        size = 6,
        hjust = 0,
        check_overlap = T
    ) +
    coord_equal() +
    theme_bw() +
    scale_color_discrete(name = "Models") +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme(
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"),
        legend.position = "right"
    )

plot_roc_test

tiff(filename = "plot/plot_roc_test.tiff", width = 10, height = 10, res = 300, units = "in", compression = "lzw")

plot_roc_test

dev.off()

# ---------------------------------------- save stack model

final_mod <- IFX_stack_fit

save(final_mod, file = "output/final_mod.Rdata")

predict(final_mod, new_data = IFX_test[1, ], type = "prob")

# ---------------------------------------- model fairness

rf_test_predictions <- IFX_test %>%
    mutate(
        predict(final_mod, new_data = ., type = "prob"),
        predict(final_mod, new_data = ., type = "class"),
        .pred_class = factor(.pred_class, levels = c("0", "1")),
        Montreal_age = factor(Montreal_age, levels = c("1", "2", "3")),
        CDAI_group = case_when(
            CDAI_before < 150 ~ "Remission",
            CDAI_before >= 150 & CDAI_before < 220 ~ "Mild",
            CDAI_before >= 220 & CDAI_before < 450 ~ "Moderate",
            TRUE ~ "Severe"
        ),
        CDAI_group = factor(CDAI_group, levels = c("Remission", "Mild", "Moderate", "Severe"))
    )

## --------------------------------------- age

metrics_by_age <- rf_test_predictions %>%
    group_by(Montreal_age) %>%
    metrics_result(truth = group, estimate = .pred_class, .pred_0) %>%
    pivot_wider(
        names_from = ".metric",
        values_from = ".estimate"
    ) %>%
    rename("group" = Montreal_age)

## --------------------------------------- gender

metrics_by_gender <- rf_test_predictions %>%
    group_by(gender) %>%
    metrics_result(truth = group, estimate = .pred_class, .pred_0) %>%
    pivot_wider(
        names_from = ".metric",
        values_from = ".estimate"
    ) %>%
    rename("group" = gender) %>%
    mutate(
        group = case_when(
            group == 1 ~ "Male",
            .default = "Female"
        )
    )

## --------------------------------------- CDAI

metrics_by_CDAI <- rf_test_predictions %>%
    group_by(CDAI_group) %>%
    metrics_result(truth = group, estimate = .pred_class, .pred_0) %>%
    pivot_wider(
        names_from = ".metric",
        values_from = ".estimate"
    ) %>%
    rename("group" = CDAI_group)

## --------------------------------------- combine and out put

subgroup_number <- bind_cols(
        `Total number` = bind_rows(
            count(rf_test_predictions, Montreal_age) %>% select(n),
            count(rf_test_predictions, gender) %>% select(n),
            count(rf_test_predictions, CDAI_group) %>% select(n)
        ),
        `Clinical response` = bind_rows(
            count(rf_test_predictions, Montreal_age, group) %>% filter(group == 0) %>% select(n),
            count(rf_test_predictions, gender, group) %>% filter(group == 0) %>% select(n),
            count(rf_test_predictions, CDAI_group, group) %>% filter(group == 0) %>% select(n),
        )
    )

names(subgroup_number) <- c("Total number", "Clinical response")

subgroup_result <- bind_rows(metrics_by_age, metrics_by_gender, metrics_by_CDAI) %>%
    mutate(across(where(is.numeric), format, digits = 3)) %>%
    select(-.estimator) %>%
    bind_cols(subgroup_number) %>%
    mutate(
        group = case_when(
            group == 1 ~ "Age ≤ 16",
            group == 2 ~ "Age (16 - 40)",
            group == 3 ~ "Age ≥ 40",
            .default = group
        )
    ) %>%
    relocate(c(`Total number`, `Clinical response`), .after = group) %>%
    rename(
        "Accuracy" = "accuracy",
        "Precision" = "precision",
        "Sensitivity" = "sensitivity",
        "Specificity" = "specificity",
        "Recall" = "recall",
        "F score" = "f_meas",
        "AUROC" = "roc_auc"
    )

write_csv(subgroup_result, file = "output/subgroup_result.csv")

# ! end
