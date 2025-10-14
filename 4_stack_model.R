# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2025-09-28 22:24:22
# | LastEditTime: 2025-09-30 13:01:23
# | FilePath: \R_scripts\4_stack_model.R
# | Description:
# | ---------------------------------------

library(tidymodels)
library(stacks)
library(future)
library(themis)

tidymodels_prefer()

set.seed(2025)

# ---------------------------------------- 

load("output/1_data_tidy.RData")

load("output/data_model.RData")

# ---------------------------------------- 

plan(multisession, workers = parallel::detectCores() - 2)

# ---------------------------------------- data split

IFX_split <- initial_split(data_tidy, prop = 0.70, strata = group)

IFX_train <- training(IFX_split)

IFX_test <- testing(IFX_split)

variable_include <- names(data_model)

IFX_recipe <- recipe(group ~ CDAI_before + ESR + Montreal_age + CRP_before + RBC, data = IFX_train) %>%
    step_normalize(all_numeric()) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_smote(group, over_ratio = 1, seed = 2025)

## ---------------------------------------- model settings

mod_plr <- logistic_reg(penalty = tune()) %>%
    set_engine("glmnet") %>%
    set_mode("classification")

mod_rf <- rand_forest(min_n = tune(), trees = tune()) %>%
    set_engine("randomForest") %>%
    set_mode("classification")

## --------------------------------------- metrics

metrics_result <- metric_set(roc_auc, accuracy, precision, sensitivity, specificity, recall, f_meas)

## ---------------------------------------- cross validaion

IFX_cv <- vfold_cv(IFX_train, strata = "group", repeats = 1, v = 10)

## ---------------------------------------- workflow

IFX_wf <- workflow_set(
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

IFX_stack <- stacks() %>%
    add_candidates(grid_result)

IFX_stack_fit <- IFX_stack %>%
    blend_predictions(
        metric = metrics_result,
        penalty = 10^(-6:-1)
    ) %>%
    fit_members()

train_pred <- predict(IFX_stack_fit, IFX_train) %>%
    bind_cols(predict(IFX_stack_fit, IFX_train, type = "prob")) %>%
    bind_cols(IFX_train %>% select(group))

test_pred <- predict(IFX_stack_fit, IFX_test) %>%
    bind_cols(predict(IFX_stack_fit, IFX_test, type = "prob")) %>%
    bind_cols(IFX_test %>% select(group))

train_metrics <- metrics_result(train_pred, truth = group, estimate = .pred_class, .pred_0) %>%
    mutate(data_type = "train")

test_metrics <- metrics_result(test_pred, truth = group, estimate = .pred_class, .pred_0) %>%
    mutate(data_type = "test")

all_metrics <- bind_rows(train_metrics, test_metrics) %>%
    pivot_wider(
        names_from = ".metric",
        values_from = ".estimate"
    ) %>%
    select(-.estimator) %>%
    rename(
        "Dataset" = "data_type",
        "Accuracy" = "accuracy",
        "Precision" = "precision",
        "Sensitivity" = "sensitivity",
        "Specificity" = "specificity",
        "Recall" = "recall",
        "F score" = "f_meas",
        "AUROC" = "roc_auc"
    ) %>%
    mutate(across(Accuracy:`F score`, ~ format(., digits = 3)))

write.csv(all_metrics, file = "output/result_stack_model.csv")

autoplot(IFX_stack_fit) +
    facet_wrap(vars(.metric), ncol = 2, scales = "free") +
    theme_bw() +
    labs(title = "混合模型成员权重")

autoplot(IFX_stack_fit, type = "members") +
    facet_wrap(vars(.metric), ncol = 2, scales = "free") +
    theme_bw() +
    labs(title = "成员模型性能比较")

autoplot(IFX_stack_fit, type = "performance") +
    theme_bw() +
    labs(title = "混合模型性能")

save(IFX_stack_fit, file = "output/IFX_stack_fit.Rdata")

# ! end

# ---------------------------------------- roc plot



