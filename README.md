<!--
 * #Author: Simplezzz
 * #Date: 2025-10-14 14:36:44
 * #LastEditTime: 2025-10-14 14:59:16
 * #FilePath: \R_scripts\README.md
 * #Description: This is a guide to use the final_mod.Rdata.
-->

We have deployed the online application at www.simplydeploy.work:3838/IFX_efficiency. However, due to server performance limitations (only 2 cores and 2GB RAM), we recommend downloading final_mod.Rdata to your local machine for batch predictions. Below are the usage instructions for final_mod.Rdata in a local environment:

```
if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}

library(pacman)

p_load(tidymodels, randomForest, glmnet, stacks, probably)

# setwd() to your path

load("final_mod.RData")

newdata <- data.frame(
    CDAI_before = 83.4,
    ESR = 23,
    CRP = 55.5,
    RBC = 3.89,
    Age_at_diagnosis = 2
) %>%
    rename(
        "CDAI_before" = "CDAI",
        "CRP_before" = "CRP",
        "Montreal_age" = "Age_at_diagnosis"
    ) # Make column names consistent with the variable names required in the model

predict(final_mod, newdata, type = "prob")

```

**R Environment & Packages**

-   `R-base`: 4.4.1
-   `tidymodels`: 1.2.0
-   `randomForest`: 4.7-1.2
-   `glmnet`: 4.1-8
-   `stacks`: 1.0.5
-   `probably`: 1.0.3