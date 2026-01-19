# | ---------------------------------------
# | Author: Simplezzz
# | Date: 2025-08-21 21:42:08
# | LastEditTime: 2025-10-17 13:15:02
# | FilePath: \R_scripts\6_deployment.R
# | Description: 
# | ---------------------------------------

library(shiny)
library(tidyverse)
library(tidymodels)
library(shinydashboard)
library(DT)
library(DALEX)
library(DALEXtra)
library(stacks)
library(probably)

load(file = "output/IFX_recipe.RData")
load(file = "output/final_mod.RData")

# 定义UI
ui <- dashboardPage(
    dashboardHeader(title = "Clinical Response Prediction Model"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Single Prediction", tabName = "single", icon = icon("user")),
            menuItem("Batch Prediction", tabName = "batch", icon = icon("table")),
            menuItem("Model Information", tabName = "model", icon = icon("info-circle"))
        )
    ),
    dashboardBody(
        tabItems(
            # Single Prediction Tab
            tabItem(
                tabName = "single",
                fluidRow(
                    box(
                        title = "Input Predictor Variables", status = "primary", solidHeader = TRUE,
                        width = 6,
                        numericInput("cdai", "CDAI", value = 83.4, min = 0, max = 100),
                        numericInput("esr", "ESR (mm/h):", value = 23, min = 0, max = 150),
                        numericInput("crp", "CRP (mg/L):", value = 55.5, min = 0, max = 200),
                        numericInput("rbc", "RBC (10^12/L):", value = 3.89, min = 1, max = 8),
                        selectInput("age", "Age at Diagnosis:",
                            choices = list(
                                `≤ 16` = "1",
                                `16-40` = "2",
                                `≥ 40` = "3"
                            ),
                            selected = "2"
                        ),
                        actionButton("predict", "Predict", class = "btn-success")
                    ),
                    box(
                        title = "Prediction Results", status = "info", solidHeader = TRUE,
                        width = 6,
                        h4("Clinical Response Prediction:"),
                        verbatimTextOutput("prediction_result"),
                        h4("Prediction Probability:"),
                        plotOutput("probability_plot", height = "200px"),
                        downloadButton("download_single", "Download Prediction Results")
                    )
                )
            ),

            # Batch Prediction Tab
            tabItem(
                tabName = "batch",
                fluidRow(
                    box(
                        title = "Upload Data", status = "primary", solidHeader = TRUE,
                        width = 12,
                        fileInput("file_upload", "Choose CSV File",
                            accept = c(".csv")
                        ),
                        helpText("File should contain the following columns: CDAI, ESR, CRP, RBC, Age_at_diagnosis"),
                        helpText("Age is encoded as '1', '2' and '3', representing '≤ 16', '16 - 40' and '≥ 40' "),
                        helpText("Age is encoded as '1', '2' and '3', representing '≤ 16', '16 - 40' and '≥ 40' "),
                        actionButton("predict_batch", "Batch Predict", class = "btn-success")
                    )
                ),
                fluidRow(
                    box(
                        title = "Prediction Results", status = "info", solidHeader = TRUE,
                        width = 12,
                        DTOutput("batch_results"),
                        downloadButton("download_batch", "Download Results")
                    )
                )
            ),

            # Model Information Tab
            tabItem(
                tabName = "model",
                fluidRow(
                    box(
                        title = "Model Features", status = "info", solidHeader = TRUE,
                        width = 12,
                        h4("Predictor Variables:"),
                        tags$ul(
                            tags$li("CDAI - Clinical Disease Activity Index"),
                            tags$li("ESR - Erythrocyte Sedimentation Rate"),
                            tags$li("CRP - C-Reactive Protein"),
                            tags$li("RBC - Red Blood Cell Count"),
                            tags$li("Age at diagnosis - Categorical variable with three levels:")
                        ),
                        tags$ul(
                            tags$li("≤16 years"),
                            tags$li("16-40 years"),
                            tags$li("≥40 years")
                        ),
                        h4("Outcome Variable:"),
                        p("Clinical Response - Clinical treatment response"),
                        h4("Model Performance Metrics:"),
                        tableOutput("model_metrics"),
                        h4("Age Category Distribution in Training Data:"),
                        plotOutput("age_distribution_plot")
                    )
                )
            )
        )
    )
)

server <- function(input, output, session) {
    # Reactive value storage
    prediction_result <- reactiveVal()
    input_data_raw <- reactiveVal()
    batch_results <- reactiveVal()

    # Single prediction
    observeEvent(input$predict, {
        # Create input data frame with categorical age
        input_data <- tibble(
            CDAI_before = input$cdai,
            ESR = input$esr,
            CRP_before = input$crp,
            RBC = input$rbc,
            Montreal_age = factor(input$age, levels = c("1", "2", "3"))
        )

        # Preprocess the data
        input_data_processed <- input_data %>%
            bake(IFX_recipe_prep, new_data = .) %>%
            mutate(
                Montreal_age = case_when(
                    Montreal_age_X1 == 1 ~ 1,
                    Montreal_age_X3 == 3 ~ 3,
                    .default = 2
                ),
                Montreal_age = as.factor(Montreal_age)
            )

        # Use model for prediction
        prediction <- tryCatch(
            {
                tibble(
                    .pred_class = predict(final_mod, new_data = input_data_processed, type = "class")$.pred_class,
                    .pred_Response = predict(final_mod, new_data = input_data_processed, type = "prob")$.pred_1,
                    .pred_Non_response = 1 - .pred_Response
                ) %>%
                    mutate(.pred_class = ifelse(.pred_class == 1, "Response", "Non-response"))
            },
            error = function(e) {
                tibble(
                    .pred_class = "Error",
                    .pred_Response = NA,
                    .pred_Non_response = NA
                )
            }
        )

        # Store prediction results
        prediction_result(prediction)
        input_data_raw(input_data)
    })

    # Display prediction results
    output$prediction_result <- renderPrint({
        req(prediction_result())
        pred <- prediction_result()
        cat("Predicted Class:", pred$.pred_class, "\n")
        cat("Response Probability:", round(pred$.pred_Response * 100, 2), "%\n")
        cat("Non-response Probability:", round(pred$.pred_Non_response * 100, 2), "%")
    })

    # Probability plot
    output$probability_plot <- renderPlot({
        req(prediction_result())
        pred <- prediction_result()

        prob_data <- data.frame(
            Category = c("Response", "Non-response"),
            Probability = c(pred$.pred_Response, pred$.pred_Non_response)
        )

        ggplot(prob_data, aes(x = Category, y = Probability, fill = Category)) +
            geom_col(alpha = 0.8) +
            geom_text(aes(label = paste0(round(Probability * 100, 1), "%")),
                vjust = -0.5, size = 5
            ) +
            scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
            ylim(0, 1) +
            theme_minimal() +
            theme(
                legend.position = "none",
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)
            )
    })

    # Batch prediction
    observeEvent(input$predict_batch, {
        req(input$file_upload)

        # Read uploaded file
        data <- read_csv(input$file_upload$datapath)

        # Validate columns' name
        required_cols <- c("CDAI", "ESR", "CRP", "RBC", "Age_at_diagnosis")
        if (!all(required_cols %in% names(data))) {
            showNotification("File is missing required columns!", type = "error")
            return()
        }

        # Validate age categories
        valid_age_categories <- c("1", "2", "3")
        invalid_ages <- setdiff(unique(data$Age_at_diagnosis), valid_age_categories)

        if (length(invalid_ages) > 0) {
            showNotification(paste("Invalid age categories found:", paste(invalid_ages, collapse = ", ")),
                type = "error"
            )
            return()
        }

        # Convert Age_at_diagnosis to factor and preprocess
        data_processed <- data %>%
            rename(
                "CDAI_before" = "CDAI",
                "CRP_before" = "CRP",
                "Montreal_age" = "Age_at_diagnosis"
            ) %>%
            mutate(Montreal_age = factor(Montreal_age, levels = c("1", "2", "3"))) %>%
            bake(IFX_recipe_prep, new_data = .) %>%
            mutate(
                Montreal_age = case_when(
                    Montreal_age_X1 == 1 ~ 1,
                    Montreal_age_X3 == 3 ~ 3,
                    .default = 2
                ),
                Montreal_age = as.factor(Montreal_age)
            )

        # Batch prediction
        predictions <- tryCatch(
            {
                tibble(
                    .pred_class = predict(final_mod, new_data = data_processed, type = "class")$.pred_class,
                    .pred_Response = predict(final_mod, new_data = data_processed, type = "prob")$.pred_1,
                    .pred_Non_response = 1 - .pred_Response
                ) %>%
                    mutate(.pred_class = ifelse(.pred_class == 1, "Response", "Non-response"))
            },
            error = function(e) {
                NULL
            }
        )

        if (!is.null(predictions)) {
            results <- bind_cols(data, predictions)
            batch_results(results)
        }
    })

    # Display batch results
    output$batch_results <- renderDT({
        req(batch_results())
        datatable(
            batch_results() %>%
                mutate(across(where(is.numeric), round, 3)),
            options = list(scrollX = TRUE, pageLength = 10)
        )
    })

    # Download single prediction results
    output$download_single <- downloadHandler(
        filename = function() {
            paste0("prediction_", Sys.Date(), ".csv")
        },
        content = function(file) {
            req(prediction_result(), input_data_raw())
            result <- bind_cols(input_data_raw(), prediction_result())
            write_csv(result, file)
        }
    )

    # Download batch results
    output$download_batch <- downloadHandler(
        filename = function() {
            paste0("batch_predictions_", Sys.Date(), ".csv")
        },
        content = function(file) {
            req(batch_results())
            write_csv(batch_results(), file)
        }
    )
    # Delete the temporary files uploaded in this session
    session$onSessionEnded(function() {
        if (!is.null(input$userFile)) {
            unlink(input$userFile$datapath, force = TRUE)
        }
    })
}

shinyApp(ui = ui, server = server)
