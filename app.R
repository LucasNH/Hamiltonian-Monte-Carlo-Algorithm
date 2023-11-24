library(shiny)
library(shinycssloaders)
options(spinner.color = "#9dabf5", spinner.color.background = "#ffffff",
        type = 6)

source("HMC_algorithm.R")

HMC <- HMC

# This is the front-end - user inputs and what are displayed.
ui <- fluidPage(
  titlePanel("Hamiltonian Monte Carlo"),
  h4("Authors from MAT332: Mark Asuncion, Simran Bilkhu, Anna Ly," +
       " Lucas Noritomi-Hartwig"),
  h4("Website is used for demonstration purposes only."),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "graph_type",
                  label = "select the type of distribution to sample from",
                  selected = "beta",
                  choices = list("Beta" = "beta", "Normal" = "norm")),
      conditionalPanel(condition = "input.graph_type == 'beta'",
        numericInput(inputId = "alpha_prime", label = "alpha parameter",
                     value = 50),
        numericInput(inputId = "beta_prime", label = "beta parameter",
                     value = 50),
      ),

      conditionalPanel(
        condition = "input.graph_type == 'norm'",
        numericInput(inputId = "mu_mean", label = "mu parameter", value = 0),
        numericInput(inputId = "sigma_sd", label = "sigma parameter",
                     value = 1),
      ),

      numericInput(inputId = "bigN", label = "The number of samples",
                   value = 5000, min = 1),
      numericInput(inputId = "current_position",
                   label = "The current position of your sample", value = 0.3),
      numericInput(inputId = "stepsize",
                   label = "How long are you stepping to the next position?",
                   value = 0.01),
      strong("Remarks"),
      p("For for beta, 0.01 leads to a reasonable result."),
      p("For normal, use 0.2 leads to a reasonable result."),
      numericInput(inputId = "num_steps",
                   label = "How many steps are you taking?", value = 2),
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      # below the textoutput is a placeholder
      # verbatimTextOutput(outputId = "HMC_results")
      withSpinner(plotOutput(outputId = "HMC_plot"), type = 6),
    )
  )
)

# This is the back-end.
server <- function(input, output) {
  HMC_values <- reactive({
    HMC(n_samples = input$bigN, U_type = input$graph_type,
        epsilon = input$stepsize, L = input$num_steps,
        current_q = input$current_position, alpha = input$alpha_prime,
        beta = input$beta_prime, mu = input$mu_mean, sigma = input$sigma_sd)
  })

  output$HMC_results <- renderPrint({HMC_values()})

  output$HMC_plot <- renderPlot({
    hist(HMC_values(), breaks = "scott", freq = FALSE, border = "#ffffff",
         xlab = "Values of the Chain",
         main = "Histogram of the Densities of the Chain")

    if("beta" == input$graph_type) {
      curve(dbeta(x, input$alpha_prime, input$beta_prime), col = "#FF6666",
            add = TRUE, lwd = 2)
    } else if ("norm" == input$graph_type) {
      curve(dnorm(x, 0, 1), col = "#6699FF", add = TRUE, lwd = 2)
    }
  })
}

shinyApp(ui = ui, server = server)
