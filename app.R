# here: want to make an R Shiny site based on the algorithm
library(shiny)
library(shinycssloaders)
# below: codes for the project

options(spinner.color="#9dabf5", spinner.color.background="#ffffff", type = 6)

ui = fluidPage(

  titlePanel("Hamiltonian Monte Carlo"),

  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "graph_type", 
                  label = "select the type of distribution to sample from",
                  selected = "beta", 
                  choices = list("Beta" = "beta",
                                 "Normal" = "norm")),
      conditionalPanel(
        condition = "input.graph_type == 'beta'",
        numericInput(inputId = "alpha_prime", label = "alpha parameter",
                     value = 50),
        numericInput(inputId = "beta_prime", label = "beta parameter",
                     value = 50),
      ),
      
      conditionalPanel(
        condition = "input.graph_type == 'norm'",
        numericInput(inputId = "mu_mean", label = "mu parameter",
                     value = 0),
        numericInput(inputId = "sigma_sd", label = "sigma parameter",
                     value = 1),
      ),
      
      numericInput(inputId = "bigN", label = "The number of samples", 
                   value = 5000, min = 1),
      numericInput(inputId = "current_position", label = "The current position of your sample", 
                   value = 0.3),
      numericInput(inputId = "stepsize", label = "How long are you stepping to the next position?",
                   value = 0.01),
      p("Remark: for normal, use 0.2 for a reasonable result."),
      numericInput(inputId = "num_steps", label = "How many steps are you taking?",
                   value = 2),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # below the textoutput is a placeholder
      #verbatimTextOutput(outputId = "HMC_results")
      withSpinner(plotOutput(outputId = "HMC_plot"), type = 6),
      
    )
  )
)

# Define server logic required to draw a histogram ----
server = function(input, output) {

  grad_U = function(q, type, alpha = 50, beta = 50){
    if(type == "beta"){
      return(-(alpha/q - beta/(1-q)))
    } else if (type == "norm"){
      return(q)
    }
  }
  
  HMC = function(U_type, epsilon, L, current_q,
                 alpha = 50, beta = 50,
                 mu = 0, sigma = 1){
    # note: parameters past current_q are inputs for the potential energies depending on the
    # distribution the user wants to utilize
    q = current_q
    p = rnorm(length(q),0,1) # independent standard normal variates
    current_p = p
    ######################################################
    # Leapfrog Algorithm 
    # Make a half step for momentum at the beginning
    p = p - epsilon * grad_U(q, U_type, alpha, beta) / 2
    # Alternate full steps for position and momentum
    for (i in 1:L){
      # Make a full step for the position
      q = q + epsilon*p
      # Make a full step for the momentum, except at end of trajectory
      if (i!=L){ p = p - epsilon * grad_U(q, U_type, alpha, beta)}
    }
    # Make a half step for momentum at the end.
    p = p - epsilon * grad_U(q, U_type, alpha, beta) / 2
    ######################################################
    # Evaluate potential and kinetic energies at start and end of trajectory
    # switches depending on the value of U
    if(U_type == "beta"){
      current_U = -dbeta(current_q, alpha, beta, log = TRUE)
      proposed_U = -dbeta(q, alpha, beta, log = TRUE)
    } else if (U_type == "norm"){
      current_U = -dnorm(current_q, mu, sigma, log = TRUE)
      proposed_U = -dnorm(q, mu, sigma, log = TRUE)
    }
    current_K = sum(current_p^2) / 2
    proposed_K = sum(p^2) / 2
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
      return (q) # accept
    }
    else{
      return (current_q) # reject
    }
  }
  
  
  HMC_values = reactive({
    cur_q = input$current_position
    
    D = numeric(length = input$bigN)
    for(i in 1:input$bigN){
      D[i] = HMC(input$graph_type, 
                 epsilon = input$stepsize, 
                 L = input$num_steps, 
                 current_q = cur_q,
                 alpha = input$alpha_prime, 
                 beta = input$beta_prime,
                 mu = input$mu_mean, 
                 sigma = input$sigma_sd)
      cur_q = D[i]
    }
    D[500:input$bigN] # this is going to represent the chain
  })
  
  output$HMC_results = renderPrint({
    HMC_values()
  })
  
  output$HMC_plot = renderPlot({
    mean(D)
    p = seq(0,1, length=1000)
    hist(HMC_values(), breaks = 'scott', freq = FALSE,
         border = "#ffffff",
         xlab = "Values of the Chain",
         main = "Histogram of the Densities of the Chain")
    n = 100

    if(input$graph_type == "beta"){
      curve(dbeta(x,input$alpha_prime, input$beta_prime), col = 'red', add = TRUE)
    } else if (input$graph_type == "norm"){
      curve(dnorm(x,0, 1), col = 'red', add = TRUE)
    }
  })
  
}

shinyApp(ui = ui, server = server)
