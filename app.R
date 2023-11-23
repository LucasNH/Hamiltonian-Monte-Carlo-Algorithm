# here: want to make an R Shiny site based on the algorithm
library(shiny)

# below: codes for the project

U = function(q){
  return(-dbeta(q, alpha_prime, beta_prime, log = TRUE))
}

grad_U = function(q){
  return(-(alpha_prime/q - beta_prime/(1-q)))
}

HMC = function(U, grad_U, epsilon, L, current_q){
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  ######################################################
  # Leapfrog Algorithm 
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L){
    # Make a full step for the position
    q = q + epsilon*p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L){ p = p - epsilon * grad_U(q)}
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  ######################################################
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
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


ui = fluidPage(

  titlePanel("Hamiltonian Monte Carlo"),

  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "alpha_prime", label = "alpha parameter",
                  min = 1, value = 50),
      numericInput(inputId = "beta_prime", label = "beta parameter",
                  min = 1, value = 50),
      numericInput(inputId = "bigN", label = "The number of samples", 
                   value = 5000, min = 1),
      numericInput(inputId = "current_position", label = "The current position of your sample", 
                   value = 0.3),
      numericInput(inputId = "stepsize", label = "How long are you stepping to the next position?",
                   value = 0.01),
      numericInput(inputId = "num_steps", label = "How many steps are you taking?",
                   value = 2),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # below the textoutput is a placeholder
      #verbatimTextOutput(outputId = "HMC_results")
      plotOutput(outputId = "HMC_plot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  HMC_values = reactive({
    cur_q = input$current_position
    
    D = numeric(length = input$bigN)
    for(i in 1:input$bigN){
      D[i] = HMC(U, grad_U, 
                 epsilon = input$stepsize, 
                 L = input$num_steps, 
                 current_q = cur_q)
      cur_q = D[i]
    }
    D[500:N] # this is going to represent the chain
  })
  
  output$HMC_results = renderPrint({
    HMC_values()
  })
  
  output$HMC_plot <- renderPlot({
    mean(D)
    p = seq(0,1, length=1000)
    hist(HMC_values(), breaks = 'scott', freq = FALSE, xlim = c(0,1))
    #curve(dbeta(x,input$alpha_prime, input$beta_prime), col = 'red', add = TRUE)
  })
  
}

shinyApp(ui = ui, server = server)
