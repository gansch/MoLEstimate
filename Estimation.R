library(shiny)
library(ggplot2)
library(DT)
library(plotly)
############################## Functions that we need later ####################################

# Function for Likelihood of Binomial distribution
llbinom.grad <- function(theta, n, y){
  #loglikelihood
  ell <- y*log(theta) + n*log(1-theta)
  #first derivative
  ellp <- y/theta - n/(1-theta)
  #seconds derivative
  ellpp <- -y/theta^2 - n/(1-theta)^2
  out<-list(l=ell,lp=ellp,lpp=ellpp)
  return(out)
}

# Function that searches for Maximum Likelihood: Iterative Method
MaxBinom <- function(theta, alpha, delta = 10^-6, n, y){
  if (20*alpha > theta) {
    stop("You must enter a bigger value for the learning rate or a smaller theta")
  }
  # define empty vectors
  thetaVec <- NULL
  LL <- NULL
  LLp <- NULL
  LLpp <- NULL
  
  # assign theta
  thetaVec[1] <- theta
  
  # calculate likelihoods
  out <- llbinom.grad(theta = thetaVec[1], n = n, y = y)
  
  # enter result of the list into the vectors
  LL[1] <- out$l
  LLp[1] <- out$lp
  LLpp[1] <- out$lpp
  
  convergence = FALSE # convergence is set to false in order to run while loop
  
  start <- 1 # starting value
  
  # as long as convergence is false, keep running the loop
  while(!convergence){
    
    thetaVec[start + 1] <- thetaVec[start] + alpha*out$lp
    
    out <- llbinom.grad(thetaVec[start + 1], n = n, y = y)
    
    LL[start + 1] <- out$l
    LLp[start + 1] <- out$lp
    LLpp[start + 1] <- out$lpp
    
    # Stop the while loop as soon as the absolute value between theta and 
    # previous theta is below convergence level
    if ((abs(thetaVec[start + 1] - thetaVec[start])) < delta){
      convergence = TRUE
    }
    
    # increase starting value
    start <- start + 1
    
  }
  
  # make dataframe of the iteration history
  # iterhistBinom <- cbind(thetaVec, LL, LLp, LLpp)
  iterhistBinom <- data.frame(Iterations = seq(from = 1, to = start), theta = thetaVec,
                              LogLikelihood = LL, FirstDerivative = LLp, SecondDerivative = LLpp)
  
  return(iterhistBinom)
}


################################### Choices and Buttons ##########################


distribution <- selectInput(inputId = "dist", label = "Choose your distribution",
                            choices = c("Binomial", "Poisson"))
thetaVal <- numericInput(inputId = "theta", label = "Choose your theta value", min = 0,
                         max = 0.99, step = 0.01, value = 0.5)
n <- numericInput(inputId = "n", label = "Choose the number of Hits", min = 1,
                  max = 20, step = 1, value = 3)
y <-numericInput(inputId = "y", label = "Choose the number of Misses", min = 1,
                 max = 20, step = 1, value = 1)
learningRate <- numericInput(inputId = "learning", label = "Choose your learning rate",
                             min = 0.005, max = 0.05, step = 0.001, value = 0.005)
convergence <- numericInput(inputId = "conv", label = "Set the Convergence Level", value = 10^-16, min = 10^-20, 
                            max = 10^-5, step= 10^-20)
#changeIter <- actionButton(inputId = "change", label ="Change iteration max value")


####################### User Interface ###################################
ui <- fluidPage(
  titlePanel("Maximum Likelihood"),
  
  sidebarLayout(
    sidebarPanel(
      distribution,
      thetaVal,
      y,
      n,
      learningRate,
      convergence
      #uiOutput("iterSlider")
      #changeIter
      
    ),
    
    mainPanel(
      plotlyOutput(outputId = "plot"),
      verbatimTextOutput(outputId = "selected"),
      DTOutput(outputId = "iterhistBinom")
    )
  
  )
  
)

###################### Server #################################
server <- function(input, output, session){
    
  # Plots
    output$plot <- renderPlotly({

      thetas <- seq(.01, .99, by=.001)
      llik <- llbinom.grad(thetas, n = input$n, y = input$y)$l
      d <- data.frame(thetas = thetas, llik = llik)
      
      df <- MaxBinom(theta = input$theta, alpha = input$learning, 
                     delta = input$conv, n = input$n, y = input$y)
      s <- input$iterhistBinom_rows_selected
      
      newdf <- data.frame(cbind(df[s, 2],df[s, 3]))
      colnames(newdf)[1] <- "theta"
      colnames(newdf)[2] <- "Loglikelihood"
      
      if (length(s)){
        
        ggplot(d, aes(x = thetas, y = llik)) + 
          geom_line() + 
          geom_point(aes(x = input$theta,
                         y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
                     colour = "blue") +
          geom_point(data = newdf, mapping = aes(x = theta, y = Loglikelihood, colour = "red"))
        
      } else {
        
        ggplot(d, aes(x = thetas, y = llik)) + 
          geom_line() + geom_point(aes(x = input$theta,
                                       y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
                                   colour = "blue")
      }
      #qplot(x = thetas, y = llbinom.grad(thetas, n = input$n, y = input$y)$l, geom = 'line',
      #      xlab = expression(theta), ylab = "log Likelihood") + 
        #scale_x_continuous(limits = c(0, .99)) + 
        #geom_point(aes(x = input$theta, y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
        #           colour = "red")

    })
      
    # Table of iterations
    output$iterhistBinom <- DT::renderDataTable(
      MaxBinom(theta = input$theta, alpha = input$learning, 
               delta = input$conv, n = input$n, y = input$y)
    )
    
    output$selected <- renderPrint({
      s <- input$iterhistBinom_rows_selected
      if(length(s)){
        cat('These iterations were selected:\n\n')
        cat(s, sep = ', ')      }
    })
    
  
}


shinyApp(ui, server)



  


