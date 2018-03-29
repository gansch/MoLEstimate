library(shiny)
library(ggplot2)


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
MaxBinom <- function(theta, alpha, delta = 10^-6, start = 1, n, y){
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


Poislik <- function(n, lam, y){
  # Likelihood
  Lik <- exp(-n*lam)*lam^(sum(y))
  
  # Loglikelihood
  #loglik <- log(exp(-n*lam)*theta^(sum(y)))
  loglik <- -n*lam + sum(y) * log(lam)
  #loglik <- log(exp(-n*lam)) + log(lam^(sum(y)))
  
  #First derivative of loglik 
  loglik1 <- n + sum(y)/lam
  
  #Second derivative lof loglik
  loglik2 <- -sum(y)*lam^(-2)
  
  #Plot
  par(mfrow = c(2, 2))
  plot(lam, Lik,type="l",xlab=expression(lam),ylab=expression(L(lam)))
  plot(lam, loglik,type="l",xlab=expression(lam),ylab=expression(L(lam)), xlim = c(0, 6),
       ylim = c(-25, 10))
  plot(lam, loglik1,type="l",xlab=expression(lam),ylab=expression(L(lam)))
  plot(lam, loglik2,type="l",xlab=expression(lam),ylab=expression(L(lam)), xlim = c(0, 6),
       ylim = c(-25, 10))
  par(mfrow = c(1, 1))
}

#Poislik(3, seq(0,6,.001), c(2, 2, 3))

################################### Choices and Buttons ##########################

distribution <- selectInput(inputId = "dist", label = "Choose your distribution",
                            choices = c("Binomial", "Poisson"))
thetaVal <- numericInput(inputId = "theta", label = "Choose your theta value", min = 0,
                         max = 0.99, step = 0.01, value = 0.5)
n <- numericInput(inputId = "n", label = "Choose your n value", min = 1,
                  max = 20, step = 1, value = 3)
y <-numericInput(inputId = "y", label = "Choose your y value", min = 1,
                 max = 20, step = 1, value = 1)
startingVal <- numericInput(inputId = "starting", label = "Choose your starting value",
                           min = 0.1, max = 1, step = 0.05, value = 0.5)
learningRate <- numericInput(inputId = "learning", label = "Choose your learning rate",
                            min = 0.005, max = 0.05, step = 0.001, value = 0.005)
iterations <- sliderInput(inputId = "iter", label = "Iteration:", value = 1, min = 1, 
                          max = 2, step= 1)
changeIter <- actionButton(inputId = "change", label ="Change iteration max value")


####################### User Interface ###################################
ui <- fluidPage(
  titlePanel("Maximum Likelihood"),
  
  sidebarLayout(
    sidebarPanel(
      distribution,
      thetaVal,
      y,
      n,
      startingVal,
      learningRate,
      iterations,
      changeIter
      
    ),
    
    mainPanel(
      plotOutput(outputId = "plot"),
      dataTableOutput(outputId = "iterhistBinom")
      #verbatimTextOutput(outputId = "iterhistBinom")
    )
  
  )
  
)

###################### Server #################################
server <- function(input, output, session){
    
  # Table of iterations
    output$iterhistBinom <- renderDataTable(
      MaxBinom(theta = input$theta, alpha = input$learning, 
               delta = 10^-16, start = input$starting, n = input$n, y = input$y)
      
    )

    observeEvent(input$change, {
      updateNumericInput(session, "iterations", value = 1, min = 1, max = nrow(output$iterhistBinom),
                         step = 1)
    })

  # Plots
    output$plot <- renderPlot({

      ## Binom Plot
      plot(curve(llbinom.grad(x, n = input$n, y = input$y)$l, from = 0, to = 2*input$theta),
           type = 'l', lwd = 2, xlab = "theta",
           ylab = "log Likelihood", bty = "l")
      points(input$theta, llbinom.grad(input$theta, n = input$n, y = input$y)$l,
             col = "red", pch = 16)

      # # Adding Tangent Line and vertical lines
      # if (input$iter[2] != 0) {
      #   segments(x0= iterhist$theta[input$iter[1]:input$iter[2]]-.3,
      #            y0 = iterhist$ell[input$iter[1]:input$iter[2]]-.3*iterhist$ellp[input$iter[1]:input$iter[2]],
      #            x1 = iterhist$theta[input$iter[1]:input$iter[2]]+.3,
      #            y1 =iterhist$ell[input$iter[1]:input$iter[2]]+.3*iterhist$ellp[input$iter[1]:input$iter[2]],
      #            col="red")
      #   segments(x0= iterhist$theta[input$iter[1]:input$iter[2]], y0=-100, y1 = iterhist$ell[input$iter[1]:input$iter[2]],
      #            col="green",
      #            lty=2)
      # }

    })

    
    

    

  
}


shinyApp(ui, server)






