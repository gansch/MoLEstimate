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

# Function for Likelihood.. of Poisson distribution
llpoisson.grad <- function(lambda, n, y){
  pell <- -n*(lambda) + sum(y)*log(lambda) # loglikelihood
  pellp <- -n + sum(y)/lambda # 1st derivative of llh
  pellpp <- -1*(sum(y)/lambda^2) # 2nd derivative of llh
  out<-list(l=pell,lp=pellp, lpp=pellpp)
  return(out)
}

###Cauchy distribution
llcauchy.grad <- function(theta, y){
  cell <-NULL
  cellp <- NULL
  cellpp <- NULL
  for (i in seq_along(theta)){
    #loglikelihood
    cell[i] <- sum(-log(pi)-log((y-theta[i])^2+1))
    #first derivative
    cellp[i] <- sum(2*(y-theta[i])/(1+(y-theta[i])^2))
    #seconds derivative
    cellpp[i] <- 2*sum((((y-theta[i])^2-1)/(((y-theta[i])^2)+1)^2))
  }
  out<-list(l=cell,lp=cellp,lpp=cellpp)
  
  return(out)
}

######## Function that searches for Maximum Likelihood: Binom
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

###### Function that searches for Maximum Likelihood: Iterative Method: POISSON
MaxPois <- function(lambda, alpha, delta = 10^-6, n, y){
  if (20*alpha > lambda) {
    stop("You must enter a bigger value for the learning rate or a smaller lambda")
  }
  # define empty vectors
  lambdaVec <- NULL
  LL <- NULL
  LLp <- NULL
  LLpp <- NULL
  
  # assign lambda
  lambdaVec[1] <- lambda
  
  # calculate likelihoods
  out <- llpoisson.grad(lambda = lambdaVec[1], n = n, y = y)
  
  # enter result of the list into the vectors
  LL[1] <- out$l
  LLp[1] <- out$lp
  LLpp[1] <- out$lpp
  
  convergence = FALSE # convergence is set to false in order to run while loop
  
  start <- 1 # starting value
  
  # as long as convergence is false, keep running the loop
  while(!convergence){
    
    lambdaVec[start + 1] <- lambdaVec[start] + alpha*out$lp
    
    out <- llpoisson.grad(lambdaVec[start + 1], n = n, y = y)
    
    LL[start + 1] <- out$l
    LLp[start + 1] <- out$lp
    LLpp[start + 1] <- out$lpp
    
    # Stop the while loop as soon as the absolute value between lambda and 
    # previous lambda is below convergence level
    if ((abs(lambdaVec[start + 1] - lambdaVec[start])) < delta){
      convergence = TRUE
    }
    
    # increase starting value
    start <- start + 1
    
  }
  
  # make dataframe of the iteration history
  iterhistPoisson <- data.frame(Iterations = seq(from = 1, to = start), lambda = lambdaVec,
                                LogLikelihood = LL, FirstDerivative = LLp, SecondDerivative = LLpp)
  
  return(iterhistPoisson)
}

#####Cauchy iteration 
MaxCauchy <- function(lambda = 1, alpha = .005, delta = 10^-6, y ){
  
  # if (10*alpha > lambda) {
  #   stop("You must enter a bigger value for the learning rate or a smaller theta")
  # }
  # define empty vectors
  thetaVec <- NULL
  LL <- NULL
  LLp <- NULL
  LLpp <- NULL
  
  # assign theta
  thetaVec[1] <- lambda
  
  # calculate likelihoods
  out <- llcauchy.grad(theta = thetaVec[1],  y = y)
  
  # enter result of the list into the vectors
  LL[1] <- out$l
  LLp[1] <- out$lp
  LLpp[1] <- out$lpp
  
  convergence = FALSE # convergence is set to false in order to run while loop
  
  start <- 1 # starting value
  
  # as long as convergence is false, keep running the loop
  while(!convergence){
    
    thetaVec[start + 1] <- thetaVec[start] + alpha*out$lp
    
    out <- llcauchy.grad(thetaVec[start + 1], y = y)
    
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
  #iterhistCauchy <- cbind(thetaVec, LL, LLp, LLpp)
  #iterhistCauchy <- data.frame(Iterations = seq(from = 1, to = start), theta = thetaVec,
  #                           LogLikelihood = LL, FirstDerivative = LLp, SecondDerivative = LLpp)
  iterhistCauchy <- cbind(thetaVec,LL,LLp,LLpp)
  #iterhistCauchy <- rbind(iterhistCauchy,interhistCauchy)
  return(iterhistCauchy)
}


################################### Choices and Buttons ##########################

########## Bionominal

distribution <- selectInput(inputId = "dist", label = "Choose your distribution",
                            choices = c("Binomial", "Poisson", "Cauchy"), selected = "Binomial")
thetaVal <- numericInput(inputId = "theta", label = "Choose your theta ($\\theta$)", min = 0,
                         max = 0.99, step = 0.01, value = 0.5)
n <- numericInput(inputId = "n", label = "Choose the number of Misses", min = 1,
                  max = 20, step = 1, value = 3)
y <-numericInput(inputId = "y", label = "Choose the number of Hits", min = 1,
                 max = 20, step = 1, value = 1)
learningRate <- numericInput(inputId = "learning", label = "Choose your learning rate",
                             min = 0.0005, max = 0.2, step = 0.001, value = 0.015)
convergence <- sliderInput(inputId = "conv", label = "Set the exponent of the convergence level (1*10^-x)", value = -6, min = -15, 
                                                        max = -2, step= 1)
deriv <- checkboxInput("deriv", label="Check box to display curve of the 1st derivative",
                       value=FALSE)

########## Poisson

lambdaVal <- numericInput(inputId = "lambda", label = "Choose your lambda ($\\lambda$)", min = 0,
                          max = 0.99, step = 0.01, value = 0.5)
observationVal <-numericInput(inputId = "observationP", label = "Choose the sum of observations ($k$)", min = 1,
                              max = 20, step = 1, value = 1)
timeVal <- numericInput(inputId = "timeP", label = "Choose numbers of time interval", min = 1,
                        max = 20, step = 1, value = 3)
convergencePois <- sliderInput(inputId = "convP", label = "Set the exponent of the convergence level (1*10^-x)", value = -6,
                               			       	min = -15, max = -2, step= 1)
learningRatePois <- numericInput(inputId = "learningP", label = "Choose your learning rate",
                                 min = 0.005, max = 0.5, step = 0.001, value = 0.015)
derivPois <- checkboxInput("derivP", label="Check box to display curve of the 1st derivative",
                           value=FALSE)

#### cauchy

lambdaValC <- numericInput(inputId = "lambdac", label = "Choose your lambda ($\\lambda$)", min = -6,
                          max = 6, step = 0.5, value = 1)

convergenceC <- sliderInput(inputId = "convc", label = "Set the exponent of the convergence level (1*10^-x)", value = -6, min = -15, 
                            max = -2, step= 1)
derivC <- checkboxInput("derivc", label="Check box to add curve of the 1st derivative",
                       value=FALSE)

learningRateC <- numericInput(inputId = "learningc", label = "Choose your learning rate",
                             min = 0.0005, max = 0.2, step = 0.001, value = 0.015)

#r <- numericInput(inputId = "r", "Pick a number of points and press button to generate.",
 #                 min=0,step=1,value=10)

#button <- actionButton("generateBt", "Generate Numbers")


####################### User Interface ###################################
ui <- fluidPage(
  titlePanel("Maximum Likelihood"),
  
  sidebarLayout(
    sidebarPanel(
      distribution,
      
      # Conditional Panel for Binomial Distribution
      conditionalPanel(
        condition = ("input.dist == 'Binomial'"),
        withMathJax(helpText("$x$ ~ $B(n, \\theta)$")),
        y,
        n,
        withMathJax(thetaVal),
        learningRate,
        convergence,
        deriv
      ),
      
      # Conditional Panel for Poisson Distribution
      conditionalPanel(
        condition = ("input.dist == 'Poisson'"),
        helpText("$x$ ~ $Po(\\lambda)$"),
        observationVal,
        timeVal,
        withMathJax(lambdaVal),
        learningRatePois,
        convergencePois,
        derivPois
      ),
      
      # conditional panel for Cauchy
      conditionalPanel(
        condition = ("input.dist == 'Cauchy'"),
        lambdaValC,
        learningRateC,
        convergenceC,
        derivC
       #r
       #button
      )
    ),
    
    mainPanel(
      plotlyOutput(outputId = "plot"),
      verbatimTextOutput(outputId = "selected"),
      DTOutput(outputId = "iterhist")
    )
    
  )
  
)

###################### Server #################################
server <- function(input, output, session){
  
  # Plots
  output$plot <- renderPlotly({
    
    # BINOMIAL PLOT    
    if (input$dist == "Binomial"){
      thetas <- seq(.01, .99, by=.001)
      llik_list <- llbinom.grad(thetas, n = input$n, y = input$y)[1:2]
      
      d <- data.frame(thetas = thetas, llik = llik_list[[1]],
                      llikp = llik_list[[2]])
      
      df <- MaxBinom(theta = input$theta, alpha = input$learning, 
                     delta = (10^(input$conv)), n = input$n, y = input$y)
      s <- input$iterhist_rows_selected
      
      newdf <- data.frame(cbind(df[s, 2], df[s, 3], df[s, 4]))
      colnames(newdf)[1] <- "theta"
      colnames(newdf)[2] <- "Loglikelihood"
      colnames(newdf)[3] <- "FirstDeriv"
      # no 1st derivative added
      if (!input$deriv){
        if (length(s)){
          
          ggplot(d, aes(x = thetas)) + 
            geom_line(aes(y = llik)) +
            geom_point(aes(x = input$theta,
                           y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
                       colour = "blue") +
           # geom_point(data = newdf[nrow(newdf),], 
            #             mapping = aes(x = theta, y = Loglikelihood, colour = "blue"), alpha=1) + 
            geom_point(data = newdf, 
                       mapping = aes(x = theta, y = Loglikelihood, colour = "red"), alpha=.3) + 
            geom_segment(data = newdf,
                         mapping = aes(x = theta - .1, xend = theta + .1,
                                       y = Loglikelihood - .1*FirstDeriv, yend = Loglikelihood + .1*FirstDeriv),
                         colour = "red") + 
            geom_segment(data = newdf,
                         mapping = aes(x = theta, xend = theta, y = min(llik_list[[1]]), yend = Loglikelihood),
                         colour = "green", alpha = .5, linetype = "dashed")
          
        } else {
          
          ggplot(d, aes(thetas)) + 
            geom_line(aes(y = llik)) + 
            geom_point(aes(x = input$theta,
                           y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
                       colour = "blue")
        }
        # 1st derivative added
      } else {
        if (length(s)){
          
          ggplot(d, aes(x = thetas)) + 
            geom_line(aes(y = llik)) +
            geom_line(aes(y = llikp), alpha = .5, linetype = "dotted") +
            geom_point(aes(x = input$theta,
                           y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
                       colour = "blue") +
            geom_point(data = newdf, 
                       mapping = aes(x = theta, y = Loglikelihood, colour = "red"), alpha=1) + 
            geom_segment(data = newdf,
                         mapping = aes(x = theta - .1, xend = theta + .1,
                                       y = Loglikelihood - .1*FirstDeriv, yend = Loglikelihood + .1*FirstDeriv),
                         colour = "red") + 
            geom_segment(data = newdf,
                         mapping = aes(x = theta, xend = theta, y = min(llik_list[[2]]), yend = Loglikelihood),
                         colour = "green", alpha = .5, linetype = "dashed")
          
        } else {
          
          ggplot(d, aes(thetas)) + 
            geom_line(aes(y = llik)) + 
            geom_line(aes(y = llikp), alpha = .6, linetype = "dotted") + 
            geom_point(aes(x = input$theta,
                           y = llbinom.grad(input$theta, n = input$n, y = input$y)$l),
                       colour = "blue")
        }
      }
      # POISSON PLOT
    } else if (input$dist == "Poisson") {
      lambdas <- seq(.01, input$observationP, by=.001)
      llik_list <- llpoisson.grad(lambdas, n = input$timeP, y = input$observationP)[1:2]
      
      d <- data.frame(lambdas = lambdas, llik = llik_list[[1]],
                      llikp = llik_list[[2]])
      
      df <- MaxPois(lambda = input$lambda, alpha = input$learningP, 
                    delta = (10^(input$convP)), n = input$timeP, y = input$observationP)
      s <- input$iterhist_rows_selected
      
      newdf <- data.frame(cbind(df[s, 2], df[s, 3], df[s, 4]))
      colnames(newdf)[1] <- "lambda"
      colnames(newdf)[2] <- "Loglikelihood"
      colnames(newdf)[3] <- "FirstDeriv"
      # no 1st derivative added
      if (!input$derivP){
        if (length(s)){
          
          ggplot(d, aes(x = lambdas)) + 
            geom_line(aes(y = llik)) +
            geom_point(aes(x = input$lambda,
                           y = llpoisson.grad(input$lambda, n = input$timeP, y = input$observationP)$l),
                       colour = "blue") +
            geom_point(data = newdf, 
                       mapping = aes(x = lambda, y = Loglikelihood, colour = "red"), alpha=1) + 
            geom_segment(data = newdf,
                         mapping = aes(x = lambda - (input$observationP/10), xend = lambda + (input$observationP/10),
                                       y = Loglikelihood - (input$observationP/10)*FirstDeriv,
                                       yend = Loglikelihood + (input$observationP/10)*FirstDeriv),
                         colour = "red") + 
            geom_segment(data = newdf,
                         mapping = aes(x = lambda, xend = lambda, y = min(llik_list[[1]]), yend = Loglikelihood),
                         colour = "green", alpha = .5, linetype = "dashed")
          
        } else {
          
          ggplot(d, aes(lambdas)) + 
            geom_line(aes(y = llik)) + 
            geom_point(aes(x = input$lambda,
                           y = llpoisson.grad(input$lambda, n = input$timeP, y = input$observationP)$l),
                       colour = "blue") 
        }
        # 1st derivative added
      } else {
        if (length(s)){
          
          ggplot(d, aes(x = lambdas)) + 
            geom_line(aes(y = llik)) +
            geom_line(aes(y = llikp), alpha = .5, linetype = "dotted") +
            geom_point(aes(x = input$lambda,
                           y = llbinom.grad(input$lambda, n = input$timeP, y = input$observationP)$l),
                       colour = "blue") +
            geom_point(data = newdf, 
                       mapping = aes(x = lambda, y = Loglikelihood, colour = "red"), alpha=1) + 
            geom_segment(data = newdf,
                         mapping = aes(x = lambda - .1, xend = lambda + (input$observationP/10),
                                       y = Loglikelihood - (input$observationP/10)*FirstDeriv,
                                       yend = Loglikelihood + (input$observationP/10)*FirstDeriv),
                         colour = "red") + 
            geom_segment(data = newdf,
                         mapping = aes(x = lambda, xend = lambda, y = min(llik_list[[1]]), yend = Loglikelihood),
                         colour = "green", alpha = .5, linetype = "dashed")
          
        } else {
          
          ggplot(d, aes(lambdas)) + 
            geom_line(aes(y = llik)) + 
            geom_line(aes(y = llikp), alpha = .6, linetype = "dotted") + 
            geom_point(aes(x = input$lambda,
                           y = llpoisson.grad(input$lambda, n = input$timeP, y = input$observationP)$l),
                       colour = "blue")
        }
      }
      
####Cauchy plot
    } else {
        thetasc <- seq(-6, 6, by=.1)
        
        ##random number generator
        ran <- c(-0.71,-0.53,3.27,5.64)
        #ran <- rcauchy(input$r)
        llik_list <- llcauchy.grad(thetasc, y = ran)[1:2]
        
        d <- data.frame(thetasc = thetasc, llik = llik_list[[1]],
                        llikp = llik_list[[2]])
        
        df <- MaxCauchy(lambda = input$lambdac, alpha = input$learningc, 
                       delta = (10^(input$convc)), y = ran)
        s <- input$iterhist_rows_selected
        
        newdf <- data.frame(cbind(df[s, 1], df[s, 2], df[s, 3]))
        colnames(newdf)[1] <- "thetac"
        colnames(newdf)[2] <- "Loglikelihood"
        colnames(newdf)[3] <- "FirstDeriv"
        # no 1st derivative added
        if (!input$derivc){
          if (length(s)){
            
            ggplot(d, aes(x = thetasc)) + 
              geom_line(aes(y = llik)) +
              geom_point(aes(x = input$lambdac,
                             y = llcauchy.grad(input$lambdac, y = ran)$l),
                         colour = "blue") +
              geom_point(data = newdf, 
                         mapping = aes(x = thetac, y = Loglikelihood, colour = "red"), alpha=1) + 
              geom_segment(data = newdf,
                           mapping = aes(x = thetac - .1, xend = thetac + .1,
                                         y = Loglikelihood - .1*FirstDeriv, yend = Loglikelihood + .1*FirstDeriv),
                           colour = "red") + 
              geom_segment(data = newdf,
                           mapping = aes(x = thetac, xend = thetac, y = min(llik_list[[1]]), yend = Loglikelihood),
                           colour = "green", alpha = .5, linetype = "dashed")
            
          } else {
            
            ggplot(d, aes(thetasc)) + 
              geom_line(aes(y = llik)) + 
              geom_point(aes(x = input$lambdac,
                             y = llcauchy.grad(input$lambdac, y = ran)$l),
                         colour = "blue")
          }
          # 1st derivative added
        } else {
          if (length(s)){
            
            ggplot(d, aes(x = thetasc)) + 
              geom_line(aes(y = llik)) +
              geom_line(aes(y = llikp), alpha = .5, linetype = "dotted") +
              geom_point(aes(x = input$lambdac,
                             y = llcauchy.grad(input$lambdac, y = ran)$l),
                         colour = "blue") +
              geom_point(data = newdf, 
                         mapping = aes(x = thetac, y = Loglikelihood, colour = "red"), alpha=1) + 
              geom_segment(data = newdf,
                           mapping = aes(x = thetac - .1, xend = thetac + .1,
                                         y = Loglikelihood - .1*FirstDeriv, yend = Loglikelihood + .1*FirstDeriv),
                           colour = "red") + 
              geom_segment(data = newdf,
                           mapping = aes(x = thetac, xend = thetac, y = min(llik_list[[1]]), yend = Loglikelihood),
                           colour = "green", alpha = .5, linetype = "dashed")
            
          } else {
            
            ggplot(d, aes(thetasc)) + 
              geom_line(aes(y = llik)) + 
              geom_line(aes(y = llikp), alpha = .6, linetype = "dotted") + 
              geom_point(aes(x = input$lambdac,
                             y = llcauchy.grad(input$lambdac, y = ran)$l),
                         colour = "blue")
          }
        }
    }
  })
  
  # Table of iterations
  output$iterhist <- DT::renderDataTable({
    if (input$dist == "Binomial"){ # if binomial selected
      MaxBinom(theta = input$theta, alpha = input$learning, 
               delta = (10^(input$conv)), n = input$n, y = input$y)
    } else if (input$dist == "Poisson"){ # or poisson
      MaxPois(lambda = input$lambda, alpha = input$learningP,
              delta = (10^(input$convP)), n = input$timeP, y = input$observationP)
    } else {
      MaxCauchy(lambda = input$lambdac, alpha = input$learningc,
                delta = (10^(input$convc)), y = ran)
    }
  })
  
  output$selected <- renderPrint({
    s <- input$iterhist_rows_selected
    if(length(s)){
      cat('These iterations were selected:\n\n')
      cat(s, sep = ', ')      }
  })
  
  
}



shinyApp(ui, server)



