library(shiny)

###Bino Demo Data Creation

llbinom.grad <- function(theta){
  ell <- 7*log(theta)+3*log(1-theta)
  ellp <- 7/theta - 3/(1-theta)
  ellpp <- -7/theta^2 - 3/(1-theta)^2
  out<-list(l=ell,lp=ellp,lpp=ellpp)
  return(out)
}

alpha=.015; delta=10^-6
th.vec=NULL; LL=NULL; LLp=NULL; LLpp=NULL
th.vec[1]=0.3
out<-llbinom.grad(th.vec[1])
LL[1]=out$l; LLp[1]=out$lp; LLpp[1]=out$lpp
convergence=FALSE
a=1
while (!convergence){
  th.vec[a+1] <- th.vec[a] + alpha*out$lp
  out<-llbinom.grad(th.vec[a+1])
  LL[a+1]=out$l; LLp[a+1]=out$lp; LLpp[a+1]=out$lpp
  if ((abs(th.vec[a+1]-th.vec[a]))<delta) {convergence=TRUE}
  a <- a+1
}
iterhist = data.frame(t=seq(1,a),theta=th.vec,ell=LL,ellp=LLp,ellpp=LLpp)
iterhist




# UI
ui <- fluidPage(
  sidebarLayout(
    
    # Input
    sidebarPanel(
      
      # Select variable for y-axis
      selectInput(inputId = "dist", 
                  label = "Distribution:",
                  choices = c("Binomial", "Poison"), 
                  selected = "Poison"),
      
      
      # Set alpha level
      sliderInput(inputId = "iter", 
                  label = "Iteration:", 
                  min = 0, max = nrow(iterhist), step= 1, 
                  value = c(1,2)
      )
    ),
    
    # Output:
    mainPanel(
      plotOutput(outputId = "plot")
      
    )
  )
)

# Define server function required to create the scatterplot-
server <- function(input, output, session) {
  
  # Create scatterplot object the plotOutput function is expecting 
  output$plot <- renderPlot({
    ### Bino Plot
    
    binocurve <- 
    poiscurve <- 
    
    plot(curve(7*log(x)+3*log(1-x)), type = 'l', lwd = 3)
    
    if (input$iter[2] != 0) {
      segments(x0= iterhist$theta[input$iter[1]:input$iter[2]]-.3, 
             y0 = iterhist$ell[input$iter[1]:input$iter[2]]-.3*iterhist$ellp[input$iter[1]:input$iter[2]],  
             x1 = iterhist$theta[input$iter[1]:input$iter[2]]+.3, 
             y1 =iterhist$ell[input$iter[1]:input$iter[2]]+.3*iterhist$ellp[input$iter[1]:input$iter[2]],
             col="red")
      segments(x0= iterhist$theta[input$iter[1]:input$iter[2]], y0=-100, y1 = iterhist$ell[input$iter[1]:input$iter[2]],
               col="green", 
               lty=2)
      } else {}
      })
      
  
  
}

# Create a Shiny app object
shinyApp(ui = ui, server = server)

###