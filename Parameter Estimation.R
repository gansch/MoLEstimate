
theta <- seq(0,6,.001)
n <- 3
y <- c(2, 2, 3)
Ltheta= exp(-n*theta)*theta^(sum(y))
plot(theta,Ltheta,type="l",xlab=expression(theta),ylab=expression(L(theta)))

Poislik <- function(n, theta, y){
  # Likelihood
  Ltheta <- exp(-n*theta)*theta^(sum(y))
  par(mfrow = c(2, 1))
  # Loglikelihood
  plot(theta,Ltheta,type="l",xlab=expression(theta),ylab=expression(L(theta)))
  Ltheta1 <- log(exp(-n*theta)*theta^(sum(y)))
  #Ltheta1 <- n*theta + sum(y) * log(theta)
  #Ltheta1 <- log(exp(-n*theta)) + log(theta^(sum(y)))
  # 
  plot(theta, Ltheta1,type="l",xlab=expression(theta),ylab=expression(L(theta)), xlim = c(0, 6),
       ylim = c(-25, 10))
}

Poislik(3, seq(0,6,.001), c(2, 2, 3))



# Taylor's method: Binomial

llbinom.grad <- function(theta){
  ell <- 7*log(theta)+3*log(1-theta)
  ellp <- 7/theta - 3/(1-theta)
  ellpp <- -7/theta^2 - 3/(1-theta)^2
  out<-list(l=ell,lp=ellp,lpp=ellpp)
  return(out)
}
alpha=.015; delta=10^-6
th.vec=NULL; LL=NULL; LLp=NULL
th.vec[1]=0.3
out<-llbinom.grad(th.vec[1])
LL[1]=out$l; LLp[1]=out$lp
convergence=FALSE
a=1
while (!convergence){
  th.vec[a+1] <- th.vec[a] + alpha*out$lp
  out<-llbinom.grad(th.vec[a+1])
  LL[a+1]=out$l; LLp[a+1]=out$lp
  if ((abs(th.vec[a+1]-th.vec[a]))<delta) {convergence=TRUE}
  a <- a+1
}
iterhist = data.frame(t=seq(1,a),theta=th.vec,ell=LL,ellp=LLp)
iterhist


# Taylor's method: Poisson 

llpoisson.grad <- function(theta){
  ell <- -3*theta + 7*log(theta)
  ellp <- 3 + 7/theta
  ellpp <- -7*theta^(-2)
  out <- list(l = ell, lp = ellp, lpp = ellpp)
  return(out)
}


alpha=.015; delta=10^-6
th.vec = NULL
LL=NULL
LLp=NULL
th.vec[1] <- 1

out <- llpoisson.grad(th.vec[1])

LL[1] <- out$l
LLp[1] <- out$lp
convergence=FALSE
a=1

while (!convergence){
  th.vec[a+1] <- th.vec[a] + alpha*out$lp
  out <- llpoisson.grad(th.vec[a+1])
  LL[a+1] <- out$l
  LLp[a+1] <- out$lp
  if ((abs(th.vec[a+1]-th.vec[a])) < delta){
    convergence=TRUE
    }
  a <- a+1
}

iterhist = data.frame(t=seq(1,a),theta=th.vec,ell=LL,ellp=LLp)
iterhist
















