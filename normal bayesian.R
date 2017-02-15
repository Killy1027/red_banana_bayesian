loglikelihood <- function(y,mu,sigma)
{
  sum(dnorm(y,mu,sigma,log=T))
}

MCMC <- function(y,iters=1000)
{
  mu_current <- mean(y)
  sigma_current <- sd(y)
  currentcandidate <- loglikelihood(y,mu_current,sigma_current)
  keep.mu <- keep.sigma <- NULL
  
  for( i in 1:iters)
  {
    mu_candidate <- rnorm(1,mu_current,10)
    candidatelikelihood <- loglikelihood(y,mu_candidate,sigma_current)
    ratio_mu <- candidatelikelihood + dnorm(mu_candidate,0,1000,log=T) - currentcandidate - dnorm(mu_current,0,1000,log=T)
    if(log(runif(1)) < ratio_mu )
    {
      currentcandidate <- candidatelikelihood
      mu_current <- mu_candidate
    }
    keep.mu[i] <- mu_current
    
    sigma_candidate <- exp(rnorm(1,log(sigma_current),1))
    candidatelikelihood <- loglikelihood(y,mu_current,sigma_candidate)
    ratio_sigma <- candidatelikelihood + dgamma(sigma_candidate,0.5,0.001,log=T) - currentcandidate - dgamma(sigma_current,0.5,0.001,log=T) 
    if(log(runif(1)) < ratio_sigma ) 
    {
      currentcandidate <- candidatelikelihood
      sigma_current <- sigma_candidate
    }
    keep.sigma[i] <- sigma_current
  }
  return(list(mu=keep.mu,sigma=keep.sigma))
}

# my range of mu
x <- seq(-5000,5000,1) 
plot(x, dnorm(x,0,1000),type="l")
# my range of sigma
x <- seq(0,10000,1)
plot(x,dgamma(x,3,0.001),type="l")

# my test
y <- 10
a <- MCMC(10,10000)

plot(a$mu,type="l")
summary(a$mu)
hist(a$mu)
plot(a$sigma,type="l")
summary(a$sigma)
hist(a$sigma)

x <- seq(-1000,1000,1) 
plot(x, dnorm(x,median(a$mu),median(a$sigma)),type="l")
abline(v=median(a$mu),col=2)

dunif(runif(1))
rnorm(1,0.1,10)
dnorm(10,3,2)
