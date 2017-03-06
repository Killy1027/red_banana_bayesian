


#======================================
#======================================
# log likelihood function
#======================================
#======================================






llike <- function(x, y, t, a, b, p){
  
  mu1 <- rep(z[sum(x[1] < new.x), sum(y[1] < new.y)], NROW(x))
  lambda[1] <- (1-p)*mu1/max(t)
  for(i in 2:length(t))
  {
    mu <- z[sum(x[i]<new.x),sum(y[i]<new.y)]
    #====================================
    # There is sth. wrong with belowing code!!!!   each element results -Inf.
    #====================================
    lambda[i] <- (1-p)*mu[i] p*sum((a*b/pi)*exp(-a*iet[1:(i-1),i] - b*ied2[1:(i-1),i])
  }
  sum(log(lambda))-(1-p)*max(t)*density.sum-p*length(t)
}


llike <- function(x, y, t, a, b, p){
  
}


#======================================
#======================================
# MCMC Sampler
#======================================
#======================================
RedBananaMCMC <- function(x, y, t,
                          iters = 1000, burn = iters/2, thin = 1,
                          prior.sd.a = 100, prior.sd.b = 100){
  
  
  #======================================
  # Initial values
  #======================================
  a <- rnorm(1, 0, prior.sd.a)
  b <- rnorm(1, 0, prior.sd.b)
  p <- runif(1)
  ll <- llike(x, y, t, a, b, p)
  
  #======================================
  # Keepers
  #======================================
  keep.pars <- matrix(0, iters, 3)
  colnames(keep.pars) <- c("a", "b", "p")
  keep.ll <- rep(0, iters)
  
  #======================================
  # START MCMC!!!
  #======================================
  for(i in 1:iters){
    for(j in 1:thin){
      
      
      #__________________________________
      # Update a
      #
      cana <- rnorm(1, a, 1)
      canll <- llike(x, y, t, cana, b, p)
      R <- sum(canll - ll) + dnorm(cana, 0, prior.sd.a) - dnorm(a, 0, prior.sd.a)
      if(log(runif(1)) < R){
        a <- cana
        ll <-  canll
      }
      
      #__________________________________
      # Update b
      #
      canb <- rnorm(1, b, 1)
      canll <- llike(x, y, t, a, canb, p)
      R <- sum(canll - ll) + dnorm(canb, 0, prior.sd.b) - dnorm(a, 0, prior.sd.b)
      if(log(runif(1)) < R){
        b <- canb
        ll <-  canll
      }
      
      #__________________________________
      # Update p
      #
      canp <- rbeta(1, p, 1-p)
      canll <- llike(x, y, t, a, b, canp)
      R <- sum(canll - ll) + dbeta(canp, p, 1-p) - dbeta(p, canp, 1-canp)
      if(log(runif(1)) < R){
        p <- canp
        ll <-  canll
      }
    

    }##end thin
    
    
    #======================================
    # Keepers
    #======================================
    keep.pars[i,] <- c(a,b,p)
    keep.ll[i] <- ll
    
    #======================================
    # Plots
    #======================================
    if(iters%%100==0){
      par(mfrow=c(3,2)){
        plot(keep.pars[1:i,1], type = "s", main = "a")
        plot(keep.pars[1:i,2], type = "s", main = "b")
        plot(keep.pars[1:i,3], type = "s", main = "p")
        plot(keep.ll[1:i], type = "s", main = "log-likelihood")
      }
    }
    
    
  }##end mcmc
  
  
  #======================================
  # Output
  #======================================
  
  return(list(pars = keep.pars,
              ll   = keep.ll))
  
  
}


