get.mu.indices <- function(xy, kern){
  
  half.grid.x <- diff(kern$x)[1]/2
  half.grid.y <- diff(kern$y)[1]/2
  new.x <- c(kern$x - half.grid.x, max(kern$x) + half.grid.x)
  new.y <- c(kern$y - half.grid.y, max(kern$y) + half.grid.y)
  
  # Calculating the corresponding position in grids
  x.index <- y.index <- NULL
  for (i in 1:NROW(xy)){
    x.index[i] <- sum(new.x < xy[i,1])
    y.index[i] <- sum(new.y < xy[i,2])
  }
  
  return(list(x = x.index, y = y.index))
}



#======================================
#======================================
# MCMC Sampler
#======================================
#======================================
RedBananaMCMC <- function(xy, t, bg,
                          iters = 1000, burn = iters/2, thin = 1,
                          prior.shape.a = 1, prior.rate.a = 10,
                          prior.shape.b = 1, prior.rate.b = 10,
                          tuner = c(0.1, 0.1, 1)
                          ){
  
  
  
  require(fields)
  
  #======================================
  # Initial values
  #======================================
  
  # Priors
  a <- rgamma(1, prior.shape.a, prior.rate.a)
  b <- rgamma(1, prior.shape.b, prior.rate.b)
  p <- runif(1)
  
  # Get interevent times
  iet <- rdist(t)
  iet[lower.tri(iet, diag = T)] <- NA
  
  # Get interevent squared distances
  ied <- rdist(xy)
  ied2 <- ied^2
  ied2[lower.tri(ied2, diag = T)] <- NA
  
  # Get background rate for each event 
  n <- length(t)
  index <- get.mu.indices(xy, bg)
  mu <- NULL
  for(i in 1:n){
    mu[i] <- bg$z[index$x[i],index$y[i]]/diff(range(t))
  }
  
  # Get integral of background rate
  x1 <- diff(range(bg$x)) + diff(bg$x)[1]
  y1 <- diff(range(bg$y)) + diff(bg$y)[1]
  intmu <- x1*y1*mean(bg$z)

  # loglikelihood function 
  trigA <- a*exp(-a*iet) #a
  trigB <- b*exp(-b*ied2) #b
  trig <- colSums(trigA*trigB, na.rm = T) #a,b
  lam1 <- (1-p)*mu #p
  lam2 <- p*trig/pi #a,b,p
  
  sumloglam <- sum(log(lam1 + lam2)) #a,b,p
  intlambda <- (1-p)*intmu + p*n #p
  
  ll <- sumloglam - intlambda
  
  # tuning 
  acc <- rep(0,3)
  att <- 0
  
  
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
      
      att <- att + 1
      
      #__________________________________
      # Update a
      #
      cana <- exp(rnorm(1, log(a), tuner[1]))
      cantrigA <- cana*exp(-cana*iet)
      cantrig <- colSums(cantrigA*trigB, na.rm = T)
      canlam2 <- p*cantrig/pi
      cansumloglam <- sum(log(lam1 + canlam2))

      Ratio.a <- cansumloglam - sumloglam + dgamma(cana, prior.shape.a, prior.rate.a, log = T) - dgamma(a, prior.shape.a, prior.rate.a, log = T)
    
      if(log(runif(1)) < Ratio.a){
        a <- cana
        trigA <- cantrigA
        trig <- cantrig
        lam2 <- canlam2
        sumloglam <-  cansumloglam
        ll <- sumloglam - intlambda
        acc[1] <- acc[1] + 1
      }
      
      #__________________________________
      # Update b
      #
      canb <- exp(rnorm(1, log(b), tuner[2]))
      cantrigB <- canb*exp(-canb*ied2)
      cantrig <- colSums(trigA*cantrigB, na.rm = T)
      canlam2 <- p*cantrig/pi
      cansumloglam <- sum(log(lam1 + canlam2))

      Ratio.b <- cansumloglam - sumloglam + dgamma(canb, prior.shape.b, prior.rate.b, log = T) - dgamma(b, prior.shape.b, prior.rate.b, log = T)
      
      if(log(runif(1)) < Ratio.b){
        b <- canb
        trigB <- cantrigB
        trig <- cantrig
        lam2 <- canlam2
        sumloglam <-  cansumloglam
        ll <- sumloglam - intlambda
        acc[2] <- acc[2] + 1
      }
      
      #__________________________________
      # Update p
      #
      canp <- rbeta(1, p/tuner[3], (1-p)/tuner[3])
      canlam1 <- (1-canp)*mu
      canlam2 <- canp*trig/pi
      cansumloglam <- sum(log(canlam1 + canlam2))
      canintlambda <- (1-canp)*intmu + canp*n
      canll <- cansumloglam - canintlambda
      
      Ratio.p <- canll - ll + dbeta(p, canp/tuner[3], (1-canp)/tuner[3], log = T) - dbeta(canp, p/tuner[3], (1-p)/tuner[3], log = T)
      
      if(log(runif(1)) < Ratio.p){
        p <- canp
        lam1 <-  canlam1
        lam2 <-  canlam2
        sumloglam <- cansumloglam
        intlambda <- canintlambda
        ll <- sumloglam - intlambda
        acc[3] <- acc[3] + 1
      }
    

    }##end thin
    
    
    #======================================
    # Keepers
    #======================================
    keep.pars[i,] <- c(a,b,p)
    keep.ll[i] <- ll
    
    
    #======================================
    # Tuning
    #======================================
    if(i < 0.75*burn & att > 50){
      tuner <- ifelse(acc/att < 0.25, 0.8*tuner, tuner)
      tuner <- ifelse(acc/att > 0.50, 1.2*tuner, tuner)
      acc <- rep(0,3)
      att <- 0
    }
    
    
    
    #======================================
    # Plots
    #======================================
    if(i%%100==0){
      par(mfrow=c(2,2))
      plot(keep.pars[1:i,1], type = "s", main = "a")
      plot(keep.pars[1:i,2], type = "s", main = "b")
      plot(keep.pars[1:i,3], type = "s", main = "p")
      plot(keep.ll[1:i], type = "s", main = "log-likelihood")
    }
    
    
  }##end mcmc
  
  
  #======================================
  # Output
  #======================================
  
  return(list(pars      = keep.pars,
              ll        = keep.ll,
              acc.ratio = acc/att,
              tuner     = tuner
              ))
  
  
}





#======================================
#======================================
#======================================
#======================================
#======================================
#======================================






if(F){
  
  #======================================
  #======================================
  # Example
  #======================================
  #======================================
  
  # Transform from lon/lat coordinates to Euclidean x/y coordinates
  require(stppResid)
  data("redbanana")
  lonlat <- SpatialPoints(coords = redbanana[,1:2])
  proj4string(lonlat) <- CRS("+proj=longlat +datum=WGS84")
  
  # Projection for Costa Rica
  xyproj <- "+proj=utm +zone=17 +ellps=WGS84"
  xy <- spTransform(lonlat, CRS(xyproj))
  xy <- coordinates(xy)
  
  # jitter locations
  xy <- xy + rnorm(788*2, 0, 3)
  
  xrange = sbox(xy)[1:2,1]
  yrange = sbox(xy)[2:3,2]
  
  require(MASS)
  h <- c(bandwidth.nrd(xy[,1]), bandwidth.nrd(xy[,2]))
  bg <- kde2d(xy[,1],xy[,2],h=h,n=400,lims=c(xrange,yrange))
  bg$z <- bg$z*NROW(xy)
  g <- RedBananaMCMC(xy, redbanana$birth, bg=bg, iters=1000, thin = 2)
  
}
