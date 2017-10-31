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
                          alpha.a = 1, alpha.b = 1, 
                          hyper.alpha = 1, hyper.beta = 10,
                          tuner = c(0.1, 0.1, 1)
                          ){
  
  
  
  require(fields)
  
  #======================================
  # Initial values
  #======================================
  
  # Priors
  beta.a <- rgamma(1, hyper.alpha , hyper.beta)
  beta.b <- rgamma(1, hyper.alpha , hyper.beta)
  a <- rgamma(1, alpha.a, beta.a)
  b <- rgamma(1, alpha.b, beta.b)
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
  keep.pars <- matrix(0, iters, 5)
  colnames(keep.pars) <- c("a", "b", "p", "beta.a", "beta.b")
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
      new.beta.a <- rgamma(1, hyper.alpha + 1 * alpha.a ,  hyper.beta + a )
      cana <- exp(rnorm(1, log(a), tuner[1]))
      cantrigA <- cana*exp(-cana*iet)
      cantrig <- colSums(cantrigA*trigB, na.rm = T)
      canlam2 <- p*cantrig/pi
      cansumloglam <- sum(log(lam1 + canlam2))

      Ratio.a <- cansumloglam - sumloglam + dgamma(cana, alpha.a, new.beta.a, log = T) - dgamma(a, alpha.a, beta.a, log = T)
    
      if(log(runif(1)) < Ratio.a){
        beta.a <- new.beta.a
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
      new.beta.b <- rgamma(1, hyper.alpha + 1 * alpha.b ,  hyper.beta + b )
      
      canb <- exp(rnorm(1, log(b), tuner[2]))
      cantrigB <- canb*exp(-canb*ied2)
      cantrig <- colSums(trigA*cantrigB, na.rm = T)
      canlam2 <- p*cantrig/pi
      cansumloglam <- sum(log(lam1 + canlam2))

      Ratio.b <- cansumloglam - sumloglam + dgamma(canb, alpha.b, new.beta.b, log = T) - dgamma(b, alpha.b, beta.b, log = T)
      
      if(log(runif(1)) < Ratio.b){
        beta.b <- new.beta.b
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
    keep.pars[i,] <- c(a,b,p,beta.a,beta.b)
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
      par(mfrow=c(3,2))
      plot(keep.pars[1:i,1], type = "s", main = "a")
      plot(keep.pars[1:i,2], type = "s", main = "b")
      plot(keep.pars[1:i,3], type = "s", main = "p")
      plot(keep.pars[1:i,4], type = "s", main = "beta.a")
      plot(keep.pars[1:i,5], type = "s", main = "beta.b")
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
  g <- RedBananaMCMC(xy, redbanana$birth, bg=bg, iters=5000, thin = 10)
  
}


if(F){
  
  # plot redbanana data onto map of Costa Rica
  library(ggmap)
  lon <- range(redbanana$longitude)
  lat <- range(redbanana$latitude)
  Costa_Rica <- get_map(location=c(lon=mean(lon),lat=mean(lat)),zoom=15)
  ggmap(Costa_Rica,extent="device")+geom_point(data=redbanana,aes(x=longitude,y=latitude),col="deeppink2")
  
  
  # for alpha
  plot(2501:5000,g$pars[2501:5000,1],type="s",main=bquote("Trace Plot of " ~ alpha),xlab="",ylab=bquote( ~ alpha))
  abline(h=mean(g$pars[2501:5000,1]),col="red")
  hist(g$pars[2501:5000,1],main=bquote("Posterior of " ~ alpha),xlab="",col="skyblue")
  abline(v=mean(g$pars[2501:5000,1]),col="red")
  # for beta
  plot(2501:5000,g$pars[2501:5000,2],type="s",main=bquote("Trace Plot of " ~ beta),xlab="",ylab=bquote( ~ beta))
  abline(h=mean(g$pars[2501:5000,2]),col="red")
  hist(g$pars[2501:5000,2],main=bquote("Posterior of " ~ beta),xlab="",col="skyblue")
  abline(v=mean(g$pars[2501:5000,2]),col="red")
  # for p
  plot(2501:5000,g$pars[2501:5000,3],type="s",main=bquote("Trace Plot of " ~ rho),xlab="",ylab=bquote( ~ rho))
  abline(h=mean(g$pars[2501:5000,3]),col="red")
  hist(g$pars[2501:5000,3],main=bquote("Posterior of " ~ rho),xlab="",col="skyblue")
  abline(v=mean(g$pars[2501:5000,3]),col="red")
  # for beta.a
  plot(2501:5000,g$pars[2501:5000,4],type="s",main=("Trace Plot of beta.a"),xlab="",ylab="beta.a")
  abline(h=mean(g$pars[2501:5000,4]),col="red")
  hist(g$pars[2501:5000,4],main=("Posterior of beta.a"),xlab="",col="skyblue")
  abline(v=mean(g$pars[2501:5000,4]),col="red")
  # for beta.b
  plot(2501:5000,g$pars[2501:5000,5],type="s",main=("Trace Plot of beta.b"),xlab="",ylab="beta.b")
  abline(h=mean(g$pars[2501:5000,5]),col="red")
  hist(g$pars[2501:5000,5],main=("Posterior of beta.b"),xlab="",col="skyblue")
  abline(v=mean(g$pars[2501:5000,5]),col="red")
  
  
  # calculating the standard deviation of parameters ( compared with MLE method)
  summary(g$pars)
  sd(g$pars[2501:5000,1])
  sd(g$pars[2501:5000,2])
  sd(g$pars[2501:5000,3])
  sd(g$pars[2501:5000,4])
  sd(g$pars[2501:5000,5])
  
  
  # creditable intetval
  quantile(g$pars[,1],c(0.025,0.975))
  quantile(g$pars[,2],c(0.025,0.975))
  quantile(g$pars[,3],c(0.025,0.975))
  quantile(g$pars[,4],c(0.025,0.975))
  quantile(g$pars[,5],c(0.025,0.975))
  
}



