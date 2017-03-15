#======================================
#======================================
# Simulate realizations from the
# ETAS model of red banana
#======================================
#======================================

sim_rb_etas <- function(theta, t1, bg){
  
  a = theta[1]; b = theta[2]; p = theta[3]
  br = p
  cat("Expected branching ratio = ", br,". ")
  if(br > 1){
    cat("The process is explosive!")
    return(0)
  }
  
  xlen <- nrow(bg$z)
  ylen <- ncol(bg$z)
  x1 <- diff(range(bg$x)) + diff(bg$x)[1]
  y1 <- diff(range(bg$y)) + diff(bg$y)[1]
  minx <- min(bg$x) - diff(bg$x)[1]/2
  miny <- min(bg$y) - diff(bg$y)[1]/2
  
  
  #======================================
  # Simulate M background events
  #======================================
  
  # M mainshocks over a span of t1 weeks
  M = rpois(1, (1-p)*x1*y1*mean(bg$z)) ##should be about (1-p)*n
  ## if p=0 (everything is background), this should be 788.
  cat("\n Number of mainshocks = ", M, ".\n")
  if(M < .5) return(99999)  ## this means no points were generated at all.
  
  #x = rep(0,M)      ## longitudes
  #y = rep(0,M)      ## latitudes
  #t = rep(0,M)  ## these will be the times of the events.
  
  prob <- as.vector(bg$z/sum(bg$z))
  samp <- sample(xlen*ylen, M, replace = FALSE, prob = prob)  
  
  yind <- ceiling(samp/ylen) ## retrieves y value of matrix (lat)
  xind <- samp%%xlen  ## retrieves x value of matrix (lon)
  xind[xind==0] <- xlen
  
  x <- (xind/xlen)*x1 + min(bg$x) ## longitudes of mainshocks
  y <- (yind/ylen)*y1 + min(bg$y) ## latitudes of mainshocks
  t <- runif(M)*t1  ## times of the mainshocks
  
  #create a vector of the type of plant: mainshock or aftershock
  type <- rep("main", M)
  
  #======================================
  # Now make aftershocks
  #======================================
  
  #vector of the number of (direct) aftershocks each mainshock will have
  after = rpois(M, br)
  
  stop1 = M
  nnext = sum(after)
  
  #set the locations and times for every aftershock
  oldx = rep(x[1:M], after)
  oldy = rep(y[1:M], after)
  oldt = rep(t[1:M], after)
  
  stop9 = 0
  while(stop9<1){
    t[(stop1+1):(stop1+nnext)] = oldt + rexp(nnext,a)
    after[(stop1+1):(stop1+nnext)] = rpois(nnext,br)
    dist1 = sqrt(rexp(nnext,b))
    
    thet1 = runif(nnext)*2*pi
    x[(stop1+1):(stop1+nnext)] = cos(thet1)*dist1 + oldx
    y[(stop1+1):(stop1+nnext)] = sin(thet1)*dist1 + oldy
    ## now get ready for next loop:
    oldt = rep(t[(stop1+1):(stop1+nnext)],after[(stop1+1):(stop1+nnext)])
    oldx = rep(x[(stop1+1):(stop1+nnext)],after[(stop1+1):(stop1+nnext)])
    oldy = rep(y[(stop1+1):(stop1+nnext)],after[(stop1+1):(stop1+nnext)])
    nextstop1 = stop1+nnext
    nnext = sum(after[(stop1+1):(stop1+nnext)])
    stop1 = nextstop1
    if(nnext < .5) stop9 = 2
  }
  
  type <- c(type, rep("after", length(t)-M))
  
  # Now get rid of aftershocks outside of S, and points with times > t1.
  keeps = which(t<t1 & (minx<x & x<(x1+minx)) & (miny<y & y<(y1+miny)))
  t <- t[keeps]
  x <- x[keeps]
  y <- y[keeps]
  type <- type[keeps]

  cat("Final number of earthquakes = ", sum(keeps), ".\n")
  

  return(list(pts   = cbind(x=x[order(t)], y=y[order(t)]),
              times = sort(t),
              type  = type[order(t)]
              ))
  
}




if(F){
  #======================================
  #======================================
  # Example
  #======================================
  #======================================
  
  theta = c(.076, .029, .577)
  theta = c(.3, 20, .96)
  sim1 <- sim_rb_etas(theta, 194.5509, bg)
  plot(cbind(sim1$lon,sim1$lat),xlim=range(bg$x),ylim=range(bg$y),pch=".")
  
  
}
