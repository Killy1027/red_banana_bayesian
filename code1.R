
# Calculating the kernel density for redbanan trees in square grids
library(MASS)
plot(xy[,1],xy[,2],xlim=range(xy[,1]),ylim=range(xy[,2]))
density <- kde2d(xy[,1],xy[,2],n=100,lims=c(range(xy[,1]),range(xy[,2])))
image(density)
z <- density$z
# To prove that the volume of density function equals 1
diff(range(xy[,2]))*diff(range(xy[,1]))/(100^2)*sum(z)
# Move those x's and y's to proper locations, crosses of grids.
half.grid.x <- diff(density$x)[1]/2
half.grid.y <- diff(density$y)[1]/2
new.x <- c(density$x - half.grid.x, density$x[100] + half.grid.x)
new.y <- c(density$y - half.grid.y, density$y[100] + half.grid.y)
# Calculating the corresponding position in grids
longitude.location <- latitude.location <- NULL
for (i in 1:nrow(xy))
{
  longitude.location[i] <- sum(xy[,1][i] < new.x)
  latitude.location[i] <- sum(xy[,2][i] < new.y)
}
# To test if vectors contain a value of 0. (impossible number)
0 %in% longitude.location
0 %in% latitude.location

density.sum <- 0
for(i in 1:nrow(xy))
{
  density.sum <- density.sum + z[longitude.location[i],latitude.location[i]]
}
density.sum


#======================================
# Get background rate for each event 
#======================================
N <- nrow(redbanana)
mu <- NULL
for(i in 1:N){
  mu[i] <- z[longitude.location[i],latitude.location[i]]/max(redbanana$birth)
}

#======================================
# Get inter-event times at each event  
#======================================
time.dist <- list()
for(i in 1:N){
  time.dist[[i]] <- t[i] - t[1:(i-1)]
}
lambda[i] <- (1-p)*mu/max(t)+p*sum((a*b/pi)*exp(-a*(iet[1:(i-1), i])-b*((x[i]-x[1:(i-1)])^2+(y[i]-y[1:(i-1)])^2)))



