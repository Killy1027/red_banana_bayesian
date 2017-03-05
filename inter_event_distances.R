#======================================
# Get interevent distances 
# (units are in meters)
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

# Calculate inter-event distances
library(fields)
ied <- rdist(xy) #distances in meters
ied2 <- ied^2 #squared distances



#======================================
# Get interevent times 
# (units are in weeks)
#======================================

iet <- rdist(redbanana$birth)
