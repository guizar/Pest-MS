# alternative highres map

library(maptools)
library(mapdata)

# Load HiRes map
map = map("worldHires", plot=F, fill = TRUE)

# Extract IDs 
IDs = sapply(strsplit(map$names, ":"), function(x) x[1])

# convert Map to SP
map = map2SpatialPolygons(map, IDs = IDs)

### Compare WorldPolyCountries VS map
setwd(wddata)
load("WorldPolyCountries.Rdata")

# look for bbox at http://bboxfinder.com/
# format xmin, xmax, ymin, ymax
ext = c(-91.45735,-87.06037,15.25057,22.27486) 
ext = c(109.643555,165.893555,-24.846565,0.263671)

# plot - map object seem to have higher res
plot(WorldPolyCountries,ylim=c(ext[3],ext[4]),xlim = c(ext[1],ext[2]))
plot(map,ylim=c(ext[3],ext[4]),xlim = c(ext[1],ext[2]))

## overlay
over = overlay(map, ALL)
summary(is.na(over))

## Result: more NAs are generated when overlaying ALL points against map (211) than ALL vs WorldPolyCountries (206)