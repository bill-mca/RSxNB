#! /usr/bin/R

library('sp')
library('raster')
library('rgeos')
library('rgdal')

#############
#GLOBAL VARS#
#############

### Projections

# The 
# area: from flat = rgeos::gArea
#       from long lat = geosphere::areaPolygon
# The former executes ~250x faster

# Cylindrical Equal Area projection for the whole globe:
cea <-  CRS("+proj=cea +lon_0=Central Meridian
            +lat_ts=Standard Parallel
            +x_0=False Easting
            +y_0=False Northing")

longlat <- CRS('+proj=longlat +datum=WGS84')

### Species Occurence Data

# The ALA banksia occurrence data:
databanks <- read.csv('donebanks.csv', row.name=1)

### Geographical Data

# Australian coast defenition:
aus <- readOGR('data/aus', 'Australia')
proj4string(aus) <- longlat
flataus <- into.cea(aus)

#The rasters of Environmental data. Used for infering niche from presence data.
rasters <- paste('data/',
                 c( "temperature", "precipitation", "seasonality", "dry", 'pH',
                   'clay'), '.asc', sep='')
rasters <- stack(rasters)

# Banksia Background:
# A set of polygons defining the regions Banksia could feasibly occur
# for the 2ndrange scatter null.
banksia.background1 <- databanks[seq(1, nrow(databanks), 4), 2:3]
print(head(banksia.background1))
banksia.background1 <- SpatialPoints(banksia.background1, longlat)
banksia.background1 <- gBuffer(into.cea(banksia.background1),
                              width=25000)
banksia.background1 <- gIntersection(banksia.background1, flataus)


banksia.background2 <- databanks[seq(2, nrow(databanks), 4), 2:3]
print(head(banksia.background2))
banksia.background2 <- SpatialPoints(banksia.background2, longlat)
banksia.background2 <- gBuffer(into.cea(banksia.background2),
                              width=25000) 
banksia.background2 <- gIntersection(banksia.background2, flataus)


banksia.background3 <- databanks[seq(3, nrow(databanks), 4), 2:3]
print(head(banksia.background3))
banksia.background3 <- SpatialPoints(banksia.background3, longlat)
banksia.background3 <- gBuffer(into.cea(banksia.background3),
                              width=25000)
banksia.background3 <- gIntersection(banksia.background3, flataus)

banksia.background4 <- databanks[seq(4, nrow(databanks), 4), 2:3]
print(head(banksia.background4))
banksia.background4 <- SpatialPoints(banksia.background4, longlat)
banksia.background4 <- gBuffer(into.cea(banksia.background4),
                              width=25000)
banksia.background4 <- gIntersection(banksia.background4, flataus)

banksia.background <- gUnion(banksia.background1, banksia.background2)
banksia.background <- gUnion(banksia.background, banksia.background3)
banksia.background <- gUnion(banksia.background, banksia.background4)

rm(banksia.background1, banksia.background2, banksia.background4,
   banksia.background3)

banksia.background <- gIntersection(flataus, gBuffer(banksia.background,
                                                     width=15000))

# min.max is a dataframe of the observed extreme min and max value occupied by
# an individual of the study group it is used to to standardise all data to the
# range of 0 (observed min) to 1 (observed max) eg:
#                min  max
#temperature    40.0  289
#precipitation 199.0 3298
#seasonality   991.0 6196 ...etc
min.max <- read.csv('data/min_max.csv', row.names=1)

#Load min.max into the standardise.env closure.
standardise.env <- env.standaridiser(min.max)

# How far from an occurence point is defined as being on the species range:
range.width <- 25000#m

# spat.niche is a data frame used to list the names of species used in the final# analysis and their empirical range size and abundance.
spat.niche <- read.csv('spat.niche.csv', row.names=1)


###########
#FUNCTIONS#
###########


extract.points <- function(i, data){
    # Extracts the occurence points for a species as specified numerically or
    # by name. works for dataframes and spatial.points.dataframes.
    if(is.numeric(i)){
        return(data[data$Species==unique(data$Species)[i],])}
    else{
        return(data[data$Species==i,])}
}
into.cea <- function(geom, wasin=CRS('+proj=longlat')){
    # Converts a geometry into cylindrical equal area projection
    if(is.null(proj4string(geom))){
        proj4string(geom) <- wasin}
    return(spTransform(geom, cea))
}

into.longlat <- function(geom, wasin=cea){
    # Converts a geometry into cylindrical equal area projection
    if(is.null(proj4string(geom))){
        proj4string(geom) <- wasin}
    return(spTransform(geom, longlat))
}

env.standardiser <- function(mm){
    #Standardises environmental variables so that values range from 0 -> 1.
    #Uses min.max which lists the extreme minimum and maximum value from the
    #empirical data. min.max must have the same columns in the same order as
    # the input data.
    return(function(data){
    for(i in 1:ncol(data)){
        low <- mm[i, 1]
        high <- mm[i, 2]
        #print(c(low, high))
        out <- (data[[i]]-low)/high
        data[[i]] <- out}
    return(data)})
}

#random colour function
rancol <- function(){rgb(runif(1), runif(1), runif(1), alpha=0.3)}

product <- function(v){
    # takes the product of all values in a vector. c(A, B, C) -> AxBxC 
    out <- 1
    for(x in v){out <- out*x}
    return(out)
}

hypervolume <- function(df){
    # Takes a dataframe of numeric values and returns a single vector. It takes
    # the product of values row wise and takes the log (base10) of the result.
    l <- nrow(df)
    out <- numeric(l)
    for(i in 1:length(out)){
        out[i] <- product(as.numeric(df[i, ])*100)}
    return(log10(out))
}

brickstack.to.raster.list <- function(brickstack){
    if((class(brickstack)!="RasterStack") & (class(brickstack)!="RasterBrick"))
	{
            print("Input must be a RasterStack or RasterBrick")
            return()
	}
    brickstack_nlayers <- nlayers(brickstack)
    brickstack_pos <- 1:brickstack_nlayers
    #raster_list=vector("list",brickstack_nlayers)
    raster_list <- mapply(
        function(brickstack,layer){raster(brickstack,layer=layer)},
        brickstack_pos,MoreArgs = list(brickstack=brickstack),
        SIMPLIFY=FALSE)
    return(raster_list)
}

