#! /usr/bin/R

library('parallel')
source('./null.header.R')

#############
#GLOBAL VARS#
#############

# The ALA banksia data for generating real banksia ranges :
banksia.data <- read.csv('./data/Banksia.data.csv', row.names=1)

# Generate polygons of each species range to randomly sample points from:
ranges <- gen.ranges(flataus, banksia.data, range.width)
# Throw out the species not being used in the final analysis:
names(ranges) <- gsub(' ', '_', names(ranges))
ranges <- ranges[names(ranges) %in% spat.niche[[1]]]

# Plot all the species' ranges on a map of Australia:
#plot(aus)
#for(range in ranges){plot(into.longlat(range), add=T, col=rancol())}

###########
#FUNCTIONS#
###########

resample.range <- function(i, rw=25000, data){
    # The resample null: conserves range position & randomises points
    pts <- extract.points(i, data)
    npts <- nrow(pts)
    range <- SpatialPoints(pts[ ,2:3], longlat)
    range <- into.cea(range)
    range <- gBuffer(range, width=rw)
    range <- gIntersection(flataus, range)
    return(spsample(range, npts, 'random'))
}

gen.ranges <- function(boundary, data, rw=25000){
    #
    ranges <- list()
    for(i in 1:length(unique(data$Species))){
        pts <- extract.points(i, data)
        range <- SpatialPoints(pts[ ,2:3], longlat)
        range <- into.cea(range)
        range <- gBuffer(range, width=rw)
        range <- gIntersection(boundary, range)
    ranges[[as.character(unique(data$Species)[i])]] <- range
    }
    print(sort(sapply(ranges, gArea)/1000000))
    return(ranges)
}

resample.null <- function(spat.niche, ranges, rasters){
    # spat.niche is a dataframe listing species names and abundances.
    # ranges is a list of polygons describing the ranges of the species listed
    # by spat.niche.
    # rasters is a raster::stack of the environmental layers of interest.
    # This function returns null niche data for each species, standardised,
    # hypervolumed and logged.
    # Ideally the parsing of the raw values extracted from the rasters would be
    # done by a seperate set of functions. 
    if(!(nrow(spat.niche) == length(ranges))){
        stop("Can't you do anything right?")}
    ffkd <- as.data.frame(matrix(ncol=dim(rasters)[3], nrow=0))
    for(i in 1:nrow(spat.niche)){
        pnts <- spsample(ranges[[i]], spat.niche[i, 2], 'random')#iter=100000)
        pnts <- into.longlat(pnts)
        ow <- as.data.frame(extract(rasters, pnts))
        ow <- standardise.env(ow)
        ow <- lapply(ow, function(x){sd(x, na.rm=TRUE)})
        ffkd <- rbind(ffkd, ow)
    }
    spat.niche[ ,'sd.h1v'] <- ffkd[ ,1]
    spat.niche[ ,'sd.h4v'] <- hypervolume(ffkd[ ,c(1, 2, 5, 6)])
    spat.niche[ ,'sd.h6v'] <- hypervolume(ffkd[ ,1:6])
    spat.niche <- cbind(spat.niche, ffkd)
    return(spat.niche)
}

doashuffle <- function(i){
    #just a wrapper that handles occasional errors thrown by the spsample method
    #would be better if they were grabbed at the source.
    outfolder <- 'shuffle_null_output/'
    tryCatch({j <- resample.null(spat.niche, ranges, rasters)
              write.csv(j, paste(outfolder, i, '_shuffle_null_raw.csv',
                                 sep=''), quote=F)},
             error=function(e){print(paste("Convergence error in ", i))}
             )
}

##########
##METHOD##
##########

system.time(mclapply(1:500, doashuffle, mc.cores=4, mc.preschedule=FALSE))

l <- list.files('./shuffle_null_output/')
l <- strsplit(l, '_')
l <- as.integer(sapply(l, '[[', 1))
failed <- (1:500)[!(1:500 %in% l)]
while(length(failed)>0){
    mclapply(failed, doashuffle, mc.cores=4, mc.preschedule=FALSE)
    l <- list.files('./scatter_null_output/')
    l <- strsplit(l, '_')
    l <- as.integer(sapply(l, '[[', 1))
    failed <- (1:501)[!(1:500 %in% l)]
}
