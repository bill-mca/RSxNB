library('sp')
library('raster')
library('rgeos')
library('rgdal')
library('ggplot2')
library('caper')

cea <-  CRS("+proj=cea +lon_0=Central Meridian
            +lat_ts=Standard Parallel
            +x_0=False Easting
            +y_0=False Northing")

longlat <- CRS('+proj=longlat +datum=WGS84')

# The real life banksia occurrence data:
databanks <- read.csv('donebanks.csv', row.name=1)

# Australian coast defenition:
aus <- readOGR('./aus', 'Australia')
proj4string(aus) <- longlat
flataus <- into.cea(aus)

rasters <- paste(c( "temperature", "precipitation", "seasonality", "dry", 'pH',
                   'clay'), '.asc', sep='')
rasters <- stack(rasters)

# min.max is a dataframe of the real min and max observed in each env layer
# for the genus:
#                min  max
#temperature    40.0  289
#precipitation 199.0 3298
#seasonality   991.0 6196 ...etc
min.max <- read.csv('min_max.csv', row.names=1)

# How far from an occurence point is defined as being on the species range:
range.width <- 25000#m

spat.niche <- read.csv('spat.niche.csv', row.names=1)

formulae <-
    c("`bufrange.cea`~`sd.h1v`", "`bufrange.cea`~`sd.h4v`", "`bufrange.cea`~`sd.h6v`")
ealumrof <- as.vector(lapply(strsplit(formulae, '~'), function(x){paste(x[2:1], collapse='~')}), mode='character')

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

standardise.env <- function(data, mm=min.max){
    for(i in 1:ncol(data)){
        low <- mm[i, 1]
        high <- mm[i, 2]
#        print(c(low, high))
        out <- (data[[i]]-low)/high
        data[[i]] <- out}
    return(data)
}

product <- function(v){
    out <- 1
    for(x in v){out <- out*x}
    return(out)
}

hypervolume <- function(df){
 l <- nrow(df)
 out <- numeric(l)
 for(i in 1:length(out)){
     out[i] <- product(as.numeric(df[i, ])*100)}
 return(log10(out))
}

grow.range <- function(size, start.position, boundary, rw=range.width){
    stopifnot(size >= pi*rw^2)
    maxtarget <- size + pi*rw^2
    mintarget <- size - ((pi*rw^2)/3)
    outr <- gIntersection(gBuffer(start.position, width=rw), boundary)
    iters <- 0
    fragif <- runif(1, 50, 500)
    cea <- get('cea', envir=.GlobalEnv)
    longlat <- get('longlat', envir=.GlobalEnv)
    while(mintarget > gArea(outr)){
        if(iters >= fragif){
            # This statement will split the range depending on how many
            # iterations it has been trying to grow for. This way the more
            # widespread / 'hemmed in' a range is, the more likely it is to
            # be split into several populations.
##            print('fraging')
            far <- gDifference(gBuffer(outr, width=7e+5), gBuffer(outr, width=1e+5))
            far <- gIntersection(boundary, far)
            trypoint <- spsample(far, 1, 'random')
            iters <- 0
        } else {if(runif(1)> 0.01){
            # ggplot2::fortify is a cheap way of finding all the vertices
            # of a complex polygon.
            pboarder <- SpatialPoints(fortify(outr)[1:2], cea)
            trypoint <- pboarder[sample(1:length(pboarder), 1), ]
        } else {trypoint <- spsample(outr, 1, type='random')}}
        # possible problem here w/ gInt going before merge
        tryrange <- gBuffer(trypoint, width=rw)
        tryrange <- gIntersection(gUnion(outr, tryrange), boundary)
        if(!gIsValid(tryrange)){
            # seems to help with issues
            tryrange <- gSimplify(tryrange, 0.01)
            if(!gIsValid(tryrange)){
                warning('Having problems with a range.')}
           }
        # The range was too small given that the loop started a new iteration
        # step. before we check if it still needs to be expanded we'll check
        # that it's not over sized just incase. (This step is paranoid as it is
        # testing for what ought to be an impossible scenario.)
        iters <- iters + 1
        if(gArea(tryrange) > maxtarget)
            {next} else {outr <- tryrange}
    }
    #plot(outr, add=T, col='')
    return(outr)
}

.rangerr <- function(Size, stps, bg=banksia.background, rw=range.width){
    new.range <- NULL
    #print(c(i, Size))
    wig <- 0
    while(is.null(new.range)){
        # Grow Range will throw error messages; see its comments.
        new.range <- tryCatch(grow.range(Size, stps, bg, rw=range.width),
                              error=function(e){new.range <- NULL})
        wig <- wig+1
        if(wig > 100){
            # Species with ranges smaller than 2000 km^2 aren't
            # workable with the grow range method.
            warning('Species ', i, ' of size ', Size,
                  ' causes an infinite loop.')
            break}
    }
    return(new.range)
}

tracer.ranges <- function(tracer.tips, banksia.background, plotit=FALSE){
    #tracer.tips is spdf with species, points, rsize, nb1, nb2. nb3
    for(i in row.names(tracer.tips@data)){
    #print(i)
        Size <- (tracer.tips@data[i, 3])*1000000
        stps <- tracer.tips[i, ]
        new.range <- .rangerr(Size, stps, banksia.background)
        if(is.null(new.range)) next
        if(plotit){
            png(paste('randranges/', i, tracer.tips@data[i,1],
                      'spat.png', sep='_'))
            plot(flataus)
            plot(banksia.background, add=T)
            plot(new.range, add=T, col=rancol())
            dev.off()
        }
        pnts <- myspsamp(new.range, tracer.tips@data[i, 2], type='hexagonal')
        pnts <- into.longlat(pnts)
        ow <- as.data.frame(extract(rasters, pnts))
        ow <- standardise.env(ow)
        ow <- lapply(ow, function(x){sd(x, na.rm=TRUE)})
        ow <- as.data.frame(ow)
        tracer.tips@data[i ,'sd.h1v'] <- ow[2] #PRECIPITATION
        tracer.tips@data[i ,'sd.h4v'] <- hypervolume(ow[c(1, 2, 5, 6)]) #w/ SOIL
        tracer.tips@data[i ,'sd.h6v'] <- hypervolume(ow[1:6])
    }
    return(tracer.tips)
}

library('parallel')
system.time(mclapply(c(1426:1500), doabit, mc.cores=3, mc.preschedule=FALSE))

dranges <- ranges[names(ranges) %in% Dryandras]
drybb <- dranges[[1]]
for(i in 2:length(dranges)){drybb <- gUnion(drybb, dranges[[i]])}
drybb <- gBuffer(drybb, width=15000)
drybb <- gIntersection(drybb, flataus)
drybb <- gSimplify(drybb, 1000, topologyPreserve=T)

library('parallel')
system.time(mclapply(c(1:250), dryandrascat, mc.cores=4, mc.preschedule=FALSE))

rancol <- function(){rgb(runif(1), runif(1), runif(1), alpha=0.15)}

png('Dryandra_ranges.png', width=480*2, height=480*2)
bb <- bbox(into.longlat(wcbb))
plot(extent(bb), col='white', axes=F, xlab='', ylab='')
plot(aus, add=T)
#plot(wcbb, add=T)
for(r in dranges[sample(1:53, 53)]){plot(into.longlat(r), col=rancol(), add=T)}
dev.off()

png('Banksia_ranges.png', width=480*2, height=480*2)
#bb <- bbox(into.longlat(wcbb))
#plot(extent(bb), col='white', axes=F, xlab='', ylab='')
plot(aus)
#plot(wcbb, add=T)
for(r in ranges){plot(into.longlat(r), col=rancol(), add=T)}
dev.off()


esri <- readOGR('./dryandra_output/890/', i)
plot(drybb)
plot(esri, add=T)
i <- i+1

sbs <- esri[esri$dec<esri$anc[1], 'species']
sbs[[1]] <- as.character(sbs[[1]])
sbs@data[ ,2:6] <- db[ , 2:6]

db <- spat.niche[0, ]
for(spec in sbs[[1]]){db <- rbind(db, spat.niche[spat.niche$species==spec, ])}
spat.niche[[1]] <- as.character(spat.niche[[1]])]

doit <- function(dsn){
dsn <- dsn
return(function(i){
    esri <- readOGR(dsn, i)
    sbs <- esri[esri$dec<esri$anc[1], 'species']
    sbs[[1]] <- as.character(sbs[[1]])
    sbs@data[ ,2:6] <- db[ , 2:6]
    j <- tracer.ranges(sbs, drybb)
    write.csv(j@data, paste(dsn, i, '_niches.csv', sep=''))
})
   }

funkys <- list(doit('./dryandra_output/400/'),
doit('./dryandra_output/890/'),
doit('./dryandra_output/2000/'))

require('parallel')
for(funk in funkys){
mclapply(1:500, funk, mc.cores=4, mc.preschedule=F)}
