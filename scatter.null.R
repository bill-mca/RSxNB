#! /usr/bin/R

library('ggplot2')
library('parallel')
source('null.header.R')

# ggplot2::fortify is an easy way of finding all the vertices of a complex
# polygon.

#############
#GLOBAL VARS#
#############

# A cruder version of the banksia.background polygon has less vertices which
# speeds up the execution of the range scatter method - even without sacrificing
# accuracy of the resulting points.
bankback <- gSimplify(banksia.background)

###########
#FUNCTIONS#
###########

myspsamp <- function(polygon, n=1, type='random'){
    # Should eliminate "iteration did not converge errors without costing speed
    tryCatch(expr={p <- spsample(polygon, n, type)},
             error={p <- spsample(polygon, n, type, iter=10000)},
             finally=return(p))
}

grow.range <- function(size, start.position, boundary, rw=25000){
    # This function takes a point in space and builds a range, of a specified
    # size (in m^2) around it. It restricts the range to occur within the
    # polygon 'boundary'. It only works well with polygons in flat coordinates.
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
        iters <- iters + 1
        if(gArea(tryrange) > maxtarget)
            {next} else {outr <- tryrange}
    }
    #plot(outr, add=T, col='')
    return(outr)
}

.rangerr <- function(Size, stps, bg, rw=25000){
    # A Wrapper function because gIntersection throws occasional errors that
    # derail grow.range. The issue can't be resolved at the gIntersection level.
    # the error is reproducible so that once a range being constructed by
    # grow.range has a fault ("self intersection") it can't be recovered.
    new.range <- NULL
    #print(c(i, Size))
    wig <- 0
    while(is.null(new.range)){
        # Grow Range will throw error messages; see its comments.
        new.range <- tryCatch(grow.range(Size, stps, bg,
                                         rw=range.width),
                              error=function(e){new.range <- NULL})
        wig <- wig+1
        if(wig > 50){
            # Species with ranges smaller than 2000 km^2 aren't
            # workable with the grow range method.
            warning('Species ', i, ' of size ', Size,
                  ' causes an infinite loop.')
            break}
    }
    return(new.range)
}

scatter.null <- function(spat.niche, boundary, plotit=FALSE){
    # Uses grow.range to build randomly placed ranges for each species and then
    # infers niche from the resulting range.
    for(i in row.names(spat.niche)){
    #print(i)
    Size <- (spat.niche[i, 3])*1000000
        stps <- myspsamp(boundary, 1, 'random')
        new.range <- .rangerr(Size, stps, boundary)
        if(is.null(new.range)) next
        if(plotit){
            png(paste('randranges/', i, spat.niche[i,1], 'spat.png', sep='_'))
            plot(flataus)
            plot(boundary, add=T)
            plot(new.range, add=T, col=rancol())
            dev.off()
        }
    pnts <- myspsamp(new.range, spat.niche[i, 2], type='random')
    pnts <- into.longlat(pnts)
    ow <- as.data.frame(extract(rasters, pnts))
    ow <- standardise.env(ow)
    ow <- lapply(ow, function(x){sd(x, na.rm=TRUE)})
    ow <- as.data.frame(ow)
    spat.niche[i ,'sd.h1v'] <- ow[1] #Temperature
    spat.niche[i ,'sd.h4v'] <- hypervolume(ow[c(1, 2, 5, 6)]) #w/ SOIL
    spat.niche[i ,'sd.h6v'] <- hypervolume(ow[1:6])
    spat.niche <- cbind(spat.niche, ow)
    }
    return(spat.niche)
}

doascatter <- function(i){
    outfolder <- './scatter_null_output/'
    nd <-suppressWarnings(scatter.null(spat.niche, banksia.background))
#    h1v.pgls <- mass.pgls(formula(formulae[1]), mcc.post, nd)
    write.csv(nd, paste(outfolder, i, '_raw_scatter_null.csv', sep=''))
}

scatter.by.spec <- function(spat.niche, boundary, rw=25000){
    # Because larger range species take much longer it is better to do all reps
    # for each species and then combine them into reps of the whole genus.
    # This code probably doesn't work yet.
    Size <- (spat.niche[1, 3])*1000000
    out <- cbind(spat.niche[0,], data.frame('temperature'=numeric(),
                 'precipitation'=numeric(), 'seasonality'=numeric(),
                 'dry'=numeric(), 'pH'=numeric(), 'clay'=numeric()))
    i <- 0
    while(i<500){
        stps <- myspsamp(boundary, 1, 'random')
        new.range <- tryCatch(grow.range(Size, stps, boundary,
                                         rw=rw),
                              error=function(e){new.range <- NULL})
        if(is.null(new.range)) next
        pnts <- myspsamp(new.range, spat.niche[i, 2], type='random')
        pnts <- into.longlat(pnts)
        ow <- as.data.frame(extract(rasters, pnts))
        # The values are standardised to the empirical data where min=0, max=1
        ow <- standardise.env(ow)
        # Then the SD breadth is taken for each condition
        ow <- lapply(ow, function(x){sd(x, na.rm=TRUE)})
        ow <- as.data.frame(ow)
        hvs <- data.frame(ow[[1]], hypervolume(ow[c(1, 2, 5, 6)]),
                    hypervolume(ow[1:6]))
        or <- cbind(spat.niche[1, 1:3], hvs, ow)
        out <- rbind(out, or)
        i <- i+1
        #Autosave
        write.csv(paste('scatter_output/', spat.niche[1,1], '.csv', sep=''))
    }
    colnames(out) <- colnames(spat.niche)
    write.csv(paste('scatter_output/', spat.niche[1,1], '.csv', sep=''))
}

##########
##METHOD##
##########

system.time(mclapply(1:500, doascatter, mc.cores=4, mc.preschedule=FALSE))

# incase it crashes part way through:
l <- list.files('./scatter_null_output/')
l <- strsplit(l, '_')
l <- as.integer(sapply(l, '[[', 1))
failed <- (1:500)[!(1:500 %in% l)]
while(length(failed)>0){
    mclapply(failed, doascatter, mc.cores=4, mc.preschedule=FALSE)
    l <- list.files('./scatter_null_output/')
    l <- strsplit(l, '_')
    l <- as.integer(sapply(l, '[[', 1))
    failed <- (1:501)[!(1:500 %in% l)]
}
