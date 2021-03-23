library('ape')
library('rgeos')
library('sp')
library('rgdal')
library('phytools')

source('../null.header.R')
source('../scatter.null.R')

bmSim <- function(ancpnt, blen, boundary, move.r=1000, move.v=100){
    # Takes a point as input and moves it randomly around the map to return
    # its decendant. Equates to a species' range shifting along its branch.
    # ancpnt = a sp point representing the centroid
    # move.r is the average displacement per step
    # move.v is the variance in the displacement per step
    if(!any(over(boundary, ancpnt)))
        stop("Ya just can't do anything right can ya?")
    # If we're taking less than one step do nothing:
    ###if(blen < grain){return(ancpnt)}
    steps <- blen
    move.r <- move.r*steps
    move.wv <- move.v*steps
    decpnt <- NULL
    tmp.decpnt <- NULL
    i <- 0
    while(is.null(decpnt)){
        r <- abs(rnorm(1, move.r, move.wv+(i*100)))
        if(r<=move.v){
            decpnt <- ancpnt} else {
        #print(paste(steps, move.r, move.wv, r))
                wf <- gDifference(gBuffer(ancpnt, width=r+move.v),
                                  gBuffer(ancpnt, width=r-move.v))
                wf <- gIntersection(wf, boundary)
                try(tmp.decpnt <- spsample(wf, 1, type='random',
                                           iter=1000000), silent=T)
                decpnt <- tmp.decpnt}
        i <- i+1
    }
    return(decpnt)
}

range.split <- function(ancpnt, boundary, target.r){
    # This method doesn't really work for splits less than ~200m so
    # if they really wanted a sympatric split:
    if(target.r < 200) {return(ancpnt@coords[c(1,1), ])}
    # Else a true allopatric split:
    bg <- gDifference(gBuffer(ancpnt, width=(target.r/2)+200),
                      gBuffer(ancpnt, width=(target.r/2)-100))
    #decendant candidates:
    deckhand <- spsample(bg, 20, 'regular')
    indxs <- which(spDists(deckhand)>target.r, arr.ind=T)
    #plot(bg)
    #plot(deckhand, add=T)
    #plot(deckhand[indxs[sample(1:nrow(indxs), 1), ], ], add=T, col='red')
    return(deckhand[indxs[sample(1:nrow(indxs), 1), ], ])
}

phylo.tracer <- function(phy, boundary, total.dispersal, total.variance,
                         allopatric.r=0, outfolder='./output/'){
    # A super closure (should really be a class). Ultimately the returned
    # function will just take a number (i) and then run a repetition of the
    # parameters. This can easily be extended to say, specify the starting
    # point by ancpnt=foo to the resulting closure.
    edges2 <- SpatialPointsDataFrame(
        coords=matrix(rep(7, 2*nrow(phy$edge)), ncol=2),
        data=data.frame("anc"=phy$edge[,1], "dec"=phy$edge[,2],
            "btime"=NA,"edgelength"=phy$edge.length),
        proj4string=boundary@proj4string
    )
    move.r <- total.dispersal/max(branching.times(phy))
    if(is.null(total.variance)){ total.variance <- total.dispersal/13}
    move.v <- total.variance/max(branching.times(phy))
    all.steps <- floor(max(branching.times(phy)))
    print(paste("Over the", all.steps, "steps of the simulation, points will",
                "move up to", total.dispersal/1000, "km."))
    ntip <- Ntip(phy)
    outfolder <- outfolder
    trace <- function(branchID, ancpnt=spsample(boundary, 1, 'random',
                       iter=100000)){
        # Hidden from the user this is the bit that uses functional recursion
        # to trace down the phylogeny spliting at each node. Keeping track of
        # where different species are located on the map.
        p <- edges2[branchID, ]
        decs <- which(edges2$anc==edges2@data[branchID, 'dec'])
        position <- bmSim(ancpnt, p@data[1, 'edgelength'], boundary, move.r, move.v)
        p@coords[1, ] <- position@coords[1, ]
        if(length(decs)==0) return(p)
         #print(p)
        dec.pnts <- range.split(position, boundary, allopatric.r)
        for(dec in decs){
            p <- rbind(p, trace(dec, position))}
        return(p)
    }
    run <- function(ancpnt=spsample(boundary, 1, 'random', iter=100000)){
        # Sets the phylogenetic tracing off from the root node.
        # Could also take arguments to override those given to the closure.
        position <- ancpnt
        decs <- which(edges2@data[ , 'anc']==edges2@data[1, 'anc'])
        dec.pnts <- range.split(position, boundary, allopatric.r)
        j <- lapply(decs, trace, ancpnt=position)
        p <- edges2[0, ]
        for(x in j){p <- rbind(p, x)}
        p[['species']] <- phy$tip[p$dec]
        return(p)
    }
    run.write <- function(i, ...){
        p <- run(...)
        dir.create(paste(outfolder, as.integer(total.dispersal/1000),
                         '/', sep=''))
        writeOGR(p, paste(outfolder, as.integer(total.dispersal/1000),
                          '/', sep=''), i,  driver="ESRI Shapefile")
    }
    return(run.write)
}

### These two clunky functions trace back over the output of phylotracer (one
### down each branch from the root node) and can plot it as lines on the map:

line.trace <-function(branchID, tracerout){
        p <- tracerout[branchID, ]
        p <- rbind(p, tracerout[tracerout@data[ , 'dec']==p$anc, ])
        p <- coordinates(p)
        l <- list(Lines(Line(p), branchID))
        if(tracerout@data[branchID, 'dec'] < 134){
            col <- fsp134[branchID]
        }else{ col <- rgb(0.3, 0.3, 0.3, alpha=0.3)}
        plot(SpatialLines(l), add=T, col=col, lwd=2)
                                        #readline('ok?')
        decs <- which(tracerout$anc==tracerout@data[branchID, 'dec'])
        for(dec in decs){
            l <- c(l, line.trace(dec, tracerout))}
                                        #print(l)
        return(l)
    }

line.trace2 <- function(tracerout, branchID=0){
    if(branchID==0){
        root.node <- tracerout@data[1, 'anc']
        p <- tracerout[tracerout@data[ , 'anc']==root.node, ]
    } else {p <- tracerout[branchID, ]}
    decs <- which(tracerout$anc %in% p$dec)
    p <- coordinates(p)
    l <- list(Lines(Line(p), branchID))
    #if(tracerout@data[branchID, 'dec'] < 134){
    #    col <- fsp134[branchID]
    #}else{ col <- rgb(0.3, 0.3, 0.3, alpha=0.3)}
    #plot(SpatialLines(l), add=T, col=col, lwd=2)
    #readline('ok?')
    
    for(dec in decs){
        l <- c(l, line.trace2(tracerout, dec))}
    return(l)
}

jj3 <- line.trace(12, j3)

# To make a flipbook of descent being over the australia map:
for(i in 1:length(jj3)){
    png(paste('animation/', i, '.png', sep=''))
    plot(flataus, xlim=bbox(SpatialLines(jj3))[1,]+100000, ylim=bbox(SpatialLines(jj3))[2,]+100000)
    plot(SpatialLines(jj3[1:i]), col=col, lwd=2, add=T)
    box()
dev.off()
}

# To plot how species moved around the map as they evolved:
plot(flataus)
plot(SpatialLines(line.trace(1, di(12))), add=T,
     col=rgb(0.5, 0.7, 1, alpha=0.6), lwd=2.0)

SpatialLines(line.trace(1, di(12)))

plot(SpatialLines(line.trace(1, di(12))), add=T, col='blue')

plot(SpatialLines(list(Lines(Line(j), 'crap'))), add=T)

dranges <- ranges[names(ranges) %in% Dryandras]
drycent <- lapply(dranges, gCentroid)
cents <- drycent[[1]]
for(x in drycent[2:length(drycent)]){cents <- gUnion(cents, x)}
j <- spDists(cents)
which(j==max(j), TRUE)

branching.times(drop.tip(phy, phy$tip[!(phy$tip %in% names(dranges[c(1, 53)])
                                        )]))

branching.times(drop.tip(phy, phy$tip[!(phy$tip %in% names(dranges))]))

drop.tip(mcc.post[[1]], mcc.post[[1]]$tip[!(mcc.post[[1]]$tip %in%
                                            names(dranges[c(1, 53)]))])

driver="ESRI Shapefile"

esri <- readOGR('dryandra_output/890/', '2')

j <- dryrun[[2]]()


dryrun <- list(phylo.tracer(drytre, drybb, 400*1000, 20000,
                           outfolder='./dryandra_output/'),
              phylo.tracer(drytre, drybb, 890*1000, 20000,
                           outfolder='./dryandra_output/'),
              phylo.tracer(drytre, drybb, 2000*1000, 800*1000,
                           outfolder='./dryandra_output/')
              )

library('parallel')

system.time({
for(func in dryrun){ 
    mclapply(1:500, func, mc.preschedule=FALSE, mc.cores=4)}
})

simd <- funcs[i]()
writeOGR(simd, paste('./drayandra_output/', i, '/' sep='') i,  driver="ESRI Shapefile")

esri[esri@data[ , 'dec']<54, ]
