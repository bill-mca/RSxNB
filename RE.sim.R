library(gstat)
library(ape)
library(sp)
library(phytools)
library(parallel)

########################
## ORIGINAL FUNCTIONS ##
########################

### These are scarcely modified from the code used for the range overlap paper:

#######-------------------------------------------------------------------------
##  plot a rectangular range defined by north, south, east & west limits:
## boundaries of the domain are determined by lims

plotRange <- function(n,s,e,w,lims=c(0,200),col="red",add=FALSE){
  if(!add) plot(lims[1]:lims[2],lims[1]:lims[2],type="n")
  arrows(e,n,w,n,length=0,col=col)
  arrows(e,s,w,s,length=0,col=col)
  arrows(e,n,e,s,length=0,col=col)
  arrows(w,n,w,s,length=0,col=col)
  }

#######-------------------------------------------------------------------------
## simulate gradual (Brownian motion) drift of range boundaries along one phylogenetic branch:
## blen is the branch length, sd is the rate parameter
## range expansion could be simulated by setting mean>0 for n & e boundaries, and mean<0 for s & w boundaries

bmSim <- function(blen,grain=0.1,mean=0,sd=0.35,n,s,e,w){
    steps <- floor(blen/grain)
    for(i in 1:steps){
        if((n-s)+(e-w)<0.02){
            n<-n
            s<-s
            e<-e
            w<-w}else{
                  n <- n+rnorm(1,mean,sd)
                  if(n<=s) n<-s+0.01
                  s <- s+rnorm(1,mean,sd)
                  if(s>=n) s<-n-0.01
                  e <- e+rnorm(1,mean,sd)
                  if(e<=w) e<-w+0.01
                  w <- w+rnorm(1,mean,sd)
                  if(w>=e) w<-e-0.01
                  }
          }
    return(c(n,s,e,w))
    }

#######-------------------------------------------------------------------------
## simulate an allopatric split of a rectangular range into two

alloSplit <- function(n,s,e,w){
    whichway <- sample(c("x","y"),1)
    if(whichway=="x"){
        xsplit <- runif(1, min=w, max=e)
        n1 <- n
        n2 <- n
        s1 <- s
        s2 <- s
        e1 <- xsplit
        e2 <- e
        w1 <- w
        w2 <- xsplit
        }
    if(whichway=="y"){
        ysplit <- runif(1, min=s, max=n)
        n1 <- ysplit
        n2 <- n
        s1 <- s
        s2 <- ysplit
        e1 <- e
        e2 <- e
        w1 <- w
        w2 <- w
        }
    return(data.frame(n1=n1,s1=s1,e1=e1,w1=w1,n2=n2,s2=s2,e2=e2,w2=w2))
    }

rangeSim <- function(phy=NULL,anc.bounds=c(150,50,150,50),evolmodel=c("BM","slowdown","labile"),plotit=TRUE,tips=20,mean=0,sd=0.35){
    require(ape)
    #if(is.null(phy))phy<-birthdeath.tree(b=0.05,d=0,time.stop=0, taxa.stop=tips, return.all.extinct=FALSE)
    if(is.null(phy)){stop('no phy provided')}
    edges <- data.frame("anc"=phy$edge[,1],"dec"=phy$edge[,2],"btime"=NA,"edgelength"=phy$edge.length,ancNth=NA,ancSth=NA,ancEast=NA,ancWest=NA,decNth=NA,decSth=NA,decEast=NA,decWest=NA)
    if(is.null(phy$nodel.label)) nl <- c((Ntip(phy)+1):(Ntip(phy)+Nnode(phy))) else{
        nl <- phy$node.label}
    # set ancestral range boundaries:
    curr.edges <- which(edges$anc==nl[1])
    edges[curr.edges[1],5:8]<-anc.bounds
    edges[curr.edges[2],5:8]<-anc.bounds
    switch(evolmodel,
    "BM"={
        spl <- alloSplit(edges$ancNth[curr.edges[1]],edges$ancSth[curr.edges[1]],edges$ancEast[curr.edges[1]],edges$ancWest[curr.edges[1]])
        edges[curr.edges[1],c(9:12)] <- bmSim(blen=edges$edgelength[curr.edges[1]],mean=mean,sd=sd,n=spl[1],s=spl[2],e=spl[3],w=spl[4])
        edges[curr.edges[2],c(9:12)] <- bmSim(blen=edges$edgelength[curr.edges[2]],mean=mean,sd=sd,n=spl[5],s=spl[6],e=spl[7],w=spl[8])
        for(i in 1:length(nl)){
            if(!is.na(edges$ancNth[which(edges$anc==nl[i])])) next
            curr.edges <- which(edges$anc==nl[i])
            edges[edges$anc==nl[i],c(5:8)] <- edges[edges$dec==nl[i],c(9:12)]
            z <- edges[edges$anc==nl[i],c(5:8)]
            spl <- alloSplit(z[1,1],z[1,2],z[1,3],z[1,4])
            edges[curr.edges[1],c(9:12)] <- bmSim(blen=edges$edgelength[curr.edges[1]],mean=mean,sd=sd,n=spl[1],s=spl[2],e=spl[3],w=spl[4])
            edges[curr.edges[2],c(9:12)] <- bmSim(blen=edges$edgelength[curr.edges[2]],mean=mean,sd=sd,n=spl[5],s=spl[6],e=spl[7],w=spl[8])
            }
    },
    "slowdown"={
        bt <- round(ltt(phy,plot=F)$times,3)
        bt<-data.frame("nd"=names(bt),bt)
        edges$btime<-bt$bt[match(edges$anc,bt$nd)]
        spl <- alloSplit(edges$ancNth[curr.edges[1]],edges$ancSth[curr.edges[1]],edges$ancEast[curr.edges[1]],edges$ancWest[curr.edges[1]])
        edges[curr.edges[1],c(9:12)] <- slowdownSim(blen=edges$edgelength[curr.edges[1]],mean=mean,sd=sd,starttime=edges$btime[curr.edges[1]],maxage=max(bt$bt),n=spl[1],s=spl[2],e=spl[3],w=spl[4])
        edges[curr.edges[2],c(9:12)] <- slowdownSim(blen=edges$edgelength[curr.edges[2]],mean=mean,sd=sd,starttime=edges$btime[curr.edges[2]],maxage=max(bt$bt),n=spl[5],s=spl[6],e=spl[7],w=spl[8])
        for(i in 1:length(nl)){
            if(!is.na(edges$ancNth[which(edges$anc==nl[i])])) next
            curr.edges <- which(edges$anc==nl[i])
            edges[edges$anc==nl[i],c(5:8)] <- edges[edges$dec==nl[i],c(9:12)]
            z <- edges[edges$anc==nl[i],c(5:8)]
            spl <- alloSplit(z[1,1],z[1,2],z[1,3],z[1,4])
            edges[curr.edges[1],c(9:12)] <- slowdownSim(blen=edges$edgelength[curr.edges[1]],mean=mean,sd=sd,starttime=edges$btime[curr.edges[1]],maxage=max(bt$bt),n=spl[1],s=spl[2],e=spl[3],w=spl[4])
            edges[curr.edges[2],c(9:12)] <- slowdownSim(blen=edges$edgelength[curr.edges[2]],mean=mean,sd=sd,starttime=edges$btime[curr.edges[2]],maxage=max(bt$bt),n=spl[5],s=spl[6],e=spl[7],w=spl[8])
            }
    },
    "labile"={
        spl <- alloSplit(edges$ancNth[curr.edges[1]],edges$ancSth[curr.edges[1]],edges$ancEast[curr.edges[1]],edges$ancWest[curr.edges[1]])
        edges[curr.edges[1],c(9:12)] <- labileSim(n=spl[1],s=spl[2],e=spl[3],w=spl[4])
        edges[curr.edges[2],c(9:12)] <- labileSim(n=spl[5],s=spl[6],e=spl[7],w=spl[8])
        for(i in 1:length(nl)){
            if(!is.na(edges$ancNth[which(edges$anc==nl[i])])) next
            curr.edges <- which(edges$anc==nl[i])
            edges[edges$anc==nl[i],c(5:8)] <- edges[edges$dec==nl[i],c(9:12)]
            z <- edges[edges$anc==nl[i],c(5:8)]
            spl <- alloSplit(z[1,1],z[1,2],z[1,3],z[1,4])
            edges[curr.edges[1],c(9:12)] <- labileSim(n=spl[1],s=spl[2],e=spl[3],w=spl[4])
            edges[curr.edges[2],c(9:12)] <- labileSim(n=spl[5],s=spl[6],e=spl[7],w=spl[8])
            }
    },
    )
    tipranges <- edges[edges$dec<=Ntip(phy),c(2,9:12)]
    names(tipranges)[c(2:5)] <- c("n","s","e","w")
    tipranges$species <- phy$tip.label
    if(tipranges$n-tipranges$s == 0 && tipranges$e-tipranges$w>0) tipranges$n<-tipranges$n+0.01    #ensures a long skinny range doesn't get a range size of zero
    if(tipranges$e-tipranges$w == 0 && tipranges$n-tipranges$s>0) tipranges$e<-tipranges$e+0.01
    tipranges$cent.lat <- (tipranges$n+tipranges$s)/2
    tipranges$cent.long <- (tipranges$e+tipranges$w)/2
    tipranges$rangesize <- round((tipranges$n-tipranges$s)*(tipranges$e-tipranges$w),3)
    # prune descendants of extinct lineages:
    tipstodrop <- tipranges$species[tipranges$rangesize==0]
    tipranges <- tipranges[!is.element(tipranges$species,tipstodrop),]
    if(Ntip(phy)-length(tipstodrop)<3){   # avoids error being returned if phy ends up with 2 tips or less
        phy<-NULL
        tipranges<-NULL
        }else{
            if(length(tipstodrop)>0) phy<- drop.tip(phy,tipstodrop)
            if(length(tipstodrop)>0) phy <- phy
            }
    if(plotit){
        plotRange(anc.bounds[1],anc.bounds[2],anc.bounds[3],anc.bounds[4],col="red",add=TRUE)
        for(i in 1:length(tipranges$dec)){
            plotRange(tipranges[i,2],tipranges[i,3],tipranges[i,4],tipranges[i,5],col="black",add=TRUE)
            }
        }
    return(list(phy=phy,tipranges=tipranges))
    }

###################
## NEW FUNCTIONS ##
###################

make.layer <- function(r, bounds=c(1, 200, 1, 200)){
    # By modifying the range parameter in the variogram model it is possible to
    # control the degree of spatial correlation. For example, by setting it at
    # 15 instead of 5 we get a random field with a coarser autocorrelation.
    xy <- expand.grid(bounds[1]:bounds[2], bounds[3]:bounds[4])
    names(xy) <- c("x","y")
    g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                     model=vgm(psill=0.025,model="Exp",range=r), nmax=20)
    y2y <- predict(g.dummy, newdata=xy, nsim=1)
    gridded(y2y) = ~x+y
    return(y2y)
}

boxtract <- function(nsew, enviroray){
    # Takes output from rangeSim (nsew) and an array of simulated environmental
    # data. It uses the range boundaries for each species (as described in nsew)
    # to infer the niche breadth of each of them.
    dimality <- dim(enviroray)[3]
    if(!is.null(dimnames(enviroray)[[3]])){
        lesnoms <- dimnames(enviroray)[[3]]} else {
            lesnoms <- paste('Var', 1:dimality, spe='')}
    env <- as.data.frame(matrix(nrow=nrow(nsew), ncol=dimality))
    colnames(env) <- lesnoms
    rangesize <- numeric(nrow(nsew))
    for(i in 1:nrow(nsew)){
        #alternative: almost cheating, grabs anything the box is partially over
        #also eliminates NAs:
        #bounds <- as.integer(c(ceiling(nsew[i, 2]), floor(nsew[1,3]),
        #                    ceiling(nsew[i, 4]), floor(nsew[1,5])))
        bounds <- as.integer(round(nsew[i, 2:5]))
        names(bounds) <- c('n','s','e','w')
        #if(any(c(bounds['n']-bounds['s'], bounds['e']-bounds['w']) <= 0)){
        #    next}
        #print(bounds)
        raw <- enviroray[bounds['e']:bounds['w'], bounds['s']:bounds['n'],
                         1:dimality, drop=F]
        rangesize[i] <- length(raw[,,1])
        raw2 <- as.data.frame(array(raw, dim=c(length(raw[,,1]), dimality)))
        env[i, ] <- lapply(raw2, function(x){sd(x, na.rm=T)})
        
        #names(env)[i] <- nsew[i, 'species']
    }
    env <- cbind(rangesize, env)
    #colnames(env)[1] <- 'rangesize'
}

sim.closure <- function(sig=0.7, ear, outfolder, ntip=150){
    # The simulation closure takes a set of parameters that will define the
    # simulation. It can then be executed multiple time to generate repetitions
    # with constant parameters. sig determines the rate at which the boundaries
    # move. ear is an array of simulated environmental data. outfolder is a
    # location where the output files will be stored. ntip is the max number of
    # tips to be generated on the yule phylogeny. However, over the course of
    # the simulation species do go extinct so not all tips will be present at
    # the conclusion of the simulation. higher sigma values cause more
    # extinctions
    sig <- sig
    ear <- ear
    outfolder <- outfolder
    ntip <- ntip
    return(function(jf){
        phy <- pbtree(n=ntip)
        raout <- rangeSim(phy, sd=sig, evolmodel='BM', plotit=F)
        env <- boxtract(raout[[2]], ear)
        rownames(env) <- raout[[2]]$species
        write.csv(env, paste(outfolder, jf, '.csv', sep=''))
        write.tree(raout[[1]], paste(outfolder, jf, '_', ntip, '.tre', sep=''))
          })
}

product <- function(v){
    out <- 1
    for(x in v){out <- out*x}
    return(out)
}

############
## METHOD ##
############

# This enviroray object is a 3D array 200 x 200 x (n grids) it can then be
# sliced to give a species' niche hypervolume.
enviroray <- array(c(make.layer(1)[[1]], make.layer(1)[[1]],
                     make.layer(1)[[1]],  make.layer(1)[[1]],
                     make.layer(5)[[1]],  make.layer(5)[[1]],
                     make.layer(5)[[1]],  make.layer(5)[[1]],
                     make.layer(25)[[1]], make.layer(25)[[1]],
                     make.layer(25)[[1]], make.layer(25)[[1]],
                     make.layer(75)[[1]], make.layer(75)[[1]],
                     make.layer(75)[[1]], make.layer(75)[[1]]),
                   dim=c(200, 200, 16),
                   dimnames=list(NULL, NULL, c(paste('l', 1:4, 'r1', sep=''),
                       paste('l', 1:4, 'r5', sep=''),
                       paste('l', 1:4, 'r25', sep=''),
                       paste('l', 1:4, 'r75', sep='')))
                )

spqr <- list(sim.closure(0.3, enviroray, 'output/sig_0.3/'),
             sim.closure(0.7, enviroray, 'output/sig_0.7/'),
             sim.closure(2.0, enviroray, 'output/sig_2.0/'))

# Run the simulation with all diferent parameters: 
for(Fun in spqr){mclapply(1:1000, Fun, mc.preschedule=F, mc.cores=4)}

# 
rownames(j) <- paste('t', rownames(boxes[[1]]), sep='')
j[rownames(j) %in% phy$tip, ]
phy2 <- drop.tip(phy, phy$tip[!(phy$tip %in% rownames(j))])

system.time({
#j <- rangeSim(mcc.post[[1]], evolmodel='BM', sd=0.7)
j <- rangeSim(pbtree(n=100, ape=T), evolmodel='BM', sd=3)
z <- j[[2]][[9]]
length(which(j[[1]]$tip %in% j[[2]]$species))
names(z) <- j[[2]]$species
phylosig(j[[1]], z, method='lambda')
})

# Plot an overlay of an environmental grid and some box ranges:
i <- 1

i <- i+1
png(paste(i, 'box_null.png', sep='_'))
#plot(l, axes=F)
l <- make.layer(35)
plot(raster(l), add=F, axes=T, legend=F, col=terrain.colors(250), width=480*2)
j <- rangeSim(phy, evolmodel='BM', sd=2.5)
dev.off()


# Plot the array of environmental data: 
png('enviroray.png', width=3000, height=3000)
par(mfrow=c(4,4), mar=c(0,0,0,0))
for(i in c(seq(1, 16, 4), seq(2, 16, 4), seq(3, 16, 4), seq(4, 16, 4))){
plot(raster(enviroray[ , ,i]), mar=c(0,0,0,0), legend=F, axes=F, box=F)
}
dev.off()
