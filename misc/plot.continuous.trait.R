#! /usr/bin/R

library('caper')
library('diversitree')
library('RColorBrewer')

cols <- list(rev(grey(1:250/251)), brewer.pal(9, 'Reds'), brewer.pal(9, 'Greens'))

###########
#FUNCTIONS#
###########

ValuedColPal <- function(x, colpal){
    # x is a numeric vector
    # colpal is a vector that R can interpret as colours
    # Creates a function that retuns a colour corresponding to a value
    # within the range of 'x'.
    # sapply(x, ValuedColPal(x, some.palette(250))) will convert the data
    # of x into colours.
    s <- seq(min(x), max(x), length.out=length(colpal))
    return(function(val){
        return(colpal[which.min(abs(s-val))])
    })
}

continuous.trait <- function(cd, wcols=NULL, colpal=rev(grey(1:250/251)), ...){
    # Extends diversitree::trait.plot to plot a contiuous trait on a circular
    # phylogeny.
    # cd is a comparative.data object from caper
    # wcols = which columns of the cd object to plot.
    # colpal is a list of colours to be used as the scale.
    # Other key word arguments can be specified and will be passsed to
    # trait.plot.
    require('diversitree')
    if(is.null(wcols)) {wcols <- 1:ncol(cd$data)}
    specdata <- cd$data[ ,wcols, drop=FALSE]
    phy <- cd$phy
    dtip <- phy$tip[!(phy$tip %in% rownames(specdata))]
    phy <- drop.tip(phy, dtip)
    #rownames(specdata) <- cd$data$species
    #colpal <- switch(palette, 'rainbow'=rainbow(250))
    if(is.character(colpal)){
        incols <- lapply(specdata, function(x){
            sapply(x, ValuedColPal(x, colpal))})
    } else {if(is.list(colpal)){
        #assert lengths
        incols <- list()
        for(i in 1:ncol(specdata)){
            x <- specdata[[i]]
            VCP <- ValuedColPal(x, colpal[[i]])
            incols[[names(specdata)[i]]] <- sapply(x, VCP)
        }
        incols <- as.data.frame(incols, stringsAsFactors=FALSE)
    } else {stop("'colpal' needs to be a vector of colours or a list of
                 vectors of colours.")}
        }
    indata <- as.data.frame(matrix(rep(1:nrow(specdata),  length(wcols)),
                            ncol=length(wcols)), row.names=rownames(specdata))
    colnames(indata) <- colnames(specdata)
    print(colnames(indata))
    trait.plot(phy, indata, incols, legend=FALSE, ...)
}
