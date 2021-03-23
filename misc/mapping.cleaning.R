
library(rgdal)
library(sp)
library(raster)
library(rgeos)

aust.basemap <- readOGR("./data/aus/", "Australia")
plot(aust.basemap)

banksia.data <- read.csv('data/Banksia.data.csv', row.names=1)

# Map each species individually:
#dir.create('./maps/')
map.data.by.species(banksia.data)

spatial.uncertainty <- function(certainty, data){
    data <- data[!is.na(data[[24]]),]
    return(data[data$Coordinate.Uncertainty.in.Metres...parsed <= certainty,])
}

# Takes a vector of record IDs and removes them from a dataset.
IDkill <- function(data, IDs){
    return(data[!(data[[1]] %in% IDs),])}

## a function to kill points of the same species and same coordinates.
## assumes data is sorted by species then location
find.successive.duplicates <- function(data){
    out <- character()
    for(i in 2:nrow(data)){
        a <- data[i-1,]
        b <- data[i,]
        # if species, latitude and longitude match:
        if(all(c(a[[22]], a[[23]]) == c(b[[22]], b[[23]]))
           & a[[15]] == b[[15]]){
            out <- append(out, as.character(b[[1]]))} #note b's ID
    }
    return(out)
}

verify.data <- function(data, fname){
    write.table(data[,c('Species...matched', 'Latitude...processed',
                        'Longitude...processed')], fname, quote=FALSE,
                sep='\t', row.names=data[[1]])
}

uzrgh <- function(speciesnum){
    #Hidden function called by cleanup
    spec <- extract.points(speciesnum, scdata) ###s3
    mapbanks(speciesnum, scdata)
    ANS <- readline("All good?\n\t")
    if(ANS == 'y') return(character())
    hamstermap(spec)
    print('Hamstermap.txt')
    bad <- numeric()
    ANS <- readline('which one?')
    while(ANS != 'x'){
        bad <- append(bad, ANS)
        ANS <- readline('which one?')
    return(spec[[1]][bad])
    }
}

# extract set of points for a specified species:
extract.points <- function(i, data){
    if(is.numeric(i)){
        return(data[data$Species...matched==unique(data$Species...matched)[i],])
    }
    else{
        return(data[data$Species...matched == i,])}
}

#unique(extract.points(7, alldata)$Species...matched)
#unique(extract.points('Banksia marginata', alldata)$Species...matched)
#extract.points(1, alldata)[[23]]

spatial.convert <- function(data){
    return(SpatialPoints(cbind(data[2],data[3])))}

#plot(aust.basemap); plot(spatial.convert(alldata), add=T)

mapbanks <- function(data, longcol=2, latcol=3){
    #hamstermap(data)
    #par(mar=c(2,2,2,2))
    spd = SpatialPoints(data[ , c(longcol, latcol)])
    if(library('rgeos', logical.return=TRUE)){
        require('rgeos')
        bbox <-spd@bbox
        #print(bbox)
        longd <- abs(bbox[1, 2] - bbox[1, 1])
        latid <- abs(bbox[2, 2] - bbox[2, 1])
        if(latid > longd){
            grow <- (latid - longd + 4)/2
            bbox <- bbox + matrix(c(-grow, -2, grow, 2), ncol=2)
        }else{if (longd > latid){
            grow <- (longd - latid + 4)/2
            bbox <- bbox + matrix(c(-2, -grow, 2, grow), ncol=2)}
        else{bbox <- bbox + matrix(c(-2, -2, 2, 2), ncol=2)}}
        bm <- gClip(aust.basemap, bbox)
    }else{bm <- aust.basemap}
    plot(bm, col='wheat', bg='lightblue', axes=T, border='black',
         las=1, lwd=0.2)
    plot(spd, add=T, axes=F, cex=0.5, pch=3, col='maroon')
    #if(is.numeric(i)){title(unique(data$Species...matched)[i])}
    #else title(i)
    title(data$Species[1])
}

# copy hamstermap.txt to www.hamstermap.com/custommap.html
# point number will correspond to cleaning.txt
hamstermap <- function(data){
    coords = data[[2:3]]
    l = length(rownames(coords))
    map <- cbind(coords, rep('numbered', l),
                 rep('black', l), 1:l)
    write.table(map, 'hamstermap.txt', row.names=F, col.names=F,
                sep='\t', quote=F)
}

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))),
              "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

## Converts numbers gathered from hamster map with uzrgh into record IDs that
## IDkill can then use to clean the data.
zzz <- function(spec.num, points){
    #Hidden function called by cleanup
    spec <- extract.points(spec.num, scdata)
    plot(aust.basemap)
    plot(spatial.convert(spec[points,]), add=TRUE, col='blue')
    title(spec[1,15])
    return(as.character(spec[[1]])[points]) #IDs
}

cleanup <- function(range){
    #interactive prompt to eliminate outliers
    out = character()
    species = unique(scdata[[15]])
    for(i in range){
        mapbanks(extract.points(i, scdata))
        prompt = paste(i, ': ', species[i], '\n\t', sep='')
        points <- as.integer(strsplit(readline(prompt), " ")[[1]])
        if(length(points) > 0){
            ids <- zzz(i, points)
            out <- append(out, ids)
        }
    }
    return(out)
}

particular.point <- function(i, p) {
    plot(aust.basemap)
    plot(banksia_records[banksia_records$binomial==unique(banksia_records$binomial)[i],][p,], axes=F,add=TRUE,pch=3,cex=0.5,col="red")
    title(paste(unique(banksia_records$binomial)[i], p))
    }

map.data.by.species<- function(data){
    for(i in 1:length(unique(data$Species))){
        tryCatch({
            svg(paste('maps/', i, '_', unique(data$Species)[i], '.svg'))
            mapbanks(extract.points(i, data))
            dev.off()},
                 error=function(e){print(paste(i, " Didn't work."))}
                 )
    }
}


##########################
## DATA CLEANING METHOD ##
##########################

alldata <- read.csv('banksia_before_cleaning.csv')
#alldata <- alldata[!is.na(alldata$Coordinate.Uncertainty.in.Metres...parsed),]

osc = length(unique(alldata$Species))

# remove data with missing lat or long (not actually needed)
alldata <- alldata[!is.na(alldata[[22]]),]
alldata <- alldata[!is.na(alldata[[23]]),]

# Finds 180 from all data all of which are in the ocean.
#sus <- alldata[alldata$Habitat.incorrect.for.species == 'true',]
alldata <- alldata[!alldata$Habitat.incorrect.for.species == 'true',]

# Finds 149 in all data and they look to be legitimate invasives.
# e.g. marginata's WA range.
# Also matches many points from adelaide botanic gardens.
#sus <- alldata[alldata$occ.Cultivated.Escapee == 'true',]
alldata <- alldata[!alldata$occ.Cultivated.Escapee == 'true',]

#sus <- alldata[alldata$Coordinates.dont.match.supplied.state == 'true',]
alldata <- alldata[!alldata$Coordinates.dont.match.supplied.state == 'true',]

# Sort data by species then latitude then longidude
sdata <- alldata[order(alldata$Species, alldata[[22]], alldata[[23]]),]

# Remove duplicate points:
sdata <- IDkill(sdata, find.successive.duplicates(sdata))

# Refine to points declared as known to less than 10km:
# eliminates ~20 species with 15000km -not run
scdata <- spatial.uncertainty(1000000, sdata)

# Before manually removing remaining bad points check that there is still
# adequate data left to work with:
print(paste('Original species count:', osc,
            'Current species count:', length(unique(scdata$Species))))

# Then use the functions above to go through by eye eliminating bad points
# such as cultivated / escaped occurences or those that are far from the
# rest of the range so that they likely represent a misidentified specimen
# or incorrect coordinates. I removed points that were more than ~100m of the
# coast but left points that didn't exceedingly agree with expert range.
# I tried to disregard the environment when deciding whether to cut or keep
# a point.


