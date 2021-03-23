#! /usr/bin/R

library('caper')
library('phytools')
library('parallel')

#############
#GLOBAL VARS#
#############

formulae <-
    c("`bufrange.cea`~`sd.h1v`", "`bufrange.cea`~`sd.h4v`",
      "`bufrange.cea`~`sd.h6v`")
formulae <- lapply(formulae, formula)

mcc.post <- read.tree('./data/mcc.post.tre')

###########
#FUNCTIONS#
###########


robust.pgls <- function(formule, phy, data){
    # Uses caper::pgls to make pgls models but feeds in a lambda value
    # calculated  by phytools::phylosig because it is faster and the caper
    # alternative occasionally throws an error with no easy solution.
    # This function also generates a caper::comparative.data object from the
    # Banksia data by applying dataset specific functions such as tips.to.keep.
    # This is so that the subspecies tip used in each iteration will be
    # randomised. Ideally it would be split into two functions.
    wyew <- as.character(c(formule[[2]], formule[[3]]))
    # Discard all infinite data:
    data <- data[is.finite(data[,wyew[2]]),]
    cd <- comparative.data(phy, data[, c('species', wyew)], 'species')
    phy <- drop.tip(phy, cd$drop[[1]])
    data <- cd$data
    wyewd <- cd$data[,wyew[2]]
    names(wyewd) <- rownames(cd$data)
    lambda <- phylosig(phy, wyewd, 'lambda')
    #print(lambda)
    return(pgls(formula(formule), cd, lambda=lambda[[1]],
                bounds=list('lambda'=c(0,3))))
}

# collapse tips to species level by randomly removing all but one subspecies:
tips.to.keep <- function(){
  rmv <- sample(c("Banksia_gardneri_var_bevidentata",
                  "Banksia_gardneri_var_gardneri",
                  "Banksia_gardneri_var_hiemalis"),1)
  rmv <- c(rmv, sample(c("Banksia_integrifolia_subsp_compar",
                         "Banksia_integrifolia_subsp_integrifolia",
                         "Banksia_integrifolia_subsp_monticola"),1))
  rmv <- c(rmv, sample(c("Banksia_laevigata_subsp_fuscolutea",
                         "Banksia_laevigata_subsp_laevigata"),1))
  rmv <- c(rmv, sample(c("Banksia_leptophylla_ ",
                         "Banksia_leptophylla_var_leptophylla",
                         "Banksia_leptophylla_var_melletica"),1))
  rmv <- c(rmv, sample(c("Banksia_meisneri_var_ascendens",
                         "Banksia_meisneri_var_meisneri"),1))
  rmv <- c(rmv, sample(c("Banksia_nivea_subsp_fuliginosa",
                         "Banksia_nivea_subsp_nivea"),1))
  rmv <- c(rmv, sample(c("Banksia_nobilis_subsp_fragrans",
                         "Banksia_nobilis_subsp_nobilis"),1))
  rmv <- c(rmv, sample(c("Banksia_nutans_var_cernuella",
                         "Banksia_nutans_var_nutans"),1))
  rmv <- c(rmv, "Banksia_serratuloides_subsp_peris")
  rmv <- c(rmv, sample(c("Banksia_sessilis_var_cordata",
                         "Banksia_sessilis_var_sessilis"),1))
  rmv <- c(rmv, sample(c("Banksia_sphaerocarpa_var_caesia",
                         "Banksia_sphaerocarpa_var_latifolia",
                         "Banksia_sphaerocarpa_var_sphaerocarpa"),1))
  rmv <- c(rmv, sample(c("Banksia_spinulosa_var_collina",
                         "Banksia_spinulosa_var_cunninghamii",
                         "Banksia_spinulosa_var_neoanglica",
                         "Banksia_spinulosa_var_spinulosa"),1))
  rmv <- c(rmv, sample(c("Banksia_splendida_subsp_macrocarpa",
                         "Banksia_splendida_subsp_splendida"),1))
  rmv <- c(rmv, sample(c("Banksia_squarrosa_subsp_argillacea",
                         "Banksia_squarrosa_subsp_squarrosa"),1))
  rmv <- c(rmv, sample(c("Banksia_subpinnatifida_var_imberbis",
                         "Banksia_subpinnatifida_var_subpinnatifida"),1))
  rmv <- c(rmv, sample(c("Banksia_tenuis_var_reptans",
                         "Banksia_tenuis_var_tenuifolia"),1))
  return(rmv)
  }

prune <- function(phy){
    # Prunes the Banksia phylogenetic tree to species level. 
    phy[[3]][grep('^Banksia_leptophylla$', phy$tip.label)] <- 'Banksia_leptophylla_ '
    phy$tip.label[grep('^Banksia_penicillata$', phy$tip.label)] <- 'Banksia_conferta'
    phy$tip.label[grep('^Banksia_integrifolia_subsp_aquilonia$', phy$tip.label)] <- 'Banksia_aquilonia'
    phy$tip.label[grep('^Banksia_proteiodes$', phy$tip.label)] <- 'Banksia_proteoides'
    phy$tip.label[grep('^Banksia_dallaneyi$', phy$tip.label)] <- 'Banksia_dallanneyi'
    tips <- phy$tip.label
    retain.o <- tips.to.keep()
    retain.n <- as.character(lapply(strsplit(tips.to.keep(), '_'),
                    function(x){return(paste(x[1:2], collapse='_'))}))
    for(i in 1:length(retain.o)){
        #print(c(tips[grep(retain.o[i], tips)], retain.n[i]))
        #tips[grep(retain.o[i], tips)] <- retain.o[i]}
        tips <- gsub(retain.o[i], retain.n[i], tips)}
    rx <- '(?<!Banksia)_.*'
    #print(tips[grep(rx, tips, perl=T)])
    phy <- drop.tip(phy, tips[grep(rx, tips, perl=T)])
    phy$tip.label <- gsub(rx, '', phy$tip.label, perl=T)
    #phy$tip.label <- tips[!grep(rx, tips, perl=T)]
    return(phy)
}

.PGLSextract <- function(pglsout){
    # A hidden function that pulls all the valuable information out of a pgls
    # model
    pgsme <- summary(pglsout)
    ce <- pgsme$coefficients
    return(c(ce[c(2,1,8,7,4,3,6,5)], pgsme$r.squared, pgsme$fstatistic[[1]]))
}

mass.pgls <- function(formule, multiphy, specdata){
    # Generates PGLS models of a single dataset for every tree in the posterior
    l <- length(multiphy)
    outfile <- data.frame(slope=rep(NA,l), intercept=rep(NA, l),
                          s.pvalue=rep(NA,l), i.pvalue=rep(NA,l),
                          s.error=rep(NA,l), i.error=rep(NA,l),
                          s.tvalue=rep(NA,l), i.tvalue=rep(NA,l),
                          rsquared=rep(NA,l), fstatistic=rep(NA,l))
    #vector('list', length(multiphy))
    for(i in 1:length(multiphy)){
        phy <- prune(multiphy[[i]])
        #print(i)
        cd <- comparative.data(phy=phy, specdata, 'species')
        corow <- pgls(formule, cd, lambda='ML')
        outfile[i, ] <- .PGLSextract(corow)
    }
    #return(out)
    return(outfile)
}

robust.mass.pgls <- function(formule, multiphy, specdata){
    # Like mass.pgls except it uses the robust method of calculating lambda
    l <- length(multiphy)
    outfile <- data.frame(slope=rep(NA,l), intercept=rep(NA, l),
                          s.pvalue=rep(NA,l), i.pvalue=rep(NA,l),
                          s.error=rep(NA,l), i.error=rep(NA,l),
                          s.tvalue=rep(NA,l), i.tvalue=rep(NA,l),
                          rsquared=rep(NA,l), fstatistic=rep(NA,l))
    #vector('list', length(multiphy))
    for(i in 1:length(multiphy)){
        phy <- prune(multiphy[[i]])
        #print(i)
        corow <- robust.pgls(formule, phy, specdata)
        outfile[i, ] <- .PGLSextract(corow)
    }
    #return(out)
    return(outfile)
}


makeapgls <- function(formule, trees){
    # A closure that wraps a formula and set of phylogenetic trees. The
    # resulting function can then be applied over a series of datasets.
    # Uses the robust method of calculating lambda.
    formule <- formule
    trees <- trees
    return(function(f){
        specdata <- read.csv(f, row.names=1)
        nd <- robust.mass.pgls(formule, trees, specdata)
        fout <- paste(as.character(formule[[3]]), '/', f, sep='')
        write.csv(nd, fout)
    })
}

### Have to be careful that the Range size data has been logged. Slopes come out
### in the thousands if this is wrong.

log.range.size <- function(f){
    nulldata <- read.csv(f, row.names=1)
    nulldata[['bufrange.cea']] <- log10(nulldata[['bufrange.cea']])
    write.csv(nulldata, f)
    return(f)
}

check.log.range.size <- function(lsf){
    # Takes a list of null output files, reads each in, and applies log10 to
    # the range size collumn if it hasn't been done already.
    l=character()
    for(f in lsf){
        nulldata <- read.csv(f, row.names=1)
        # check if it has already been logged by reference to B. marginata:
        if((nulldata[66, 3]==338389.264417969))
            l <- c(l, f)
    }
    if(length(l)==0){warn('It seems all the files are already logged or not in',
                 ' the familiar format.')} 
    return(unlist(mclapply(l, log.range.size)))
}

##########
##METHOD##
##########

# Navigate into the directory with all the raw data:
setwd('shuffle_null')

# List all the raw data files:
lsf <- list.files('./', pattern='*.csv')

# Create directories to dump all the PGLS model data into. One directory per
# hypervolume:
if(!file.exists('./sd.h1v/')){dir.create('./sd.h1v/')}
if(!file.exists('./sd.h4v/')){dir.create('./sd.h4v/')}
if(!file.exists('./sd.h6v/')){dir.create('./sd.h6v/')}

# Check that range size has been logged: 
check.log.range.size(lsf)

#6D first then 4D and 1D last:
pgs <- list(makeapgls(formulae[[3]], mcc.post),
            makeapgls(formulae[[2]], mcc.post),
            makeapgls(formulae[[1]], mcc.post))

# Actually create the model (takes ages to execute).
for(func in pgs){mclapply(lsf, func, mc.cores=4, mc.preschedule=FALSE)}

# To resume mid-way through  you can just shorten the list of files and formulae
# E.g. if it crashed mid way through the 4D hypervolume:
lsf <- list.files('.', pattern='*.csv')
lsf2 <- list.files('sd.h4v/', pattern='*.csv')
lsf <- lsf[!(lsf %in% lsf2)]
for(func in pgs[2:3]){mclapply(lsf, func, mc.cores=4, mc.preschedule=FALSE)}
