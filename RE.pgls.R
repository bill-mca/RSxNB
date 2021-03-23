library('caper')
library('phytools')
library('parallel')

infolder <- '../RE.sim/sig_2.0/'

lscf <- list.files(infolder, pattern='*.csv', full.names=T)
lstf <- list.files(infolder, pattern='*.tre', full.names=T)

product <- function(v){
    out <- 1
    for(x in v){out <- out*x}
    return(out)
}

make.hvfunc <- function(...){
    vars  <- c("l1r1", "l2r1", "l3r1", "l4r1", "l1r5", "l2r5", 
               "l3r5", "l4r5", "l1r25", "l2r25", "l3r25", "l4r25", "l1r75", 
               "l2r75", "l3r75", "l4r75")
    print(paste(vars[c(...)]))
    vars <- vars[c(...)]
    return(function(data){return(log10(product(data[,vars]*100)))})
}


#The hypervolumes i played with for the thesis:
hvfuncs <- list(
'uni75.h4v'=make.hvfunc(13:16),
#"l1r75" "l2r75" "l3r75" "l4r75"
'mix.h4v'=make.hvfunc(5, 9, 10, 13),
#"l1r5"  "l1r25" "l2r25" "l1r75"
'uni25.h4v'=make.hvfunc(9:12),
#"l1r25" "l2r25" "l3r25" "l4r25"
'mix.h6v'=make.hvfunc(1, 5, 9, 10 , 13, 14),
#"l1r1"  "l1r5"  "l1r25" "l2r25" "l1r75" "l2r75"
    'mix2.h6v'=make.hvfunc(2, 6, 11, 12 , 15, 16),
'mix.h8v'=make.hvfunc(9:16),
    'uni5.h4v'=make.hvfunc(5:8),
#    'uni25.h1v'=make.hvfunc(9),
#    'uni75.h1v'=make.hvfunc(13),
#    'uni5.h1v'=make.hvfunc(5),
    'uni1.h4v'=make.hvfunc(1:4)
    )

# I only continued with the uni (all one AC level) HVs
hvfuncs <- hvfuncs[c(1,3,7,8)]

formulae <- paste('rangesize ~ `', names(hvfuncs), '`', sep='')
formulae <- lapply(formulae, formula)
    
formulae <- c('rangesize ~ `uni75`', 'rangesize ~ `mix.h4v`',
              'rangesize ~ `uni25`',  'rangesize ~ `mix.h6v`')
formulae <- lapply(formulae, formula)

.PGLSextract <- function(pglsout){
    pgsme <- summary(pglsout)
    ce <- pgsme$coefficients
    return(c(ce[c(2,1,8,7,4,3,6,5)], pgsme$r.squared, pgsme$fstatistic[[1]]))
}

pglser <- function(infolder, formule){
    infolder <- infolder
    formule <- formule
    hvfunc <- hvfuncs[[as.character(formule[[3]])]]
    #take formula
    #make hyper volume.
    #pgls as in rp
    #ntips, lambda, extracted pgls *1000 in a csv
    return(function(n){
        data <- read.csv(paste(infolder, n, '.csv', sep=''), row.names=1)
        phy <- paste(infolder, n, '_150.tre', sep='')
        #print(phy)
        phy <- read.tree(phy)
    if(is.null(data$species)){
        data$species <- rownames(data)}
    wyew <- as.character(c(formule[[2]], formule[[3]]))
    # Discard all infinite data:
    data[[wyew[2]]] <- hvfunc(data)
    data <- data[!is.na(data[,wyew[2]]), ]
    data <- data[is.finite(data[,wyew[2]]),]
    data[['rangesize']] <- log10(data[['rangesize']])
    NTips <- nrow(data)
    #print(NTips)
    cd <- comparative.data(phy, data[, c('species', wyew)], 'species')
    phy <- drop.tip(phy, cd$drop[[1]])
    data <- cd$data
    wyewd <- cd$data[,wyew[2]]
    names(wyewd) <- rownames(cd$data)
    lambda <- phylosig(phy, wyewd, 'lambda')
    #print(lambda)
    geez <- pgls(formula(formule), cd, lambda=lambda[[1]],
                 bounds=list('lambda'=c(0,3)))
    return(c(NTips, lambda[[1]], .PGLSextract(geez)))
    })
}

parse.out <- function(lapout){
    lapout <- t(as.data.frame(lapout))
    colnames(lapout) <- outcols
    row.names(lapout) <- NULL
    return(lapout)
}


outcols <- c('tips', 'lambda', 'slope', 'intercept', 's.pvalue', 'i.pvalue',
             's.error', 'i.error', 's.tvalue', 'i.tvalue', 'rsquared',
             'fstatistic')

spqr <- pglser('sig_2.0/', formulae[[2]], mix.h4v)
cb <- mclapply(1:1000, spqr, mc.cores=4, mc.preschedule=F)
write.csv(cb, '2.0_mix.h4v_sun26.csv')

spqr <- pglser('sig_0.7/', formulae[[4]], mix.h6v)
cb <- mclapply(1:1000, spqr, mc.cores=4, mc.preschedule=F)
write.csv(cb, '0.7_mix.h6v_.csv')

spqr <- pglser('sig_0.7/', formulae[[1]], uni75.h4v)
cb <- mclapply(1:1000, spqr, mc.cores=4, mc.preschedule=F)
write.csv(cb, '0.7_uni75_tues21.csv')

system.time({
spqr <- pglser('sig_0.7/', formulae[[3]], uni25.h4v)
cb <- mclapply(1:1000, spqr, mc.cores=4, mc.preschedule=F)
write.csv(cb, '0.7_uni25_tues21.csv')})

evo <- c('0.7', '0.3') #'2.0'
for(e in evo){
    infolder <- paste('sig_', e, '/', sep='')
    for(formule in formulae){
        spqr <- pglser(infolder, formule)
        cb <- mclapply(1:1000, spqr, mc.cores=4, mc.preschedule=F)
        bad <- which(sapply(cb, class)=='try-error')
        cb[bad] <- rep(NA, 12)
        cb <- t(as.data.frame(cb))
        rownames(cb) <- NULL
        colnames(cb) <- outcols
        write.csv(cb, paste('pgls/', e, '/',
                            as.character(formule[[3]]), sep=''))
    }
}

lsf <- list.files('pgls/', full.names=T, recursive=T)
check <- function(f){
    j <- read.csv(f, row.names=1)
    return(length(which(is.na(j[])))/12)}
lapply(lsf, check)

#
