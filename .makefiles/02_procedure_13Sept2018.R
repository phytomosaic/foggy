######################################################################
# phylogenetically-corrected trait syndrome responses to environment
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 13 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)


# require(devtools)
# devtools::install_github('phytomosaic/foggy')
require(ade4)
require(ecole)
require(vegan)
require(mgcv)
require(picante)

### data load
load('./data/veg.rda', verbose=T)
# load('./data/invert.rda', verbose=T)
# load('./data/pillar.rda', verbose=T)

### transformations
spe <- data.frame(genlogtrans(spe))
env <- data.frame(decostand(scale(env, center=F), 'range'))
tra <- data.frame(decostand(tra, 'range'))

### phylogenetic correction of traits
ptra <- phylo_corr(phy, tra)

### make CWM matrices
cwm  <- data.frame(makecwm(spe, tra)) # usual
pcwm <- data.frame(makecwm(spe, ptra))# phylo-corr

### NMDS ordination function
`ordfn` <- function(comm, distance, k, ...){
     metaMDS(comm=comm, distance=distance, k=k, try=100, trymax=200,
             autotransform=F, trace=0, maxit=500, weakties=TRUE)
}

###
m1 <- ordfn(spe, 'bray', 2)
m2 <- ordfn(cwm, 'altGower', 2)
m3 <- ordfn(pcwm,'altGower', 2)
set_par(3)
plot(m1)
plot(m2)
plot(m3)

set_par(2)

plot(m3, display='si', cex=standardize(pcwm$spikiness)+1)
plot(m3, display='si', cex=standardize( cwm$spikiness)+1)


###   END   ########################################################
