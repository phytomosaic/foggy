######################################################################
# Simulate species in communities using coenoflex
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 24 Jan 2017
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

rm(list=ls())
pkg <- c('vegan', 'coenoflex', 'ecole')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=T)
rm(pkg, has)

### simulate community by rejection sampling: iterate til connected
`simfn` <- function(G=2, tol, cmp, ...){
     environment(coenoflex) <- environment()
     i <- 1
     while( i < 999 ){
          spe <- coenoflex(numgrd=G, numplt=100, numspc=100,
                           grdtyp=rep('e',G), grdlen=rep(100,G),
                           width=rep(tol,G), variab=rep(60,G),
                           grdprd=rep(0,G), alphad=rep(1.2,G),
                           pdist='g',sdist='r', skew=3.0,aacorr=0.0,
                           cmpasy=cmp, cmpphy=0.0,
                           maxtot=100, noise=1.0, slack=0.10,
                           autlin='irm(1,2)')$veg
          if( all(rowSums(spe)!=0) &
              sum(unique(distconnected(vegdist(spe),1,F)))==1){
               break
          }
          i <- i + 1
     }
     if(i>1) cat(sprintf('\n%d iters of rejection sampling',i))
     if(i==999) cat(', still disconnected, try increasing iters')
     dimnames(spe)[[1]] <- paste0('SU', seq(1:100))
     dimnames(spe)[[2]] <- paste0('sp', seq(1:100))
     spe <- log10(spe+1)
     spe
}

###  simulate data, varying both niche widths and competition:
set.seed(19)             # for reproducibility
R    <- 10               # replicates
tols <- seq(20,65,5)     # vary niche tolerances
cmps <- seq(1,19,2)      # vary competition asymmetry
grd  <- expand.grid(tols=tols, cmps=cmps)      # all combinations
egrd <- grd[rep(seq_len(nrow(grd)), each=R),]  # rep ea combn R times
spesim <- mapply(FUN=function(G,tol,cmp)simfn(G,tol,cmp),
                 tol=egrd$tols, cmp=egrd$cmps,
                 MoreArgs=(list(G=2)), SIMPLIFY='array', USE.NAMES=F)
idsim  <- expand.grid(x=seq(0,100,len=10), y=seq(0,100,len=10))
head(spesim[,,1])   # peek at first of 1000 species datasets
head(idsim)         # 2-dimensional environmental gradients

### diversity attributes of simulated datasets:
d <- t(apply(spesim, 3, mx_diversity))
r <- colvec(as.factor(egrd$tols))
dimnames(d)[[2]] <- c('gamma','alpha','beta','halfchanges',
                         'dbi','prop0','propnoshare','N')
pairs(d, pch=20, cex=0.7, col=r, lower.panel=NULL, las=2)

###   end generate simulated community data   #######################
