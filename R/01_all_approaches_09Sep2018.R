######################################################################
# Procrustes/nonlinear gradient approach
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 09 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

require(ade4)
require(ecole)
require(vegan)
require(mgcv)

### data load
load('./data/veg.rda')
# load('./data/invert.rda')
d   <- veg
xy  <- d$xy
spe <- d$spe
env <- d$env
tra <- d$tra
phy <- d$phy

### split species/traits into two (to simulate 'plant vs lichen')
###    (break is at deepest phylogenetic node)
i1   <- 1:22
i2   <- 23:56
spe1 <- mx_dropzero(spe[,i1])
spe2 <- mx_dropzero(spe[,i2])
tra1 <- tra[rownames(tra) %in% names(spe1),]
tra2 <- tra[rownames(tra) %in% names(spe2),]

### zero-adjustment for empty SUs
spe1 <- cbind(spe1, dummy=rep(1, nrow(spe1)))
spe2 <- cbind(spe2, dummy=rep(1, nrow(spe2)))
tra1 <- rbind(tra1, dummy=colMeans(tra))
tra2 <- rbind(tra2, dummy=colMeans(tra))

### transformations
spe1 <- data.frame(genlogtrans(spe1))
spe2 <- data.frame(genlogtrans(spe2))
env  <- data.frame(decostand(scale(env, center=F), 'range'))
tra1 <- data.frame(decostand(tra1, 'range'))
tra2 <- data.frame(decostand(tra2, 'range'))
cwm1 <- data.frame(makecwm(spe1, tra1))
cwm2 <- data.frame(makecwm(spe2, tra2))

### NMDS ordinations
`ordfn` <- function(comm, distance, k, ...){
     metaMDS(comm=comm, distance=distance, k=k, try=100, trymax=200,
             autotransform=F, trace=0, maxit=500, weakties=TRUE)
}
m1 <- ordfn(spe1, 'bray', 2)
m2 <- ordfn(spe2, 'bray', 2)
m3 <- ordfn(cwm1, 'altGower', 2)
m4 <- ordfn(cwm2, 'altGower', 2)
set_par(4)
plot(m1)
plot(m2)
plot(m3)
plot(m4)

### procrustes
p1 <- protest(m1,m2) # plant vs lichen species
p2 <- protest(m3,m4) # plant vs lichen CWM traits
p1$t0  # plant vs lichen species
p2$t0  # plant vs lichen CWM traits

### fit GAMs to regress each enviro variable on NMDS scores
`gamfit` <- function(ord, env, ...){
     ndim <- ord$ndim
     scr  <- scores(ord)
     nm   <- dimnames(env)[[2]]
     nenv <- length(nm)
     stopifnot(identical(rownames(scr), rownames(env)))
     xn <- rep(NA, ndim)
     ml <- vector('list', ndim)
     st <- matrix(NA, nrow=nenv, ncol=ndim+1)
     dimnames(st)[[1]] <- nm
     dimnames(st)[[2]] <- c(dimnames(scr)[[2]],'Adj_R2')
     for(i in 1:ndim){
          xn[i] <- paste0('s(scr[,',i,'])')
     }
     right <- paste(xn, collapse='+')
     for(i in 1:nenv){
          left <- paste0('env[,',i,'] ~ ')
          fmla <- as.formula(paste(left, right))
          ml[[i]] <- mgcv::gam(fmla)
          ss      <- summary(ml[[i]])
          st[i,1:ndim] <- as.numeric(sprintf('%.3f',round(ss$s.pv,3)))
          st[i,ndim+1] <- as.numeric(sprintf('%.3f',round(ss$r.sq,3)))
     }
     out <- list(mods=ml, sumtab=st)
     class(out) <- 'gamfit'
     out
}
print.gamfit <- function(x, ...){
     print(x[['sumtab']])
}
summary.gamfit <- function(object, ...){
     object[['sumtab']]
}
# envfit(m1, env, perm=999)    # LINEAR correlations
gamfit(m1, env)  # lichen species fit to env
gamfit(m2, env)  #  plant species fit to env
gamfit(m3, env)  # lichen traits CWM fit to env
gamfit(m4, env)  #  plant traits CWM fit to env
### TODO: need row/col permutation tests...

### fourthcorner
(f1 <- fourthcorner(env, spe1, tra1, nrepet=999, modeltype=6))
(f2 <- fourthcorner(env, spe2, tra2, nrepet=999, modeltype=6))

### RLQ and fourthcorner.rlq
`rlqfn` <- function(spe, env, tra, ndim=2, ...){
     ca     <-       dudi.coa(spe, scannf=F, nf=ndim)
     pc_env <- dudi.hillsmith(env, scannf=F, nf=ndim, row.w=ca$lw)
     pc_tra <- dudi.hillsmith(tra, scannf=F, nf=ndim, row.w=ca$cw)
     rlqres <- rlq(pc_env, ca, pc_tra, scannf=F, nf=ndim)
     out <- list(rlqres=rlqres, ca=ca, pc_env=pc_env, pc_tra=pc_tra)
     class(out) <- 'rlqres'
     out
}
`fc_envfit` <- function(x, ...){
     if(!class(x)[1] == 'rlqres') stop('must be class `rlqres`')
     rlqres <- x[[1]]
     ca     <- x[[2]]
     pc_env <- x[[3]]
     pc_tra <- x[[4]]
     fc   <- fourthcorner.rlq(rlqres, type='R.axes')
     nm   <- fc$colnames.R
     nenv <- length(nm)
     tt   <- cbind(r_obsvd  = fc$tabD2$obs,
                   adj_pval = fc$tabD2$adj.pvalue)
     out  <- cbind(tt[1:nenv,], tt[(nenv+1):(nenv*2),])
     dimnames(out)[[1]] <- nm
     out
}
r1 <- rlqfn(spe1, env, tra1)
r2 <- rlqfn(spe2, env, tra2)
fc_envfit(r1)
fc_envfit(r2)

###   END   ########################################################
