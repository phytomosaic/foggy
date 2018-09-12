######################################################################
# Procrustes/nonlinear gradient approach
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 09 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

require(ade4)
require(ecole)
require(vegan)
require(mgcv)

### data load
load('./data/veg.rda', verbose=T)
# load('./data/invert.rda')
# d   <- veg
# xy  <- d$xy
# spe <- d$spe
# env <- d$env
# tra <- d$tra
# phy <- d$phy

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

### plot species WA and range
`plot_wa` <- function(spe, mod, xexp=1.6, yexp=1, pick=NULL, lcol, ...){
     tmp <- spe
     tmp[tmp>0]  <- 1
     tmp[tmp<=0] <- NA
     max1 <- sapply(tmp * mod$points[,1], max, na.rm=TRUE)
     min1 <- sapply(tmp * mod$points[,1], min, na.rm=TRUE)
     max2 <- sapply(tmp * mod$points[,2], max, na.rm=TRUE)
     min2 <- sapply(tmp * mod$points[,2], min, na.rm=TRUE)
     cent <- wascores(mod$points, spe1, expand=F)
     c1   <- cent[,1]
     c2   <- cent[,2]
     if(missing(lcol)) lcol <- '#00000050'
     if(is.null(pick)){
          buff <- 0.1
          plot(c1, c2, xlim=c(min(min1),max(max1))*(1+buff),
               ylim=c(min(min2),max(max2))*(1+buff),
               xlab='NMDS1', ylab='NMDS2', ...)
          for(i in 1:NROW(c1) ){
               lines(c(min1[i],max1[i]), rep(c2[i],ea=2), col=lcol)
               lines(rep(c1[i],ea=2), c(min2[i],max2[i]), col=lcol)
          }
     } else {
          ## single axis
          xlab <- c('NMDS1','NMDS2')[pick]
          cc   <- list(c1,c2)[[pick]]
          minn <- list(min1,min2)[[pick]]
          maxx <- list(max1,max2)[[pick]]
          o <- order(cc)
          op <- par(mfrow = c(1, 1),
                    mar = c(4.1 * yexp, 5 * xexp, 0, 0) + 0.1,
                    oma = c(0, 0, 0, 0), font = 2)
          plot(cc[o], 1:length(cc),
               xlim=c(min(minn),max(maxx))*(1+buff),
               xlab=xlab, ylab='', axes=F, ...)
          for (i in 1:length(cc)) {
               lines(c(minn[o][i],maxx[o][i]), rep(i,ea=2), col=lcol)
          }
          axis(1, at=pretty(cc), label=pretty(cc), las=3, cex.axis=0.6)
          axis(2, at=1:length(cc), label=names(cc)[o], las=1,
               cex.axis=0.6, tick=F)
     }
}
set_par(1)
plot_wa(spe1, m1, cex=0.7)
plot_wa(spe1, m1, pick=1, cex=0.7)
plot_wa(spe1, m1, pick=2, cex=0.7)
plot_wa(cwm1, m3, pick=1, cex=0.7)
plot_wa(cwm1, m3, pick=2, cex=0.7)
plot_wa(cwm2, m4, pick=1, cex=0.7)
plot_wa(cwm2, m4, pick=2, cex=0.7)


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
     # for(i in 1:ndim){
     #      xn[i] <- paste0('s(scr[,',i,'])')
     # }
     # right <- paste(xn, collapse='+')
     for(i in 1:ndim){
          xn[i] <- paste0('scr[,',i,']')
     }
     right <- paste0('s(',paste(xn,collapse=','),')')
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

g1 <- gamfit(m1, env)  # lichen species fit to env

g1
summary(g1$mods[[6]])

gam.check(g1$mods[[6]])


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



###   GAM sensitivity analysis   ###################################
###      add noise to each NMDS axis in turn, calc MAE
x0   <- gamfit() # baseline iteration
nrep <- 100    # average the sensitivity over 100 iterations
`f` <- function(perc=5.0, ...){ # convenience function
     xx   <- gamfit()
     # xx <- litetvi(spe, noise(id$mwmt,perc=perc), m_mwmt, na.rm=T)
     c(mae(x0[,1],xx[,1], stdz=T),
       mae(x0[,2],xx[,2], stdz=T)
       # mae(x0[,3],xx[,3], stdz=T)
     )
}
# TIME WARN, 5 min for 100 reps, about 3 s per rep
Q <- matrix(nrow=nrep,ncol=3) # initiate sensitivity matrix
for(i in 1:nrep){
     cat('round',i,'of',nrep,'\n')
     Q[i,] <- f(perc=5.0)     # add 5% uncertainty, get Q values
}
colMeans(Q, na.rm=T)          # summary of mean sensitivity at 5%
# TIME WARN, about 3 s per rep
Q <- matrix(nrow=nrep,ncol=3) # initiate sensitivity matrix
for(i in 1:nrep){
     cat('round',i,'of',nrep,'\n')
     Q[i,] <- f(perc=10.0)    # add 10% uncertainty, get Q values
}
colMeans(Q, na.rm=T)          # summary of mean sensitivity at 10%
#
rm(Q, f, nrep, noise)
###   end sensitivity analysis   #####################################












###   END   ########################################################
