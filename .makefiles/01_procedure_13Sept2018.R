######################################################################
# Procedure: phylogenetically-corrected trait syndromes in relation to
#     environment; compare to fourthcorner and RLQ
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 13 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

devtools::install_github('phytomosaic/foggy')
require(foggy)

### data load
data(veg)
spe <- veg$spe
env <- veg$env
tra <- veg$tra
phy <- veg$phy

### transformations
spe <- data.frame(genlogtrans(spe))
env <- data.frame(decostand(scale(env, center=F), 'range'))
tra <- data.frame(decostand(tra, 'range'))
ptra <- phylo_corr(phy, tra) # phylogenetic correction of traits

### make CWM matrices
cwm  <- data.frame(makecwm(spe, tra)) # usual
pcwm <- data.frame(makecwm(spe, ptra))# phylo-corr

### check effects of phylogenetic correction
ecole::set_par(NCOL(tra))
for(i in 1:NCOL(tra)){
     plot(tra[,i], ptra[,i], pch=16, cex=0.7, col='#00000050',
          xlab=dimnames(tra)[[2]][i], ylab=dimnames(ptra)[[2]][i])
}
for(i in 1:NCOL(cwm)){
     plot(cwm[,i], pcwm[,i], pch=16, cex=0.7, col='#00000050',
          xlab=dimnames(cwm)[[2]][i], ylab=dimnames(pcwm)[[2]][i])
}

### first select dimensionality of NMDS ordinations
screeNMDS <- function(spe, distance, kk=5, ...) {
     stress <- rep(NA, kk)
     for (i in 1:kk) {
          cat(i, 'of', kk, 'dimensions\n')
          m <- ordfn(spe=spe, distance=distance, k=i, ...)
          stress[i] <- m$stress
     }
     plot(1:kk, stress, main='', xlab='Dimension', ylab='Stress',
          ylim=c(0, max(stress)*1.05), pch=16, las=1, bty='l')
     lines(1:kk, stress)
     abline(0.20, 0, col='red', lty = 2)
     data.matrix(stress)
}
screeNMDS(spe, distance='bray', kk=5)         # 2 dimensions
screeNMDS(cwm, distance='altGower', kk=5)     # 2 dimensions
screeNMDS(pcwm, distance='altGower', kk=5)    # 2 dimensions

### NMDS ordinations
m1 <- ordfn(spe, 'bray', 2)
m2 <- ordfn(cwm, 'altGower', 2)
m3 <- ordfn(pcwm,'altGower', 2)

### GAM fit of enviro to NMDS
?gamfit
# envfit(m1, env, perm=999)   # for LINEAR correlations
g1 <- gamfit(m1, env)         # lichen species fit to env
g2 <- gamfit(m2, env)         # lichen traits CWM fit to env
g3 <- gamfit(m3, env)         # lichen traits PCWM fit to env

### plot GAM gradient regressions
set_par(6)
plot(m1)
plot(m2)
plot(m3)
plot(g1)  # species compositions have nonlinear relation to clay
plot(g2)  # but trait syndromes have nearly linear!
plot(g3)  # as do phylo-corr trait syndromes!

### plot GAM gradient regressions
set_par(12)
par(mfrow=c(4,3))
for (i in 1:4){
     plot(g1, i)
     plot(g2, i)
     plot(g3, i)
}
for (i in 5:8){
     plot(g1, i)
     plot(g2, i)
     plot(g3, i)
}
for (i in 9:11){
     plot(g1, i)
     plot(g2, i)
     plot(g3, i)
}

### TODO: relate trait syndromes to axes

### fourth-corner analysis
(f1 <- fourthcorner(env, spe, tra, nrepet=999, modeltype=6))
(f2 <- fourthcorner(env, spe, ptra, nrepet=999, modeltype=6))
plot(f1$tabD2$obs, f2$tabD2$obs); abline(a=0,b=1,h=0,v=0)

### RLQ and fourthcorner.rlq
`rlqfn` <- function(...){
     `f` <- function(spe, env, tra, ndim=2, ...){
          ca     <-       dudi.coa(spe, scannf=F, nf=ndim)
          pc_env <- dudi.hillsmith(env, scannf=F, nf=ndim, row.w=ca$lw)
          pc_tra <- dudi.hillsmith(tra, scannf=F, nf=ndim, row.w=ca$cw)
          rlqres <- rlq(pc_env, ca, pc_tra, scannf=F, nf=ndim)
          out <- list(rlqres=rlqres, ca=ca, pc_env=pc_env, pc_tra=pc_tra)
          class(out) <- 'rlqres'
          out
     }
     f(...)
}
`fc_envfit` <- function(...){
     `f` <- function(x, ...){
          if(!class(x)[1] == 'rlqres') stop('must be class `rlqres`')
          environment(fourthcorner.rlq) <- environment()
          rlqres <- x[['rlqres']]
          ca     <- x[['ca']]
          pc_env <- x[['pc_env']]
          pc_tra <- x[['pc_tra']]
          eval(fc <- fourthcorner.rlq(rlqres,type='R.axes'),
               parent.frame())
          nm   <- fc$colnames.R
          nenv <- length(nm)
          tt   <- cbind(r_obsvd  = fc$tabD2$obs,
                        adj_pval = fc$tabD2$adj.pvalue)
          out  <- cbind(tt[1:nenv,], tt[(nenv+1):(nenv*2),])
          dimnames(out)[[1]] <- nm
          out
     }
     f(...)
}
r1 <- rlqfn(spe, env, tra)
fc_envfit(x=r1)

###   END   ########################################################
