######################################################################
# Compare plants and lichens (simulated data)
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 18 Nov 2018
#   CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

# GITHUB_PAT <- '01591b7a42837969a5e9360bb61089d5dc1822a0'
# devtools::install_github('phytomosaic/foggy', auth_token=GITHUB_PAT)

require(foggy)

### load data
data(simdata)
spe <- simdata[[1]]$spe
env <- simdata[[1]]$env
tra <- simdata[[1]]$tra
phy <- simdata[[1]]$phy

### split species/traits into clades ('plant vs lichen')
i1   <- 1:25
i2   <- 26:50
spe1 <- mx_dropzero(spe[,i1])
spe2 <- mx_dropzero(spe[,i2])
tra1 <- tra[rownames(tra) %in% names(spe1),]
tra2 <- tra[rownames(tra) %in% names(spe2),]
names(tra2) <- c('Y1','Y2','Y3') # incommensurate traits
phy1 <- drop.tip(phy, names(spe2))
phy2 <- drop.tip(phy, names(spe1))

### visualize
set_par(2)
plot(phy1, cex=0.6, no.margin=TRUE)
plot(phy2, cex=0.6, no.margin=TRUE)

### zero-adjustment for empty SUs (if needed)
if (!all(mx_valid(spe1), mx_valid(spe2))){
     cat('adding dummy column')
     spe1 <- cbind(spe1, dummy=rep(1, nrow(spe1)))
     spe2 <- cbind(spe2, dummy=rep(1, nrow(spe2)))
}
if (!all(mx_valid(tra1), mx_valid(tra2))){
     cat('adding dummy column')
     tra1 <- rbind(tra1, dummy=colMeans(tra))
     tra2 <- rbind(tra2, dummy=colMeans(tra))
}

### transformations
spe1 <- data.frame(genlogtrans(spe1))
spe2 <- data.frame(genlogtrans(spe2))
env  <- data.frame(decostand(scale(env, center=F), 'range'))
tra1 <- data.frame(decostand(tra1, 'range'))
tra2 <- data.frame(decostand(tra2, 'range'))

### detect phylogenetic signal for traits
sapply(tra1, FUN=function(j){
     names(j) <- rownames(tra1)
     round(picante::phylosignal(j, multi2di(phy1)),3)})
sapply(tra2, FUN=function(j){
     names(j) <- rownames(tra2)
     round(picante::phylosignal(j, multi2di(phy2)),3)})

### phylogenetic correction of traits (optional, given phylo signal)
ptra1 <- phylo_corr(phy1, tra1)
ptra2 <- phylo_corr(phy2, tra2)
ptra1 <- decostand(ptra1, 'range')
ptra2 <- decostand(ptra2, 'range')

### see effects of phylogenetic correction
set_par(NCOL(tra1))
for(i in 1:NCOL(tra1)){
     plot(tra1[,i], ptra1[,i], pch=16, cex=0.8, col='#00000050',
          xlab=dimnames(tra1)[[2]][i], ylab=dimnames(ptra1)[[2]][i])
}
set_par(NCOL(tra2))
for(i in 1:NCOL(tra2)){
     plot(tra2[,i], ptra2[,i], pch=16, cex=0.8, col='#00000050',
          xlab=dimnames(tra2)[[2]][i], ylab=dimnames(ptra2)[[2]][i])
}

### community-weighted means of phylo-corrected traits (PCWM),
###     interpreted here as phylo-corrected trait syndromes
pcwm1 <- makecwm(spe1, ptra1)
pcwm2 <- makecwm(spe2, ptra2)

### NMDS ordination of phylo-corrected trait syndromes
(m1 <- ordfn(pcwm1, 'altGower', 2))
(m2 <- ordfn(pcwm2, 'altGower', 2))

### procrustes
p <- protest(m1,m2,perm=999,symm=T) # plant v lichen trait syndromes
p$t0 ; p$signif
m2$points <- p$Yrot # replace old scores with procrustes aligned ones

### plot configurations after procrustes alignment
set_par(2)
plot(m1, type='t')
plot(m2, type='t')

### GAM gradient regressions give *nonlinear* trait-environment fit
(g1 <- gamfit(m1, env))
(g2 <- gamfit(m2, env))

### plot gradient regressions for each enviro variable
set_par(NCOL(env)*2)
for (i in 1:NCOL(env)){ plot(g1, i, lcol='#FF000080', lwd=1) }
for (i in 1:NCOL(env)){ plot(g2, i, lcol='#FF000080', lwd=1) }

# H1: Trait convergence of non-analogous communities will be greater
# at the extremes of the environmental stress gradients (because
# extremes more strongly select for similar trait syndromes???).
# Convergence is measured as 1) Procrustes residuals or 2) overlap of
# env surfaces in trait syndrome space.

### H1: Procrustes disagreement (residuals) higher at env extremes?
resid <- residuals(p)
set_par(2)
plot(env$env1, resid, pch=16) ; abline(lm(resid ~ env$env1))
plot(env$env2, resid, pch=16) ; abline(lm(resid ~ env$env2))
# greater disagreement of plant vs lichen at extremes of env2 (but not env1)

### H1: enviro deviation surfaces greatest at enviro extremes?
y1 <- fitted(g1)
y2 <- fitted(g2)
dev <- decostand(y1, 'range') - decostand(y2, 'range')
(gdev <- gamfit(m1, dev))   # deviation surface
set_par(6)
for(j in 1:NCOL(env)){
     plot(g1,   j, lcol='#FF000080', lwd=1)
     plot(g2,   j, lcol='#FF000080', lwd=1)
     plot(gdev, j, lcol='#FF000080', lwd=1)
}

# H2: Suites of traits (e.g. water retention traits, photoprotection
# traits, etc) will show similar responses to environmental gradients
# in both vascular plants and lichens.

# H2: overlay gradient surfaces per (phylo-corrected) trait
gt1 <- gamfit(m1, pcwm1)
gt2 <- gamfit(m2, pcwm2)
yt1 <- fitted(gt1)
yt2 <- fitted(gt2)
devt <- decostand(yt1, 'range') - decostand(yt2, 'range')
(gdevt <- gamfit(m1, devt))   # deviation surface
set_par(NCOL(pcwm1)*3)
for (j in 1:NCOL(pcwm1)){
     plot(gt1, j, lcol='#FF000080', lwd=1)
}
for (j in 1:NCOL(pcwm1)){
     plot(gt2, j, lcol='#FF000080', lwd=1)
}
for (j in 1:NCOL(pcwm1)){
     plot(gdevt, j, lcol='#FF000080', lwd=1)
}


# H3: Resource delivery form will affect taxonomic diversity of clades
# differentially.
# Lichen species diversity increase with fog (altitude), no change rain (latitude)
# Plant species diversity no change with fog (altitude), increase with rain (latitude)
s1 <- rowSums((spe1>0)*1, na.rm=T)
s2 <- rowSums((spe2>0)*1, na.rm=T)
set_par(NCOL(env))
for (j in 1:NCOL(env)){
     plot(env[,j], s1, col='#00000070', pch=16, ylim=range(c(s1,s2)))
     abline(lm(s1~env[,j]), col='#00000070')
     points(env[,j], s2, col='#FF000070', pch=16)
     abline(lm(s2~env[,j]), col='#FF000070')
}


# H4: Resource delivery form will affect phylogenetic diversity and
# phylogenetic richness of clades differentially.
# Lichen P increase with fog (altitude), no change rain (latitude)
# Plant P no change with fog (altitude), increase with rain (latitude)
pd1 <- cbind(psd(spe1,phy1,F)) ; pairs(pd1, upper.panel=NULL)
pd2 <- cbind(psd(spe2,phy2,F)) ; pairs(pd2, upper.panel=NULL)
`plot_pd` <- function(y1, y2, x, pick=1, ...){
     y1   <- y1[,pick]
     y2   <- y2[,pick]
     yrng <- range(c(y1,y2))
     plot(x, y1, col='#00000070', pch=16, ylim=yrng, ...)
     abline(lm(y1~x), col='#00000070')
     points(x, y2, col='#FF000070', pch=16)
     abline(lm(y2~x), col='#FF000070')
}
set_par(10); par(mfcol=c(2,5))
for (jj in 1:NCOL(pd1)){
     for (j in 1:NCOL(env)){
          plot_pd(pd1, pd2, env[,j], pick=jj,
                  xlab=dimnames(env)[[2]][j],
                  ylab=dimnames(pd1)[[2]][jj])
     }
}


### TODO:
###   further Hypotheses: verify performance against existing methods

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
### fourthcorner
(f1 <- fourthcorner(env, spe1, tra1, nrepet=999, modeltype=6))
(f2 <- fourthcorner(env, spe2, tra2, nrepet=999, modeltype=6))
## RLQ and fourthcorner.rlq
r1 <- rlqfn(spe1, env, tra1)
r2 <- rlqfn(spe2, env, tra2)
fc_envfit(r1)
fc_envfit(r2)

####   END   ######################################################
