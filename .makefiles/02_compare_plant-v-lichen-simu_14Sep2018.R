######################################################################
# Compare plants and lichens (simulated data)
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 18 Nov 2018
#   CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

# GITHUB_PAT <- '01591b7a42837969a5e9360bb61089d5dc1822a0'
# devtools::install_github('phytomosaic/foggy', auth_token=GITHUB_PAT)

require(foggy)
rm(list=ls())

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

### see effects of phylogenetic correction
set_par(NCOL(tra1))
for(i in 1:NCOL(tra1)){
     plot(tra1[,i], ptra1[,i], pch=16, cex=0.7, col='#00000050',
          xlab=dimnames(tra1)[[2]][i], ylab=dimnames(ptra1)[[2]][i])
}
set_par(NCOL(tra2))
for(i in 1:NCOL(tra2)){
     plot(tra2[,i], ptra2[,i], pch=16, cex=0.7, col='#00000050',
          xlab=dimnames(tra2)[[2]][i], ylab=dimnames(ptra2)[[2]][i])
}

### community-weighted means of phylo-corrected traits (PCWM),
###     interpreted here as phylo-corrected trait syndromes
pcwm1 <- makecwm(spe1, ptra1)
pcwm2 <- makecwm(spe2, ptra2)

### NMDS ordination of phylo-corrected trait syndromes
m1 <- ordfn(pcwm1, 'altGower', 2)
m2 <- ordfn(pcwm2, 'altGower', 2)

### GAM gradient regressions give trait-environment fit
(g1 <- gamfit(m1, env))               # *nonlinear* fit
(g2 <- gamfit(m2, env))               # *nonlinear* fit

### plot gradient regressions for each enviro variable
set_par((NCOL(env)+1)*2)
plot(m1, type='t')
for (i in 1:NCOL(env)){
     plot(g1, i, lcol='#FF000080', lwd=1)
}
plot(m2, type='t')
for (i in 1:NCOL(env)){
     plot(g2, i, lcol='#FF000080', lwd=1)
}

### deviation surfaces
fit1 <- fitted(g1)
fit2 <- fitted(g2)
dev  <- ecole::standardize(y1) - ecole::standardize(y2) # calc devn
(gdev <- gamfit(m1, dev))               # *nonlinear* fit
plot(gdev)

### procrustes
p <- protest(m1,m2) # plant vs lichen trait syndromes
p$t0
env$resid <- residuals(p)
plot(env$env1, env$resid)
plot(env$env2, env$resid)

# ### fourthcorner
# (f1 <- fourthcorner(env, spe1, tra1, nrepet=999, modeltype=6))
# (f2 <- fourthcorner(env, spe2, tra2, nrepet=999, modeltype=6))
#
# ## RLQ and fourthcorner.rlq
# r1 <- rlqfn(spe1, env, tra1)
# r2 <- rlqfn(spe2, env, tra2)
# fc_envfit(r1)
# fc_envfit(r2)













