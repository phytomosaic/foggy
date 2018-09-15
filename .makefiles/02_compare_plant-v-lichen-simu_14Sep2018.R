######################################################################
# Simulate data to compare plants and lichens
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 14 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

require(foggy)

set.seed(33)
n   <- 25      # n species per group
po  <- rcoal(2, tip.label=paste0('po',1:2))      # basal clade
pl  <- rcoal(n, tip.label=paste0('lichen',1:n))  # lichens
pp  <- rcoal(n, tip.label=paste0('plant',1:n))   # plants
phy <- bind.tree(bind.tree(po, pl, 2), pp, 1)    # joint phylogeny

# simulate trait evolution for plants and lichens separately
tra <- rbind(replicate(4, rTraitCont(pl, sigma=1, root.value=10)),
             replicate(4, rTraitCont(pp, sigma=10, root.value=33)))
tra <- sapply(data.frame(tra), standardize)
tra <- data.frame(tra, row.names=phy$tip.label)
head(tra)

### phylogenetic correction of traits
ptra <- phylo_corr(phy, tra)

### see effects of phylogenetic correction
set_par(NCOL(tra))
par(mfcol=c(2,4))
for(i in 1:NCOL(tra)){
     plot(tra[,i], ptra[,i], pch=16, cex=0.7,
          xlab=dimnames(tra)[[2]][i], ylab=dimnames(ptra)[[2]][i],
          col=colvec(as.factor(c(rep('lichen',n),rep('plant',n)))))
     color.plot.phylo(phy,data.frame(lab=phy$tip.label, tra),
                      colnames(tra)[i],'lab',
                      leg.cex=0.5, col.names=viridisLite::viridis(12))
}

##############################################################



### basic transformations
spe <- data.frame(genlogtrans(spe))
env <- data.frame(decostand(scale(env, center=F), 'range'))
tra <- data.frame(decostand(tra, 'range'))

### phylogenetic correction of traits
ptra <- phylo_corr(phy, tra)

### community-weighted means of phylo-corrected traits (PCWM),
###     interpreted here as phylo-corrected trait syndromes
pcwm <- data.frame(makecwm(spe, ptra))

### NMDS ordination of phylo-corrected trait syndromes
m <- ordfn(pcwm,'altGower', 2)

### GAM gradient regressions give trait-environment fit
vegan::envfit(m, env, perm=999)    # *linear* fit may be unrealistic
g <- gamfit(m, env)                # *nonlinear* fit

### plot gradient regressions for each enviro variable;
###     trait syndrome responses to enviro are commonly nonlinear!
set_par(NCOL(env))
plot(m) # NMDS ordination
for (i in 1:NCOL(env)){
     plot(g, i, lcol='#FF000080', lwd=1)
}

### can also overlay phylo-corrected traits in the same space
gt <- gamfit(m, pcwm)
set_par(NCOL(pcwm))
for (i in 1:NCOL(pcwm)){
     plot(gt, i, lcol='#FF000080', lwd=1)
}





