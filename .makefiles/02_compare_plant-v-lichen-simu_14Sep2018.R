######################################################################
# Simulate data to compare plants and lichens
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 14 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

# GITHUB_PAT <- '01591b7a42837969a5e9360bb61089d5dc1822a0'
# devtools::install_github('phytomosaic/foggy', auth_token=GITHUB_PAT)

require(foggy)

### load data
data(simdata)
spe <- simdata[[1]]$spe
env <- simdata[[1]]$env
tra <- simdata[[1]]$tra
phy <- simdata[[1]]$phy

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
for (i in 1:NCOL(env)){
     plot(g, i, lcol='#FF000080', lwd=1)
}

### can also overlay phylo-corrected traits in the same space
gt <- gamfit(m, pcwm)
set_par(NCOL(pcwm))
for (i in 1:NCOL(pcwm)){
     plot(gt, i, lcol='#FF000080', lwd=1)
}


