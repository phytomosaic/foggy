# foggy

Gradient regression to integrate phylogeny, species, and trait responses to environment.


## Contributors

Daniel Stanton  
Reinaldo Vargas-Castillo  
Peter Nelson  
Robert Smith


## Motivation

Analysis repository to faciliate work on species, traits, and phylogenies of plants and lichens in relation to Pacific fog gradients.  Development of multivariate analysis tools. 
Gradient regression to integrate phylogeny, species, and trait responses to environment.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/foggy')
devtools::install_github('phytomosaic/ecole')
```


## Load data

```r
require(foggy)
?veg
data(veg)
d   <- veg
xy  <- d$xy    # spatial
spe <- d$spe   # species
env <- d$env   # environment
tra <- d$tra   # traits
phy <- d$phy   # phylogeny
```

## Visualize

```r
# spatial
plot(xy, pch=19, col='#00000050')

#species
plot_heatmap(spe, xord=FALSE, logbase=10)

# environment
plot_heatmap(sapply(env, function(x)100+scale(x)), xord=FALSE)

# traits
plot_heatmap(sapply(tra, function(x)100+scale(x)), xord=FALSE)

# phylogeny
plot(phy, cex=0.6, no.margin=TRUE)

```

### Introducing gradient regression

![Gradient regression](https://github.com/phytomosaic/foggy/blob/master/.makefiles/workflow.png)

```r
### basic transformations
spe <- data.frame(genlogtrans(spe))
env <- data.frame(decostand(scale(env, center=F), 'range'))
tra <- data.frame(decostand(tra, 'range'))

### detect phylogenetic signal for traits
K <- sapply(tra, FUN=function(j){
     names(j) <- rownames(tra)
     round(picante::phylosignal(j, multi2di(phy)),6)})
as.matrix(K) # Blomberg's K statistic for continuous traits

### phylogenetic correction of traits (given phylo signal)
ptra <- phylo_corr(phy, tra)

### see effects of phylogenetic correction
set_par(NCOL(tra))
for(i in 1:NCOL(tra)){
     plot(tra[,i], ptra[,i], pch=16, cex=0.7, col='#00000050',
          xlab=dimnames(tra)[[2]][i], ylab=dimnames(ptra)[[2]][i])
}

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

#### end ####
```




