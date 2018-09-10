# foggy

Species and trait responses to fog gradients.


## Contributors

Daniel Stanton  
Reinaldo Vargas-Castillo  
Peter Nelson  
Robert Smith


## Motivation

Analysis repository to faciliate work on species, traits, and phylogenies of plants and lichens in relation to Pacific fog gradients.  Development of multivariate analysis tools.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/foggy')
devtools::install_github('phytomosaic/ecole')
```


## Load data

Data source:  
Pavoine, S., Vela, E., Gachet, S., de Bélair, G. and Bonsall, M. B. 2011. Linking patterns in phylogeny, traits, abiotic variables and space: a novel approach to linking environmental filtering and plant community assembly. Journal of Ecology 99:165–175. doi:10.1111/j.1365-2745.2010.01743.x

```r
require(foggy)
require(ade4)
require(ecole)
require(picante)

data(veg)
d <- veg

# data(invert)
# d <- invert

# spatial data
xy  <- d$xy

# species data
spe <- d$spe

# environmental data
env <- d$env

# traits data
tra <- d$tra

# phylogenetic data
phy <- d$phy
```

## Visualize

```r
### load('./data/veg.rda')
### load('./data/invert.rda')

# spatial
plot(xy, pch=19, col='#00000050')
plot(xy, pch=19, col=ecole::colvec(spe$bolboschoenus_maritimus, alpha=1))

#species
ecole::plot_heatmap(spe, xord=FALSE, logbase=10)

# environment
ecole::plot_heatmap(sapply(env, function(x)100+scale(x)), xord=FALSE)

# traits
ecole::plot_heatmap(sapply(tra, function(x)100+scale(x)), xord=FALSE)

# phylogeny
plot(phy, cex=0.6, no.margin=TRUE)

```
