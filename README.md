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

?mafragh
data(mafragh)
d <- mafragh

# spatial data
xy  <- d$xy

# environmental data
env <- d$env
colnames(env) <- ecole::clean_text(colnames(env), lower=TRUE)

# species data
spe <- d$flo
colnames(spe) <- ecole::clean_text(d$spenames$scientific, lower=TRUE)

# traits data
tra <- d$traits
tra <- data.frame(tra[[1]],tra[[2]],tra[[3]],tra[[4]])

# phylogenetic data
phy <- d$tre
cat(paste0(phy), file = 'phy.tre', sep = '\n')
phy <- read.tree('phy.tre')
unlink('phy.tre')
```

## Visualize

```r
# spatial
plot(xy)

# environment
ecole::plot_heatmap(sapply(env, function(x)100+scale(x)), xord=FALSE)

#species
ecole::plot_heatmap(spe, xord=FALSE)

# traits
ecole::plot_heatmap(sapply(tra, function(x)100+scale(x)), xord=FALSE)

# phylogenetic
plot(phy, cex=0.8, no.margin=TRUE)
```
