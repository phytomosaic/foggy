# foggy

Species and trait responses to fog gradients.


## Contributors

Daniel Stanton
Reinaldo Vargas-Castillo
Peter Nelson
Robert Smith


## Motivation

Analysis repository created to faciliate work on species, traits, and phylogenies in relation to Pacific fog gradients.  Development of multivariate tools.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/foggy')
require(foggy)
```


## Usage
```r
data(pavoine)

# spatial data
xy  <- d$mxy

# environmental data
env <- d$env

# species data
spe <- d$flo

# traits data
tra <- d$traits

# phylogenetic data
phy <- d$phy
```
