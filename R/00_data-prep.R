######################################################################
# Clean and save complete datasets: macroloire and mafragh
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 10 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

require(picante)
require(ade4)
require(ecole)

?macroloire
data(macroloire) # macroinvertebrates of Loire River, France
d <- macroloire

# labels
lab <- d$labels
lab$lab <- dimnames(lab)[[1]]
lab$nm  <- ecole::clean_text(as.character(lab$latin),lower=TRUE)

# species
spe <- d$fau
dimnames(spe)[[1]] <- lab$nm[match(dimnames(spe)[[1]], lab$lab)]
spe <- data.frame(t(spe))

# environment
env <- d$envir
env <- env[,2:4]
dimnames(env)[[2]] <- ecole::clean_text(dimnames(env)[[2]], TRUE)

# traits
tra <- d$traits
dimnames(tra)[[2]] <- ecole::clean_text(dimnames(tra)[[2]], TRUE)
dimnames(tra)[[1]] <- lab$nm[match(dimnames(tra)[[1]], lab$lab)]
head(tra)

# phylogeny (from a taxon classification)
tax <- d$taxo
dimnames(tax)[[1]] <- lab$nm[match(dimnames(tax)[[1]], lab$lab)]
tax <- tax[match(dimnames(spe)[[2]], dimnames(tax)[[1]]),]
phy <- taxo2phylog(tax)
rm(tax)
phy <- phy$tre
cat(paste0(phy), file = 'phy_tmp.tre', sep = '\n')
phy <- read.tree('phy_tmp.tre')
unlink('phy_tmp.tre')
phy <- ape::compute.brlen(phy, power=0.75)
phy <- ape::multi2di(phy)

### check congruence of names
`names_match` <- function(spe, env, tra, phy, ...){
     stopifnot(class(phy)=='phylo')
     nspe <- dimnames(spe)
     nenv <- dimnames(env)
     ntra <- dimnames(tra)
     nphy <- phy$tip.label
     n1 <- identical(nspe[[1]], nenv[[1]]) # site names
     n2 <- identical(nspe[[2]], ntra[[1]]) # species names
     n3 <- identical(nspe[[2]], nphy)      # species names
     all(n1,n2,n3)
}

# check
names_match(spe, env, tra, phy)

# save to file
save(spe, env, tra, phy, file = './data/invert.rda')

################################

?mafragh
data(mafragh)  # Mafragh, Algeria plant communities
d <- mafragh

# spatial
xy  <- d$xy

# labels
lab <- d$spenames
lab$lab <- dimnames(lab)[[1]]
lab$nm  <- ecole::clean_text(as.character(lab$scientific),lower=TRUE)
lab$nm[lab$nm == 'alisma_plantago'] <- 'alisma_plantago_aquatica'

# species
spe <- d$flo
dimnames(spe)[[2]] <- lab$nm[match(dimnames(spe)[[2]], lab$lab)]

# environment
env <- d$env
colnames(env) <- ecole::clean_text(colnames(env), lower=TRUE)

# traits
tra <- d$traits
tra <- data.frame(tra[[1]],tra[[2]],tra[[3]],tra[[4]])
colnames(tra) <- ecole::clean_text(colnames(tra), lower=TRUE)
tra <- tra[,!names(tra)%in%c('seasonnal','succulence')]

# phylogeny
phy <- d$tre
cat(paste0(phy), file = 'phy_tmp.tre', sep = '\n')
phy <- read.tree('phy_tmp.tre')
unlink('phy_tmp.tre')

# check
names_match(spe, env, tra, phy)

# save to file
save(xy, spe, env, tra, phy, file = './data/veg.rda')

####    END    ######################################################
