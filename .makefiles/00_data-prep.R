######################################################################
# Clean and save complete datasets: macroloire, mafragh, pillar
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 12 Sep 2018
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

require(picante)
require(ade4)
require(ecole)

#####################################################################
#####################################################################

## macroinvertebrates of Loire River, France

?macroloire
data(macroloire)
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
invert <- list(spe, env, tra, phy)
names(invert) <- c('spe', 'env', 'tra', 'phy')
save(invert, file = './data/invert.rda')

#####################################################################
#####################################################################

## Mafragh, Algeria plant communities

?mafragh
data(mafragh)
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
veg <- list(xy, spe, env, tra, phy)
names(veg) <- c('xy', 'spe', 'env', 'tra', 'phy')
save(veg, file = './data/veg.rda')

#####################################################################
#####################################################################

### Pillar and Duarte (2010) grassland plants data:
# Data source:  https://github.com/SrivastavaLab/syncsa/

spe <- read.csv('C:/Users/Rob/Desktop/spe.csv')[,-1]
env <- read.csv('C:/Users/Rob/Desktop/env.csv')[,-1]
tra <- read.csv('C:/Users/Rob/Desktop/tra.csv')
Dp  <- read.csv('C:/Users/Rob/Desktop/Dp.csv')
dimnames(tra)[[1]] <- as.character(tra[,1])
tra <- tra[,-1]
dimnames(Dp)[[1]] <- as.character(Dp[,1])
Dp  <- Dp[,-1]
env <- data.frame(nitrogen=env)
cl  <- hclust(as.dist(Dp))
phy <- hclust2phylog(cl)$tre
cat(paste0(phy), file = 'phy_tmp.tre', sep = '\n')
phy <- read.tree('phy_tmp.tre')
unlink('phy_tmp.tre')
phy  <- ape::compute.brlen(phy, power=.7)
nphy <- phy$tip.label
spe  <- spe[,nphy]
tra  <- tra[nphy,]
Dp   <- Dp[nphy,nphy]

# check
names_match(spe, env, tra, phy)

# save to file
pillar <- list(spe, env, tra, phy, Dp)
names(pillar) <- c('spe', 'env', 'tra', 'phy', 'Dp')
save(pillar, file='./data/pillar.rda')

#####################################################################
#####################################################################

### simulated data, for plants and lichens

`simfn` <- function(...){

     require(ecolottery)
     require(foggy)

     # parameter setup
     time0 <- Sys.time()
     n     <- 25    # n species per plants/lichens
     nspe  <- n*2   # n species
     ntra  <- 3     # n traits
     nsite <- 33    # n sites
     nenv  <- 4     # n enviro
     sd_li <- 1     # sd of trait values
     sd_pl <- 2     # sd of trait values
     mn_li <- 10    # mean of trait values
     mn_pl <- 15    # mean of trait values

     # simulate phylogeny, including both lichens and plants
     po  <- rcoal(2, tip.label=paste0('po',1:2))      # basal clade
     pl  <- rcoal(n, tip.label=paste0('lichen',1:n))  # lichens
     pp  <- rcoal(n, tip.label=paste0('plant',1:n))   # plants
     phy <- bind.tree(bind.tree(po, pl, 2), pp, 1)    # both

     # simulate a species pool
     J     <- 100        # n individuals in local community
     Jpool <- nspe*100   # n individuals in regional pool
     pool  <- data.frame(ind=1:Jpool,
                         sp=rep(phy$tip.label,100),
                         tra=rep(NA,Jpool))

     # simulate trait evolution for plants and lichens
     tra <- rbind(
          replicate(ntra,rTraitCont(pl,sigma=sd_li,root.value=mn_li)),
          replicate(ntra,rTraitCont(pp,sigma=sd_pl,root.value=mn_pl)))
     tra <- sapply(data.frame(tra), standardize)
     tra <- data.frame(tra, row.names=phy$tip.label)
     pool$tra <- tra[pool$sp,]

     # species abundances from filtering on multiple traits
     spe <- matrix(NA, ncol=nspe, nrow=nsite,
                   dimnames=list(NULL,dimnames(tra)[[1]]))
     `filt_gauss` <- function(topt, x, sigma){
          exp(-(x - topt)^2/(2*sigma^2))
     }
     `f` <- function(x) {
          filt_gauss(0.5, x[1], 0.6)*
               filt_gauss(0.25, x[2], 0.3)*
               filt_gauss(0.75, x[3], 0.3)
     }
     for(i in 1:nsite){
          r   <- coalesc(J, m=1, pool=pool, filt=f, traits=tra)
          com <- abund(r)$com
          spe[i,] <- com$ab[match(dimnames(spe)[[2]],com$sp)]
     }
     spe[is.na(spe)] <- 0
     if(any(zero <- which(colSums(spe, na.rm=T)==0))){
          cat('adding trace values to zero-sum species\n')
          rnd <- sample.int(nsite, length(zero))
          for(i in 1:length(zero)){
               spe[rnd[i], zero[i]] <- min(spe, na.rm=T)
          }
     }
     spe <- as.data.frame(spe)

     # simulate correlated enviro variables (in a quirky way)
     scr <- prcomp(spe)$x[,3:4] # 3rd and 4th PCs
     env <- scr + (runif(nsite*2, 0, diff(range(scr))*0.1) *
                        sample(c(1,-1), nsite*2, repl=T))
     dimnames(env)[[2]] <- c('env1','env2')
     env <- as.data.frame(env)

     # result
     cat('time elapsed:', Sys.time()-time0, '\n')
     list(spe=spe, env=env, tra=tra, phy=phy)
}

# replicate several simulations
nrep    <- 5   # number of replicates (increase for production...)
simdata <- lapply(1:nrep, FUN=simfn)

# save to file
save(simdata, file = './data/simdata.rda')

####    END    ######################################################
