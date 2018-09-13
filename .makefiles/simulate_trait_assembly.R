## Nathan Kraft
## nkraft@umd.edu

## Trait-based community assembly simulation functions in R source:
#http://life.umd.edu/biology/kraftlab/Code_files/community_assembly.R

## Supplemental code from Kraft, N. J. B. and D. D. Ackerly. 2010.
## Functional trait and phylogenetic tests of community assembly
## across spatial scales in an Amazonian forest. Ecological Monographs
## 80:401-422.

## Functions here are based on the EVELYN community assembly
## model(which was originally written in JAVA) from: Kraft, N. J. B.,
## W. K. Cornwell, C. O. Webb, and D. D. Ackerly. 2007. Trait
## evolution, community assembly, and the phylogenetic structure of
## ecological communities. American Naturalist 170:271-283.

## for the original verbal description of the competition algorithm
## used here and general inspiration see: Colwell, R. K., and D. W.
## Winkler. 1984. A null model for null models in biogeography., Pages
## 344-359 in D. R. Strong, D. S. Simberloff, L. G. Abele, and A. B.
## Thistle, eds. Ecological communities: conceptual issues and the
## evidence. Princeton, NJ, Princeton University Press.

## the key functions are: random_assembly(), compete_until() and
## filter_until() Everything else supports these three assembly
## functions.


########################
### GENERAL FUNCTIONS ##
########################

#check for duplicate species names in a community vector
validate_community <- function(community=community){
     spnames <- community[,1]
     if(length(spnames)>length(unique(spnames))){
          print("warning- duplicate species names")
          return(FALSE)
     }
     return(TRUE)
}

#remove a species from the community (general function used in all models)
zap <- function(target, community=community){
     #target can be either a species name (letters) or a column number(numeric).
     #Can also be a vector of names. Defaults to looking for a vector called
     #community with two columns (first for name and second for trait value)
     #currently no error is given if some targets match species in the community
     #and others do not- only matching species removed

     if(!is.numeric(target)){
          target<-which(target==community[,1])

          if(length(target)<1){
               print("target does not match anyone in community- can not zap")
               return(community)
          }
     }
     community<-community[-target,]
     return(community)
}

####################
## RANDOM ASSEMBLY #
####################

## Randomly removes species, weighted by abundance, from the species pool until
## a final_richness value is reached
random_assembly <- function(pool=pool, final_richness, abund=NULL){
     size     <- nrow(pool)
     n_victim <- size-final_richness
     if(n_victim<1){
          print("pool has equal or lesser richness to final size- can't do random assembly")
          return(pool)
     }

     if(is.null(abund)){
          alive <- sample(1:size, final_richness)
          return(pool[alive,])
     }
     merge(pool, abund)->poolPlus
     alive <- sample((1:nrow(poolPlus)), final_richness, prob=poolPlus$abund)
     return(poolPlus[alive,c(1,2)])
}


######################
## COMPETITION MODEL #
######################

#find the most similar taxa within the community:
find_most_similar_taxa <- function(community=community, speciesNamesCol1=TRUE){
     if(speciesNamesCol1!=TRUE){
          print("first col needs to be species names for now")
          return(community)
     }
     ##somewhat of an awkward fix here:
     labels <-community[,1]
     traits <-community[,-1]
     m  <- dist(traits)
     m2 <- as.matrix(m)
     m3 <- which(m2==min(m), arr.ind=TRUE)
     a  <- as.vector(m3[,1])
     b  <- as.vector(m3[,2])
     ##will contain multiples but this is good- more weight to a "sandwiched" taxa being removed
     taxa_indices <- c(a[which(a!=b)], b[which(a!=b)])
     return(labels[taxa_indices])
}

## find most similar species and remove it from a community
one_round_competition <- function(community=community){
     threatened_list<-find_most_similar_taxa(community)
     n<-length(threatened_list)
     if(n<2){
          print("couldn't find most similar taxa for competition")
          return(community)
     }
     return(zap(sample(threatened_list, 1), community=community))
}

## remove species via competition from the community until a specified number of
## taxa are left (nfinal)
compete_until <- function(nfinal, community=community){
     nrow(community)->start
     if(nfinal>start){
          print("community is already smaller than target- can't compete")
          return(community)
     }
     if(nfinal<4){
          print("target for competition is too small- can't compete")
          return(community)
     }
     to_kill <- start - nfinal
     for(i in 1:to_kill){
          community <- one_round_competition(community)
     }
     if(nrow(community)!=nfinal){
          print('error- competition ended with incorrect number of species')
     }
     return(community)
}

######################
## Habitat Filtering #
######################

#Used to identify which species in a community is farthest from the trait optima used for habitat filtering- farthest species are removed first.

find_farthest_from_optima <- function(optima=optima,
                                      community=community,
                                      speciesNamesCol1=TRUE){
     if(speciesNamesCol1!=TRUE){
          print("first col needs to be species names for now")
          return(community)
     }

     if( (length(optima)!=(ncol(community)-1))){
          print("error- filtering optima is not right dimension for community")
          return(community)
     }

     start  <- paste(community[,1])
     labels <-c(start,"optima")
     traits <-rbind(community[,-1], optima)

     dist(traits)->m
     as.matrix(m)->m2
     ##get last row, which is the distance of each taxa from optimum
     m2[,nrow(m2)]->comp
     which(comp==max(comp))->index
     return(labels[index])
}


##identify farthest species from optima and remove it
one_round_filtering_optima <- function(community=community,
                                       optima=optima){
     threatened_list <- find_farthest_from_optima(optima, community)
     n <- length(threatened_list)
     if(n<1){
          print("couldn't find farthest taxa for filtering")
          return(community)
     }
     return(zap(sample(threatened_list, 1), community=community))
}


##run filtering model until only nfinal number of taxa are left in
## the community
filter_until <- function(nfinal, community=community, optima=optima){
     nrow(community)->start
     if(nfinal>start){
          print("community is already smaller than target- can't filter")
          return(community)
     }

     if(nfinal<4){
          print("target for filtering is too small- can't filter")
          return(community)
     }

     start-nfinal->to_kill

     for(i in 1:to_kill){
          one_round_filtering_optima(community, optima)->community

     }

     if(nrow(community)!=nfinal){
          print('error- filtering ended with incorrect number of species')
     }
     return(community)

}

#################################################################
#################################################################

######################
## Example/ tutorial #
######################

## Make sure that all of the functions in the sections above have been run in the terminal before running the code below.

# Create dummy communities with 1 or 2 traits:
sp   <- c("A","B","C", "D", "E", "F")
tr1  <- c(1,2,9,10,9, 13)
tr2  <- c(1,1,100,1000, 10000, 1000000)

## test communities, with one and two traits
test <- data.frame(sp, tr1)
test2<- data.frame(sp, tr1, tr2)

## for competition, use the compete_until() function
## It needs species names and at least one trait- it will
## use as many traits as you give it.
## First argument specifies how many taxa you want to end up with.

compete_until(4, test) # compete until 4 spp remain
## ties are broken randomly, so multiple calls give different communities:
compete_until(4, test)
compete_until(4, test)
compete_until(4, test)


## you need to end up with at least 4 taxa (you could change this in the code)
## I imposed this for trait based tests, where you need n=4 for some metrics
compete_until(3, test) # throws error, returns original data

## Works with multiple traits too:
compete_until(4, test2)


## for filtering use filter_until(). You need to specify what the optimum trait
## value is- in the AmNat paper we used the ancestral trait value for the clade
## or the most extreme (largest) derived trait values observed in the data. In
## the current Ecological Monographs paper we randomly place the optima in the
## environment.
filter_until(4, test, optima=c(0) )

## for multiple traits, the optima needs to be of the same number of dimensions:
filter_until(4, test2, optima=c(0) ) # throws error
filter_until(4, test2, optima=c(0,0) ) # works
## again, ties are broken randomly so repeat calls can give different outcomes
