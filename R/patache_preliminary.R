######################################################################
# Preliminary Patache data
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 24 Jan 2017
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

rm(list=ls())
pkg <- c('vegan', 'ecole')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=T)
rm(pkg, has)

### helper function to standardize 0-1 an entire matrix
`stdz_matrix` <- function(m, ...){
     v <- as.vector(decostand(as.vector(m), 'range', na.rm=T))
     matrix(v, nrow=nrow(m), ncol=ncol(m))
}

### load data from local
setwd('~/data/')
spe <- read.csv('spe_27Jul2017.csv', header=T, stringsAsFactors=F,
                na.strings=c('','-9999','NA'), strip.white=T,
                row.names=1)
colnames(spe) <- tolower(colnames(spe))

id <- read.csv('id_27Jul2017.csv', header=T, stringsAsFactors=F,
               na.strings=c('','-9999','NA'), strip.white=T,
               row.names=1)
colnames(id) <- tolower(colnames(id))
identical(row.names(spe), row.names(id))

# the elevation-richness relationship:
plot(id$elev, id$n_spp, xlim=c(425, 875))
boxplot(id$n_spp~id$elev, boxwex=0.2, las=1)

# arbitrarily partition half the species
s1 <- spe[,seq(1, ncol(spe), 2),]
s2 <- spe[,seq(1, ncol(spe), 2)+1,]
`rm_empty` <- function(s){
     if (any(which(colSums(s)==0))){ # check, remove zero-sum cols
          s <- s[ , !( colSums(s) == 0 ) ]
          print('removed zero-sum COLUMNS in spe data')
     }else{print('no zero-sum columns found in spe data')}
     if( any(which(rowSums(s)==0)) ){ # check, remove zero-sum rows
          s <- s[ !( rowSums(s) == 0 ) , ]
          print('removed zero-sum ROWS in spe data')
     }else{print('no zero-sum rows found in spe data')}
     return(s)
}
s1 <- rm_empty(s1) ; s2 <- rm_empty(s2)
e1 <- e2 <- id
keeps <- intersect(row.names(s1),row.names(s2))
keeps <- keeps[!keeps%in%c('500_A_1','450_A_5','500_A_4','500_A_15',
                           '850_A_5','850_A_7')]
s1 <- s1[row.names(s1)%in%keeps,]; s2 <- s2[row.names(s2)%in%keeps,]
e1 <- e1[row.names(e1)%in%keeps,]; e2 <- e2[row.names(e2)%in%keeps,]

### ordinations
D1 <- vegdist(s1); D2 <- vegdist(s2)
nms1 <- metaMDS(D1, k=2)
nms2 <- metaMDS(D2, k=2)
nms1$points <- standardize(nms1$points)  # *(-1)
nms2$points <- standardize(protest(nms1, nms2, symm=T, perm=0)$Yrot)

# tiff('~/fig/patache.tif',wid=8,hei=6,bg='transparent')
par(mfrow=c(2,3), las=1, bty='l')
m1 <- ordisurf(nms1~r_median,e1,col=4,main='Lichens',
               xlim=c(0,1),ylim=c(0,1))
m2 <- ordisurf(nms2~r_median,e2,col=4,main='Plants',
               xlim=c(0,1),ylim=c(0,1))
m3 <- stdz_matrix(m1$grid$z) - stdz_matrix(m2$grid$z) # deviation
contour(m3, col=4, main='Deviation surface',
        xlim=c(0,1),ylim=c(0,1))
persp_heat(m1$grid$z,    theta=5,phi=33,r=9,axes=F,expand=.3)
persp_heat(m2$grid$z,    theta=5,phi=33,r=9,axes=F,expand=.3)
persp_heat(m3,pal='redblue',theta=5,phi=33,r=9,axes=F,expand=.3)
# dev.off()

###   end Atacama preliminary data   #################################
