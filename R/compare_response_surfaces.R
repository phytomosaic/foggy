######################################################################
# Comparing two (nonlinear, multidimensional) response surfaces
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 24 Jan 2017
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

###   begin preamble   ###############################################
rm(list=ls())
pkg <- c('vegan', 'coenoflex', 'RColorBrewer', 'ecole')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=T)
rm(pkg, has)
###   end preamble   #################################################

###   begin define functions   #######################################

### simulate community by rejection sampling: iterate until connected
`simfn` <- function(G=2, tol, cmp, ...){
     environment(coenoflex) <- environment()
     i <- 1
     while( i < 999 ){
          spe <- coenoflex(numgrd=G, numplt=100, numspc=100,
                           grdtyp=rep('e',G), grdlen=rep(100,G),
                           width=rep(tol,G), variab=rep(60,G),
                           grdprd=rep(0,G), alphad=rep(1.2,G),
                           pdist='g',sdist='r', skew=3.0,aacorr=0.0,
                           cmpasy=cmp, cmpphy=0.0,
                           maxtot=100, noise=1.0, slack=0.10,
                           autlin='irm(1,2)')$veg
          if( all(rowSums(spe)!=0) &
              sum(unique(distconnected(vegdist(spe),1,F)))==1){
               break
          }
          i <- i + 1
     }
     if(i>1) cat(sprintf('\n%d iters of rejection sampling',i))
     if(i==999) cat(', still disconnected, try increasing iters')
     dimnames(spe)[[1]] <- paste0('SU', seq(1:100))
     dimnames(spe)[[2]] <- paste0('sp', seq(1:100))
     spe <- log10(spe+1)
     spe
}
# helper function to standardize 0-1 an entire matrix
`stdz_matrix` <- function(m, ...){
     v <- as.vector(decostand(as.vector(m), 'range', na.rm=T))
     matrix(v, nrow=nrow(m), ncol=ncol(m))
}
# helper function to add heatmap colors to a perspective plot
`persp_heat` <- function(z, pal='heat', ...){
     nrz <- nrow(z) ; ncz <- ncol(z)
     zfacet <- z[-1,-1] + z[-1,-ncz] + z[-nrz,-1] + z[-nrz,-ncz]
     if(pal=='RdBu'){
          cc <- colorRampPalette(brewer.pal(7, 'RdBu'))(111)[
               as.numeric(cut(zfacet,breaks=111))]
     }
     if(pal=='heat'){
          cc <- heat.colors(100)[cut(zfacet, 100)]
     }
     # if(pal=='div0'){
     #      cc <- colorRampPalette(brewer.pal(7, 'RdBu'))(111)[
     #           as.numeric(cut(zfacet,breaks=seq(-1,1,len=111) ))]
     # }
     persp(z, col=cc, ...)
}
###   end define functions   ##########################################

###   begin generate simulated community data   ########################
#  simulate data, varying both niche widths and competition:
set.seed(19)             # for reproducibility
R    <- 10               # replicates
tols <- seq(20,65,5)     # vary niche tolerances
cmps <- seq(1,19,2)      # vary competition asymmetry
grd  <- expand.grid(tols=tols, cmps=cmps)      # all combinations
egrd <- grd[rep(seq_len(nrow(grd)), each=R),]  # rep ea combn R times
spesim <- mapply(FUN=function(G,tol,cmp)simfn(G,tol,cmp),
                 tol=egrd$tols, cmp=egrd$cmps,
                 MoreArgs=(list(G=2)), SIMPLIFY='array', USE.NAMES=F)
idsim  <- expand.grid(x=seq(0,100,len=10), y=seq(0,100,len=10))
str(idsim)
str(spesim)
head(spesim[,,1])   # peek at first of 1000 species datasets
head(idsim)         # 2-dimensional environmental gradients
# # diversity attributes of simulated datasets:
# div1 <- c(apply(spesim, 3, mx_diversity))
# r <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(length(tols))[
#      as.factor(egrd$tols)]
# pairs(div1[], pch=20, cex=0.7, col=r, lower.panel=NULL, las=2)
###   end generate simulated community data   ########################

###   begin deviation surfaces   #####################################
### 1d example
m1 <- function(x, xo, ym)   -ym * (x-xo)/x
m2 <- function(x, xo, ym, k)-ym * (1-exp(-log(2)*(x-xo)/k))
x <- c(0.52,1.21,1.45,1.64,1.89,2.14,2.47,3.20,4.47,5.31,6.48)
ym <- c(0.00,0.35,0.41,0.49,0.58,0.61,0.71,0.83,0.98,1.03,1.06)
xo <- 1:11
k <- 1.9
y1 <- standardize(m1(x, xo, ym))         # normalize 0-1
y2 <- rev(standardize(m2(x, xo, ym, k))) # normalize 0-1
dev <- y1-y2 # deviations measure how much two curves differ at ea pt
# pdf('~/responsesurf_compare.pdf',wid=6,hei=6,bg='transparent')
par(mfrow=c(1,1), oma=c(0,0,0,0))
plot(x, y1, xlim=c(0,8), ylim=c(-1,1), las=1, bty='l',
     xlab='Ordination axis 1', ylab='Normalized value') ; lines(x, y1)
points(x, y2, col=2) ; lines(x, y2, col=2)
points(x, dev, col=3) ; lines(x, dev, col=3) ; abline(0,0)
legend('bottomright', leg=c('Trait A', 'Enviro A', 'Deviation'),
       pch=1, col=c(1,2,3), lty=1, bty='n')
text(3,-.9,paste('Rho:',formatC(cor(y1,y2,method='spe'),dig=2,form='f')),pos=4)
# dev.off()


### 2d example, borrow NMS from vegan and use GAM surfaces
data(varespec); data(varechem) ; set.seed(47)
nms <- monoMDS(vegdist(varespec))
nms$points[,1] <- nms$points[,1]*(-1) # simple reflection
# # # pdf('C:/Users/Rob/Desktop/deviationsurfaces.pdf',
# # #     wid=12,hei=8,bg='transparent')
# # tiff('C:/Users/Rob/Desktop/deviationsurfaces.tif',
# #      wid=12,hei=8,units='in',compr='lzw+p',res=300)
# # png('C:/Users/Rob/Desktop/deviationsurfaces.png',
# #      wid=12,hei=8,units='in',bg='transparent',res=300)
# jpeg('C:/Users/Rob/Desktop/deviationsurfaces.jpg',
#      wid=12,hei=8,units='in',bg='transparent',res=300)
par(mfrow=c(2,3), las=1, bty='l')
m1 <- ordisurf(nms~Diphcomp, varespec,col=4,main='')
title(main='Lichens', cex.main=3)
m2 <- ordisurf(nms~Cladarbu, varespec,col=4,main='')
title(main='Plants', cex.main=3)
y1 <- calibrate(m1)  # get fitted values from GAM
y2 <- calibrate(m2)  # get fitted values from GAM
dev <- standardize(y1)-standardize(y2) # calc deviation surface values
m3 <- ordisurf(nms ~ dev, col=4, main='')
title(main='Deviation surface', cex.main=3)
persp_heat(m1$grid$z,theta=5,phi=33,r=9,axes=F,expand=.3)
persp_heat(m2$grid$z,theta=5,phi=33,r=9,axes=F,expand=.3)
persp_heat(m3$grid$z,pal='RdBu',theta=5,phi=33,r=9,axes=F,expand=.3)
# dev.off()
# ###   end deviation surfaces   #######################################


###   begin use Atacama preliminary data   ###########################
setwd('C:/Users/Rob/Documents/_prj/6_foggy/data/')
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
     if( any(which(colSums(s)==0)) ){ # check, remove zero-sum cols
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
D1 <- vegdist(s1); D2 <- vegdist(s2)
nms1 <- metaMDS(D1, k=2)
nms2 <- metaMDS(D2, k=2)
nms1$points <- standardize(nms1$points)  # *(-1)
nms2$points <- standardize(protest(nms1, nms2, symm=T, perm=0)$Yrot)
# pdf('patache.pdf',wid=8,hei=6,bg='transparent')
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
persp_heat(m3,pal='RdBu',theta=5,phi=33,r=9,axes=F,expand=.3)
# dev.off()
###   end Atacama preliminary data   #################################




###   END   ###
