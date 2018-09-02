######################################################################
# Comparing two (nonlinear, multidimensional) response surfaces
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 24 Jan 2017
## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

###   begin deviation surfaces   #####################################

### 1d example
m1 <- function(x, xo, ym)   - ym * (x-xo)/x
m2 <- function(x, xo, ym, k)- ym * (1-exp(-log(2)*(x-xo)/k))
x <- c(0.52,1.21,1.45,1.64,1.89,2.14,2.47,3.20,4.47,5.31,6.48)
ym <- c(0.00,0.35,0.41,0.49,0.58,0.61,0.71,0.83,0.98,1.03,1.06)
xo <- 1:11
k <- 1.9
y1 <- ecole::standardize(m1(x, xo, ym))         # normalize 0-1
y2 <- rev(ecole::standardize(m2(x, xo, ym, k))) # normalize 0-1
dev <- y1-y2 # devns measure how much two curves differ at ea pt

# tiff('~/fig/responsesurf_compare.tif',
#      wid=6,hei=6,units='in',compr='lzw+p',res=500)
par(mfrow=c(1,1), oma=c(0,0,0,0))
plot(x, y1, xlim=c(0,8), ylim=c(-1,1), las=1, bty='l',
     xlab='Ordination axis 1', ylab='Normalized value')
lines(x, y1)
points(x, y2, col=2) ; lines(x, y2, col=2)
points(x, dev, col=3) ; lines(x, dev, col=3) ; abline(0,0)
legend('bottomright', leg=c('Trait A', 'Enviro A', 'Deviation'),
       pch=1, col=c(1,2,3), lty=1, bty='n')
text(3,-.9, paste('Rho:',
                 formatC(cor(y1,y2,method='spe'),dig=2,form='f')),
     pos=4)
# dev.off()


### 2d example, borrow NMS from vegan and use GAM surfaces
data(varespec); data(varechem) ; set.seed(47)
nms <- vegan::monoMDS(vegan::vegdist(varespec))
nms$points[,1] <- nms$points[,1]*(-1) # simple reflection

# tiff('~/fig/deviationsurfaces.tif',
#      wid=12,hei=8,units='in',compr='lzw+p',res=500)
par(mfrow=c(2,3), las=1, bty='l')
m1 <- vegan::ordisurf(nms~Diphcomp, varespec,col=4,main='')
title(main='Lichens', cex.main=3)
m2 <- vegan::ordisurf(nms~Cladarbu, varespec,col=4,main='')
title(main='Plants', cex.main=3)
y1 <- calibrate(m1)  # get fitted values from GAM
y2 <- calibrate(m2)  # get fitted values from GAM
dev <- ecole::standardize(y1) - ecole::standardize(y2) # calc devn
m3 <- ordisurf(nms ~ dev, col=4, main='')
title(main='Deviation surface', cex.main=3)
persp_heat(m1$grid$z,theta=5,phi=33,r=9,axes=F,expand=.3)
persp_heat(m2$grid$z,theta=5,phi=33,r=9,axes=F,expand=.3)
persp_heat(m3$grid$z,pal='redblue',theta=5,phi=33,r=9,axes=F,expand=.3)
# dev.off()

###   end deviation surfaces   ######################################
