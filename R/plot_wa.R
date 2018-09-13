#' @title Plot centroids and ranges
#'
#' @description Plot species centroids and ranges in ordinations.
#'
#' @param spe  species abundance matrix used in \code{mod}
#'
#' @param mod  ordination model from \code{vegan}
#'
#' @param xexp,yexp	 numeric, expansion factor making room for x,y
#'     labels
#'
#' @param pick  numeric, single ordination axis to plot (default =
#'     \code{NULL} will plot first 2 axes as a biplot)
#'
#' @param lcol  color for range lines
#'
#' @param buff  numeric, expansion factor giving buffer around range
#'     lines
#'
#' @param ...  additional arguments passed to
#'     \code{\link[vegan]{ordiplot}}
#'
#' @return
#' Plots to device.
#'
#' @details
#' Plot species centroids and ranges in ordination spaces or single
#'     ordination axes.
#'
#' @examples
#' data(veg) # load('./data/veg.rda', verbose=T)
#' spe <- veg$spe
#' m <- vegan::metaMDS(vegan::vegdist(spe), trace=0)
#' plot_wa(spe, m, cex=0.7, las=1)
#'
#' @seealso
#' \code{\link[vegan]{ordiplot}}
#'
#' @export
#' @rdname plot_wa
`plot_wa` <- function(spe, mod, xexp=1.6, yexp=1, pick=NULL, lcol,
                      buff=0.1, ...){
     tmp <- spe
     tmp[tmp>0]  <- 1
     tmp[tmp<=0] <- NA
     max1 <- sapply(tmp * mod$points[,1], max, na.rm=TRUE)
     min1 <- sapply(tmp * mod$points[,1], min, na.rm=TRUE)
     max2 <- sapply(tmp * mod$points[,2], max, na.rm=TRUE)
     min2 <- sapply(tmp * mod$points[,2], min, na.rm=TRUE)
     cent <- vegan::wascores(mod$points, spe, expand=F)
     c1   <- cent[,1]
     c2   <- cent[,2]
     if(missing(lcol)) lcol <- '#00000050'
     if(is.null(pick)){
          ## two axis biplot
          plot(c1, c2, xlim=c(min(min1),max(max1))*(1+buff),
               ylim=c(min(min2),max(max2))*(1+buff),
               xlab='NMDS1', ylab='NMDS2', ...)
          for(i in 1:NROW(c1) ){
               lines(c(min1[i],max1[i]), rep(c2[i],ea=2), col=lcol)
               lines(rep(c1[i],ea=2), c(min2[i],max2[i]), col=lcol)
          }
     } else {
          ## single axis
          xlab <- c('NMDS1','NMDS2')[pick]
          cc   <- list(c1,c2)[[pick]]
          minn <- list(min1,min2)[[pick]]
          maxx <- list(max1,max2)[[pick]]
          o <- order(cc)
          op <- par(mfrow = c(1, 1),
                    mar = c(4.1 * yexp, 5 * xexp, 0, 0) + 0.1,
                    oma = c(0, 0, 0, 0), font = 2)
          plot(cc[o], 1:length(cc),
               xlim=c(min(minn),max(maxx))*(1+buff),
               xlab=xlab, ylab='', axes=F, ...)
          for (i in 1:length(cc)) {
               lines(c(minn[o][i],maxx[o][i]), rep(i,ea=2), col=lcol)
          }
          axis(1, at=pretty(cc), label=pretty(cc), las=3, cex.axis=0.6)
          axis(2, at=1:length(cc), label=names(cc)[o], las=1,
               cex.axis=0.6, tick=F)
     }
}
