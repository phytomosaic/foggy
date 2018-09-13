#' @title GAM for environment-ordination fit
#'
#' @description Generalized Additive Model (GAM) to regress each
#'     environmental variable on NMDS ordination scores.
#'
#' @param ord ordination model from \code{vegan}
#'
#' @param env environmental matrix with variables to evaluate
#'
#' @param x,object object of class \code{gamfit}
#'
#' @param pick numeric, the column number of environmental variable to
#'     plot as a nonlinear regression surface in the ordination space
#'
#' @param pcol,lcol point and line colors
#'
#' @param pcex,lwd point size and line width
#'
#' @param title logical, add variable name as title? (default=TRUE)
#'
#' @param ...  additional arguments passed to
#'     \code{\link[mgcv]{gam}}
#'
#' @return
#' Object of class \code{gamfit}.
#'
#' @details
#' Nonlinear regression of each environmental variable in turn on the
#'     site scores resulting from NMDS ordination.  Returns p-values
#'     for the joint smooth term, and adjusted-R2 expressing the
#'     strength of (nonlinear) relationship between each environmental
#'     variable and the ordination. Print, summary and plot methods
#'     exist. Sensitivity analysis is yet to do.
#'
#' @examples
#' data(veg) # load('./data/veg.rda', verbose=T)
#' spe <- veg$spe
#' env <- veg$env
#' m <- vegan::metaMDS(vegan::vegdist(spe), trace=0)
#' g <- gamfit(m, env)
#' g
#' plot(g, pick=10)
#'
#' @seealso
#' \code{\link[mgcv]{gam}}
#'
#' @export
#' @rdname gamfit
`gamfit` <- function(ord, env, ...){
     if (!any(match(class(ord),'metaMDS')))
          stop('`ord` must be of class `metaMDS`')
     scr   <- vegan::scores(ord)
     scrnm <- dimnames(scr)[[2]]
     envnm <- dimnames(env)[[2]]
     ndim  <- ord$ndim
     nenv  <- length(envnm)
     stopifnot(identical(dimnames(scr)[[1]], dimnames(env)[[1]]))
     xx    <- data.frame(env, scr)
     # dimnames(xx)[[2]] <- c(envnm, scrnm)
     xn <- rep(NA, ndim)
     ml <- vector('list', ndim)
     st <- matrix(NA, nrow=nenv, ncol=2)
     dimnames(st)[[1]] <- envnm
     dimnames(st)[[2]] <- c('pval','adj_r2')
     for(i in 1:ndim){
          xn[i] <- scrnm[i]
     }
     right <- paste0('s(',paste(xn,collapse=','),')')
     for(i in 1:nenv){
          left    <- paste0(envnm[i], ' ~ ')
          fmla    <- as.formula(paste(left, right))
          ml[[i]] <- mgcv::gam(fmla, data=xx, ...)
          ss      <- summary(ml[[i]])
          st[i,1] <- as.numeric(sprintf('%.3f',round(ss$s.pv,3)))
          st[i,2] <- as.numeric(sprintf('%.3f',round(ss$r.sq,3)))
     }
     out <- list(mods=ml, sumtab=st)
     class(out) <- 'gamfit'
     out
}
#' @export
#' @rdname gamfit
`print.gamfit` <- function(x, ...){
     print(x[['sumtab']])
}
#' @export
#' @rdname gamfit
`summary.gamfit` <- function(object, ...){
     object[['sumtab']]
}
#' @export
#' @rdname gamfit
`plot.gamfit` <- function(x, pick=1, pcol, lcol, pcex, lwd,
                          title=TRUE, ...){
     x <- x[['mods']][[pick]] # currently best w 2 smooth predictors
     m <- x$model
     xx <- m[,2:NCOL(m)]
     if(missing(lcol)) lcol <- '#00000095'
     if(missing(pcol)) pcol <- '#00000080'
     if(missing(pcex)) pcex <-  0.5
     if(missing(lwd))  lwd  <-  2
     # currently works but throws warning for 'type':
     if(title) main <- dimnames(m)[[2]][1] else main <- ''
     plot(x, se=F, col=lcol, lwd=lwd, main=main, type='n', ...)
     points(xx, pch=16, col=pcol, cex=pcex, ...)
}

# ###   GAM sensitivity analysis   ###################################
# ###      add noise to each NMDS axis in turn, calc MAE
#
# load('./data/veg.rda', verbose=T)
# m    <- vegan::metaMDS(vegan::vegdist(spe), trace=0)
# x0   <- gamfit(m, env) # baseline iteration
# nrep <- 100    # average the sensitivity over 100 iterations
# `f` <- function(perc=5.0, ...){ # convenience function
#      xx   <- gamfit(m, noise(env))
#      c(mae(x0[,1],xx[,1], stdz=T),
#        mae(x0[,2],xx[,2], stdz=T)
#      )
# }
# # TIME WARN
# Q <- matrix(nrow=nrep,ncol=3) # initiate sensitivity matrix
# for(i in 1:nrep){
#      cat('round',i,'of',nrep,'\n')
#      Q[i,] <- f(perc=5.0)     # add 5% uncertainty, get Q values
# }
# colMeans(Q, na.rm=T)          # summary of mean sensitivity at 5%
# # TIME WARN
# Q <- matrix(nrow=nrep,ncol=3) # initiate sensitivity matrix
# for(i in 1:nrep){
#      cat('round',i,'of',nrep,'\n')
#      Q[i,] <- f(perc=10.0)    # add 10% uncertainty, get Q values
# }
# colMeans(Q, na.rm=T)          # summary of mean sensitivity at 10%
# #
# rm(Q, f, nrep, noise)
# ###   end sensitivity analysis   #####################################
