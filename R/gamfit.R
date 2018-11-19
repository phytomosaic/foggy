#' @title GAM for environment-ordination fit
#'
#' @description Generalized Additive Model (GAM) to regress one or
#'     many environmental variables on NMDS ordination scores.
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
#' @param nrep number of repetitions for sensitivity analysis
#'     (default=99)
#'
#' @param perc numeric, percent noise to add to env variables for
#'     sensitivity analysis
#'
#' @param r object of class \code{gamfit_sens}, from sensitivity
#'     analysis
#'
#' @param stat one of \code{c('r2','pval')} to plot from sensitivity
#'     analysis
#'
#' @param pltype one of \code{c('heat','joy')} for a heatmap or
#'     joyplot, respectively
#'
#' @param ...  additional arguments passed to
#'     \code{\link[mgcv]{gam}}
#'
#' @return
#' Object of class \code{gamfit} or \code{gamfit_sens} for those
#'     functions respectively, or plots to device.
#'
#' @details
#' Nonlinear regression of each environmental variable in turn on the
#'     site scores resulting from NMDS ordination.  Returns p-values
#'     for the joint smooth term, and adjusted-R2 expressing the
#'     strength of (nonlinear) relationship between each environmental
#'     variable and the ordination. Print, summary and plot methods
#'     exist. Sensitivity analysis currently exists but needs
#'     refinement.
#'
#' @examples
#' data(veg)
#' spe <- veg$spe
#' env <- veg$env
#' m <- vegan::metaMDS(vegan::vegdist(spe), trace=0)
#' g <- gamfit(m, env)
#' g
#' plot(g, pick=10)
#'
#' # sensitivity analysis
#' r <- gamfit_sens(g, env)
#' plot(r, pltype='heat')
#' plot(r, pltype='joy')
#' (MAE <- summary(r, g))
#' plot(MAE, xaxt='n', las=1)
#' axis(1, 1:NCOL(env), names(env), cex.axis=0.7)
#' sapply(data.frame(r[,,1]<0.05),sum) / 99 # proportion 'significant'
#'
#' @seealso
#' \code{\link[mgcv]{gam}} for the internal fitting function,
#'     \code{\link[vegan]{ordisurf}} for a similar procedure that only
#'     admits variables one-at-a-time, and \code{\link[vegan]{envfit}}
#'     for a linear alternative.  Plotting methods for sensitivity
#'     analysis follow \code{\link[ecole]{plot_joy}} and
#'     \code{\link[ecole]{plot_heatmap}}.
#'
#' @export
#' @rdname gamfit
#'
`gamfit` <- function(ord, env, ...){
     # if (!inherits(ord, c('metaMDS','procrustes'))){
     #      stop('`ord` must be of class `metaMDS` or `procrustes`')
     # }
     if (inherits(ord, 'metaMDS')){
          scr  <- vegan::scores(ord)
          ndim <- ord$ndim
     } else {
          scr  <- as.data.frame(ord)
          ndim <- (dim(ord)[[2]])
     }
     scrnm <- dimnames(scr)[[2]]
     envnm <- dimnames(env)[[2]]
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
`fitted.gamfit` <- function(object, ...){
     ml <- object$mods
     nc <- length(ml)
     nr <- length(ml[[1]]$fitted.values)
     out<- data.frame(matrix(NA, nrow=nr, ncol=nc))
     cn <- dimnames(object$sumtab)[[1]]
     for(j in 1:nc){
          out[[j]] <- ml[[j]]$fitted.values
     }
     dimnames(out)[[2]] <- cn
     out
}
#' @export
#' @rdname gamfit
`plot.gamfit` <- function(x, pick=1, pcol, lcol, pcex, lwd,
                          title=TRUE, ...){
     x <- x[['mods']][[pick]] # currently best w 2 smooth predictors
     m <- x$model
     xx <- m[,2:NCOL(m)]
     if(missing(lcol)) lcol <- '#FF000080'
     if(missing(pcol)) pcol <- '#00000080'
     if(missing(pcex)) pcex <-  0.5
     if(missing(lwd))  lwd  <-  1
     if(title) main <- dimnames(m)[[2]][1] else main <- ''
     plot(x, se=F, col=lcol, lwd=lwd, main=main, ...)
     points(xx, pch=16, col=pcol, cex=pcex, ...)
}

###   sensitivity functions   #####
#' @export
#' @rdname gamfit
`gamfit_sens` <- function(x, env, nrep=99, perc=5, ...){
     x <- x$sumtab
     # add noise then fit GAM
     `f` <- function(...){
          `anon` <- function(perc=perc, ...){
               nenv   <- env
               nenv[] <- sapply(nenv, noise, perc=perc)
               nenv[] <- sapply(nenv, standardize)
               xx     <- gamfit(m, nenv, ...)$sumtab
          }
          anon(...)
     }
     # iterate
     r <- array(NA, dim=c(dim(x),nrep)) # initialize
     for(i in 1:nrep){
               cat('round',i,'of',nrep,'\n')
          r[,,i] <- f(perc=perc)
     }
     r <- aperm(r, c(3,1,2) ) # rearrange
     dimnames(r)[[1]] <- 1:nrep              # rows are reps
     dimnames(r)[[2]] <- dimnames(env)[[2]]  # columns are env
     dimnames(r)[[3]] <- c('pval','adjr2')
     class(r) <- 'gamfit_sens'
     r
}
#' @export
#' @rdname gamfit
`plot.gamfit_sens` <- function(r, stat='r2', pltype='joy', ...){
     stat <- pmatch(stat, c('r2','pval'))
     pltype <- pmatch(pltype, c('heat','joy'))
     xlab <- c('Adj-R2', 'p-value')[stat]
     if (stat==1){
          if(pltype==1){
               plot_joy(r[,,2], ypad=1.01, xlab=xlab)
          } else {
               plot_heatmap(r[,,2], xord=F)
          }
     }
     if (stat==2){
          if(pltype==1){
               plot_joy(r[,,1], ypad=1.01, xlab=xlab)
          } else {
               plot_heatmap(r[,,1]<0.05, xord=F)
          }
     }
}
#' @export
#' @rdname gamfit
`summary.gamfit_sens` <- function(r, x, ...){
     x <- x$sumtab
     apply(r[,,2], 2, mean)                # mean fit to env
     devn <- abs(sweep(r[,,2], 2, x[,2])) # abs devns from observed
     MAE  <- apply(devn, 2, mean)          # unstandardized ! MAE
     MAE
}
### unexported
`noise` <- function(z, perc=0.0001, unif=TRUE){
     ### add random uncertainty, bounded by some % of the range
     if(is.vector(z)) wasvec <- TRUE else wasvec <- FALSE
     if(!is.matrix(z)) z <- as.matrix(z)
     nr  <- dim(z)[1]
     rng <- diff(range(z, na.rm=T))*perc/100
     n   <- prod(dim(z))
     if(unif){
          zr <- matrix(runif(n,rng*(-1),rng),nrow=nr)
     }else{
          zr <- matrix(rnorm(n,0,rng),nrow=nr)
     }
     out <- zr + z
     if(wasvec) as.vector(out) else out
}
