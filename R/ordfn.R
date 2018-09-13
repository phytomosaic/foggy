#' @title Nonmetric multidimensional scaling helper
#'
#' @description Light wrapper around \code{metaMDS} with some
#'     pre-specified arguments.
#'
#' @param spe species abundance matrix
#'
#' @param distance dissimilarity index, per \code{vegdist}
#'
#' @param k number of ordination axes sought
#'
#' @param ... additional arguments passed to \code{metaMDS}
#'
#' @return
#' Object of class \code{metaMDS}.
#'
#' @details
#' Very light wrapper around \code{metaMDS}.
#'
#' @examples
#' data(veg) # load('./data/veg.rda', verbose=T)
#' spe <- veg$spe
#' (m <- ordfn(spe, dist='bray', k=2))
#' plot(m)
#'
#' @seealso
#' \code{metaMDS}
#'
#' @export
#' @rdname ordfn
`ordfn` <- function(spe, distance, k, ...){
     vegan::metaMDS(comm=spe, distance=distance, k=k,
                    try=100, trymax=200, autotransform=F, trace=0,
                    maxit=500, weakties=TRUE, ...)
}
