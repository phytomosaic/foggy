#' @title Add heatmap colors to a perspective plot
#'
#' @description Helper function to add heatmap colors to a perspective
#'     plot.
#'
#' @param z  matrix of values
#'
#' @param pal character string, one of \code{'heat'} for heat colors,
#'     or \code{'redblue'} for diverging red-blue palette
#'
#' @param ...  additional arguments passed to \code{persp}
#'
#' @return
#' \code{persp} plot.
#'
#' @details
#' Plots a perspective plot with facets colored by z-value.
#'
#' @examples
#' # from the persp documentation:
#' x <- seq(-10, 10, length= 30)
#' y <- x
#' f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
#' z <- outer(x, y, f)
#' persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = 'lightblue')
#' persp_heat(z, theta = 30, phi = 30, expand = 0.5)
#' persp_heat(z, pal='redblue', theta = 30, phi = 30, expand = 0.5)
#'
#' @seealso
#' \code{\link[vegan]{ordiplot}}
#'
#' @export
#' @rdname persp_heat
### helper function to add heatmap colors to a perspective plot
`persp_heat` <- function(z, pal='heat', ...){
     nrz <- nrow(z) ; ncz <- ncol(z)
     zfacet <- z[-1,-1] + z[-1,-ncz] + z[-nrz,-1] + z[-nrz,-ncz]
     if(pal=='redblue'){
          cc <- colorRampPalette(c('red','transparent','blue'))(111)[
               as.numeric(cut(zfacet,breaks=111))]
     }
     if(pal=='heat'){
          cc <- heat.colors(100)[cut(zfacet, 100)]
     }
     persp(z, col=cc, ...)
}
