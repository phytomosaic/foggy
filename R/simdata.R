#' @name simdata
#' @title Simulated data for lichens and plants
#' @aliases simdata
#' @docType data
#' @description
#' This data set gives information about species, environment, traits,
#'      and phylogeny for simulated replicate communities of lichens
#'      and plants.
#'
#' @format A list of lists each containing:\cr
#'     - \code{spe} 33 observations of 50 species,\cr
#'     - \code{env} 33 observations of 2 enviro variables,\cr
#'     - \code{tra} 3 traits for 50 species,\cr
#'     - \code{phy} phylogeny for 50 species.
#'
#' @details
#' Simulated with R package \code{ecolottery} per Munoz et al. (2018).
#'
#' @references
#' Munoz, F., M. Grenié, P. Denelle, A. Taudière, F. Laroche,
#'      C. Tucker, and C. Violle. 2018. ecolottery: Simulating and
#'      assessing community assembly with environmental filtering and
#'      neutral dynamics in R. Methods in Ecology and Evolution
#'      9:693–703.
#'
#' @examples
#' data(simdata)
#' x <- simdata[[1]]
#' plot_heatmap(x$spe, xord=FALSE)
#' plot_heatmap(x$env, xord=FALSE)
#' plot_heatmap(x$tra, xord=FALSE)
#' plot(x$phy, cex=0.7)
#'
#' @keywords datasets
"simdata"
