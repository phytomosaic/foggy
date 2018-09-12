#' @name pillar
#' @title Macroinvertebrates of Loire River data
#' @aliases pillar
#' @docType data
#' @description
#' This data set gives information about plant species,
#'     environment, traits, and phylogeny from Eldorado do Sul,
#'     Brazil.
#'
#' @format A list containing:\cr
#'     - \code{spe} 14 observations of 81 plant species,\cr
#'     - \code{env} 14 observations of 1 env variable (nitrogen),\cr
#'     - \code{tra} 8 traits for 81 species,\cr
#'     - \code{phy} phylogeny for 81 species, from \code{Dp},\cr
#'     - \code{Dp} phylogenetic distances for 81 species.
#'
#' @details
#' From the original paper: "an experiment evaluating the effect of
#'     N-fertilizer and grazing levels on natural grassland, located
#'     in Eldorado do Sul, Brazil...  Fourteen experimental plots were
#'     subjected during 5 years to limited combinations of
#'     N-fertilizer (0, 30, 100, 170 and 200 kg N ha-1 year-1)...
#'
#' @source Originally from Pillar and Duarte (2010), but here sourced
#'     from https://github.com/SrivastavaLab/syncsa.
#'
#' @references
#' Pillar, V. D., and L. d. S. Duarte. 2010. A framework for
#'     metacommunity analysis of phylogenetic structure. Ecology
#'     Letters 13:587–596.
#'
#' @examples
#' data(pillar)
#' plot_heatmap(spe, logbase=10)
#'
#' @keywords datasets
"pillar"
