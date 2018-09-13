#' @name invert
#' @title Macroinvertebrates of Loire River data
#' @aliases invert
#' @docType data
#' @description
#' This data set gives information about invertebrate species,
#'     environment, traits, and phylogeny.
#'
#' @format A list containing:\cr
#'     - \code{spe} 38 observations of 40 invertebrate species,\cr
#'     - \code{env} 38 observations of 3 aquatic enviro variables,\cr
#'     - \code{tra} 11 traits for 40 species,\cr
#'     - \code{phy} phylogeny for 40 species.
#'
#' @details
#' See \code{?macroloire} in \code{ade4}.
#'
#' @source Ivol et al. (1997).
#'
#' @references
#' Ivol, J. M., Guinand, B., Richoux, P. and Tachet, H. 1997.
#'      Longitudinal changes in Trichoptera and Coleoptera assemblages
#'      and environmental conditions in the Loire River (France).
#'      Archiv for Hydrobiologie 138: 525–557.
#'
#' @examples
#' \dontrun{?macroloire} # in package 'ade4'
#' data(invert)
#'
#' @keywords datasets
"invert"
