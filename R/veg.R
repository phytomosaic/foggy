#' @name veg
#' @title Mafragh, Algeria vegetation data
#' @aliases veg
#' @docType data
#' @description
#' This data set gives information about spatial coordinates, species,
#'     environment, traits, and phylogeny.
#'
#' @format A list containing:\cr
#'     - \code{xy}  97 observations of 2 spatial coordinates,\cr
#'     - \code{spe} 97 observations of 56 plant species,\cr
#'     - \code{env} 347 observations of 11 soil enviro variables,\cr
#'     - \code{tra} 12 traits for 56 plant species,\cr
#'     - \code{phy} phylogeny for 56 plant species.
#'
#' @details
#' See \code{?mafragh} in \code{ade4}.
#'
#' @source Pavoine et al. (2011).
#'
#' @references
#' Pavoine, S., Vela, E., Gachet, S., de Bélair, G. and Bonsall, M. B.
#'      2011. Linking patterns in phylogeny, traits, abiotic variables
#'      and space: a novel approach to linking environmental filtering
#'      and plant community assembly. Journal of Ecology 99:165–175.
#'      doi:10.1111/j.1365-2745.2010.01743.x
#'
#' @examples
#' \dontrun{?mafragh} # in package 'ade4'
#' data(veg)
#'
#' @keywords datasets
"veg"
