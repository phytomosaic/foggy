#' @title Phylogenetic correction of traits
#'
#' @description Remove potential phylogenetic signal from a set of
#'     continuous or discrete traits.
#'
#' @param phy phylogenetic tree of class \code{'phylo'}
#'
#' @param tra matrix or data.frame, rows=species and cols=traits
#'
#' @param ... additional arguments (currently ignored)
#'
#' @return
#' Data frame of same dimensions as \code{tra}, but with
#'     phylogenetically-corrected values.
#'
#' @details
#' Removes phylogenetic signal from trait values, by applying
#'     correction factors obtained from Cholesky decomposition of the
#'     phylogenetic variance-covariance matrix (Butler et al. 2000).
#'     Current implementation follows Eklöf and Stouffer (2016) in
#'     permitting both continuous and categorical traits.
#'
#' @examples
#' data(veg) # load('./data/veg.rda', verbose=T)
#' phy <- veg$phy
#' tra <- veg$tra
#' ptra <- phylo_corr(phy, tra)
#' ecole::set_par(8)
#' for(i in 1:8){
#'      plot(tra[,i], ptra[,i], pch=16, cex=0.7, col='#00000050',
#'      xlab=dimnames(tra)[[2]][i], ylab=dimnames(ptra)[[2]][i])
#' }
#'
#' @references
#' Butler, M.A., T.W. Schoener, and J.B. Losos. 2000. The relationship
#'     between sexual size dimorphism and habitat use in Greater
#'     Antillean Anolis lizards. Evolution 54:259-272.\cr
#' Eklöf, A., and D.B. Stouffer. 2016. The phylogenetic component of
#'     food web structure and intervality. Theoretical Ecology
#'     9:107–115.\cr
#'
#' @seealso
#' Code adapted from: https://github.com/stoufferlab/phyloint/
#'
#' @export
#' @rdname phylo_corr
`phylo_corr` <- function(phy, tra, ...){

     # checks
     if (class(phy) != 'phylo') stop('phy must be of class `phylo`')
     if (!is.data.frame(tra)) tra <- as.data.frame(tra)
     rn <- dimnames(tra)[[1]] # species names
     cn <- dimnames(tra)[[2]] # traits names
     ntra <- dim(tra)[[2]]    # n of traits
     if (!identical(phy$tip.label, rn)) stop('species name mismatch')

     # calculate G matrix based on the phylogeny
     G <- ape::vcv.phylo(phy)

     # phylogenetic 'correction factor' per Butler et al. (2000)
     corfac <- chol(solve(G))

     # initialize U for 'corrected' trait values
     U <- as.data.frame(
          matrix(NA, nrow=dim(tra)[[1]], ncol=dim(tra)[[2]],
                 dimnames=list(rn,cn)))

     # clean up factors, and correct single-level factors
     tra <- droplevels(tra) # drop unused factor levels
     tra[,sapply(tra,function(x)nlevels(x)==1)] <- 1

     # corrections per individual trait
     for(j in 1:ntra){
          M <- model.matrix(as.formula(paste0("~0+",cn[j])),data=tra)
          corrtra <- data.frame(corfac %*% M)
          U[,j]   <- corrtra[rn, ,drop=FALSE]
     }
     return(U)
}
