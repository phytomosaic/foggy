###  SYNCSA tutorial   ############################################

require(SYNCSA)
require(foggy)

### usage
# data(flona)
# (res <- syncsa(comm   = flona$community,
#                traits = flona$traits,
#                phylo  = flona$phylo,
#                envir  = flona$environment,
#                ro.method = 'procrustes',
#                permu  = 99))

### data load
data(veg)
spe <- veg$spe
env <- veg$env
tra <- veg$tra
phy <- veg$phy

### SYNCSA
(res <- syncsa(comm   = spe,
               traits = tra,
               phylo  = Dp,
               envir  = env,
               ro.method = 'procrustes',
               permu  = 999))
names(res$matrices)

PP <- matrix.p(spe, Dp)
# matrix.w = Stdzed community matrix, rows = SUs, cols = species.
# matrix.q = Stdzed matrix containing degree of belonging of species.
# matrix.P = Phylogeny-weighted species composition matrix.

XX <- matrix.x(spe, tra) # ! forces gower or euc
# matrix.w = Stdzed community mx, rows = SUs, cols = species.
# matrix.u = Stdzed matrix containing degree of belonging of species.
# matrix.X = Trait-weighted species composition mx.

TT <- matrix.t(spe, tra)
# matrix.w = Standardized community mx, where rows are communities and columns species. Row totals (communities) = 1.
# matrix.b = Matrix of traits, exactly the same data input.
# matrix.T = Matrix containing trait averages at community level. If Scale = TRUE the matrix T is standardized within the traits.

QQ <- belonging(Dp, standardize = TRUE)
identical(QQ,  PP$`matrix.q`)           # TRUE, phylo degree of belonging
QQQ <- belonging(as.matrix(FD::gowdis(cwm)), standardize = TRUE)
identical(QQQ,  XX$`matrix.u`)          # FALSE (bc forces gower/euc)
identical(spe, PP$`matrix.w`)           # FALSE (bc row-stdzed)
identical(XX$`matrix.w`, PP$`matrix.w`) # TRUE, Stdzed comm mx
identical(XX$`matrix.u`, PP$`matrix.q`) # FALSE (bc trait vs phylo)

names(res$matrices)
mantel( dist(res$matrices$E), dist(res$matrices$T))
protest(dist(res$matrices$E), dist(res$matrices$T))
