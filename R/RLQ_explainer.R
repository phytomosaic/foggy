# By: Nathan Lemoine
# https://www.r-bloggers.com/r-for-ecologists-rlq-analysis-semi-explained/

# RLQ analysis is a method by which one can uncover how the
# environment filters certain species traits. For example, you can
# determine whether a particular environment selects for species with
# rapid growth rates, high reproductive output, or whatever trait you
# choose to measure. It accomplishes this by, more or less, linking a
# description of the environment to species traits by measurements of
# species abundances.
# You start with three data tables: The R matrix
# is a site x environment table: sites are rows and columns are
# environmental descriptors. The L matrix is a site x species table,
# where rows are sites and columns are abundances of specific species.
# The Q matrix is a species x trait table, where rows are species and
# columns are biological traits of those species. What RLQ analysis
# does, simply, is makes a new matrix that I’ll call V, which is a
# environment x traits matrix, and you then perform your standard
# PCA-esque eigendecomposition on that. Although not technically
# correct, RLQ analysis is simply thought of as nothing more than PCA
# on the matrix V. The vast majority of the work is actually in
# constructing V. The ‘ade4′ package in R can do RLQ analyses. First,
# you do a principle components analysis on both R and Q and a
# correspondance analysis on L. You then pass these analyses to the
# rlq() function. To figure out how RLQ works, I took apart the rlq()
# function, then several secondary functions called by rlq(). It turns
# out, the PCA and CA on the R, L, and Q matrices aren’t actually
# used, you can do the whole thing by hand without those preliminary
# analyses.
# results verified with rlq() function

require(ade4)
require(ecole)
#### READ IN THE DATA AND CLEAN IT UP ####
data(mafragh)
d <- mafragh
env <- d$env
colnames(env) <- ecole::clean_text(colnames(env), lower=TRUE)
spe <- d$flo
colnames(spe) <- ecole::clean_text(d$spenames$scientific, lower=TRUE)
tra <- d$traits
tra <- data.frame(tra[[1]],tra[[2]],tra[[3]],tra[[4]])

# RLQ analysis operates on the matrix RLQ, which is calculated as
# R‘ D_site L D_spe Q, where D_site and D_spe are diagonal matrices of
# the row and column weights from the species matrix. As shown below,
# this is the same as R’ P Q, where P is the centered probability mx.

# First, we have a site by species matrix, N, of raw abundances
N <- spe

# Convert N to relativized species matrix P, where p_ij =
# n_ij / n++, where n++ is the total n of individuals in matrix N.
P <- N/sum(N)

# Now divide each observation by its row weight (p_i+) and column
# weight (p_j+). The row and column weights are simply the sum of the
# observations in a row divided by the sum of the matrix (and
# similarly for columns). This gives p_ij = p_ij/(p_i+ p+j+)
row.w <- apply(P, 1, function(x) sum(x)/sum(P))
col.w <- apply(P, 2, function(x) sum(x)/sum(P))
P <- sweep(P, 1, row.w, '/')
P <- sweep(P, 2, col.w, '/')

# Next, subtract 1 from each observation, giving p_ij =
# p_ij(p_i+ p_j+) – 1, which equals (p_ij – p_i+p_j+)/(p_i+p_j+).
P <- P-1

# This IS the chi-distance matrix used in correspondance analysis.
# You can verify this by checking the table from dudi.coa function.

# However, we only want the centered matrix, we need to remove the
# weights in the denominator. Create diagonal matrices D_site and
# D_spe of the row and column weights respectively. Then pre- and
# post-multiply the matrix P. This will yield a matrix L, where
# l_ij = p_ij – p_i+ p_j+
D_site <- diag(row.w)
D_spe  <- diag(col.w)
L      <- D_site %*% as.matrix(P) %*% D_spe

# Now make the R’LQ matrix. First, center and standardize the columns
# of R and Q. The center is taken as the WEIGHTED average where the
# weights are the row weights (for the environment matrix) and species
#  weights (for the species matrix).
# Calculate the weighted average for each trait and site
traAvg  <- apply(tra, 2, function(x) sum(x*col.w)/sum(col.w))
envAvg  <- apply(env, 2, function(x) sum(x*row.w)/sum(row.w))
traCent <- sweep(tra, 2, traAvg)
envCent <- sweep(env, 2, envAvg)

# Calculate the weighted standard deviation. Since the values are now
# in deviations from the mean, the weighted variance is
# sum(x^2w) / sum(weights), and the std deviation is the sqrt of this.
traSD <- apply(traCent, 2, function(x)sqrt(sum(x^2*col.w)/sum(col.w)))
envSD <- apply(envCent, 2, function(x)sqrt(sum(x^2*row.w)/sum(row.w)))
traitScale <- sweep(traCent, 2, traSD, '/')
envScale   <- sweep(envCent, 2, envSD, '/')
R <- as.matrix(envScale)
Q <- as.matrix(traitScale)

# Next, V is just the R’ L Q product. This is actually the correlation
# matrix between the environment traits and the species traits,
# mediated by species abundances.
V <- t(R) %*% L %*% Q
round(V, 3)
# This is identical to the matrix operated on by the rlq() command.
# You can check this with by examining table ($tab) returned by rlq()
# Next, get the cross-product matrix because the correlation matrix is
# not guaranteed to be either square or symmetric:
Z <- crossprod(V, V)

# The rest is a standard PCA-like eigen decomposition of V.
eigVals <- eigen(Z)$values
eigVecs <- eigen(Z)$vectors
sum(eigVals)

## THE SPECIES TRAIT LOADINGS ON EACH AXIS ARE THE EIGENVECTORS
traitLoad <- data.frame(eigVecs)
rownames(traitLoad) <- colnames(V)
colnames(traitLoad) <- paste('Axis', 1:length(eigVals))

## THE ENVIRONMENTAL TRAIT SCORES ARE CALCULATED EXACTLY AS IN PCA
envscr <- data.frame(V %*% eigVecs)
names(envscr) <- paste('Axis', 1:length(eigVals))

plot(envscr)
plot(envscr[,1:2])

###   END   #########################################################
