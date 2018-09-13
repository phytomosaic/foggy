# ######################################################################
# # Multiple approaches: RLQ, fourthcorner, and gradient regression
# #  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 12 Sep 2018
# ## CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
#
# require(foggy)
#
# ### data load
# data(veg)
# spe <- veg$spe
# env <- veg$env
# tra <- veg$tra
# phy <- veg$phy
#
# ### split species/traits into two (to simulate 'plant vs lichen')
# ###    (break is at deepest phylogenetic node)
# i1   <- 1:22
# i2   <- 23:56
# spe1 <- mx_dropzero(spe[,i1])
# spe2 <- mx_dropzero(spe[,i2])
# tra1 <- tra[rownames(tra) %in% names(spe1),]
# tra2 <- tra[rownames(tra) %in% names(spe2),]
#
# ### zero-adjustment for empty SUs
# spe1 <- cbind(spe1, dummy=rep(1, nrow(spe1)))
# spe2 <- cbind(spe2, dummy=rep(1, nrow(spe2)))
# tra1 <- rbind(tra1, dummy=colMeans(tra))
# tra2 <- rbind(tra2, dummy=colMeans(tra))
#
# ### transformations
# spe1 <- data.frame(genlogtrans(spe1))
# spe2 <- data.frame(genlogtrans(spe2))
# env  <- data.frame(decostand(scale(env, center=F), 'range'))
# tra1 <- data.frame(decostand(tra1, 'range'))
# tra2 <- data.frame(decostand(tra2, 'range'))
# cwm1 <- data.frame(makecwm(spe1, tra1))
# cwm2 <- data.frame(makecwm(spe2, tra2))
#
# ### NMDS ordinations
# m1 <- ordfn(spe1, 'bray', 2)
# m2 <- ordfn(spe2, 'bray', 2)
# m3 <- ordfn(cwm1, 'altGower', 2)
# m4 <- ordfn(cwm2, 'altGower', 2)
# set_par(4)
# plot(m1)
# plot(m2)
# plot(m3)
# plot(m4)
#
# ### procrustes
# p1 <- protest(m1,m2) # plant vs lichen species
# p2 <- protest(m3,m4) # plant vs lichen CWM traits
# p1$t0  # plant vs lichen species
# p2$t0  # plant vs lichen CWM traits
#
# ### fit GAMs to regress each enviro variable on NMDS scores
# # envfit(m1, env, perm=999)    # LINEAR correlations
# gamfit(m1, env)  # lichen species fit to env
# (g1 <- gamfit(m1, env))  # lichen species fit to env
# gamfit(m2, env)  #  plant species fit to env
# gamfit(m3, env)  # lichen traits CWM fit to env
# gamfit(m4, env)  #  plant traits CWM fit to env
#
# ### fourthcorner
# (f1 <- fourthcorner(env, spe1, tra1, nrepet=999, modeltype=6))
# (f2 <- fourthcorner(env, spe2, tra2, nrepet=999, modeltype=6))
#
# ### RLQ and fourthcorner.rlq
# # r1 <- rlqfn(spe1, env, tra1)
# # r2 <- rlqfn(spe2, env, tra2)
# # fc_envfit(r1)
# # fc_envfit(r2)
#
#
# ###   END   ########################################################
