
# This code builds the IPM for three different assumptions of below-dam penalty on early mortality
# (none, +50%, +100%) and a set of combinations of decreases in harvest mortality 
# and dam passage mortality of above-dam spawners (= background mortality of smaller than
# average above-dam spawners). 

# The relevant sections in the paper are:
# - Methods: 2.3.4
# - Results: 3.3
# - Figures: Figures 5 a) & S1.5

# NOTE: 
# Due to the large number of high-dimensional matrices created and analysed here,
# this code takes a substantial amount of time to run. 

#------
# SETUP
#------

## Execute code for IPM building blocks
source('HunderTroutIPM_IPM_Building.R')

## Set blackbox parameters
m0 <- -log(0.082) 
Mu.mj <- -log(0.353)
Mu.mdam <- -log(0.75)
Mu.mOs <- 0.62
beta2.mj <- -0.2
beta2.mdam <- 0.05
beta2.mOs <- -0.7


#-----------------------
# BASELINE IPMs - LAMBDA
#-----------------------

## ORIGINAL
IPM1.orig <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.orig <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.orig <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM


#-----------------------------------------------
# CONTINUOUS MITIGATION SCENARIOS - CALCULATIONS
#-----------------------------------------------

## Write function to calculate lambda and store in a data frame
pop.analysis <- function(t, sex, orig, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0, dPenalty.j, l.limit, u.limit, mH.factor, mO1size.factor){
	
	# Build IPM
	IPM <- build.IPM(t, sex, orig, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0, dPenalty.j, l.limit, u.limit, mH.factor, mO1size.factor)$IPM
	
	# Eigen analysis
	output <- wvlambda(IPM)
	
	# Asymptotic population growth rate	
	lambda <- output$lambda
	
	# Collate results
	results <- data.frame(mH.factor = mH.factor, mO1size.factor = mO1size.factor, lambda = lambda)
	
	return(results)
	
}

## Make different perturbation combinations
mH.fac <- seq(1, 0, by = -1/50)
mO.fac <- seq(1, 0, by = -1/50)

pert.params <- expand.grid(mH.fac, mO.fac)

## Calculate lambda for all possible combinations
test1.X <- do.call("rbind", sapply(1:nrow(pert.params), FUN = function(z) pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=pert.params[z,1], mO1size.factor=pert.params[z,2]), simplify = FALSE))

res1.mat <- matrix(test1.X$lambda, nrow = length(mH.fac), ncol = length(mO.fac), dimnames = list(mH.fac, mO.fac))

test2.X <- do.call("rbind", sapply(1:nrow(pert.params), FUN = function(z) pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=pert.params[z,1], mO1size.factor=pert.params[z,2]), simplify = FALSE))

res2.mat <- matrix(test2.X$lambda, nrow = length(mH.fac), ncol = length(mO.fac), dimnames = list(mH.fac, mO.fac))

test3.X <- do.call("rbind", sapply(1:nrow(pert.params), FUN = function(z) pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=pert.params[z,1], mO1size.factor=pert.params[z,2]), simplify = FALSE))

res3.mat <- matrix(test3.X$lambda, nrow = length(mH.fac), ncol = length(mO.fac), dimnames = list(mH.fac, mO.fac))


#################################
# UNCHANGED SMOLT DAM MORTALITY #
#################################

library(viridis)

# Prepare vectors for labelling
perc.dec <- c('0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%')
prop.ind <- rev(seq(1, 51, 5))

# Set plotting limits
llim <- 0.706
ulim <- 1.154

# Plot
par(mar=c(5.1,5.1,4.1,2.1))
image.plot(res1.mat, col = magma(100), zlim = c(llim, ulim), main = expression('Asymptotic growth rate'~lambda~'(no penalty)'), xlab = 'Decrease in harvest mortality', ylab = '', xaxt= "n", yaxt= "n", cex.lab = 1.3, cex.main = 1.5, axis.args = list(cex.axis = 1.2))
title (ylab = 'Decrease in dam mortality (small spawners)', line=3.8, cex.lab=1.3)
axis(1, at = mH.fac[prop.ind], labels = perc.dec, cex.axis = 1.2)
axis(2, at = mO.fac[prop.ind], labels = perc.dec, las = 1, cex.axis = 1.2)
contour(res1.mat, add = TRUE, levels = c(0.8, 0.9, 1.0, 1.1), col = 'white', labcex = 1, method = 'edge', lty = c(2,2,1,2), crt = 90, drawlabels = F)

par(mar=c(5.1,5.1,4.1,2.1))
image.plot(res2.mat, col = magma(100), zlim = c(llim, ulim),  main = expression('Asymptotic growth rate'~lambda~'(50% penalty)'), xlab = 'Decrease in harvest mortality', ylab = '', xaxt= "n", yaxt= "n", cex.lab = 1.3, cex.main = 1.5, axis.args = list(cex.axis = 1.2))
title (ylab = 'Decrease in dam mortality (small spawners)', line=3.8, cex.lab=1.3)
axis(1, at = mH.fac[prop.ind], labels = perc.dec, cex.axis = 1.2)
axis(2, at = mO.fac[prop.ind], labels = perc.dec, las = 1, cex.axis = 1.2)

par(mar=c(5.1,5.1,4.1,2.1))
image.plot(res3.mat, col = magma(100), zlim = c(llim, ulim),  main = expression('Asymptotic growth rate'~lambda~'(100% penalty)'), xlab = 'Decrease in harvest mortality', ylab = '', xaxt= "n", yaxt= "n", cex.lab = 1.3, cex.main = 1.5, axis.args = list(cex.axis = 1.2))
title (ylab = 'Decrease in dam mortality (small spawners)', line=3.8, cex.lab=1.3)
axis(1, at = mH.fac[prop.ind], labels = perc.dec, cex.axis = 1.2)
axis(2, at = mO.fac[prop.ind], labels = perc.dec, las = 1, cex.axis = 1.2)
contour(res3.mat, add = TRUE, levels = c(0.8, 0.9, 1.0, 1.1), col = 'white', labcex = 1, method = 'edge', lty = c(2,2,1,2), crt = 90, drawlabels = F)

