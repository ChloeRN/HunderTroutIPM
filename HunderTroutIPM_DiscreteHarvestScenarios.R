# This code builds the IPM for three different assumptions of below-dam penalty on early mortality
# (none, +50%, +100%) and analyses it / projects it forward in time for different sets of harvest 
# rules. 

# The relevant sections in the paper are:
# - Methods: 2.3.4
# - Results: 3.3
# - Figures: Figures 4 and S1.4

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


#---------------------------------------
# DISCRETE MITIGATION SCENARIOS - LAMBDA
#---------------------------------------

## ORIGINAL
IPM1.orig <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.orig <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.orig <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM


## A) NO HARVEST
IPM1.A <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=0, mO1size.factor=1)$IPM

IPM2.A <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=0, mO1size.factor=1)$IPM

IPM3.A <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=0, mO1size.factor=1)$IPM


## B) NO HARVEST < 500mm
IPM1.B <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=500, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.B <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=500, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.B <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=500, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

## C) NO HARVEST > 700mm
IPM1.C <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=700, mH.factor=1, mO1size.factor=1)$IPM

IPM2.C <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=700, mH.factor=1, mO1size.factor=1)$IPM

IPM3.C <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=700, mH.factor=1, mO1size.factor=1)$IPM


## D) HARVEST ONLY 500-700mm
IPM1.D <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=500, u.limit=700, mH.factor=1, mO1size.factor=1)$IPM

IPM2.D <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=500, u.limit=700, mH.factor=1, mO1size.factor=1)$IPM

IPM3.D <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=500, u.limit=700, mH.factor=1, mO1size.factor=1)$IPM


## E) NO DAM MORTALITY OF SMOLT
IPM1.E <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam=0, Mu.mOs, beta2.mj, beta2.mdam=0, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.E <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam=0, Mu.mOs, beta2.mj, beta2.mdam=0, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.E <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam=0, Mu.mOs, beta2.mj, beta2.mdam=0, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

## LAMBDAS
wvlambda(IPM1.orig)$lambda
wvlambda(IPM2.orig)$lambda
wvlambda(IPM3.orig)$lambda
wvlambda(IPM1.A)$lambda
wvlambda(IPM2.A)$lambda
wvlambda(IPM3.A)$lambda
wvlambda(IPM1.B)$lambda
wvlambda(IPM2.B)$lambda
wvlambda(IPM3.B)$lambda
wvlambda(IPM1.C)$lambda
wvlambda(IPM2.C)$lambda
wvlambda(IPM3.C)$lambda
wvlambda(IPM1.D)$lambda
wvlambda(IPM2.D)$lambda
wvlambda(IPM3.D)$lambda
wvlambda(IPM1.E)$lambda
wvlambda(IPM2.E)$lambda
wvlambda(IPM3.E)$lambda

#--------------------------------------------
# DISCRETE MITIGATION SCENARIOS - PROJECTIONS
#--------------------------------------------

## Build IPM for stocked fish
IPM1.stock <- build.IPM(t=test.year, sex=1, orig=2, m0, Mu.mj=-log(0), Mu.mdam, Mu.mOs, beta2.mj=0, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.stock <- build.IPM(t=test.year, sex=1, orig=2, m0, Mu.mj=-log(0), Mu.mdam, Mu.mOs, beta2.mj=0, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.stock <- build.IPM(t=test.year, sex=1, orig=2, m0, Mu.mj=-log(0), Mu.mdam, Mu.mOs, beta2.mj=0, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

## Function for manual projection (with stocking & mitigation measures)
IPM.project3 <- function(IPM.wild, IPM.stock, IPM.mit, initSSD, steps, stock){
	
	wild.pop <- stock.pop <- matrix(NA, nrow = length(initSSD), ncol = steps*2)
	wild.pop[,1] <- initSSD
	stock.pop[,1] <- stock
	
	for(t in 2:steps){
	  
	  # Project stocked population
	  stock.pop[,t] <- IPM.stock %*% stock.pop[,t-1] + stock
	  
	  # Determine offspring produced by stocked population
	  stock.off <- rep(0, nrow(wild.pop))
	  stock.off[1:600] <- stock.pop[1:600,t]
	  
	  # Project wild population & add offspring of stocked fish
	  wild.pop[,t] <- IPM.wild %*% wild.pop[,t-1] + stock.off
	}

	for(t in (steps+1):(steps*2)){
	  
	  # Project stocked population
	  stock.pop[,t] <- IPM.stock %*% stock.pop[,t-1]
	  
	  # Determine offspring produced by stocked population
	  stock.off <- rep(0, nrow(wild.pop))
	  stock.off[1:600] <- stock.pop[1:600,t]
	  
	  # Project wild population & add offspring of stocked fish
	  wild.pop[,t] <- IPM.mit %*% wild.pop[,t-1] + stock.off
	}
	
	# Remove juveniles from stocked population
	stock.pop[1:600,] <- 0
	
	return(list(wild.pop = wild.pop, stock.pop = stock.pop, tot.pop = wild.pop+stock.pop))
}

## Running projection for different harvest scenarios
set.seed <- 12
test.stock <- stock.fun(test.year, Mu.mdam, beta2.mdam)
sSSD <- rep(1, 1800)

run.orig <- IPM.project3(IPM1.orig, IPM1.stock, IPM1.orig, sSSD, 200, test.stock)
run.A <- IPM.project3(IPM1.orig, IPM1.stock, IPM1.A, sSSD, 200, test.stock)
run.B <- IPM.project3(IPM1.orig, IPM1.stock, IPM1.B, sSSD, 200, test.stock)
run.C <- IPM.project3(IPM1.orig, IPM1.stock, IPM1.C, sSSD, 200, test.stock)
run.D <- IPM.project3(IPM1.orig, IPM1.stock, IPM1.D, sSSD, 200, test.stock)

run2.orig <- IPM.project3(IPM2.orig, IPM2.stock, IPM2.orig, sSSD, 200, test.stock)
run2.A <- IPM.project3(IPM2.orig, IPM2.stock, IPM2.A, sSSD, 200, test.stock)
run2.B <- IPM.project3(IPM2.orig, IPM2.stock, IPM2.B, sSSD, 200, test.stock)
run2.C <- IPM.project3(IPM2.orig, IPM2.stock, IPM2.C, sSSD, 200, test.stock)
run2.D <- IPM.project3(IPM2.orig, IPM2.stock, IPM2.D, sSSD, 200, test.stock)

run3.orig <- IPM.project3(IPM3.orig, IPM3.stock, IPM3.orig, sSSD, 200, test.stock)
run3.A <- IPM.project3(IPM3.orig, IPM3.stock, IPM3.A, sSSD, 200, test.stock)
run3.B <- IPM.project3(IPM3.orig, IPM3.stock, IPM3.B, sSSD, 200, test.stock)
run3.C <- IPM.project3(IPM3.orig, IPM3.stock, IPM3.C, sSSD, 200, test.stock)
run3.D <- IPM.project3(IPM3.orig, IPM3.stock, IPM3.D, sSSD, 200, test.stock)

## Plot projections
par(mar=c(5.1,4.1,3.1,2.1))
plot(log(colSums(run.orig$tot.pop[,150:270])), type = 'l', ylab = 'Log population size (age 1+)', xlab = 'Number of years', lwd = 2, ylim = c(0, 20), cex.axis = 1.2, cex.lab = 1.3)
lines(c(51:120), log(colSums(run.A$tot.pop[,201:270])), col = 'black', lwd = 2, lty = 3)
lines(c(51:120), log(colSums(run.B$tot.pop[,201:270])), col = 'black', lwd = 2, lty = 2)
lines(c(51:120), log(colSums(run.C$tot.pop[,201:270])), col = '#F3626F', lwd = 2)
lines(c(51:120), log(colSums(run.D$tot.pop[,201:270])), col = '#F3626F', lwd = 2, lty = 2)
abline(v = 51, col = 'grey', lwd = 1.5, lty = 2)
legend(0, 8, legend = c('none', 'No harvest', 'No harvest < 500mm', 'No harvest > 700mm', 'Harvest only for 500-700mm'), col = c(1, 1, 1, '#F3626F', '#F3626F'), lty = c(1,3,2,1,2), lwd = c(2,2,2,2,2), bty = 'n', title = expression(bold('Mitigation measure')), title.adj = 0.1, cex = 1.2)


par(mfrow = c(2,1),mar=c(5.1,4.1,3.1,2.1))
plot(log(colSums(run2.orig$tot.pop[,150:270])), type = 'l', ylab = 'Log population size (age 1+)', xlab = 'Number of years', lwd = 2, ylim = c(-5, 20), cex.axis = 1.2, cex.lab = 1.3)
lines(c(51:120), log(colSums(run2.A$tot.pop[,201:270])), col = 'black', lwd = 2)
lines(c(51:120), log(colSums(run2.B$tot.pop[,201:270])), col = 'black', lwd = 2, lty = 2)
lines(c(51:120), log(colSums(run2.C$tot.pop[,201:270])), col = '#F3626F', lwd = 2)
lines(c(51:120), log(colSums(run2.D$tot.pop[,201:270])), col = '#F3626F', lwd = 2, lty = 2)
abline(v = 51, col = 'grey', lwd = 1.5, lty = 2)
legend(0, 5, legend = c('none', 'No harvest', 'No harvest < 500mm', 'No harvest > 700mm', 'Harvest only for 500-700mm'), col = c(1, 1, 1, '#F3626F', '#F3626F'), lty = c(1,3,2,1,2), lwd = c(2,2,2,2,2), bty = 'n', title = expression(bold('Mitigation measure')), title.adj = 0.1, cex = 1.2)
title('a) Below-dam penalty: +50% early mortality', adj = 0)

plot(log(colSums(run3.orig$tot.pop[,150:270])), type = 'l', ylab = 'Log population size (age 1+)', xlab = 'Number of years', lwd = 2, ylim = c(-5, 20), cex.axis = 1.2, cex.lab = 1.3)
lines(c(51:120), log(colSums(run3.A$tot.pop[,201:270])), col = 'black', lwd = 2)
lines(c(51:120), log(colSums(run3.B$tot.pop[,201:270])), col = 'black', lwd = 2, lty = 2)
lines(c(51:120), log(colSums(run3.C$tot.pop[,201:270])), col = '#F3626F', lwd = 2)
lines(c(51:120), log(colSums(run3.D$tot.pop[,201:270])), col = '#F3626F', lwd = 2, lty = 2)
abline(v = 51, col = 'grey', lwd = 1.5, lty = 2)
title('b) Below-dam penalty: +100% early mortality', adj = 0)


