
# This code builds the IPM for three different assumptions of below-dam penalty on early mortality
# (none, +50%, +100%) and runs an elasticity analysis of equilibrium population size under stocking.

# The relevant sections in the paper are:
# - Methods: 2.3.3
# - Results: 3.2
# - Figures: Figure S1.3 a)

################
#### SETUP #####
################
 
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

## Draw a stocking vector for testing
set.seed <- 12
test.stock <- stock.fun(test.year, Mu.mdam, beta2.mdam)

## Re-define survival functions to allow for perturbation of mortality hazard rates

Sa.pred <- function(size, dam, t, orig, min.size, max.size, spawn.year, fdy){ # Added fdy here (has to be a vector of 2)
	
	# Calculate the vector of survival probabilities (whole size range)
	mTOT <- (mH.2yr.pred(size, t, orig, l.limit = Lx, u.limit = Ux, mH.factor = 1)*fdy[2]) + (mO.2yr.pred(size, dam, t, orig, mO1size.factor = 1)*fdy[1])	
	Sa <- exp(-mTOT)

	# Calculate minimum & maximum size survival
	mTOT.min <- (mH.2yr.pred(min.size, t, orig, l.limit = Lx, u.limit = Ux, mH.factor = 1)*fdy[2]) + (mO.2yr.pred(min.size, dam, t, orig, mO1size.factor = 1)*fdy[1])
	mTOT.max <- (mH.2yr.pred(max.size, t, orig, l.limit = Lx, u.limit = Ux, mH.factor = 1)*fdy[2]) + (mO.2yr.pred(max.size, dam, t, orig, mO1size.factor = 1)*fdy[1])
	Sa.min <- exp(-mTOT.min)
	Sa.max <- exp(-mTOT.max)
	
	# Replace survival probabilites for sizes beyond thresholds
	Sa[which(size < min.size)] <- Sa.min
	Sa[which(size > max.size)] <- Sa.max
	
	if(spawn.year==1){Sa.final <- Sa}
	if(spawn.year!=1){Sa.final <- rep(1, length(Sa))}
	return(Sa.final)
}

Ss.pred <- function(size, t, orig, Mu.mOs, beta2.mOs, fdy){ # Added fdy (should be a vector of 2)
	
	mOs <- mOs.pred(size, Mu.mOs, beta2.mOs)*fdy[1]
	mHs <- (mH.2yr.pred(size, t, orig, l.limit = Lx, u.limit = Ux, mH.factor = 1)*fdy[2])/2
	mTOT <- mHs+ mOs	
	Ss <- exp(-mTOT)
	return(Ss)
}

Sj.pred <- function(size, Mu.mj, beta2.mj, dPenalty.j, fdy){ # Added fdy (= 1 number)
	
	size.std <- (size-157)/50
	
	log.mj <- log(Mu.mj) + beta2.mj*size.std
	mj <- exp(log.mj)*dPenalty.j*fdy
	return(exp(-mj))
}

S0.pred <- function(m0, dPenalty.0, fdy){ # Added dy (should be 1 number)
	
	m0.new <- m0*dPenalty.0*fdy
	return(exp(-m0.new))
}

Sdam.pred <- function(size, Mu.mdam, beta2.mdam, fdy){ # Added dy (should be 1 number)
	
	size.std <- (size-250)/50
	
	log.mdam <- log(Mu.mdam) + beta2.mdam*size.std
	mdam <- exp(log.mdam)*fdy
	return(exp(-mdam))
}


## Re-define kernel functions to allow allow for perturbation of mortality hazard rates

#----------------------------------------------
# a) Juvenile to Juvenile (upriver & downriver)
#----------------------------------------------

K.JJ.fun  <-function(t, dam, Mu.mj, beta2.mj, dPenalty.j, fdy){ # Added dam argument here
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
    if(dam == 1){
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, 1, fdy[which(VR.names=='mj.u')])
  	}else{
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, dPenalty.j, fdy[which(VR.names=='mj.d')])
  } 
  Tvec <- 1-(pS.pred(x, t))

  for(i in 1:n){
    Smat[,i] <- Svec[i]*Tvec[i]*dist.grR.pred(x[i], x, t, 1)*dx
  }
  return(Smat)
}


#----------------------------------------------
# b) Juvenile to Subadult (upriver & downriver)
#----------------------------------------------

K.JS.fun  <-function(t, orig, dam, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j, fdy){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  if(dam == 1){
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, 1, fdy[which(VR.names=='mj.u')])*Sdam.pred(x, Mu.mdam, beta2.mdam, fdy[which(VR.names=='mdam')])
  }else{
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, dPenalty.j, fdy[which(VR.names=='mj.d')])
  }  
  Tvec <- pS.pred(x, t)
  
  for (i in 1:n){
    Smat[,i] <- Svec[i]*Tvec[i]*dist.grL.pred(x[i], x, t, 1, 0, orig, 0)*dx

  }
  return(Smat)
}


#------------------------
# c) Subadult to Subadult
#------------------------

K.SS.fun  <-function(t, sex, orig, Mu.mOs, beta2.mOs, fdy){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  Svec <- Ss.pred(x, t, orig, Mu.mOs, beta2.mOs, fdy[which(VR.names%in%c('mOs','mH'))])
  Tvec <- 1-(pM.pred(x, t, sex, orig))
  
  for (i in 1:n){
    Smat[,i] <- Svec[i]*Tvec[i]*dist.grL.pred(x[i], x, t, 1, 0, orig, 0)*dx
  }
  return(Smat)
}


#---------------------------------------------
# d) Subadult to Spawner (upriver & downriver)
#---------------------------------------------

K.SAs.fun  <-function(t, dam, sex, orig, Mu.mOs, beta2.mOs, fdy){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  Svec <- Ss.pred(x, t, orig, Mu.mOs, beta2.mOs, fdy[which(VR.names%in%c('mOs','mH'))])
  Tvec1 <- pM.pred(x, t, sex, orig)
  
  if(dam == 1){
  	Tvec2 <- pL.pred(x, t, orig)
  }else{
  	Tvec2 <- 1-pL.pred(x, t, orig)
  }
  
  # 1. Survive, mature, and grow (given current size, with reproduction cost)
  for (i in 1:n){
    Smat[,i] <- Svec[i]*Tvec1[i]*dist.grL.pred(x[i], x, t, 1, 1, orig, 1)
  }
  
  # 2. Ascend the ladder or not (given next size)
  for(j in 1:n){
		 Smat[j,] <- Smat[j,]*Tvec2[j]*dx
  }
  return(Smat)
}


#-----------------------------------------------------
# e) Spawner to Post-spawner (upriver & downriver)
#-----------------------------------------------------

K.AsAn.fun  <-function(t, dam, sex, orig, min.size, max.size, fdy){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  Tvec <- rep(1, n) # 1 because the transition is deterministic (and independent of the dam)
  
  if(dam == 1){
  	Svec <- Sa.pred(x, dam, t, orig, min.size, max.size, spawn.year = 1, fdy[which(VR.names%in%c('mO1','mH'))])
  }else{
  	 Svec <- Sa.pred(x, dam, t, orig, min.size, max.size, spawn.year = 1, fdy[which(VR.names%in%c('mO2','mH'))])
  }

  # Deterministic transition, survive, and grow (without reproduction cost)
  for (i in 1:n){
    Smat[,i] <- Tvec[i]*Svec[i]*dist.grL.pred(x[i], x, t, 1, 0, orig, 1)*dx
  }
  return(Smat)
}


#-----------------------------------------------------
# f) Post-spawner to Spawner (upriver & downriver)
#-----------------------------------------------------

K.AnAs.fun  <-function(t, dam, sex, orig, min.size, max.size, fdy){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  if(dam == 1){
  	Tvec <- pL.pred(x, t, orig)
  }else{
  	Tvec <- 1-pL.pred(x, t, orig)
  }
  
  Svec <- Sa.pred(x, dam, t, orig, min.size, max.size, spawn.year = 0, 0) # No fdy needed here because this always returns 1
  
  # 1. Survive and grow (given current size, with reproduction cost)
  for (i in 1:n){
    Smat[,i] <- Svec[i]*dist.grL.pred(x[i], x, t, 1, 1, orig, 1)
  }
  
  # 2. Ascend the ladder or not (given next size)
  for(j in 1:n){
    Smat[j,] <- Smat[j,]*Tvec[j]*dx
  }

  return(Smat)
}

#--------------------------------------------------
# g) Pre-spawners to Juvenile (upriver & downriver)
#--------------------------------------------------

K.AsJ.fun  <-function(t, m0, dam, orig, dPenalty.0, fdy){
  
  # Define matrix
  Bmat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  if(dam == 1){
  	Bvec <- fec.pred(x)*S0.pred(m0, 1, fdy[which(VR.names=='m0.u')])*0.5
  }else{
  	Bvec <- fec.pred(x)*S0.pred(m0, dPenalty.0, fdy[which(VR.names=='m0.d')])*0.5
  }
  
  # Reproduction, egg survival, growth to 1-year size
  for (i in 1:n){
    Bmat[,i] <- Bvec[i]*dist.sizeOff.pred(x, t, 1) * dx
  }
  return(Bmat)
}


#########################################
#### HAZARD RATE ELASTICITY ANALYSIS ####
#########################################

#-------------------
# PERTURBATION SETUP
#-------------------

## Names of all vital rates to include in analysis
VR.names <- c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH')

## Perturbation factor
DY <- 1e-5

## Collection of vectors for relative perturbation
fdy <- matrix(1, nrow = length(VR.names), ncol = length(VR.names))
diag(fdy) <- DY+1
fdy <- cbind(rep(1, length(VR.names)), fdy)

rownames(fdy) <- VR.names
colnames(fdy) <- c('orig', VR.names)

# --> Each column (or row) of this matrix corresponds to a small absolute perturbation in a given mortality hazard rate

## Separation into wild and stocked perturbation matrices
fdy.wild <- fdy
fdy.stock <- fdy

#-----------------------------------
# FUNCTION FOR POPULATION PROJECTION
#-----------------------------------

IPM.project2 <- function(IPM.wild, IPM.stock, initSSD, steps, stock){
	
	wild.pop <- stock.pop <- matrix(NA, nrow = length(initSSD), ncol = steps)
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
	
	# Remove juveniles from stocked population
	stock.pop[1:600,] <- 0
	
	return(list(wild.pop = wild.pop, stock.pop = stock.pop, tot.pop = wild.pop+stock.pop))
}


#----------------------------------
# FUNCTION FOR EQUILIBRIUM ANALYSIS
#----------------------------------

# The following function builds the IPM for a specific hazard rate perturbation given by the vector dy and calculates the following metrics at equilibrium (with stocking):
# 1) Population size
# 2) Proportions of the population in each stage
# 3) Mean size of individuals in each stage
# 4) Size variability (SD) of individuals in each stage

pop.analysis <- function(t, sex, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0, dPenalty.j, sSSD, N, test.stock, fdy.wild, fdy.stock){
	
	## Build separate kernels (wild & stocked)
	wK.JJ.u <- K.JJ.fun(t, 1, Mu.mj, beta2.mj, dPenalty.j, fdy.wild)
	wK.JJ.d <- K.JJ.fun(t, 0, Mu.mj, beta2.mj, dPenalty.j, fdy.wild)
	wK.JS.u <- K.JS.fun(t, orig=1, 1, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j, fdy.wild)
	wK.JS.d <- K.JS.fun(t, orig=1, 0, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j, fdy.wild)
	wK.SS <- K.SS.fun(t, sex, orig=1, Mu.mOs, beta2.mOs, fdy.wild)
	wK.SAs.u <- K.SAs.fun(t, 1, sex, orig=1, Mu.mOs, beta2.mOs, fdy.wild)
	wK.SAs.d <- K.SAs.fun(t, 0, sex, orig=1, Mu.mOs, beta2.mOs, fdy.wild)
	wK.AsJ.u <- K.AsJ.fun(t, m0, 1, orig=1, dPenalty.0, fdy.wild)
	wK.AsJ.d <- K.AsJ.fun(t, m0, 0, orig=1, dPenalty.0, fdy.wild)
	wK.AsAn.u <- K.AsAn.fun(t, 1, sex, orig=1, Lx, Ux, fdy.wild)
	wK.AsAn.d <- K.AsAn.fun(t, 0, sex, orig=1, Lx, Ux, fdy.wild)
	wK.AnAs.u <- K.AnAs.fun(t, 1, sex, orig=1, Lx, Ux, fdy.wild)
	wK.AnAs.d <- K.AnAs.fun(t, 0, sex, orig=1, Lx, Ux, fdy.wild)

	sK.JJ.u <- K.JJ.fun(t, 1, Mu.mj=-log(0), beta2.mj=0, dPenalty.j, fdy.stock)
	sK.JJ.d <- K.JJ.fun(t, 0, Mu.mj=-log(0), beta2.mj=0, dPenalty.j, fdy.stock)
	sK.JS.u <- K.JS.fun(t, orig=2, 1, Mu.mj=-log(0), beta2.mj=0, Mu.mdam, beta2.mdam, dPenalty.j, fdy.stock)
	sK.JS.d <- K.JS.fun(t, orig=2, 0, Mu.mj=-log(0), beta2.mj=0, Mu.mdam, beta2.mdam, dPenalty.j, fdy.stock)
	sK.SS <- K.SS.fun(t, sex, orig=2, Mu.mOs, beta2.mOs, fdy.stock)
	sK.SAs.u <- K.SAs.fun(t, 1, sex, orig=2, Mu.mOs, beta2.mOs, fdy.stock)
	sK.SAs.d <- K.SAs.fun(t, 0, sex, orig=2, Mu.mOs, beta2.mOs, fdy.stock)
	sK.AsJ.u <- K.AsJ.fun(t, m0, 1, orig=2, dPenalty.0, fdy.stock)
	sK.AsJ.d <- K.AsJ.fun(t, m0, 0, orig=2, dPenalty.0, fdy.stock)
	sK.AsAn.u <- K.AsAn.fun(t, 1, sex, orig=2, Lx, Ux, fdy.stock)
	sK.AsAn.d <- K.AsAn.fun(t, 0, sex, orig=2, Lx, Ux, fdy.stock)
	sK.AnAs.u <- K.AnAs.fun(t, 1, sex, orig=2, Lx, Ux, fdy.stock)
	sK.AnAs.d <- K.AnAs.fun(t, 0, sex, orig=2, Lx, Ux, fdy.stock)
	
		
	## Build mega-matrix kernel (wild & stocked)
	IPM.wild <- trout.megamatrix(wK.JJ.u, wK.JJ.d, wK.JS.u, wK.JS.d, wK.SS, wK.SAs.u, wK.SAs.d, wK.AsAn.u, wK.AsAn.d, wK.AsJ.u, wK.AsJ.d, wK.AnAs.u, wK.AnAs.d)$IPM

	IPM.stock <- trout.megamatrix(sK.JJ.u, sK.JJ.d, sK.JS.u, sK.JS.d, sK.SS, sK.SAs.u, sK.SAs.d, sK.AsAn.u, sK.AsAn.d, sK.AsJ.u, sK.AsJ.d, sK.AnAs.u, sK.AnAs.d)$IPM
	
	## Project population for N years
	test <- IPM.project2(IPM.wild, IPM.stock, sSSD, N, test.stock)$tot.pop
	
	## Asymptotic population size	
	eq.size <- sum(test[,N])
	
	## Proportion in each stage
	stage.prop <- c(sum(test[1:300,N]), sum(test[301:600,N]), sum(test[601:900,N]), sum(test[901:1200,N]), sum(test[1201:1500,N]), sum(test[1501:1800,N]))/eq.size
	
	## Mean size in each stage
	# = sum(frequency[size]*size/sum(frequency))
	stage.size.mean <- c(sum(test[1:300,N]*x/sum(test[1:300,N])), sum(test[301:600,N]*x/sum(test[301:600,N])), sum(test[601:900,N]*x/sum(test[601:900,N])), sum(test[901:1200,N]*x/sum(test[901:1200,N])), sum(test[1201:1500,N]*x/sum(test[1201:1500,N])), sum(test[1501:1800,N]*x/sum(test[1501:1800,N])))

	## Size variation in each stage
	# = sum(frequency[size]*size^2/sum(frequency)) - size.mean
	stage.size.var <- c(sum(test[1:300,N]*x^2/sum(test[1:300,N])), sum(test[301:600,N]*x^2/sum(test[301:600,N])), sum(test[601:900,N]*x^2/sum(test[601:900,N])), sum(test[901:1200,N]*x^2/sum(test[901:1200,N])), sum(test[1201:1500,N]*x^2/sum(test[1201:1500,N])), sum(test[1501:1800,N]*x^2/sum(test[1501:1800,N]))) - stage.size.mean^2
stage.size.SD <- sqrt(stage.size.var)
	
	## Collate results
	results <- data.frame(stage = c('Juvenile (u)', 'Juvenile (d)', 'Subadult', 'Spawner (u)', 'Spawner (d)', 'Post-spawner'), eq.size = rep(eq.size,6), stage.prop = stage.prop, mean.size = stage.size.mean, SD.size = stage.size.SD)
	
	# Add perturbed vital rate
	results$perturbation <- ifelse(sum(fdy.wild, fdy.stock) > length(fdy.wild)*2, names(fdy.wild)[which(fdy.wild!=1)], 'orig')
	
	return(results)	
}

## Set parameters for projection
sSSD <- rep(1, 1800)
N <- 200

## Run population analysis (for 3 different penalties)
orig.data1 <- pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,1], fdy.stock = fdy.stock[,1])

pert.data1 <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,z], fdy.stock = fdy.stock[,z]), simplify = FALSE))

orig.data2 <- pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,1], fdy.stock = fdy.stock[,1])

pert.data2 <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,z], fdy.stock = fdy.stock[,z]), simplify = FALSE))

orig.data3 <- pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,1], fdy.stock = fdy.stock[,1])

pert.data3 <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,z], fdy.stock = fdy.stock[,z]), simplify = FALSE))


#--------------------
# ELASTICITY ANALYSIS
#--------------------

elas.fun <- function(orig.data, pert.data, VR, dy){
	
	orig.data$eq.size.s <- 1
	
	# Take subset pertaining to perturbation of interest
	data <- subset(pert.data, perturbation == VR)
	data$eq.size.s <- data$eq.size/orig.data$eq.size[1]
	
	# Calculate sensitivity
	elas <- (data[,c(3:5,7)] - orig.data[,c(3:5,7)])/(orig.data[,c(3:5,7)]*dy)
	
	# Collate results
	elas$stage <- data$stage
	elas$perturbation <- data$perturbation
	return(elas)
}

test1 <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data1, pert.data1, VR = VR.names[z], DY), simplify = FALSE))
test1$Penalty <- 'none'

test2 <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data2, pert.data2, VR = VR.names[z], DY), simplify = FALSE))
test2$Penalty <- '+50% early mortality'

test3 <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data3, pert.data3, VR = VR.names[z], DY), simplify = FALSE))
test3$Penalty <- '+100% early mortality'

test <- rbind(test1, test2, test3)
test$Penalty <- factor(test$Penalty, levels = c('none', '+50% early mortality', '+100% early mortality'))
test$perturbation <- factor(test$perturbation, levels = c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH'))

#-----------------
# PLOTTING RESULTS
#-----------------

## Writing parameter names
par.label = expression(m['0,u'], m['0,d'], m['j,u'], m['j,d'], m[dam], m[s]^O, m['a,u']^O, m['a,d']^O, m^H)
 
## Function to re-format labels
scaleFUN <- function(x) sprintf("%.1f", x)

ggplot(subset(test, stage == 'Subadult'), aes(x = perturbation, y = eq.size.s, group = Penalty)) + geom_bar(aes(fill = Penalty), stat = 'identity', position = 'dodge') + ylab('Elasticity of population size') + ggtitle('a) Population with stocking') +
scale_fill_manual(name = 'Below-dam penalty', values = c('#73C76D', '#008087', 'black')) +
scale_x_discrete(labels = par.label) + scale_y_continuous(labels = scaleFUN) +
theme_bw() + theme(panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.35,0.2), plot.title = element_text(face = 'bold', hjust = -0.2), axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = 'bold'), panel.grid.minor = element_blank())


#----------------------------------------
# ANALYSES FOR WILD VS. STOCKED MORTALITY
#----------------------------------------

fdy.none <- matrix(1, nrow = length(VR.names), ncol = length(VR.names))
fdy.none <- cbind(rep(1, length(VR.names)), fdy.none)
rownames(fdy.none) <- VR.names
colnames(fdy.none) <- c('orig', VR.names)

## Only perturbing wild
pert.data1A <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,z], fdy.stock = fdy.none[,z]), simplify = FALSE))

pert.data2A <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,z], fdy.stock = fdy.none[,z]), simplify = FALSE))

pert.data3A <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.wild[,z], fdy.stock = fdy.none[,z]), simplify = FALSE))

test1A <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data1, pert.data1A, VR = VR.names[z], DY), simplify = FALSE))
test1A$Penalty <- 'no penalty'

test2A <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data2, pert.data2A, VR = VR.names[z], DY), simplify = FALSE))
test2A$Penalty <- '+50% early mortality penalty'

test3A <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data3, pert.data3A, VR = VR.names[z], DY), simplify = FALSE))
test3A$Penalty <- '+100% early mortality penalty'

testA <- rbind(test1A, test2A, test3A)
testA$Penalty <- factor(testA$Penalty, levels = c('no penalty', '+50% early mortality penalty', '+100% early mortality penalty'))
testA$perturbation <- factor(test$perturbation, levels = c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH'))
testA$model <- 'Wild'


## Only perturbing stocked
pert.data1B <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.none[,z], fdy.stock = fdy.stock[,z]), simplify = FALSE))
pert.data1B$perturbation <- rep(c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH'), each = 6)

pert.data2B <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.none[,z], fdy.stock = fdy.stock[,z]), simplify = FALSE))
pert.data2B$perturbation <- rep(c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH'), each = 6)

pert.data3B <- do.call("rbind", sapply(2:ncol(fdy.wild), FUN = function(z) pop.analysis(test.year, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, sSSD, N, test.stock, fdy.wild = fdy.none[,z], fdy.stock = fdy.stock[,z]), simplify = FALSE))
pert.data3B$perturbation <- rep(c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH'), each = 6)

test1B <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data1, pert.data1B, VR = VR.names[z], DY), simplify = FALSE))
test1B$Penalty <- 'no penalty'

test2B <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data2, pert.data2B, VR = VR.names[z], DY), simplify = FALSE))
test2B$Penalty <- '+50% early mortality penalty'

test3B <- do.call("rbind", sapply(1:length(VR.names), FUN = function(z) elas.fun(orig.data3, pert.data3B, VR = VR.names[z], DY), simplify = FALSE))
test3B$Penalty <- '+100% early mortality penalty'

testB <- rbind(test1B, test2B, test3B)
testB$Penalty <- factor(testB$Penalty, levels = c('no penalty', '+50% early mortality penalty', '+100% early mortality penalty'))
testB$perturbation <- factor(test$perturbation, levels = c('m0.u', 'm0.d', 'mj.u', 'mj.d', 'mdam', 'mOs', 'mO1', 'mO2', 'mH'))
testB$model <- 'Stocked'

## Plotting combined
test.comb <- rbind(testA, testB)

ggplot(subset(test.comb, stage == 'Subadult'), aes(x = perturbation, y = eq.size.s, group = model)) + geom_bar(aes(fill = model), stat = 'identity', position = 'dodge') + ylab('Elasticity of equilibrium population size') + facet_wrap(~Penalty, ncol = 1, scales = 'free_y') + 
scale_fill_manual(name = 'Origin', values = c('#ca0020', 'black')) +
scale_x_discrete(labels = par.label) + scale_y_continuous(labels = scaleFUN) +
theme_bw() + theme(panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor = element_blank())


