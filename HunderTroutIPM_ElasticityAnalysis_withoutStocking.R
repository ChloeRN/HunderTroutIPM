
# This code builds the IPM for three different assumptions of below-dam penalty on early mortality
# (none, +50%, +100%) and runs an elasticity analysis of population growth rate without stocking.

# The relevant sections in the paper are:
# - Methods: 2.3.3
# - Results: 3.2
# - Figures: Figure S1.3 b)

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


#---------------------------------
# FUNCTION FOR POPULATION ANALYSIS
#---------------------------------

# The following function builds the IPM for a specific hazard rate perturbation given by the vector dy and calculates the following:
# 1) Asymptotic population growth rate (lambda)
# 2) Proportions of the population in each stage
# 3) Mean size of individuals in each stage
# 4) Size variability (SD) of individuals in each stage

pop.analysis <- function(t, sex, orig, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0, dPenalty.j, fdy){
	
	## Build separate kernels
	K.JJ.u <- K.JJ.fun(t, 1, Mu.mj, beta2.mj, dPenalty.j, fdy)
	K.JJ.d <- K.JJ.fun(t, 0, Mu.mj, beta2.mj, dPenalty.j, fdy)
	K.JS.u <- K.JS.fun(t, orig, 1, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j, fdy)
	K.JS.d <- K.JS.fun(t, orig, 0, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j, fdy)
	K.SS <- K.SS.fun(t, sex, orig, Mu.mOs, beta2.mOs, fdy)
	K.SAs.u <- K.SAs.fun(t, 1, sex, orig, Mu.mOs, beta2.mOs, fdy)
	K.SAs.d <- K.SAs.fun(t, 0, sex, orig, Mu.mOs, beta2.mOs, fdy)
	K.AsJ.u <- K.AsJ.fun(t, m0, 1, orig, dPenalty.0, fdy)
	K.AsJ.d <- K.AsJ.fun(t, m0, 0, orig, dPenalty.0, fdy)
	K.AsAn.u <- K.AsAn.fun(t, 1, sex, orig, Lx, Ux, fdy)
	K.AsAn.d <- K.AsAn.fun(t, 0, sex, orig, Lx, Ux, fdy)
	K.AnAs.u <- K.AnAs.fun(t, 1, sex, orig, Lx, Ux, fdy)
	K.AnAs.d <- K.AnAs.fun(t, 0, sex, orig, Lx, Ux, fdy)
	
	## Build mega-matrix kernel
	IPM <- trout.megamatrix(K.JJ.u, K.JJ.d, K.JS.u, K.JS.d, K.SS, K.SAs.u, K.SAs.d, K.AsAn.u, K.AsAn.d, K.AsJ.u, K.AsJ.d, K.AnAs.u, K.AnAs.d)
	
	## Eigen analysis
	output <- wvlambda(IPM$IPM)
	
	## Asymptotic population growth rate	
	lambda <- output$lambda
	
	## Proportion in each stage
	stage.prop <- c(sum(output$w[1:300]), sum(output$w[301:600]), sum(output$w[601:900]), sum(output$w[901:1200]), sum(output$w[1201:1500]), sum(output$w[1501:1800]))
	
	## Mean size in each stage
	# = sum(frequency[size]*size/sum(frequency))
	# (formula from IPM book pp. 27, and verified with dummy data)
	stage.size.mean <- c(sum(output$w[1:300]*x/sum(output$w[1:300])), sum(output$w[301:600]*x/sum(output$w[301:600])), sum(output$w[601:900]*x/sum(output$w[601:900])), sum(output$w[901:1200]*x/sum(output$w[901:1200])), sum(output$w[1201:1500]*x/sum(output$w[1201:1500])), sum(output$w[1501:1800]*x/sum(output$w[1501:1800])))
	
	## Size variation in each stage
	# = sum(frequency[size]*size^2/sum(frequency)) - size.mean
	# (formula from IPM book pp. 39, but failed to verify with dummy data)
	stage.size.var <- c(sum(output$w[1:300]*x^2/sum(output$w[1:300])), sum(output$w[301:600]*x^2/sum(output$w[301:600])), sum(output$w[601:900]*x^2/sum(output$w[601:900])), sum(output$w[901:1200]*x^2/sum(output$w[901:1200])), sum(output$w[1201:1500]*x^2/sum(output$w[1201:1500])), sum(output$w[1501:1800]*x^2/sum(output$w[1501:1800]))) - stage.size.mean^2
	stage.size.SD <- sqrt(stage.size.var)
	
	## Collate results
	results <- data.frame(stage = c('Juvenile (u)', 'Juvenile (d)', 'Subadult', 'Pre-spawner', 'Post-spawner (u)', 'Post-spawner (d)'), lambda = rep(lambda,6), stage.prop = stage.prop, mean.size = stage.size.mean, SD.size = stage.size.SD)
	
	# Add perturbed vital rate
	results$perturbation <- ifelse(sum(fdy) > length(fdy), names(fdy)[which(fdy!=1)], 'orig')
	
	return(results)	
}

## Run population analysis (for 3 different penalties)
orig.data1 <- pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, fdy = fdy[,1])

pert.data1 <- do.call("rbind", sapply(2:ncol(fdy), FUN = function(z) pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, fdy = fdy[,z]), simplify = FALSE))

orig.data2 <- pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, fdy = fdy[,1])

pert.data2 <- do.call("rbind", sapply(2:ncol(fdy), FUN = function(z) pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, fdy = fdy[,z]), simplify = FALSE))

orig.data3 <- pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, fdy = fdy[,1])

pert.data3 <- do.call("rbind", sapply(2:ncol(fdy), FUN = function(z) pop.analysis(test.year, 1, 1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, fdy = fdy[,z]), simplify = FALSE))

#--------------------
# ELASTICITY ANALYSIS
#--------------------

elas.fun <- function(orig.data, pert.data, VR, dy){
	
	# Take subset pertaining to perturbation of interest
	data <- subset(pert.data, perturbation == VR)
	
	# Calculate sensitivity
	elas <- (data[,2:5] - orig.data[,2:5])/(orig.data[,2:5]*dy)
	
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

ggplot(subset(test, stage == 'Subadult'), aes(x = perturbation, y = lambda, group = Penalty)) + geom_bar(aes(fill = Penalty), stat = 'identity', position = 'dodge') + ylab(expression('Elasticity of '~lambda)) + ggtitle('b) Population without stocking') +
scale_fill_manual(name = 'Below-dam penalty', values = c('#73C76D', '#008087', 'black')) +
scale_x_discrete(labels = par.label) + scale_y_continuous(breaks = c(0, -0.1, -0.2), labels = scaleFUN) +
theme_bw() + theme(panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.65,0.2), plot.title = element_text(face = 'bold', hjust = -0.2), axis.title = element_text(size = 13), axis.text = element_text(size = 10), legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = 'bold'), panel.grid.minor = element_blank())
