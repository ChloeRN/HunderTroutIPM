#********************************#
#*** DAM IPM FOR HUNDER TROUT ***#
#********************************#

## General Notes:
# This IPM is based on the life cycle that does distinguishes between spawners & juveniles upriver and downriver of the dam ("Trout_LifeCycle_IPM.pdf")
# It assumes that for mature adults, all mortality for spawning + non-spawning year happens before the spawning-non-spawning transition (splitting mortality hazard rate 100:0 over two years)
# The census is placed just after mature trout do / do not ascend the ladder
# The functions specified here already contain functionality for mortality perturbation analyses.

# This code is an alteration of "HunderTroutIPM_IPM_Building.R" that also allows for perturbations of
# below-dam spawner background mortality (mO2)


## Load packages
library(fields)
library(Matrix)
library(ggplot2)

## General settings for ignoring year and individual variation
chosenYear <- 1991 
test.year <- no.years <- chosenYear - 1951
no.inds <- 1

# -------------------------------------------------------------------------------

##################
#### 1) SETUP ####
##################

#--------------------------
#INITIALIZE THE MODEL 
#--------------------------

## Lower and upper limits for dynamic trait x (body length, mmm)
Lx <- 0
Ux <- 1300
 
## Number of "mesh points" in x
n <- 300
 
## Vector of dynamic trait values
x <- seq(Lx, Ux, length=n)

## Bin size
dx <- x[2]-x[1]

## Bin width 
h <- (Ux-Lx)/n

#-----------------------------------------------
#SET RANDOM EFFECTS AND ENVIRONMENTAL COVARIATES
#-----------------------------------------------

## General settings for ignoring year and individual variation
chosenYear <- 1991 
test.year <- no.years <- chosenYear - 1951
no.inds <- 1

## Setting individual random effect levels to 0
epsilon.grR.i <- rnorm(no.inds, 0, 0)
epsilon.grL.i <- rnorm(no.inds, 0, 0)
epsilon.muI.i <- rnorm(no.inds, 0, 0)

## Setting individual random effect levels to 0
epsilon.grR.i <- rnorm(no.inds, 0, 0)
epsilon.grL.i <- rnorm(no.inds, 0, 0)
epsilon.muI.i <- rnorm(no.inds, 0, 0)

## Setting year random effect levels to 0
epsilon.grR <- rnorm(no.years, 0, 0)
epsilon.grL <- rnorm(no.years, 0, 0)
epsilon.mH <- rnorm(no.years, 0, 0)
epsilon.mO <- rnorm(no.years, 0, 0)
epsilon.pS <- rnorm(no.years, 0, 0)
epsilon.pM <- rnorm(no.years, 0, 0)
epsilon.pL <- rnorm(no.years, 0, 0)

## Setting discharge covariate to 0
discF.std <- rep(0, no.years)
discS.std <- rep(0, no.years)


############################################
#### 2) VITAL RATE PREDICTION FUNCTIONS ####
############################################

#-------------
# RIVER GROWTH 
#-------------

## Parameter values 
mu0 <- 1.548 # Average length at hatching
h0 <- 69.206 # Baseline river growth increment (time = 0, year = 1951)
betaYR <- -0.224 # Time trend in river growth (times = 1:51, years = 1952:2002)
varR.i <- 59.521 # Variance of among-individual differences in river growth
varR.t <- 10.375 # Variance of among-year differences in river growth
varR.R <- 188.074 # Residual variance in river growth

## Prediction function for average river growth
mean.grR.pred <- function(size, t, i){
	
	# Calculate increment (linear)
	inc <- h0 + betaYR*t + epsilon.grR[t] + epsilon.grR.i[i]
	
	# Add increment to size
	size.next <- size + inc
	return(size.next)
}

## Prediction function for river growth distribution
dist.grR.pred <- function(size, size.next, t, i){
  
  # Calculate growth and set variance parameter
  mu <- mean.grR.pred(size, t, i)
  var <- varR.R
  
  # Make a switch to account for potential negative growth 
  z <- ifelse(size.next-size > 0, size.next, 0)
  
  # Calculate densities 
  y <- dnorm(z, mean = mu, sd = sqrt(var))
  
  # Scale and return densities
  if(sum(y*dx)==0){
    return (c(rep(0,n-1),1/dx))
  }
  else return(y/sum(y*dx))
  y/sum(y*dx)
}

## Sampling random effect levels
#epsilon.grR <- rnorm(no.years, 0, sqrt(varR.t))
#epsilon.grR.i <- rnorm(no.inds, 0, sqrt(varR.i))


#------------
# LAKE GROWTH 
#------------

## Parameter values (1st value = wild, 2nd value = stocked)
k0 <- c(0.155, 0.139) # Baseline lake growth rate
mu_inf0 <- c(1295.570, 1435.445) # Average asymptotic length
betaS <- -1.770 # Effect of spawning status on log(k)
varL.i.k <- 0.024 # Variance of among-individual differences in log(k)
varL.i.mu <- 0 # Variance of among-individual differences in asymptotic size
varL.t <- 0.014 # Variance of among-year differences in log(k)
varL.R <-  396.822 # Residual variance in lake growth increment

#NOTE: It was not possible to estimate varL.i.mu reliably with the updated data

## New growth reduction factor for mature fish
redFactor <- 0.63

## Prediction function for river growth
mean.grL.pred <- function(size, t, i, spawn.year, orig, mature){ 
	
	# Calculate lake growth rate and asymptotic size
    k <- exp(log(k0[orig]) + betaS*spawn.year + epsilon.grL[t] + epsilon.grL.i[i]) 
    mu_inf <- mu_inf0[orig] + epsilon.muI.i[i]
        
	# Calculate increment (von Bertalanffy)
	if(mature == 0){
		inc <- (mu_inf - size)*(1-exp(-k))
	}else{
		inc <- (mu_inf - size)*(1-exp(-k))*redFactor
	}
	
	# Include residual variation in increment and add to size
	#size.next <- size + rnorm(1, inc, sqrt(varL.R))
	size.next <- size + inc
	return(size.next)
}

## Prediction function for lake growth distribution
dist.grL.pred <- function(size, size.next, t, i, spawn.year, orig, mature){ 
  
  # Calculate growth and set variance parameter
  mu <- mean.grL.pred(size, t, i, spawn.year, orig, mature)
  var <- varL.R
  
  # Make a switch to account for potential negative growth 
  #z <- ifelse(size.next-size > 0, size.next, 0)
  z <- size.next
  
  # Calculate densities 
  y <- dnorm(z, mean = mu, sd = sqrt(var))
  
  # Scale and return densities
  if(sum(y*dx)==0){
    return (c(rep(0,n-1),1/dx))
  }
  else return(y/sum(y*dx))
  y/sum(y*dx)
}

# -------------------------------------------------------------------------------

#---------------
# ADULT SURVIVAL 
#---------------

## Parameter values (1st value = wild, 2nd value = stocked)
Mu.mH <- c(1.272, 1.258) # Median harvest mortality hazard rate
Mu.mO1 <- c(0.104, 0.103) # Median background mortality hazard rate (above-dam)
Mu.mO2 <- c(0.158, 0.157) # Median background mortality hazard rate (below-dam)
beta2.mH <- -0.345 # Linear size effect on log(mH)
beta4.mH <- -0.153 # Quadratic size effect on log(mH)
beta1.mO1 <- 0.639 # Linear discharge effect on log(mO1)
beta1.mO2 <- 0.063 # Linear discharge effect on log(mO2)
beta2.mO1 <- -0.315 # Linear size effect on log(mO1)
beta2.mO2 <- 0.781 # Linear size effect on log(mO2)
sigma.mH <- 0.135 # SD of random year variation in log(mH)
sigma.mO <- 1.322 # SD of random year variation in log(mO1) and log(mO2)

## Prediction function for harvest mortality (2 years)
mH.2yr.pred <- function(size, t, orig, l.limit, u.limit, mH.factor){
	
		
	# Scale length
	size.std <- (size - 669.1941)/108.6021
	
	# Make prediction 
	y.mH <- log(Mu.mH[orig]) + beta2.mH*size.std + beta4.mH*(size.std^2) + epsilon.mH[t]
	mH <- exp(y.mH)*mH.factor
	
	# Set mH to 0 below/above limits
	mH[which(size < l.limit)] <- 0
	mH[which(size > u.limit)] <- 0
	
	return(mH)
}

## Prediction function for background mortality (2 years)
mO.2yr.pred <- function(size, dam, t, orig, mO2.factor){

	# Scale length
	size.std <- (size - 669.1941)/108.6021
	
	# Make prediction (based on spawning location)
	if(dam == 1){
		y.mO <- log(Mu.mO1[orig]) + beta1.mO1*discF.std[t] + beta2.mO1*size.std + epsilon.mO[t]
		y.mO[which(size.std < 0)] <- log(Mu.mO1[orig]) + beta1.mO1*discF.std[t] + beta2.mO1*size.std[which(size.std < 0)] + epsilon.mO[t]
		
		mO <- exp(y.mO)
	}
	if(dam == 0){
		y.mO <- log(Mu.mO2[orig]) + beta1.mO2*discF.std[t] + beta2.mO2*size.std + epsilon.mO[t]
		mO <- exp(y.mO)*mO2.factor
	}
	return(mO)	
}


## Prediction function for annual adult survival

Sa.pred <- function(size, dam, t, orig, min.size, max.size, spawn.year, l.limit, u.limit, mH.factor, mO2.factor){ 
		
	# Calculate the vector of survival probabilities (whole size range)
	mTOT <- mH.2yr.pred(size, t, orig, l.limit, u.limit, mH.factor) + mO.2yr.pred(size, dam, t, orig, mO2.factor)	
	Sa <- exp(-mTOT)

	# Calculate minimum & maximum size survival
	mTOT.min <- (mH.2yr.pred(min.size, t, orig, l.limit, u.limit, mH.factor) + mO.2yr.pred(min.size, dam, t, orig, mO2.factor))
	mTOT.max <- (mH.2yr.pred(max.size, t, orig, l.limit, u.limit, mH.factor) + mO.2yr.pred(max.size, dam, t, orig, mO2.factor))
	Sa.min <- exp(-mTOT.min)
	Sa.max <- exp(-mTOT.max)
	
	# Replace survival probabilites for sizes beyond thresholds
	Sa[which(size < min.size)] <- Sa.min
	Sa[which(size > max.size)] <- Sa.max
	
	if(spawn.year==1){Sa.final <- Sa}
	if(spawn.year!=1){Sa.final <- rep(1, length(Sa))}
	return(Sa.final)
	# EDITED: Now returns either 2-year survival or 1 (depending on spawning state)
}

## Sampling random effect levels
#epsilon.mH <- rnorm(no.years, 0, sigma.mH)
#epsilon.mO <- rnorm(no.years, 0, sigma.mO)


#------------------
# SUBADULT SURVIVAL 
#------------------

## Prediction function for subadult background mortality  
mOs.pred <- function(size, Mu.mOs, beta2.mOs){
	
	size.std <- (size-437)/150 # 437 chosen as median (halfway between length at smolting/maturation based on growth data), 150 chosen as scaling factor
	
	log.mOs <- log(Mu.mOs) + beta2.mOs*size.std
	mOs <- exp(log.mOs) 

	
	return(mOs)
}

 
## Prediction function for annual subadult survival
Ss.pred <- function(size, t, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor){
	
	mOs <- mOs.pred(size, Mu.mOs, beta2.mOs)
	mTOT <- (mH.2yr.pred(size, t, orig, l.limit, u.limit, mH.factor)/2) + mOs	
	Ss <- exp(-mTOT)
	return(Ss)
}


#------------------
# JUVENILE SURVIVAL
#------------------


Sj.pred <- function(size, Mu.mj, beta2.mj, dPenalty.j){
	
	size.std <- (size-157)/50 # 157 chosen as median (halway between length at age 1 and length at smolting)
	
	log.mj <- log(Mu.mj) + beta2.mj*size.std
	mj <- exp(log.mj)*dPenalty.j
	return(exp(-mj))
}


#---------
# SMOLTING
#---------

## Parameter values 
mu.pS <- -2.482 # Intercept of logit(pS)
beta2.pS <- 2.935 # Size effect on logit(pS)
betaY2.pS <- 0.294 # Size:time effect on logit(pS)
var.pS <- 0.205 # Among-year variance in logit(pS)

## Prediction function for smolting probability
pS.pred <- function(size, t){
	
	# Calculate relevant year from t (t=0 corresponds to 1951)
	year <- t + 1951
	
	# Scale length and year
	size.std <- (size - 164.5047)/78.71436
	year.std <- (year - 1977.933)/12.48213
	
	# Make prediction 
	y.pS <- mu.pS + beta2.pS*size.std + betaY2.pS*size.std*year.std + epsilon.pS[t]
	pS <- plogis(y.pS)
	
	return(pS)
}

## Sampling random effect levels
#epsilon.pS <- rnorm(no.years, 0, sqrt(var.pS))


#------------
# MATURATION
#------------

## Parameter values (1st dimension = origin (1 = wild, 2 = stocked), 2nd dimension = sex (1 = female, 2 = male))
betaMale <- 0.334
betaWild <- -0.068
int0 <- -1.165

betaMale.size <- -0.944
betaWild.size <- -0.380
slp0.size <- 2.556

mu.pM <- cbind(c(int0+betaWild, int0), c(int0+betaWild+betaMale, int0+betaMale)) # Intercept of logit(pS)
beta2.pM <- cbind(c(slp0.size+betaWild.size, slp0.size), c(slp0.size+betaWild.size+betaMale.size, slp0.size+betaMale.size)) # Size effect on logit(pS)
betaY.pM <- 0.413 # Time trend in logit(pS)
var.pM <- 0.547 # Among-year variance in logit(pM)

## Prediction function for smolting probability
pM.pred <- function(size, t, sex, orig){
	
	# Calculate relevant year from t (t=0 corresponds to 1951)
	year <- t + 1951
	
	# Scale length and year
	size.std <- (size - 511.4021)/128.2904
	year.std <- (year - 1984.539)/11.64529
	
	# Make prediction 
	y.pM <- mu.pM[orig,sex] + beta2.pM[orig,sex]*size.std + betaY.pM*year.std + epsilon.pM[t]
	pM <- plogis(y.pM)
	
	return(pM)
}

## Sampling random effect levels
#epsilon.pM <- rnorm(no.years, 0, sqrt(var.pM))


#-------------
# LADDER USAGE
#-------------

# Parameter values (1st value = wild, 2nd value = stocked)
Mu.pL <- c(0.532, 0.476) # Average ladder usage probability
beta1.pL <- 0.271 # Linear discharge effect on logit(p)
beta2.pL <- 0.155 # Linear size effect on logit(p)
beta3.pL <- -0.281 # Interactive discharge:size effect on logit(p)
beta4.pL <- -0.622 # Quadratic size effect on logit(p)
sigma.pL <- 0.383 # SD of random year variation in logit(p)

# Prediction function for annual ladder usage probability
pL.pred <- function(size, t, orig){
	
	# Scale length
	size.std <- (size - 669.1941)/108.6021
	
	# Make prediction 
	y.pL <- qlogis(Mu.pL[orig]) + beta1.pL*discS.std[t] + beta2.pL*size.std + beta3.pL*discS.std[t]*size.std + beta4.pL*(size.std^2) + epsilon.pL[t]
	pL <- plogis(y.pL)
	
	return(pL)
}

## Sampling random effect levels
#epsilon.pL <- rnorm(no.years, 0, sigma.pL)


#-----------
# FECUNDITY
#-----------

## Parameter values
mean.log.fec <- -7.0218 # Intercept on the log scale
size.eff.fec <- 2.3422 # Effect of log(size) on log(fecundity)

## Prediction function for fecundity
fec.pred <- function(size){
	
	y.fec <- mean.log.fec + size.eff.fec*log(size)
	fec <- exp(y.fec)
	return(fec)
}


#---------------
# EARLY SURVIVAL 
#---------------

S0.pred <- function(m0, dPenalty.0){ # Added dy (should be 1 number)
	
	m0.new <- m0*dPenalty.0
	return(exp(-m0.new))
}


#----------------
# OFFSPRING SIZE 
#----------------

# Prediction function for average offspring size (Option 2 - Growth model) 
mean.sizeOff.pred <- function(t, i){
	
	# Calculate increment (linear)
	inc <- h0 + betaYR*t + epsilon.grR[t] + epsilon.grR.i[i]
	
	# Add increment to size at hatching
	size.next <- mu0 + rnorm(1, inc, 0)
	return(size.next)
}

## Prediction function for river growth distribution
dist.sizeOff.pred <- function(size.next, t, i){
  
  # Calculate growth and set variance parameter
  mu <- mean.sizeOff.pred(t, i)
  var <- varR.R
    
  # Calculate densities 
  y <- dnorm(size.next, mean = mu, sd = sqrt(var))
  
  # Scale and return densities
  if(sum(y*dx)==0){
    return (c(rep(0,n-1),1/dx))
  }
  else return(y/sum(y*dx))
  y/sum(y*dx)
}


#----------------------
# DAM SURVIVAL (SMOLTS)
#----------------------

Sdam.pred <- function(size, Mu.mdam, beta2.mdam){ 
	
	size.std <- (size-250)/50
	
	log.mdam <- log(Mu.mdam) + beta2.mdam*size.std
	mdam <- exp(log.mdam)
	return(exp(-mdam))
}


###########################################################
### 3) CREATE KERNELS FOR ALL STAGE-MATRIX TRANSITIONS ####
###########################################################

#----------------------------------------------
# a) Juvenile to Juvenile (upriver & downriver)
#----------------------------------------------

K.JJ.fun  <-function(t, dam, Mu.mj, beta2.mj, dPenalty.j){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
    if(dam == 1){
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, 1)
  	}else{
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, dPenalty.j)
  } 
  Tvec <- 1-(pS.pred(x, t))

  # Survive, do not smolt, and grow (river growth)
  for(i in 1:n){
    Smat[,i] <- Svec[i]*Tvec[i]*dist.grR.pred(x[i], x, t, 1)*dx
  }
  return(Smat)
}


#----------------------------------------------
# b) Juvenile to Subadult (upriver & downriver)
#----------------------------------------------

K.JS.fun  <-function(t, orig, dam, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  if(dam == 1){
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, 1)*Sdam.pred(x, Mu.mdam, beta2.mdam)
  }else{
  	Svec <- Sj.pred(x, Mu.mj, beta2.mj, dPenalty.j)
  }  
  Tvec <- pS.pred(x, t)
  
  # Survive, smolt, and grow (lake growth, without reproduction cost)
  for (i in 1:n){
    Smat[,i] <- Svec[i]*Tvec[i]*dist.grL.pred(x[i], x, t, 1, 0, orig, 0)*dx

  }
  return(Smat)
}


#------------------------
# c) Subadult to Subadult
#------------------------

K.SS.fun  <-function(t, sex, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  Svec <- Ss.pred(x, t, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor)
  Tvec <- 1-(pM.pred(x, t, sex, orig))
  
  # Survive, do not mature, and grow (lake growth, without reproduction cost)
  for (i in 1:n){
    Smat[,i] <- Svec[i]*Tvec[i]*dist.grL.pred(x[i], x, t, 1, 0, orig, 0)*dx
  }
  return(Smat)
}


#---------------------------------------------
# d) Subadult to Spawner (upriver & downriver)
#---------------------------------------------

K.SAs.fun  <-function(t, dam, sex, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  Svec <- Ss.pred(x, t, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor)
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

K.AsAn.fun  <-function(t, dam, sex, orig, min.size, max.size, l.limit, u.limit, mH.factor, mO2.factor){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  Tvec <- rep(1, n) # 1 because the transition is deterministic (and independent of the dam)
  Svec <- Sa.pred(x, dam, t, orig, min.size, max.size, spawn.year = 1, l.limit, u.limit, mH.factor, mO2.factor) 

  # Deterministic transition, survive, and grow (without reproduction cost)
  for (i in 1:n){
    Smat[,i] <- Tvec[i]*Svec[i]*dist.grL.pred(x[i], x, t, 1, 0, orig, 1)*dx
  }
  return(Smat)
}


#-----------------------------------------------------
# f) Post-spawner to Spawner (upriver & downriver)
#-----------------------------------------------------

K.AnAs.fun  <-function(t, dam, sex, orig, min.size, max.size, l.limit, u.limit, mH.factor, mO2.factor){
  
  # Define matrix
  Smat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  if(dam == 1){
  	Tvec <- pL.pred(x, t, orig)
  }else{
  	Tvec <- 1-pL.pred(x, t, orig)
  }
  
  Svec <- Sa.pred(x, dam, t, orig, min.size, max.size, spawn.year = 0, l.limit, u.limit, mH.factor, mO2.factor) # Survival is 1 here
  
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



#----------------------------------------------
# g) Spawners to Juvenile (upriver & downriver)
#----------------------------------------------

K.AsJ.fun  <-function(t, m0, dam, orig, dPenalty.0){
  
  # Define matrix
  Bmat <- matrix(0, n, n)
  
  # Fill matrices using vital rate functions
  if(dam == 1){
  	Bvec <- fec.pred(x)*S0.pred(m0, 1)*0.5
  }else{
  	Bvec <- fec.pred(x)*S0.pred(m0, dPenalty.0)*0.5
  }
  
  # Reproduction, egg survival, growth to 1-year size
  for (i in 1:n){
    Bmat[,i] <- Bvec[i]*dist.sizeOff.pred(x, t, 1) * dx
  }
  return(Bmat)
}


##################################
#### 4) BUILD MEGA-MATRIX IPM ####
##################################

## Defining transition matrix 
# All entries in the transition matrix are 1, because the size-specific transition probabilities are already included in the Kernels
stg.mat <- matrix(c(1, 0, 0, 1, 0, 0,
                    0, 1, 0, 0, 1, 0,
                    1, 1, 1, 0, 0, 0,
                    0, 0, 1, 0, 0, 1,
                    0, 0, 1, 0, 0, 1,
                    0, 0, 0, 1, 1, 0), nrow=6, ncol=6, byrow=TRUE)
stg.mat # Changes to columns 3-6

# Create an additional component matrix filled with only 0's
K0 <- matrix(0, nrow=n, ncol=n)
K.sp <- as(K0, 'sparseMatrix')

## Fill the kernels into the transition matrix 
trout.megamatrix = function(K.JJ.u, K.JJ.d, K.JS.u, K.JS.d, K.SS, K.SAs.u, K.SAs.d, K.AsAn.u, K.AsAn.d, K.AsJ.u, K.AsJ.d, K.AnAs.u, K.AnAs.d){ # Kernels updated
   
   # Complete IPM kernel - updated
   J.u <- cbind(K.JJ.u, K0, K0, K.AsJ.u, K0, K0)
   J.d <- cbind(K0, K.JJ.d, K0, K0, K.AsJ.d, K0) 
   S <- cbind(K.JS.u, K.JS.d, K.SS, K0, K0, K0)
   As.u <- cbind(K0, K0, K.SAs.u, K0, K0, K.AnAs.u) 
   As.d <- cbind(K0, K0, K.SAs.d, K0, K0, K.AnAs.d)
   An <- cbind(K0, K0, K0, K.AsAn.u, K.AsAn.d, K0)
   
   mega.K <- rbind(J.u, J.d, S, As.u, As.d, An)
   
   # Pmatrix only
   J2.u <- cbind(K.JJ.u, K0, K0, K0, K0, K0)
   J2.d <- cbind(K0, K.JJ.d, K0, K0, K0, K0)
   mega.P <- rbind(J2.u, J2.d, S, As.u, As.d, An)
   
   # Complete IPM kernel with NA's for empty spots - updated
   y <- matrix(NA, nrow=n, ncol=n)

   J3.u <- cbind(K.JJ.u, y, y, K.AsJ.u, y, y)
   J3.d <- cbind(y, K.JJ.d, y, y, K.AsJ.d, y) 
   S3 <- cbind(K.JS.u, K.JS.d, K.SS, y, y, y)
   As3.u <- cbind(y, y, K.SAs.u, y, y, K.AnAs.u) 
   As3.d <- cbind(y, y, K.SAs.d, y, y, K.AnAs.d)
   An3 <- cbind(y, y, y, K.AsAn.u, K.AsAn.d, y)
      
   mega.K.na <- rbind(J3.u, J3.d, S3, As3.u, As3.d, An3)
   
   # Full Kernel without 0's in the possible transitions (for plotting) - updated
   x <- 1e-200
   
   Jx.u <- cbind(K.JJ.u+x, K0, K0, K.AsJ.u+x, K0, K0)
   Jx.d <- cbind(K0, K.JJ.d+x, K0, K0, K.AsJ.d+x, K0) 
   Sx <- cbind(K.JS.u+x, K.JS.d+x, K.SS+x, K0, K0, K0)
   Asx.u <- cbind(K0, K0, K.SAs.u+x, K0, K0, K.AnAs.u+x) 
   Asx.d <- cbind(K0, K0, K.SAs.d+x, K0, K0, K.AnAs.d+x)
   Anx <- cbind(K0, K0, K0, K.AsAn.u+x, K.AsAn.d+x, K0)
   
   mega.Kx <- rbind(Jx.u, Jx.d, Sx, Asx.u, Asx.d, Anx)

   return(list(IPM = mega.K, Pmatrix = mega.P, IPMna = mega.K.na, IPMx = mega.Kx))
 }


#-------------------------------------
# FUNCTION FOR BUILDING DIFFERENT IPMs
#-------------------------------------

# The following function builds the IPM for a specific size range and binwidth

build.IPM <- function(t, sex, orig, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0, dPenalty.j, l.limit, u.limit, mH.factor, mO2.factor){

	## Build separate kernels
	K.JJ.u <- K.JJ.fun(t, 1, Mu.mj, beta2.mj, dPenalty.j)
	K.JJ.d <- K.JJ.fun(t, 0, Mu.mj, beta2.mj, dPenalty.j)
	K.JS.u <- K.JS.fun(t, orig, 1, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j)
	K.JS.d <- K.JS.fun(t, orig, 0, Mu.mj, beta2.mj, Mu.mdam, beta2.mdam, dPenalty.j)
	K.SS <- K.SS.fun(t, sex, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor)
	K.SAs.u <- K.SAs.fun(t, 1, sex, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor) # new kernel (+ dam argument)
	K.SAs.d <- K.SAs.fun(t, 0, sex, orig, Mu.mOs, beta2.mOs, l.limit, u.limit, mH.factor) # new kernel (+ dam argument)
	K.AsJ.u <- K.AsJ.fun(t, m0, 1, orig, dPenalty.0)
	K.AsJ.d <- K.AsJ.fun(t, m0, 0, orig, dPenalty.0)
	K.AsAn.u <- K.AsAn.fun(t, 1, sex, orig, Lx, Ux, l.limit, u.limit, mH.factor, mO2.factor)
	K.AsAn.d <- K.AsAn.fun(t, 0, sex, orig, Lx, Ux, l.limit, u.limit, mH.factor, mO2.factor)
	K.AnAs.u <- K.AnAs.fun(t, 1, sex, orig, Lx, Ux, l.limit, u.limit, mH.factor, mO2.factor) # new kernel (+ dam argument)
	K.AnAs.d <- K.AnAs.fun(t, 0, sex, orig, Lx, Ux, l.limit, u.limit, mH.factor, mO2.factor) # new kernel (+ dam argument)
	
	## Build mega-matrix kernel
	IPM <- trout.megamatrix(K.JJ.u, K.JJ.d, K.JS.u, K.JS.d, K.SS, K.SAs.u, K.SAs.d, K.AsAn.u, K.AsAn.d, K.AsJ.u, K.AsJ.d, K.AnAs.u, K.AnAs.d) # Kernels updated
	
	return(IPM)	
}

#----------------------------
# FUNCTION FOR EIGEN ANALYSIS
#----------------------------

wvlambda <- function(Kmat){ # Function returning lambda, stable structure w and reproductive values v scaled so that sum(v*w)=1
  ev <- eigen(Kmat)
  tev <- eigen(t(Kmat))
  lmax <- which.max(Re(ev$values))
  W <- ev$vectors
  V <- tev$vectors
  w <- as.matrix(abs(Re(W[, lmax]))/sum(abs(Re(W[, lmax]))))
  w <- w/(sum(w))
  v <- as.matrix(abs(Re(V[, lmax])))
  v <- v/sum(w*v )
  v <- ifelse(w*v <= 0, 0, v)
  return(list("lambda"=max(Re(ev$values)),"w"=w,"v"=v))
}


##############################
#### QUANTIFYING STOCKING ####
##############################

## Set parameters (from stocking data)
stock.no.above <- 11908.84
stock.no.below <- 15110.59

stock.size.mean <- 201.902
stock.size.sd <- 26.66747

## Make a function for a calculating size-distribution of stocked fish
stock.fun  <- function(t, Mu.mdam, beta2.mdam){
  
  # Define initial sizes of hatchery smolt released above and below the dam
  init.size.above <- rnorm(stock.no.above/2, stock.size.mean, stock.size.sd)
  init.size.below <- rnorm(stock.no.below/2, stock.size.mean, stock.size.sd)
  
  # Above-dam releases: determine survivors of dam passage
  survP <- rbinom(length(init.size.above), 1, Sdam.pred(init.size.above, Mu.mdam, beta2.mdam))
  
  # Add above- and below-dam releases together
  init.size <- c(init.size.below, init.size.above[which(survP==1)])
  
  # Let survivors have their first year of lake growth
  end.size <- mean.grL.pred(init.size, t, 1, 0, 2, 0) + rnorm(length(init.size), mean = 0, sd = sqrt(varL.R))
  
  # Sort and count individuals by bins
  stock.distS <- unname(table(cut(end.size, breaks = x)))
  
  # Write complete stocking vector (all life stages)
  stock.dist <- c(rep(0,length(x)*2), c(stock.distS,0), rep(0,length(x)*3))
  return(stock.dist)
}