
# This code builds the IPM for three different assumptions of below-dam penalty on early mortality
# (none, +50%, +100%) and projects it forward in time with and without stocking.
# It further calculates equilibrium population size, population growth rates, and size distributions. 

# The relevant sections in the paper are:
# - Methods: 2.3.1 & 2.3.2
# - Results: 3.1
# - Figures: Figures 2 & S1.1

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


## Draw a stocking vector for testing
set.seed <- 12
test.stock <- stock.fun(test.year, Mu.mdam, beta2.mdam)


#---------------------------------------------------------
# IPM PROJECTION WITH STOCKING AND DIFFERENT DAM PENALTIES
#---------------------------------------------------------

## Build IPMs with different below-dam penalties - wild fish
IPM1.wild <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.wild <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.wild <- build.IPM(t=test.year, sex=1, orig=1, m0, Mu.mj, Mu.mdam, Mu.mOs, beta2.mj, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

## Build IPMs with different below-dam penalties - stocked fish
# --> Vital rates of stocked fish apply
# --> Juvenile survival = 0, because natural offspring of stocked fish join the wild population (not the stocked)
IPM1.stock <- build.IPM(t=test.year, sex=1, orig=2, m0, Mu.mj=-log(0), Mu.mdam, Mu.mOs, beta2.mj=0, beta2.mdam, beta2.mOs, dPenalty.0=1, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM2.stock <- build.IPM(t=test.year, sex=1, orig=2, m0, Mu.mj=-log(0), Mu.mdam, Mu.mOs, beta2.mj=0, beta2.mdam, beta2.mOs, dPenalty.0=1.5, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM

IPM3.stock <- build.IPM(t=test.year, sex=1, orig=2, m0, Mu.mj=-log(0), Mu.mdam, Mu.mOs, beta2.mj=0, beta2.mdam, beta2.mOs, dPenalty.0=2, dPenalty.j=1, l.limit=Lx, u.limit=Ux, mH.factor=1, mO1size.factor=1)$IPM


## Function for manual projection (with stocking & explicit vital rates)
IPM.project2 <- function(IPM.wild, IPM.stock, initSSD, steps, stock){
	
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
	  wild.pop[,t] <- IPM.wild %*% wild.pop[,t-1] + stock.off
	}
	
	# Remove juveniles from stocked population
	stock.pop[1:600,] <- 0
	
	return(list(wild.pop = wild.pop, stock.pop = stock.pop))
}

## Testing - from stable size-by-stage distribution
sSSD <- rep(1, 1800)
test1 <- IPM.project2(IPM1.wild, IPM1.stock, sSSD, 200, test.stock)
test2 <- IPM.project2(IPM2.wild, IPM2.stock, sSSD, 200, test.stock)
test3 <- IPM.project2(IPM3.wild, IPM3.stock, sSSD, 200, test.stock)

test1.all <- test1$wild.pop + test1$stock.pop
test2.all <- test2$wild.pop + test2$stock.pop
test3.all <- test3$wild.pop + test3$stock.pop

sum(test1.all[1:300,200]) 
sum(test1.all[301:600,200]) 
sum(test1.all[601:900,200])
sum(test1.all[901:1200,200]) 
sum(test1.all[1201:1500,200]) 
sum(test1.all[1501:1800,200])

## Plot projection
par(mar=c(5.1,4.1,3.1,2.1))
plot(colSums(test1.all[,150:270]), type = 'l', ylab = 'Total population size (age 1+)', xlab = 'Number of years', lwd = 2, cex.axis = 1.2, cex.lab = 1.3, col = '#73C76D')
lines(colSums(test2.all[,150:270]), col = '#008087', lwd = 2)
lines(colSums(test3.all[,150:270]), col = 'black', lwd = 2)
abline(v = 51, col = 'grey', lwd = 1.5, lty = 2)
legend(70, 140000, legend = c('none', '+50% early mortality', '+100% early mortality'), col = c('#73C76D', '#008087', 1), lty = c(1,1,1), lwd = c(2,2,2), bty = 'n', title = expression(bold('Below-dam penalty')), title.adj = 0.1, cex = 1.2)
dev.off()

## Population growth rates and equilibrium population sizes
wvlambda(IPM1.wild)$lambda 
wvlambda(IPM2.wild)$lambda
wvlambda(IPM3.wild)$lambda 

sum(test1.all[,200]) 
sum(test2.all[,200]) 
sum(test3.all[,200]) 

#----------------------------
# SIZE-BY-STAGE DISTRIBUTIONS
#----------------------------
library(ggplot2)

data.stock <- data.frame(Stage = rep(c('Juveniles', 'Subadults', 'Spawners', 'Post-spawners'), each = 300), Size = rep(x, 4), Number = c(test1.all[1:300,150]+test1.all[301:600,150], test1.all[601:900,150], test1.all[901:1200,150] + test1.all[1201:1500,150], test1.all[1501:1800,150]))
data.stock$Stage <- factor(data.stock$Stage, levels = c("Juveniles","Subadults","Spawners","Post-spawners"))

dat_text <- data.frame(
  label = c(paste('Total =',round(sum(test1.all[1:600,150]))), paste('Total =',round(sum(test1.all[601:900,150]))),paste('Total =',round(sum(test1.all[901:1200,150]))), paste('Total =',round(sum(test1.all[1201:1800,150])))),
  Stage   = c('Juveniles', 'Subadults', 'Spawners', 'Post-spawners'),
  x     = 1100,
  y     = c(25500, 850, 16, 5.5)
)

## Plot size-by-stage distributions for an equilibrium population with ant without stocking
ggplot(data.stock, aes(x = Size, y = Number)) + geom_density(color = '#00C0E2', fill = '#00C0E2', stat = 'identity', alpha = 0.5) + facet_wrap(~Stage, scales = 'free_y') + ylab('Number of individuals of each size') + xlab('Body length (mm)') + theme_bw() + theme(panel.grid.minor = element_blank()) + geom_text(dat_text, mapping = aes(x = x, y = y, label = label), col = '#00C0E2')
dev.off()

ggplot(data.stock, aes(x = Size, y = Number)) + geom_density(aes(color = Stage, fill = Stage), stat = 'identity', alpha = 0.5)

aSSD <- wvlambda(IPM1.wild)$w*sum(test1.all[,150])
data.nostock <- data.frame(Stage = rep(c('Juveniles', 'Subadults', 'Spawners', 'Post-spawners'), each = 300), Size = rep(x, 4), Number = c(aSSD[1:300]+aSSD[301:600], aSSD[601:900], aSSD[901:1200] + aSSD[1201:1500], aSSD[1501:1800]))
data.nostock$Stage <- factor(data.nostock$Stage, levels = c("Juveniles","Subadults","Spawners","Post-spawners"))

data <- rbind(data.stock, data.nostock)
data$Model <- rep(c(' With stocking     ', ' Without stocking     '), each = 1200)

ggplot(data, aes(x = Size, y = Number, group = Model)) + geom_density(aes(color = Model, fill = Model), stat = 'identity', alpha = 0.3) + facet_wrap(~Stage, scales = 'free_y') + ylab('') + xlab('Body length (mm)') + ylab('Scaled number of individuals of each size') + 
scale_fill_manual(values = c('#00C0E2', '#ca0020')) + 
scale_color_manual(values = c('#00C0E2', '#ca0020')) + 
theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = 'top', legend.title = element_blank())


#--------------------------
# PERCENTAGES IN EACH STAGE
#--------------------------

### NO PENALTY
## Stocked population
sum(test1.all[1:300,200])/sum(test1.all[,200]) 
sum(test1.all[301:600,200])/sum(test1.all[,200]) 
sum(test1.all[601:900,200])/sum(test1.all[,200]) 
sum(test1.all[901:1200,200])/sum(test1.all[,200]) 
sum(test1.all[1201:1500,200])/sum(test1.all[,200]) 
sum(test1.all[1501:1800,200])/sum(test1.all[,200]) 

# Juveniles
(sum(test1.all[1:300,200])+sum(test1.all[301:600,200]))/sum(test1.all[,200]) 

# Spawners
(sum(test1.all[901:1200,200])+sum(test1.all[1201:1500,200]))/sum(test1.all[,200]) 

## Wild population
SSD <- wvlambda(IPM1.wild)$w
sum(SSD[1:300]) 
sum(SSD[301:600]) 
sum(SSD[601:900]) 
sum(SSD[901:1200]) 
sum(SSD[1201:1500]) 
sum(SSD[1501:1800]) 

# Juveniles
sum(SSD[1:300])+sum(SSD[301:600]) 

# Spawners
sum(SSD[901:1200])+sum(SSD[1201:1500]) 



### 50% PENALTY
## Stocked population
sum(test2.all[1:300,200])/sum(test2.all[,200]) 
sum(test2.all[301:600,200])/sum(test2.all[,200]) 
sum(test2.all[601:900,200])/sum(test2.all[,200]) 
sum(test2.all[901:1200,200])/sum(test2.all[,200]) 
sum(test2.all[1201:1500,200])/sum(test2.all[,200]) 
sum(test2.all[1501:1800,200])/sum(test2.all[,200]) 

# Juveniles
(sum(test2.all[1:300,200])+sum(test2.all[301:600,200]))/sum(test2.all[,200]) 

# Spawners
(sum(test2.all[901:1200,200])+sum(test2.all[1201:1500,200]))/sum(test2.all[,200]) 

## Wild population
SSD <- wvlambda(IPM2.wild)$w
sum(SSD[1:300]) 
sum(SSD[301:600]) 
sum(SSD[601:900]) 
sum(SSD[901:1200]) 
sum(SSD[1201:1500]) 
sum(SSD[1501:1800]) 

# Juveniles
sum(SSD[1:300])+sum(SSD[301:600]) 

# Spawners
sum(SSD[901:1200])+sum(SSD[1201:1500]) 



### 100% PENALTY
## Stocked population
sum(test3.all[1:300,200])/sum(test3.all[,200]) 
sum(test3.all[301:600,200])/sum(test3.all[,200]) 
sum(test3.all[601:900,200])/sum(test3.all[,200]) 
sum(test3.all[901:1200,200])/sum(test3.all[,200]) 
sum(test3.all[1201:1500,200])/sum(test3.all[,200]) 
sum(test3.all[1501:1800,200])/sum(test3.all[,200]) 

# Juveniles
(sum(test3.all[1:300,200])+sum(test3.all[301:600,200]))/sum(test3.all[,200]) 

# Spawners
(sum(test3.all[901:1200,200])+sum(test3.all[1201:1500,200]))/sum(test3.all[,200]) 

## Wild population
SSD <- wvlambda(IPM3.wild)$w
sum(SSD[1:300]) 
sum(SSD[301:600]) 
sum(SSD[601:900]) 
sum(SSD[901:1200]) 
sum(SSD[1201:1500]) 
sum(SSD[1501:1800]) 

# Juveniles
sum(SSD[1:300])+sum(SSD[301:600]) 

# Spawners
sum(SSD[901:1200])+sum(SSD[1201:1500]) 

