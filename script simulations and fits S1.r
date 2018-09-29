setwd('C:/Users/louvrier/Documents/these/projets/misid/analyses/simulations/gamma=0.1/psi1=0.1/b=0.8')

# important quantities
R = 100 # number of sites
J = 3 # number of replicate surveys or visits
K = 5 # number of years or seasons
psi1 = 0.1 # occupancy probability in first year / season
valp10 = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30) # 
p11 = 0.4 #
delta = 0.8 # b
epsilon = 0.1 # extinction
gamma = 0.1 # colonization



for(p10 in valp10) {
#simulation
########################################################
# simule data from an occupancy model				   #
# wirh false positives 								   #
########################################################
simulfp = function(R, J, K, psi1, p10, p11, delta, epsilon, gamma){
# pre-allocate memory
muZ <- z <- array(dim = c(R, K)) # Expected and realized occurrence
y <- array(NA, dim = c(R, J, K)) # Detection histories

# STATE PROCESS
# occupancy states at once for all sites
# first year
z[,1] <- rbinom(R, 1, psi1) # Initial occupancy state
# following years
for(i in 1:R){ # Loop over sites
	for(k in 2:K){ # Loop over years
		muZ[k] <- z[i, k-1]*(1-epsilon) + (1-z[i, k-1])*gamma # Prob for occ.
		z[i,k] <- rbinom(1, 1, muZ[k])
		}
}

# OBSERVATION PROCESS
# detections/non-detections
   for(i in 1:R){ # loop sites
      for(k in 1:K){ # loop years
      	  for(j in 1:J){ # loop replicated surveys
				# if site unoccupied
      	  		if (z[i,k] == 0){
				   # first, detection (first step in GEPAT)
					y_temp <- rbinom(1,1,p10)
      	  		   # second, assignment (second step in GEPAT)
      	  		    if (y_temp == 0) y[i,j,k] <- 0
      	  		    if (y_temp == 1) y[i,j,k] <- 2
      	  		    } # end if
				# if site occupied
      	  		if (z[i,k] == 1){
				   # first, detection (first step in GEPAT)
					y_temp <- rbinom(1,1,p11)
      	  		   # second, assignment (second step in GEPAT)
      	  		    if (y_temp == 0) y[i,j,k] <- 0
      	  		    if (y_temp == 1) y[i,j,k] <- (1-rbinom(1,1,delta)) + 1
      	  		    } # end if
      	  		  } # end loop replicated surveys
      	  	} # end loop years
     } # end loop sites

# format data 
yy <- matrix(y, R, J*K)

}

# fonction to avoid the explosion of log(0)
logprot <- function(v){
eps <- 2.2204e-016
u <- log(eps) * (1+vector(length=length(v)))
index <- (v>eps)
u[index] <- log(v[index])
u
}


#maximum likelihood function for data with uncertainty
maxlikfp <- function(yy, J, K){

##########################################################################
# likelihood of a dynamic occupancy model accounting for false positive  # 
# formulated as HMM							                             #
##########################################################################

# fonction to estimate -log(likelihood) of a dynamic occupancy model accounting for false positives
devcolextfp <- function(b,data,eff,garb,nh,primary,secondary)
{
# b = vector of parameters on real scale
# data = site histories
# eff = nb of sites with that particular history
# garb = initial states
# nh = nb of sites


#---------------------------------------------
# logit link
psi <- 1/(1+exp(-b[1]))
gamma <- 1/(1+exp(-b[2]))
epsilon <- 1/(1+exp(-b[3]))
p10 <- 1/(1+exp(-b[4]))
p11 <- 1/(1+exp(-b[5]))
delta <- 1/(1+exp(-b[6]))
#---------------------------------------------

#---------------------------------------------
# useful quantities
K <- length(primary)
J2 <- length(secondary) 
J <- J2/K
N <- J * K
#---------------------------------------------

#---------------------------------------------
# initial states matrix
PI <- c(1-psi,psi)

# Transition matrix
Aprimary <- matrix(c((1-gamma),gamma,epsilon,(1-epsilon)),nrow=2,ncol=2,byrow=T)
Asecondary <- matrix(c(1,0,0,1),nrow=2,ncol=2,byrow=T)
A <- array(NA,c(2,2,N))
for (j in primary) A[1:2,1:2,j] <- Aprimary
for (k in secondary[-primary]) A[1:2,1:2,k] <- Asecondary
		
# Probabilities of events (lines) conditionnal to states (columns)
B <- matrix(c(1-p10,1-p11,0,delta*p11,p10,(1-delta)*p11),nrow=3,ncol=2,byrow=T)
#---------------------------------------------

#---------------------------------------------
# calculating -log(likelihood) 
l <- 0
for (i in 1:nh) # loop on sites
{
	oe <- garb[i] + 1 # initial events
	evennt <- data[,i] + 1 # non detections : 1, detections : 2
	ALPHA <- PI*B[oe,]
	for (j in 2:N){ # loop on time
		ALPHA <- (ALPHA %*% A[1:2,1:2,j-1])*B[evennt[j],]
		}
      	l <- l + logprot(sum(ALPHA))*eff[i]
   }
l <- -l
l
#---------------------------------------------

}

##########################################################################
# fitting the model by maximizing the likelihood					     #
##########################################################################
 
nh <- dim(yy)[1] # nb of sites
eff <- rep(1,nh) # number of sites with the specific detection history
data <- t(yy)  
garb <- yy[,1] # initial states

# primary and secondary occasions
primary <- seq(J,J*K,by=J)
secondary <- 1:(J*K)



# because we suspect local minima, we fit the model
# several times using different inits
binit <- runif(6)
tmpmin <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin2 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin3 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- c(0.1,0.1,0.1,p10,0.4,0.8)
tmpmin4 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin5 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin6 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin7 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin8 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- c(0.1,0.1,0.1,p10,0.4,0.8)
tmpmin9 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(6)
tmpmin10 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))

# pick lowest deviance
global.dev <- which.min(c(tmpmin$value,tmpmin2$value,tmpmin3$value,tmpmin4$value,tmpmin5$value,tmpmin6$value,tmpmin7$value,tmpmin8$value, tmpmin9$value, tmpmin10$value))
if (global.dev==1) xx <- tmpmin$par
if (global.dev==2) xx <- tmpmin2$par
if (global.dev==3) xx <- tmpmin3$par
if (global.dev==4) xx <- tmpmin4$par
if (global.dev==5) xx <- tmpmin5$par
if (global.dev==6) xx <- tmpmin6$par
if (global.dev==7) xx <- tmpmin7$par
if (global.dev==8) xx <- tmpmin8$par
if (global.dev==9) xx <- tmpmin9$par
if (global.dev==10) xx <- tmpmin10$par


# back-transorms estimates
psi <- exp(xx[1])/(1+exp(xx[1]))
gamma <- exp(xx[2])/(1+exp(xx[2]))
epsilon <- exp(xx[3])/(1+exp(xx[3]))
p10 <- exp(xx[4])/(1+exp(xx[4]))
p11 <- exp(xx[5])/(1+exp(xx[5]))
delta <- exp(xx[6])/(1+exp(xx[6]))

c(psi, gamma, epsilon, p10, p11, delta)
}


#maximum likelihood function for data WITHOUT uncertainty
maxlik <- function(yy, J, K) {

##########################################################################
# likelihood of a dynamic occupancy model without false positive 		 # 
# formulated as HMM							                             #
##########################################################################

# fonction to estimate -log(likelihood) of a dynamic occupancy model without false positives
devcolext <- function(b,data,eff,garb,nh,primary,secondary)
{
# b = vector of parameters on real scale
# data = site histories
# eff = nb of sites with that particular history
# garb = initial states
# nh = nb of sites

#---------------------------------------------
# logit link
psi <- 1/(1+exp(-b[1]))
gamma <- 1/(1+exp(-b[2]))
epsilon <- 1/(1+exp(-b[3]))
p <- 1/(1+exp(-b[4]))
#---------------------------------------------

#---------------------------------------------
# useful quantities
K <- length(primary)
J2 <- length(secondary) 
J <- J2/K
N <- J * K
#---------------------------------------------

#---------------------------------------------
# initial states
PI <- c(1-psi,psi)

#Transition matrix
Aprimary <- matrix(c((1-gamma),gamma,epsilon,(1-epsilon)),nrow=2,ncol=2,byrow=T)
Asecondary <- matrix(c(1,0,0,1),nrow=2,ncol=2,byrow=T)
A <- array(NA,c(2,2,N))
for (j in primary) A[1:2,1:2,j] <- Aprimary
for (k in secondary[-primary]) A[1:2,1:2,k] <- Asecondary
		
# Events probabiltiy matrix
B <- matrix(c(1,1-p,0,p),nrow=2,ncol=2,byrow=T)
#---------------------------------------------

#---------------------------------------------
# we calulate the log(likelihood)
l <- 0
for (i in 1:nh) # loop on sites
{
	oe <- garb[i] + 1 #initial event
	evennt <- data[,i] + 1 # non detections : 1 and detections : 2
	ALPHA <- PI*B[oe,]
	for (j in 2:N){ # loop on time
		ALPHA <- (ALPHA %*% A[1:2,1:2,j-1])*B[evennt[j],]
		}
      	l <- l + logprot(sum(ALPHA))*eff[i]
   }
l <- -l
l
#---------------------------------------------

}

##########################################################################
# fitting the model by maximizing the likelihood					     #
##########################################################################
 
nh <- dim(yy)[1] # nb of sites
eff <- rep(1,nh) # number of sites with the particular deteciton history
data <- t(yy) # 
garb <- yy[,1] # initial states

# primary and secondary occasions
primary <- seq(J,J*K,by=J)
secondary <- 1:(J*K)

# because we suspect local minima, we fit the model
# several times using different inits
binit <- runif(4)
tmpmin <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin2 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin3 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- c(0.1,0.1,0.1,0.2)
tmpmin4 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin5 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin6 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin7 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin8 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- c(0.1,0.1,0.1,0.2)
tmpmin9 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin10 <- optim(binit,devcolext,NULL,hessian=FALSE,data,eff,garb,nh,primary,secondary,method="BFGS",control=list(trace=1, REPORT=1))

# pick lowest deviance
global.dev <- which.min(c(tmpmin$value,tmpmin2$value,tmpmin3$value,tmpmin4$value,tmpmin5$value,tmpmin6$value,tmpmin7$value,tmpmin8$value, tmpmin9$value, tmpmin10$value))
if (global.dev==1) xx <- tmpmin$par
if (global.dev==2) xx <- tmpmin2$par
if (global.dev==3) xx <- tmpmin3$par
if (global.dev==4) xx <- tmpmin4$par
if (global.dev==5) xx <- tmpmin5$par
if (global.dev==6) xx <- tmpmin6$par
if (global.dev==7) xx <- tmpmin7$par
if (global.dev==8) xx <- tmpmin8$par
if (global.dev==9) xx <- tmpmin9$par
if (global.dev==10) xx <- tmpmin10$par

# back-transorme the estimations 
psi <- exp(xx[1])/(1+exp(xx[1]))
gamma <- exp(xx[2])/(1+exp(xx[2]))
epsilon <- exp(xx[3])/(1+exp(xx[3]))
p <- exp(xx[4])/(1+exp(xx[4]))

c(psi, gamma, epsilon, p)

}


#############################################################################
#						the loop											#
#############################################################################

#number of simulations
H = 500

#vectors with initial simulations
simul12_original <- vector("list",H) #with the 1s: certains and 2s: uncertains 
simul1_original <- vector("list",H) #only the 1s

#vectors with simulations
estim_incert_original <- vector("list", H) #estimation of parameters accounting for misidentification
estim_cert_original <- vector("list", H) #estimation of parameters without misidentification
for (i in 1:H)
{
   simul12_original[[i]] <- simulfp(R, J, K, psi1, p10, p11, delta, epsilon, gamma)
   estim_incert_original[[i]] <- maxlikfp(simul12_original[[i]], J, K)
   simul1_original[[i]] <- simul12_original[[i]]
   simul1_original[[i]][simul1_original[[i]] > 1] <- 0
   estim_cert_original[[i]] <- maxlik(simul1_original[[i]], J, K)
   
}



save(simul12_original, file = paste0('simul12_', as.character(p10), '.RData'))
save(simul1_original, file = paste0('simul1_', as.character(p10),'.RData'))
save(estim_incert_original, file = paste0('estim_incert_original_', as.character(p10),'.RData'))
save(estim_cert_original, file = paste0('estim_cert_original_', as.character(p10),'.Rdata'))

###########################################################
#		accounting for outliers							  #
###########################################################
#psi, gamma, epsilon appear in this order 
replace <- function(O) {


##check for the outliers
x <- matrix(0, nrow = 3, ncol = 1)
y <- matrix(0, nrow = 3, ncol = 1)
z <- matrix(0, nrow = 3, ncol = 1)

for(i in 1:H)
{
x[i] <- O[[i]][1]
y[i] <- O[[i]][2]
z[i] <- O[[i]][3]
}

#by checkin x we can see that x above 0.20 are the outliers
outl <- which(x[]>0.25 | x[]<0.01 | y[]>0.25 | y[]<0.01 | z[]>0.25 | z[]<0.01  )

N <- length(outl)

#we simulate again for data that are going to replace
#number of simulations
H = N

#vectors of simulations
simul12_extr <- vector("list",H) 
simul1_extr <- vector("list",H) 

#vectors of estimations
estim_incert_extr <- vector("list", H)  #estimation of parameters accounting for misidentification
estim_cert_extr <- vector("list", H) #estimation of parameters without misidentification

for (i in 1:H)
{
   simul12_extr[[i]] <- simulfp(R, J, K, psi1, p10, p11, delta, epsilon, gamma)
   estim_incert_extr[[i]] <- maxlikfp(simul12_extr[[i]], J, K)
   simul1_extr[[i]] <- simul12_extr[[i]]
   simul1_extr[[i]][simul1_extr[[i]] > 1] <- 0
   estim_cert_extr[[i]] <- maxlik(simul1_extr[[i]], J, K)
   
}

extr <- vector("list", 5)

extr[[1]] <- outl
extr[[2]] <- simul12_extr
extr[[3]] <- simul1_extr 
extr[[4]] <- estim_incert_extr
extr[[5]] <- estim_cert_extr

return(extr)

}


#for psi 
extr1 <- replace(estim_incert_original)

#final dataset 
simul12_final <- simul12_original
simul1_final <- simul1_original
estim_incert_final <- estim_incert_original
estim_cert_final <- estim_cert_original

#we replace the outliers
simul12_final[extr1[[1]]] <- extr1[[2]]
simul1_final[extr1[[1]]] <- extr1[[3]]
estim_incert_final[extr1[[1]]] <- extr1[[4]]
estim_cert_final[extr1[[1]]] <- extr1[[5]]


save(simul12_final, file = paste0('simul12_final_', as.character(p10),'.RData'))
save(simul1_final, file = paste0('simul1_final_', as.character(p10),'.RData'))
save(estim_incert_final, file = paste0('estim_incert_final_', as.character(p10),'.RData'))
save(estim_cert_final, file = paste0('estim_cert_final_', as.character(p10),'.Rdata'))


##########################################################
#					BIAIS	AND  MSE					 #
##########################################################

#extraction
H=500

xincert <- matrix(0, nrow = H, ncol =3) #matrix of the 3 ecological parameters estimates from the MCU for the 500 simulations
xcert <- matrix(0, nrow = H, ncol =3) #matrix of the 3 ecological parameters estimates from the MC for the 500 simulations
diff_incert <- matrix(0, nrow = H, ncol = 3) #matrix of difference between the estimates of the 3 parameters of MCU and their true value for the 500 simulations
diff_cert <- matrix(0, nrow = H, ncol = 3) #matrix of difference between the estimates of the 3 parameters of MC and their true value for the 500 simulations
param_incert = c(psi1, gamma, epsilon) #true values 
param_cert = c(psi1, gamma, epsilon) #true values


#########################################################################################################################################################################################################
### 			EXTRACTION

#MCU
for(i in 1:H){#loop on simulations
	for(j in 1:length(param_incert)){ #loop on parameters
xincert[i,j] <- estim_incert_final[[i]][j] #parameters estimates
diff_incert[i,j] <- (xincert[i,j]-param_incert[j]) #differences between estimates and true values
	}
}

#MC
for(i in 1:H){#loop on simulations
	for(j in 1:length(param_cert)){ #loop on parameters
xcert[i,j]<-estim_cert_final[[i]][j] #paremeters estimates
diff_cert[i,j] <- (xcert[i,j]-param_cert[j]) #differences between estimates and true values
	}
}


#calculating bias and MSE for MCU
for(i in 1:length(param_incert)) { #loop on parameters
Biais_incert[[as.character(p10)]][i] <- ((mean(xincert[,i])-param_incert[i])/param_incert[i])*100 #bias
MSE_incert[[as.character(p10)]][i] <- mean((diff_incert[,i])^2) #MSE
}

#calculating bias and MSE for MC
for(i in 1:length(param_cert)) { #loop on parameters
Biais_cert[[as.character(p10)]][i] <- (mean(xcert[,i])-param_cert[i])/param_cert[i]*100 #bias
MSE_cert[[as.character(p10)]][i] <- mean((diff_cert[,i])^2) #MSE 
}

save(Biais_incert, file ='Biais_incert.RData')
save(Biais_cert, file = 'Biais_cert.RData')
save(MSE_incert, file = 'MSE_incert.RData')
save(MSE_cert, file = 'MSE_cert.RData')