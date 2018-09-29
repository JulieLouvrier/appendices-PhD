#----- STEP 1: simulate data from a static occupancy model with misid and heterogeneity


SIMUL_FIT = function(C){
pi = C[1] # first value of vector C, pi is the proportion of sites A
pA11 = p11 # true positive detection probability for sites A
pB11 = p11 # true positive detection probability for sites B
pA10 = C[2] # second value of vector C, false positive detection probability for sites A
pB10 = C[3] # third value of vector C, false positive detection probability for sites B

# ns = number of sites
# J = number of surveys
# psi1 = initial occupancy probability here we are dealing with a static occupancy model so psi1 = psi 
# p10 = probability of false positive detection
# p11 = probability of false positive detection
# delta = probability to classify a true positive detection as certain

#function to simulate ns sites with species misidentification
#here we are using the following function which we are going to run twice with pi*ns sites A and (1-pi)*ns sites B and combine both data sets 

# simul function
simul_occ = function(ns, J, psi1, p10, p11,delta){

# pre-allocate memory
y <- array(NA, dim = c(ns, J)) # Detection histories

# STATE PROCESS
# occupancy states at once for all sites
z <- rbinom(ns, 1, psi1) 

# OBSERVATION PROCESS 
# detections/non-detections
   for(i in 1:ns){ # loop sites
      	  for(j in 1:J){ # loop replicated surveys
				# if site unoccupied
      	  		if (z[i] == 0){
				   # first, false positive detection 
					y_temp <- rbinom(1,1,p10)
      	  		   # second, assignment 
      	  		    if (y_temp == 0) y[i,j] <- 0
      	  		    if (y_temp == 1) y[i,j] <- 2
      	  		    } # end if
				# if site occupied
      	  		if (z[i] == 1){
				   # first, true positive detection (first step in GEPAT)
					y_temp <- rbinom(1,1,p11)
      	  		   # second, assignment (second step in GEPAT)
      	  		    if (y_temp == 0) y[i,j] <- 0
      	  		    if (y_temp == 1) y[i,j] <- (1-rbinom(1,1,delta)) + 1
      	  		    } # end if
      	  		  } # end loop replicated surveys
     } # end loop sites

# format data 
yy <- matrix(y, ns, J)

}


#----- STEP 2: fit a static occupancy model with misid and heterogeneity to data
#maxlikhet function for data with uncertainty and accounting for heterogeneity
maxlikhet <- function(data) {

####################################################################################
# likelihood function for data with uncertainty  accouting for heterogeneity       # 
# formulated as HMM										                           #
####################################################################################

nh = nrow(data) # nb sites
J = ncol(data) # nb surveys
eff = rep(1,nh) # nb sites with this particular history
garb = data[,1] # initial states
data = t(data) # transpose


# Deviance of static occupancy model with false positives and heterogeneity
dev_occufphet <- function(b,data,eff,garb,nh,nocc)
{
# b = vector of parameters on real scale
# data = site histories
# eff = nb of sites with that particular history
# garb = initial states
# nh = nb of sites

#---------------------------------------------
# logit link
pi <- 1/(1+exp(-b[1]))
psi1 <- 1/(1+exp(-b[2]))
pA10 <- 1/(1+exp(-b[3]))
pA11 <- 1/(1+exp(-b[4]))
pB10 <- 1/(1+exp(-b[5]))
pB11 <- 1/(1+exp(-b[6]))
delta <- 1/(1+exp(-b[7]))
#---------------------------------------------

#---------------------------------------------
# initial states matrices
PI1 <- c(pi,1-pi) # allocation of class A and B
PI2 <- matrix(c(1-psi1,0,psi1,0,0,1-psi1,0,psi1),nrow=2,ncol=4,byrow=T) # occupancy probability in sites of class A and sites of class B
PI <- PI1 %*% PI2 #initial states matrix

# transition matrix, static model so epsilon = gamma = 0
A <- diag(4)

# observation matrix
# first, detections probabilities
B1 <- matrix(c(
1-pA10,0,pA10, # false positives for sites A
1-pB10,0,pB10, # false positives for sites B
1-pA11,pA11,0, # true positives for sites A
1-pB11,pB11,0),nrow=4,ncol=3,byrow=T) # true positives for sites B

# second, classification into certain or uncertain
B2 <- matrix(c(
1,0,0,
0,delta,1-delta, #probability of classifying a true positive detection into certain
0,0,1),nrow=3,ncol=3,byrow=T)
B <- t(B1 %*% B2) # detection matrix
#---------------------------------------------
logprot <- function(v){
# avoid explosion of log(0)
eps <- 2.2204e-016
u <- log(eps) * (1+vector(length=length(v)))
index <- (v>eps)
u[index] <- log(v[index])
u
}
#---------------------------------------------
# calculate -log(lik) and AIC
k= 7
l <- 0
for (i in 1:nh) # loop on sites
{
	oe <- garb[i] + 1 # first obs
	evennt <- data[,i] + 1 # non-det/det -> 1/2
	ALPHA <- PI*B[oe,]
	for (j in 2:nocc){ # loop on time
		ALPHA <- (ALPHA %*% A)*B[evennt[j],]
		}
      	l <- l + logprot(sum(ALPHA))*eff[i]
   }
l <- -l
l

#---------------------------------------------
}



# because we suspect local minima, we fit the model
# several times using different inits
set.seed(1)
# run 1 with inits = values that were used to simulate the data
binit <- c(pi,psi1,pA10,pA11,pB10,pB11,delta)
tmpmin1 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 2 with random inits
binit <- runif(7)
tmpmin2 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 3 with random inits
binit <- runif(7)
tmpmin3 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 4 with random inits
binit <- runif(7)
tmpmin4 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 5 with random inits
binit <- runif(7)
tmpmin5 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 6 with random inits
binit <- runif(7)
tmpmin6 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 7 with random inits
binit <- runif(7)
tmpmin7 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 8 with random inits
binit <- runif(7)
tmpmin8 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 9 with random inits
binit <- runif(7)
tmpmin9 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
# run 10 with random inits
binit <- runif(7)
tmpmin10 = optim(binit, dev_occufphet,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))

# pick lowest deviance
res = c(tmpmin1$value, tmpmin2$value, tmpmin3$value, tmpmin4$value, tmpmin5$value, tmpmin6$value, tmpmin7$value, tmpmin8$value, tmpmin9$value,tmpmin10$value)
ind = which.min(res)


# display results from model with lowest deviance
# we don't know in which order the estimates are produced, 
# we order them in order to make sure each parameters corresponds to its estimates 
name=paste('tmpmin',ind,sep='')
list(estim = cbind(c(true=sort(c(pi,psi1,pA10,pA11,pB10,pB11,delta)),estimates=sort(round(1/(1+exp(-eval(parse(text=name))$par)),3)))), AIC = 2*(eval(parse(text=name))$value)+(2*7), l=(eval(parse(text=name))$value))
}


#----- STEP 3: fit a static occupancy model with misid and without heterogeneity to data
# maximum likelihood function for data with uncertainty not accouting for heterogeneity
maxlikfp <- function(yy, J){

#######################################################################################
# likelihood function for data with uncertainty not accouting for heterogeneity       # 
# formulated as HMM										                              #
#######################################################################################

# Deviance of static occupancy model with false positives only
devcolextfp <- function(b,data,eff,garb,nh,nocc)
{

#---------------------------------------------
# logit link
psi1 <- 1/(1+exp(-b[1])) # initial occupancy probability
p10 <- 1/(1+exp(-b[2])) # probability of false positives
p11 <- 1/(1+exp(-b[3])) # probability of true positives
delta <- 1/(1+exp(-b[4])) # probability of classifying a true positive detection as certain
#---------------------------------------------

#---------------------------------------------
# useful parameters
J <- length(nocc) # J = number of surveys
#---------------------------------------------

#---------------------------------------------
# matrix of initial states
PI <- c(1-psi1,psi1) #initial occupancy probability

# transition matrix, static model so epsilon = gamma = 0
A <- diag(2)
		
# Probabilités des événements (lignes) conditionnellement aux états (colonnes)
B <- matrix(c(1-p10,1-p11,0,delta*p11,p10,(1-delta)*p11),nrow=3,ncol=2,byrow=T)
#---------------------------------------------
logprot <- function(v){
# avoid explosion of log(0)
eps <- 2.2204e-016
u <- log(eps) * (1+vector(length=length(v)))
index <- (v>eps)
u[index] <- log(v[index])
u
}
#---------------------------------------------
# on forme la vraisemblance, ou plutot la log(vraisemblance) 
# on linearise car c'est plus stable numeriquement de maximiser une somme qu'un produit
# en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
l <- 0
k=4
for (i in 1:nh) # loop on sites
{
	oe <- garb[i] + 1 # initial events
	evennt <- data[,i] + 1 # non detections become 1 and detections become 2
	ALPHA <- PI*B[oe,]
	for (j in 2:nocc){ # loop on surveys
		ALPHA <- (ALPHA %*% A)*B[evennt[j],]
		}
      	l <- l + logprot(sum(ALPHA))*eff[i]
   }
l <- -l
l

#---------------------------------------------

}

##########################################################################
# fitting static occupancy model by maximizing the likelihood		     #
##########################################################################
 
nh <- dim(yy)[1] # number of sites
J = ncol(yy) # number of surveys
eff <- rep(1,nh) # nb sites with this particular history
data <- t(yy) # transposes
garb <- yy[,1] # initial states

# several times using different inits
set.seed(1)

# optimisation with 10 sets of initial values  
binit <- runif(4)
tmpmin1 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin2 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin3 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- c(psi1,p10, p11, delta)
tmpmin4 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin5 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin6 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin7 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin8 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))
binit <- runif(4)
tmpmin9 <- optim(binit,devcolextfp,NULL,hessian=FALSE,data,eff,garb,nh,J,method="BFGS",control=list(trace=1, REPORT=1))

# parametres estimes selon le run avec la deviance la plus petite
global.dev <- which.min(c(tmpmin1$value,tmpmin2$value,tmpmin3$value,tmpmin4$value,tmpmin5$value,tmpmin6$value,tmpmin7$value,tmpmin8$value, tmpmin9$value))
if (global.dev==1) xx <- tmpmin1$par
if (global.dev==2) xx <- tmpmin2$par
if (global.dev==3) xx <- tmpmin3$par
if (global.dev==4) xx <- tmpmin4$par
if (global.dev==5) xx <- tmpmin5$par
if (global.dev==6) xx <- tmpmin6$par
if (global.dev==7) xx <- tmpmin7$par
if (global.dev==8) xx <- tmpmin8$par
if (global.dev==9) xx <- tmpmin9$par


ind = global.dev 
name=paste('tmpmin',ind,sep='')

# display results from model with lowest deviance
list(estim = cbind(c(true=sort(c(psi1,p10,p11,delta)),estimates=sort(round(1/(1+exp(-eval(parse(text=name))$par)),3)))),AIC = (2*(eval(parse(text=name))$value)+ (2*4)), l=(eval(parse(text=name))$value))


}


###########################################################################
#### SIMULATIONS 
###########################################################################


#simulations vectors
simul12_A = vector("list",H) #with certain (1) and uncertain (2) data
simul12_B = vector("list",H) #with certain (1) and uncertain (2) data
simul12_AB = vector("list",H) #with certain (1) and uncertain (2) data


#estimations vectors
estim12_AB = vector("list", H) # estimations of parameters from the model accounting for misidentification and heterogeneity
estim12_C = vector("list", H) # estimations of parameters from the model accounting for misidentification only
delta_AIC = matrix(0,ncol = 1, nrow = H) # difference of AIC
l_AB = matrix(0,ncol = 1, nrow = H) # likelihood of the model accounting for misidentification and heterogeneity
l_C = matrix(0,ncol = 1, nrow = H) # likelihood of the model accounting for misidentification only
mean_DAIC = matrix(0, ncol = 1, nrow = 1) #to calculate the mena difference of AICs for the H simulations/fits



for (i in 1:H)
{
simul12_A[[i]] <- simul_occ((round(pi*R)), J, psi1, pA10, pA11, delta) #simulation of pi*R sites A
simul12_B[[i]] <- simul_occ((round((1-pi)*R)), J, psi1, pB10, pB11, delta) # simulation of (1-pi)*R sites B
# mix both
simul12_AB[[i]] <- rbind(simul12_A[[i]],simul12_B[[i]])
}

#maxlik
for (i in 1:H)
{
estim12_AB[[i]]<- maxlikhet(simul12_AB[[i]]) # fitting the model accounting for misidentification and heterogeneity on the dataset
}

for (i in 1:H){
estim12_C[[i]]<- maxlikfp(simul12_AB[[i]], J) # fitting the model accounting for misidentification only on the dataset
}

for (i in 1:H){
l_C[i] = estim12_C[[i]]$l # likelihood of MMO
l_AB[i] = estim12_AB[[i]]$l # likelihood of MMH
delta_AIC[i] = estim12_C[[i]]$AIC - estim12_AB[[i]]$AIC #delta_AIC
mean_DAIC = mean(delta_AIC) #mean of delta_AICs
}

#display results
return(list(simul12_AB = simul12_AB, estim12_AB = estim12_AB,estim12_C = estim12_C, delta_AIC = delta_AIC, mean_DAIC = mean_DAIC))
}



