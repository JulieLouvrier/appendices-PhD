##3 different scenarios
# pi = 0.2, pB10 = 0.7 and pA10 = 0.1
# pi = 0.5, pB10 = 0.8 and pA10 = 0.1
# pi = 0.8, pB10 = 0.9 and pA10 = 0.1

# muA = PA10
# muB = PB10

# heterogeneity parameter
# mu = (pi*muA)+((1-pi)*muB)

# variance between components
# sigm_car = (pi*(muA-mu)^2) + ((1-pi)*(muB-mu)^2)

## heterogeneity coefficient
# nu = sigm_car/(mu*(1-mu))

#MMH: model dealing with misidentification and heterogeneity
#MMO: model dealing with misidentification only
source('associatedfunctions.r') # calling the function to simulate data with heterogeneity and misidentification and fit them with MMMH and MMO

####################################################################################################################
psi1 = 0.8 # initial occupancy probability
delta = 0.7 # probability to classify a true detection as certain
p11=pA11=pB11 = 0.5 # true positive detection
pi=c(0.2,0.5,0.8) # proportion of site A 
pA10 = c(0.1,0.1,0.1) # false positive probability on sites of class A
pB10 = c(0.7,0.8,0.9) # false positive probability on sites of class B
p10 = (pi*pA10) + ((1-pi)*pB10)
comb = matrix(0,nrow = 3, ncol = 3) # combination of three parameters
####################################################################################################################
##########################################################################
# heterogeneity formula											         # 
##########################################################################

# allocating parameters values
comb[,1] <-pi
comb[,2] <- pA10
comb[,3] <- pB10

#function to calculate the heterogeneity oarameter
heter = function(M){
mu = (M[1]*M[2])+((1-M[1])*M[3])
sigm_car = (M[1]*(M[2]-mu)^2) + ((1-M[1])*(M[3]-mu)^2)
nu = sigm_car/(mu*(1-mu))

return(c(mu = mu, sigm_car = sigm_car, nu = nu))
}


calc_heter <- apply(comb,1,heter) # application of function
row.names(calc_heter) <- c('mu','sigm_car','nu')

##########################################################################
# SIMULATION													         # 
##########################################################################
R = 100 # number of sites
J = 10 # number of replicate surveys or visits
p11 = 0.5 # probability of true detections
# allocating parameters values
comb[,1] <-pi
comb[,2] <- pA10
comb[,3] <- pB10

#number of simulations
H = 200

Simul <- apply(comb, 1, SIMUL_FIT) # application of function 
save(Simul, file = 'Simul.RData')

##########################################################################
# BIAS AND MSE													         # 
##########################################################################

#to estimate bias and mse from the Simul.RData file: 
#Simul[[1]]$estim12_AB[[1]]$estim[c(1:7)] are the true values set for the 7 parameters for the MMH (pi, psi1, pA11, pB11, pA10, pB10, delta) for the first scenario (Simul[[1]])
#Simul_3eme_essai[[1]]$estim12_C[[1]]$estim[c(1:4)] are the values for the 4 parameters we set for the MMO (psi1, p11, p10, delta)
#then for each simulation (for i n 1:H) extract the value of the estimates from both models (Simul[[1]]$estim12_AB[[i]]$estim[j+7] for MMH and Simul[[1]]$estim12_C[[i]]$estim[j+4] for MMO)
#also calculate the difference of value between estimates and true values
#from then calculate the bias and MSE

#estimations of bias and MSE for the first scneario
#extraction
H=200

#preparing matrices
xAB_1 <- matrix(0, nrow = H, ncol =7) #matrix of estimates the 7 parameters of MMH for the 200 simulations
xC_1 <- matrix(0, nrow = H, ncol =4) #matrix of estimates the 4 parameters of MMO for the 200 simulations
diff_AB_1 <- matrix(0, nrow = H, ncol = 7) #matrix of difference between the estimates of the 7 parameters of MMH and their true value for the 200 simulations
diff_C_1 <- matrix(0, nrow = H, ncol = 4) #matrix of difference between the estimates of the 4 parameters of MMO and their true value for the 200 simulations
param_AB_1 = c(Simul[[1]]$estim12_AB[[1]]$estim[c(1:7)]) #true values set for the 7 parameters for the MMH (pi, psi1, pA11, pB11, pA10, pB10, delta) for the first scenario (Simul[[1]])
param_C_1 = c(Simul[[1]]$estim12_C[[1]]$estim[c(1:4)]) # values for the 4 parameters we set for the MMO (psi1, p11, p10, delta)


#########################################################################################################################################################################################################
### 			EXTRACTION

#MMH
for(i in 1:H){ #loop on simulations
	for(j in 1:length(param_AB_1)){ #loop on parameters
xAB_1[i,j] <- Simul[[1]]$estim12_AB[[i]]$estim[j+7] #parameter estimates
diff_AB_1[i,j] <- (xAB_1[i,j]-param_AB_1[j]) #difference 

	}
}

#MMO
for(i in 1:H){#loop on simulations
	for(j in 1:length(param_C_1)){ #loop on parameters
xC_1[i,j]<-Simul[[1]]$estim12_C[[i]]$estim[j+4] #parameter estimates
diff_C_1[i,j] <- (xC_1[i,j]-param_C_1[j]) #difference
	}
}

#preparing matrices
Biais_1_AB <- matrix(0, nrow = 7, ncol =1) #bias for the 7 parameters of the MMH
MSE_1_AB <- matrix(0, nrow = 7, ncol =1) #MSE for MMH
Biais_1_C <- matrix(0, nrow = 4, ncol =1) #bias for MMO
MSE_1_C <- matrix(0, nrow = 4, ncol =1) #MSE for MMO


#calculating bias and MSE for MMH
for(i in 1:length(param_AB_1)) { #loop on parameters
Biais_1_AB[i,] <- ((mean(xAB_1[,i])-param_AB_1[i])/param_AB_1[i])*100 #Bias
MSE_1_AB[i,] <- sqrt(mean((diff_AB_1[,i])^2) #MSE
}

#calculaing bias and MSE for MMO
for(i in 1:length(param_C_1)) {
Biais_1_C[i,] <- (mean(xC_1[,i])-param_C_1[i])/param_C_1[i]*100
MSE_1_C[i,] <- sqrt(mean((diff_C_1[,i])^2))
}

#combininb bias and MSE
BMAB_1 <- cbind(Biais_1_AB,MSE_1_AB) #for MMH
BMC_1 <- cbind(Biais_1_C,MSE_1_C) #for MMO

rownames(BMAB_1) <- c('pA10','pi','pB11','pA11','pB10','delta','psi1')
rownames(BMC_1) <- c('p11','p10','delta','psi1')

#here we were interested in the bias on psi1, after it is easy to extract this value

