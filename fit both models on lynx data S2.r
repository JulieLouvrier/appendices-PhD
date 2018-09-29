setwd('C:/Users/louvrier/Documents/these/projets/misid/analyses/modèle TOTAL avec incertitude et covariables/analyses Marc/nouvel essai avec données marc et mes c3')

###### (1) Reading in the detection and the covariate data and doing some data management

# Read in the data and put it into 3D array
ldat <- read.table(file = "lynx_data_xmarc_161110.txt", header = T)
str(ldat) # lynx data
ldat_o <- ldat[order(ldat$ID),]

aggregate(ldat_o$country, by = list(ldat_o$cname), FUN = mean)
tapply(ldat_o$country, ldat_o$cname, mean)
    # Austria      France     Germany       Italy    Slovenia Switzerland 

          # 4           1           0           2           5           3 

# But also want country in analysis to add up z's
table(ldat_o$cname)
    # Austria      France     Germany       Italy    Slovenia Switzerland 
        # 607         449         103         575          82         367 

table(as.numeric(ldat_o$country))
  # 0   1   2   3   4   5 
# 103 449 575 367 607  82 


# Select detection data   # Data from 1995 until 2013 (19 years)
yy <- as.matrix(ldat_o[,19:75])
str(yy)
table(yy)         # Detection frequencies per occasion !
yy[yy > 1] <- 1   # Turn to binary detection/nondetection data



# Define function to put data into 3D array
reformat.fn <- function(D2.data = yy, nrep = 3, nyear = 19) {
# This functions reformats detection/nondetection data from a 
# balanced 2D format into a balanced 3D format for nyear = number of years
# nrep is the number of reps per year

nsite <- nrow(yy)			# nsites
y <- array(NA, dim = c(nsite, nrep, nyear))

for(i in 1:nyear){
   y[,,i] <- D2.data[,(((i-1)*nrep+1):(i*nrep))]
}
y				# Simple output
}

# Execute function
y <- reformat.fn(D2.data = yy)
str(y)
 # int [1:2183, 1:3, 1:21] NA NA NA NA NA NA NA NA NA NA ...

# Sum check
sum(yy, na.rm = T)   ;   sum(y, na.rm = T)     # yields 2188 both times



### Grab, standardize and bundle the covariates for the analysis
elev <- ldat_o[,'elev']
forest <- ldat_o[,'forest']
humdens <- ldat_o[,'humdens'] /  100           # Express density as per 1km2
(cov.mn <- apply(cbind(elev, forest, humdens), 2, mean) )
(cov.sd <- apply(cbind(elev, forest, humdens), 2, sd) )

# Put elevation (linear and squared), forest and humdens (linear and squared) into a matrix
# Old version with all normalised
# str(cov <- cbind(elev = (elev - cov.mn[1]) / cov.sd[1], elev2 = ((elev - cov.mn[1]) / cov.sd[1])^2, forest = (forest - cov.mn[2]) / cov.sd[2], humdens = (humdens - cov.mn[3]) / cov.sd[3], humdens2 = ((humdens - cov.mn[3]) / cov.sd[3])^2) )

# New version with all normalised except for humdens, which is special
humdens <- (ldat_o[,'humdens'] /  100) / max(ldat_o[,'humdens'] /  100) -  0.03030104 # Centered on mean
str(cov <- cbind(elev = (elev - cov.mn[1]) / cov.sd[1], elev2 = ((elev - cov.mn[1]) / cov.sd[1])^2, forest = (forest - cov.mn[2]) / cov.sd[2], humdens = humdens, humdens2 = humdens^2) )


### Bundle covariate for distance to nearest release site
distrel.orig <- as.matrix(ldat_o[,79:99]) / 1000    # Express in units of 100 km
(mn.distrel <- mean(distrel.orig))
(sd.distrel <- sd(distrel.orig))
distrel <- (distrel.orig - mn.distrel) / sd.distrel
matplot(1994:2014, t(distrel.orig), type = 'l', lty = 1, xlab = 'Year', ylab = 'Distance to nearest release site')    # Funny pic


### Covariate for season (= occasion or Nov/Dec, Jan/Feb, Mar/Apr)
season <- rep(1:3, 21)


### Effort covariate: presence of a 'network' (nw) of observers, 1 (little) -- 3 (much)
nw <- cbind(ldat_o$searcheffort99, ldat_o$searcheffort99, ldat_o$searcheffort99, 
            ldat_o$searcheffort99, ldat_o$searcheffort99, ldat_o$searcheffort99,  
			ldat_o$searcheffort04, ldat_o$searcheffort04,ldat_o$searcheffort04, 
			ldat_o$searcheffort04, ldat_o$searcheffort04,
			ldat_o$searcheffort09, ldat_o$searcheffort09, ldat_o$searcheffort09,
			ldat_o$searcheffort09, ldat_o$searcheffort09,
			ldat_o$searcheffort14, ldat_o$searcheffort14, ldat_o$searcheffort14, 
			ldat_o$searcheffort14, ldat_o$searcheffort14)
nw <- nw + 1     # Factors must be coded 1 upwards !!



### Autocovariate and the neighbourhood information for the spatial part of the model

# Observed (detection-naive) value of autocovariate

# Observed values of the presence/absence matrix
(obsZ <- apply(y, c(1, 3), max, na.rm = TRUE))
obsZ[obsZ == "-Inf"] <- NA
dimnames(obsZ) <- list(ldat_o$cname, as.character(1994:2014))


# neighbourhood info

load('autocov.RData')

datju12 <- ldat_o
yju12 <- y
obsZju <- obsZ



#----------------------------------------------------------------------------------
#MODEL WITHOUT MISID
#----------------------------------------------------------------------------------


{

#### ----- Bundle data ---------
country <- datju12$country
country[country ==0] <- 6       # Germany is 6 now
# Frankreich = 1, Italien = 2, Schweiz = 3, Österreich = 4, Slowenien = 5, Germany = 6
str(bdata <- list(y = yju12, nsite = dim(yju12)[1], nrep = dim(yju12)[2], nyear = dim(yju12)[3], 
        country = country, nw = nw, cov = cov, distrel = distrel, D = autocov))

zst <- obsZju         # Inits for z for this data set


#### --------- Specify model in BUGS language ----------
# ------------------------------------------------------
sink("Dynocc9_nomisID.txt")
cat("
model {

### Specify priors
alpha.lpsi1 <- logit(alpha.psi1)
alpha.lphi <- logit(alpha.phi)
alpha.lgamma <- logit(alpha.gamma)
alpha.psi1 ~ dbeta(1, 1)
alpha.phi ~ dbeta(1, 1)
alpha.gamma ~ dbeta(1, 1)
for(i in 1:6){                      # Country-specific intercepts for p only (fix)
  alpha.lp[i] <- logit(alpha.p[i])
  alpha.p[i] ~ dbeta(1, 1)
}

beta.lpsi1.f ~ dnorm(0, 0.1) 		#effect of forest on psi1
beta.lphi.f ~ dnorm(0, 0.1)			#effect of forest on phi
beta.lgamma.elev ~ dnorm(0, 0.1)	#effect of elevation on gamma
beta.lgamma.hdens ~ dnorm(0, 0.1)	#effect of hdens on gamma
beta.lp.elev ~ dnorm(0, 0.1)		#effect of elev on detection
beta.lp.f ~ dnorm(0, 0.1)			#effect of forest on detection


beta.lpsi1.rd ~ dnorm(0, 0.1)        # Regression coefs for distance to release site
beta.lgamma.rd ~ dnorm(0, 0.1)

betaDphi ~ dnorm(0, 0.1)             # Regression coefficients of autocovariate
betaDgamma ~ dnorm(0, 0.1)

lp.season[1] <- 0                    # Effects of season on detection
for(i in 2:3){
  lp.season[i] ~ dnorm(0, 0.1)
}

lp.net[1] <- 0                       # Network effects in p
lp.net[2] ~ dnorm(0, 0.1)
lp.net[3] ~ dnorm(0, 0.1)

# Random year-specific effects
for (k in 1:(nyear-1)){
  lphi.year[k] ~ dnorm(0, tau.lphi) # means are country-specific intercepts
  lgamma.year[k] ~ dnorm(0, tau.lgamma)
  lp.year[k] ~ dnorm(0, tau.lp)
}
lp.year[nyear] ~ dnorm(0, tau.lp)

tau.lphi <- pow(sd.lphi, -2)
tau.lgamma <- pow(sd.lgamma, -2)
tau.lp <- pow(sd.lp, -2)
tau.eps <- pow(sd.eps, -2)
sd.lphi ~ dnorm(0, 0.2)I(0.01,)      # curve(dnorm(x, 0, sqrt(1/0.2)), 0, 10)
sd.lgamma ~ dnorm(0, 0.2)I(0.01,)    # Truncated effectively at 0
sd.lp ~ dnorm(0, 0.2)I(0.01,)
sd.eps ~ dnorm(0, 0.2)I(0.01,)

### Ecological submodel
for (i in 1:nsite){
  z[i,1] ~ dbern(psi1[i])
  psi1[i] <- 1 / (1 + exp(-lpsi1[i]))
#   lpsi1[i] <- min(200, max(-200, lpsi1.tmp[i]))
  lpsi1[i] <- alpha.lpsi1 + beta.lpsi1.f * cov[i,3] + beta.lpsi1.rd * distrel[i,1]
   
  for (k in 2:nyear){
    z[i,k] ~ dbern(muZ[i,k])
    muZ[i,k]<- z[i,k-1] * phi[i, k-1] + (1-z[i,k-1]) * gamma[i, k-1]
    phi[i, k-1] <- 1 / (1 + exp(-lphi[i, k-1]))
    lphi[i, k-1] <- alpha.lphi + lphi.year[k-1] + 
	beta.lphi.f * cov[i,3] + betaDphi * D[i,k-1]
    gamma[i, k-1] <- 1 / (1 + exp(-lgamma[i, k-1]))
    lgamma[i, k-1] <- alpha.lgamma + lgamma.year[k-1] + 
	beta.lgamma.elev * cov[i,1] + beta.lgamma.hdens * cov[i,4] + beta.lgamma.rd * distrel[i,k-1] + betaDgamma * D[i,k-1]
  } #k
} #i

### Observation model (also generate new data set for use in PPD checking of GOF)
for (i in 1:nsite){
  for (k in 1:nyear){
    for (j in 1:nrep){
      y[i,j,k] ~ dbern(muy[i,j,k])
      muy[i,j,k] <- z[i,k] * p[i,j,k]
      p[i,j,k] <- 1 / (1 + exp(-lp[i,j,k]))
#       lp[i,j,k] <- min(200, max(-200, lp.tmp[i,j,k]))
      lp[i,j,k] <- alpha.lp[country[i]] + lp.year[k] + lp.net[nw[i,k]] + lp.season[j] 
         	  + beta.lp.elev * cov[i,1] + beta.lp.f * cov[i,3] + eps[i,k]
      y.new[i,j,k] ~ dbern(muy[i,j,k])	# Generate replicate data set under same model
    } #j
    eps[i,k] ~ dnorm(0, tau.eps)            # random site.year effect in detection
  } #k
} #i

### Posterior predictive check of GOF based on detection frequency per year
for (i in 1:nsite){
  for (k in 1:nyear){
    sum.y[i,k] <- sum(y[i,,k])	                                      # Observed actual data
    eval[i,k] <- max(0.01, sum(muy[i,,k]))                            # Expected data (sum of muy)
    E[i,k] <- pow((sum.y[i,k] - eval[i,k]),2) / (eval[i,k] + 0.01)    # Chi2 discrepancy measure actual data set
    sum.y.new[i,k] <- max(0.01, sum(y.new[i,,k]))                     # Observed replicate data
    E.new[i,k] <- pow((sum.y.new[i,k] - eval[i,k]),2) / (eval[i,k] + 0.01)# Chi2 discrepancy measure 'perfect' data set
    # detfreq.ratio[i,k] <- sum.y[i,k] / sum.y.new[i,k]               # Ratio of detection frequencies observed, expected ----- DONT NEED NOW ---------
    # df.diff[i,k] <- sum.y[i,k] - sum.y.new[i,k]                     # Difference ..... (same)  ---- DONT REALLY NEED THIS ----------
    } #k
} #i
fit <- sum(E[,])                    # Sum up to overall discrepancy measure
fit.new <- sum(E.new[,])

### Number of occupied quads per country and year
for (k in 1:nyear){
  for (i in 1:nsite){
    tmpF[i,k] <- equals(country[i], 1) * z[i,k]      
    tmpI[i,k] <- equals(country[i], 2) * z[i,k]
    tmpCH[i,k] <- equals(country[i], 3) * z[i,k]
    tmpAu[i,k] <- equals(country[i], 4) * z[i,k]
    tmpSlo[i,k] <- equals(country[i], 5) * z[i,k]
    tmpGer[i,k] <- equals(country[i], 6) * z[i,k]
  }
  n.occ.F[k] <- sum(tmpF[,k])       # France
  n.occ.I[k] <- sum(tmpI[,k])       # Italy
  n.occ.CH[k] <- sum(tmpCH[,k])     # Switzerland
  n.occ.Au[k] <- sum(tmpAu[,k])     # Austria
  n.occ.Slo[k] <- sum(tmpSlo[,k])   # Slovenia
  n.occ.Ger[k] <- sum(tmpGer[,k])   # Germany
}
}
",fill = TRUE)
sink()


# beta.lpsi1.f ~ dnorm(0, 0.1) 		#effect of forest on psi1
# beta.lphi.f ~ dnorm(0, 0.1)			#effect of forest on phi
# beta.lgamma.elev ~ dnorm(0, 0.1)	#effect of elevation on gamma
# beta.lgamma.hdens ~ dnorm(0, 0.1)	#effect of hdens on gamma
# beta.lp.elev ~ dnorm(0, 0.1)		#effect of elev on detection
# beta.lp.f ~ dnorm(0, 0.1)			#effect of forest on detection


#### -------------- Initial values -----------
inits <- function(){ list(z = zst, 
alpha.psi1 = runif(1), beta.lpsi1.f = rnorm(1), beta.lpsi1.rd = rnorm(1),
alpha.phi = runif(1), lphi.year = rnorm(18), beta.lphi.f = rnorm(1), 
alpha.gamma = runif(1), lgamma.year = rnorm(18), beta.lgamma.elev = rnorm(1), beta.lgamma.hdens = rnorm(1),  beta.lgamma.rd = rnorm(1),
alpha.p = runif(6, 0, 1), lp.year = rnorm(19), beta.lp.elev = rnorm(1), beta.lp.f = rnorm(1), lp.season = c(NA, rnorm(2)),
lp.net = c(NA, rnorm(2)), 
sd.lphi = runif(1, 0.1, 1), sd.lgamma = runif(1, 0.1, 1), sd.lp = runif(1, 0.1, 1), sd.eps = runif(1, 0.1, 1),
betaDphi = rnorm(1), betaDgamma = rnorm(1))}



#### ----- Parameters monitored ---------    (including z)
params <- c("alpha.lpsi1", "beta.lpsi1.f", "beta.lpsi1.rd", 
"alpha.lphi", "lphi.year", "beta.lphi.f", 
"alpha.lgamma", "lgamma.year", "beta.lgamma.elev", "beta.lgamma.hdens","beta.lgamma.rd ",
"betaDphi", "betaDgamma", 
"alpha.lp", "lp.year", "beta.lp.elev", "beta.lp.f", "lp.season", "lp.net", 
"sd.lphi", "sd.lgamma", "sd.lp", "sd.eps", 
"n.occ.F", "n.occ.I", "n.occ.CH", "n.occ.Au", "n.occ.Slo", "n.occ.Ger", 
"fit", "fit.new", "z")

#2 chaines 20000 itérations, 5000 burn in
# na <- 5  ;  ni <- 6   ;   nt <- 1   ;   nb <- 3   ;   nc <- 1   # Super-short test

na <- 1000  ;  ni <- 20000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 2   


library(jagsUI)
#### ---------- Call JAGS from R, check convergence and summarize posteriors -----------
out9ju_nomisID <- jags(bdata, inits, params, "Dynocc9_nomisID.txt", parallel = TRUE,
         n.chains = nc, n.adapt = na, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(3,3), mar = c(5,3,2,2))  ;  traceplot(out9ju_nomisID)

save(out9ju_nomisID, file = "out9ju_nomisID.RData")
}





#----------------------------------------------------------------------------------
#MODEL WITH MISID
#----------------------------------------------------------------------------------

datju <- read.csv2('lynx_data_xjulie_160212.csv',h=T)
datju_o <- datju[order(datju$ID),]
yyju3 <- as.matrix(datju_o[,c(138:194)])
yyju3[which(yyju3== '999')] <- 0
yyju3[yyju3 > 1] <- 1  
yyju3[which(is.na(yyju3))] <- 0 
table(yyju3)
yyju3[which(yyju3 == 1)] <- 2

yyju123 <- yy+yyju3
yyju123[which(yyju123==3 )] <- 1

# Execute function
yju123 <- reformat.fn(D2.data = yyju123)
str(yju123)
 # int [1:2183, 1:3, 1:21] NA NA NA NA NA NA NA NA NA NA ...

# Sum check
sum(yyju123, na.rm = T)   ;   sum(yju123, na.rm = T)     # yields 4272 both times (2188+(2*1042))

save(yju123, file = 'yju123.RData')
yju123 <- yju123 +1 


#----------------------------------------------------------------------------------
#MODEL WITH MISID
#----------------------------------------------------------------------------------


{

#### ----- Bundle data ---------
country <- datju12$country
country[country ==0] <- 6       # Germany is 6 now
# Frankreich = 1, Italien = 2, Schweiz = 3, Österreich = 4, Slowenien = 5, Germany = 6
str(bdata <- list(y = yju123, nsite = dim(yju12)[1], nrep = dim(yju12)[2], nyear = dim(yju12)[3], 
        country = country, nw = nw, cov = cov, distrel = distrel, D = autocov))

zst <- obsZju         # Inits for z for this data set




#### --------- Specify model in BUGS language ----------
# ------------------------------------------------------
sink("Dynocc9_misID.txt")
cat("
model {

### Specify priors
alpha.lpsi1 <- logit(alpha.psi1)
alpha.lphi <- logit(alpha.phi)
alpha.lgamma <- logit(alpha.gamma)
alpha.lb <- logit(alpha.b)
alpha.psi1 ~ dbeta(1, 1)
alpha.phi ~ dbeta(1, 1)
alpha.gamma ~ dbeta(1, 1)
alpha.b ~ dbeta(1,1)

for(i in 1:6){                      # Country-specific intercepts for p11 only (fix)
  alpha.lp11[i] <- logit(alpha.p11[i])
  alpha.p11[i] ~ dbeta(1, 1)
  alpha.lp10[i] <- logit(alpha.p10[i])
  alpha.p10[i] ~ dbeta(1, 1)
}

beta.lpsi1.f ~ dnorm(0, 0.1) 		#effect of forest on psi1
beta.lphi.f ~ dnorm(0, 0.1)			#effect of forest on phi
beta.lgamma.elev ~ dnorm(0, 0.1)	#effect of elevation on gamma
beta.lgamma.hdens ~ dnorm(0, 0.1)	#effect of hdens on gamma
beta.lp11.elev ~ dnorm(0, 0.1)		#effect of elev on detection
beta.lp10.elev ~ dnorm(0, 0.1)		#effect of elev on detection
beta.lp11.f ~ dnorm(0, 0.1)			#effect of forest on detection
beta.lp10.f ~ dnorm(0, 0.1)			#effect of forest on detection


beta.lpsi1.rd ~ dnorm(0, 0.1)        # Regression coefs for distance to release site
beta.lgamma.rd ~ dnorm(0, 0.1)

betaDphi ~ dnorm(0, 0.1)             # Regression coefficients of autocovariate
betaDgamma ~ dnorm(0, 0.1)

beta.trendyearp10 ~ dnorm(0, 0.1)	# Regression coefficients of year on p10

lp11.season[1] <- 0                    # Effects of season on p11
for(i in 2:3){
  lp11.season[i] ~ dnorm(0, 0.1)
}

lp11.net[1] <- 0                       # Network effects in p11
lp11.net[2] ~ dnorm(0, 0.1)
lp11.net[3] ~ dnorm(0, 0.1)
 
lp10.season[1] <- 0                    # Effects of season on p10
for(i in 2:3){
  lp10.season[i] ~ dnorm(0, 0.1)
}

lp10.net[1] <- 0                       # Network effects in p10
lp10.net[2] ~ dnorm(0, 0.1)
lp10.net[3] ~ dnorm(0, 0.1)


# Random year-specific effects
for (k in 1:(nyear-1)){
  lphi.year[k] ~ dnorm(0, tau.lphi) # means are country-specific intercepts
  lgamma.year[k] ~ dnorm(0, tau.lgamma)
  lp11.year[k] ~ dnorm(0, tau.lp11)
  lp10.year[k] ~ dnorm(0, tau.lp10)

}
lp11.year[nyear] ~ dnorm(0, tau.lp11)
lp10.year[nyear] ~ dnorm(0, tau.lp10)

tau.lphi <- pow(sd.lphi, -2)
tau.lgamma <- pow(sd.lgamma, -2)
tau.lp11 <- pow(sd.lp11, -2)
tau.lp10 <- pow(sd.lp10, -2)
tau.eps11 <- pow(sd.eps11, -2)
tau.eps10 <- pow(sd.eps10, -2)
sd.lphi ~ dnorm(0, 0.2)I(0.01,)      # curve(dnorm(x, 0, sqrt(1/0.2)), 0, 10)
sd.lgamma ~ dnorm(0, 0.2)I(0.01,)    # Truncated effectively at 0
sd.lp11 ~ dnorm(0, 0.2)I(0.01,)
sd.lp10 ~ dnorm(0, 0.2)I(0.01,)
sd.eps11 ~ dnorm(0, 0.2)I(0.01,)
sd.eps10 ~ dnorm(0, 0.2)I(0.01,)


### Ecological submodel
for (i in 1:nsite){
  z[i,1] ~ dbern(psi1[i])
  psi1[i] <- 1 / (1 + exp(-lpsi1[i]))
#   lpsi1[i] <- min(200, max(-200, lpsi1.tmp[i]))
  lpsi1[i] <- alpha.lpsi1 + beta.lpsi1.f * cov[i,3] + beta.lpsi1.rd * distrel[i,1]
   
  for (k in 2:nyear){
    z[i,k] ~ dbern(muZ[i,k])
    muZ[i,k]<- z[i,k-1] * phi[i, k-1] + (1-z[i,k-1]) * gamma[i, k-1]
    phi[i, k-1] <- 1 / (1 + exp(-lphi[i, k-1]))
    lphi[i, k-1] <- alpha.lphi + lphi.year[k-1] + 
	beta.lphi.f * cov[i,3] + betaDphi * D[i,k-1]
    gamma[i, k-1] <- 1 / (1 + exp(-lgamma[i, k-1]))
    lgamma[i, k-1] <- alpha.lgamma + lgamma.year[k-1] + 
	beta.lgamma.elev * cov[i,1] + beta.lgamma.hdens * cov[i,4] + beta.lgamma.rd * distrel[i,k-1] + betaDgamma * D[i,k-1]
  } #k
} #i

#DELTA IS CONSTANT
b <-  1 / (1 + exp(-lb))
lb <- alpha.lb

### Observation model (also generate new data set for use in PPD checking of GOF)
for (i in 1:nsite){
  for (k in 1:nyear){
    for (j in 1:nrep){
	  
      p11[i,j,k] <- 1 / (1 + exp(-lp11[i,j,k]))
#     lp11[i,j,k] <- min(200, max(-200, lp11.tmp[i,j,k]))
      lp11[i,j,k] <- alpha.lp11[country[i]] + lp11.year[k] + lp11.net[nw[i,k]] + lp11.season[j] + 
         	  beta.lp11.elev * cov[i,1] + beta.lp11.f * cov[i,3] + eps11[i,k]
			  

	  p10[i,j,k] <- 1 / (1 + exp(-lp10[i,j,k]))
#     lp10[i,j,k] <- min(200, max(-200, lp10.tmp[i,j,k]))
      lp10[i,j,k] <- alpha.lp10[country[i]] + lp10.year[k] + lp10.net[nw[i,k]] + lp10.season[j] + 
         	  beta.lp10.elev * cov[i,1] + beta.lp10.f * cov[i,3] + beta.trendyearp10 * (k-1) + eps10[i,k]
			  

	muy[i,j,k,1]<- (1-z[i,k])*(1-p10[i,j,k])+z[i,k]*(1-p11[i,j,k])  ### P(O=0)
    muy[i,j,k,2]<- z[i,k]*p11[i,j,k]*b                           	### P(O=1)
    muy[i,j,k,3]<- z[i,k]*(1-b)*p11[i,j,k]+(1-z[i,k])*p10[i,j,k] 	### P(O=2)
        
	
	y[i,j,k] ~ dcat(muy[i,j,k,])
	y.new[i,j,k] ~ dcat(muy[i,j,k,])	# Generate replicate data set under same model
    } #j
    eps11[i,k] ~ dnorm(0, tau.eps11)            # random site.year effect in detection
	eps10[i,k] ~ dnorm(0, tau.eps10)
  } #k
} #i

### Posterior predictive check of GOF based on detection frequency per year
for (i in 1:nsite){
  for (k in 1:nyear){
    sum.y[i,k] <- sum(y[i,,k])	                                      # Observed actual data
    eval[i,k] <- max(0.01, sum(muy[i,,k,]))                            # Expected data (sum of muy)
    E[i,k] <- pow((sum.y[i,k] - eval[i,k]),2) / (eval[i,k] + 0.01)    # Chi2 discrepancy measure actual data set
    sum.y.new[i,k] <- max(0.01, sum(y.new[i,,k]))                     # Observed replicate data
    E.new[i,k] <- pow((sum.y.new[i,k] - eval[i,k]),2) / (eval[i,k] + 0.01)# Chi2 discrepancy measure 'perfect' data set
    # detfreq.ratio[i,k] <- sum.y[i,k] / sum.y.new[i,k]               # Ratio of detection frequencies observed, expected ----- DONT NEED NOW ---------
    # df.diff[i,k] <- sum.y[i,k] - sum.y.new[i,k]                     # Difference ..... (same)  ---- DONT REALLY NEED THIS ----------
    } #k
} #i
fit <- sum(E[,])                    # Sum up to overall discrepancy measure
fit.new <- sum(E.new[,])

### Number of occupied quads per country and year
for (k in 1:nyear){
  for (i in 1:nsite){
    tmpF[i,k] <- equals(country[i], 1) * z[i,k]      
    tmpI[i,k] <- equals(country[i], 2) * z[i,k]
    tmpCH[i,k] <- equals(country[i], 3) * z[i,k]
    tmpAu[i,k] <- equals(country[i], 4) * z[i,k]
    tmpSlo[i,k] <- equals(country[i], 5) * z[i,k]
    tmpGer[i,k] <- equals(country[i], 6) * z[i,k]
  }
  n.occ.F[k] <- sum(tmpF[,k])       # France
  n.occ.I[k] <- sum(tmpI[,k])       # Italy
  n.occ.CH[k] <- sum(tmpCH[,k])     # Switzerland
  n.occ.Au[k] <- sum(tmpAu[,k])     # Austria
  n.occ.Slo[k] <- sum(tmpSlo[,k])   # Slovenia
  n.occ.Ger[k] <- sum(tmpGer[,k])   # Germany
}
}
",fill = TRUE)
sink()




#### -------------- Initial values -----------
inits <- function(){ list(z = zst, 
alpha.psi1 = runif(1), beta.lpsi1.f = rnorm(1), beta.lpsi1.rd = rnorm(1),
alpha.phi = runif(1), lphi.year = rnorm(18), beta.lphi.f = rnorm(1), 
alpha.gamma = runif(1), lgamma.year = rnorm(18), beta.lgamma.elev = rnorm(1), beta.lgamma.hdens = rnorm(1),  beta.lgamma.rd = rnorm(1),
alpha.p11 = runif(6, 0, 1), 
lp11.year = rnorm(19), beta.lp11.elev = rnorm(1), beta.lp11.f = rnorm(1), lp11.season = c(NA, rnorm(2)),
lp11.net = c(NA, rnorm(2)), 
alpha.p10 = runif(6, 0, 1), 
lp10.year = rnorm(19), beta.lp10.elev = rnorm(1), beta.lp10.f = rnorm(1), lp10.season = c(NA, rnorm(2)),
lp10.net = c(NA, rnorm(2)), beta.trendyearp10 = rnorm(1),
alpha.b = runif(1),
sd.lphi = runif(1, 0.1, 1), sd.lgamma = runif(1, 0.1, 1), sd.lp11 = runif(1, 0.1, 1), sd.eps11 = runif(1, 0.1, 1),
betaDphi = rnorm(1), betaDgamma = rnorm(1))}



#### ----- Parameters monitored ---------    (including z)
params <- c("alpha.lpsi1", "beta.lpsi1.f", "beta.lpsi1.rd", 
"alpha.lphi", "lphi.year", "beta.lphi.f", 
"alpha.lgamma", "lgamma.year", "beta.lgamma.elev", "beta.lgamma.hdens","beta.lgamma.rd ",
"betaDphi", "betaDgamma", 
"alpha.lp11", "lp11.year", "beta.lp11.elev", "beta.lp11.f", "lp11.season", "lp11.net", 
"alpha.lp10", "lp10.year", "beta.lp10.elev", "beta.lp10.f", "lp10.season", "lp10.net", "beta.trendyearp10",
"alpha.lb",
"sd.lphi", "sd.lgamma", "sd.lp11", "sd.eps11", "sd.eps10",
"n.occ.F", "n.occ.I", "n.occ.CH", "n.occ.Au", "n.occ.Slo", "n.occ.Ger", 
"fit", "fit.new", "z")



#2 chaines 20000 itérations, 5000 burn in

na <- 1000  ;  ni <- 20000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 2   

# na <- 5  ;  ni <- 3   ;   nt <- 1   ;   nb <- 1   ;   nc <- 1   # Super-short test


library(jagsUI)
#### ---------- Call JAGS from R, check convergence and summarize posteriors -----------
out9ju_misID <- jags(bdata, inits, params, "Dynocc9_misID.txt", parallel = TRUE,
         n.chains = nc, n.adapt = na, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(3,3), mar = c(5,3,2,2))  ;  traceplot(out9ju_misID)

save(out9ju_misID, file = "out9ju_misID.RData")
}
