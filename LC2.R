## MORTALITY RATES DATA PREPARATION ##
ageLimit <- 90																# Constants
ageGroups <- ageLimit + 1
timePeriods <- 98
forcHorizon <- 52
mrRaw <- read.csv(file.choose(), header=FALSE) 								# Import dataset
colnames(mrRaw) <- c("Year", "Age", "Females", "Males", "Combined")		 	# Name columns
mrReducedAge <- subset(mrRaw, Age >= 0 & Age <= ageLimit)					# Reduce to age interval
mrVec <- t(mrReducedAge[-c(1:4)])											# Remove all columns except "Combined"
mrMatrix <- matrix(mrVec, nrow = ageGroups, ncol = timePeriods)				# Fold into matrix
colnames(mrMatrix) <- 1921:2018												# Assign appropriate column names
rownames(mrMatrix) <- 0:ageLimit 											# Assign appropriate row names
mode(mrMatrix) = "numeric"													# Make matrix numeric

## 2ND ORDER LEE-CARTER METHOD ##
## Fitting procedure ##
a_x = rowMeans(log(mrMatrix))												# Create vector of means
M = log(mrMatrix) - a_x														# Define matrix M (centralize)
svdM <- svd(M)																# Perform SVD on M (sing.vec. = 1)
b_x1 <- -svdM$u[,1]															# Define b_x parameters
b_x2 <- -svdM$u[,2]																
k_t1 <- -svdM$d[1]*svdM$v[,1]												# Define k_t parameters
k_t2 <- -svdM$d[2]*svdM$v[,2]											
Mhat2 <- matrix(nrow = ageGroups, ncol = timePeriods)						# Initiate M-hat matrix
for (i in 1:timePeriods){													# Fit data
	Mhat2[,i] <- a_x + b_x1*k_t1[i] + b_x2*k_t2[i]
}
colnames(Mhat2) <- 1921:2018												# Appropriate naming
rownames(Mhat2) <- 0:ageLimit
## Forecasting procedure ##
dhat1 <- (k_t1[timePeriods] - k_t1[1])/(timePeriods - 1)					# Define drift parameters
dhat2 <- (k_t2[timePeriods] - k_t2[1])/(timePeriods - 1)	
kForc1_t <- matrix(nrow = forcHorizon, ncol = 1)							# Initiate forcast k_t vector
kForc2_t <- matrix(nrow = forcHorizon, ncol = 1)
for (i in 1:forcHorizon) {													# Calculate forecasted k_t
	kForc1_t[i] = k_t1[timePeriods] + dhat1*i								# Random walk
	kForc2_t[i] = k_t2[timePeriods] + dhat2*i
}
Mtilde2 <- matrix(nrow = ageGroups, ncol = forcHorizon)						# Initiate M-tilde matrix
for (i in 1:forcHorizon) {													# Fit data
	Mtilde2[,i] <- a_x + b_x1*kForc1_t[i] + b_x2*kForc2_t[i]
}
colnames(Mtilde2) <- 2019:2070												# Appropriate naming
rownames(Mtilde2) <- 0:ageLimit	
Mstar2 <- cbind(Mhat2,Mtilde2)												# Concatenate to complete model
## Proportion variance explained procedure ##
num2 <- matrix(nrow = timePeriods, ncol = 1)								# Initiate temp vectors
den2 <- matrix(nrow = timePeriods, ncol = 1)	
etaSqrd2 <- matrix(nrow = ageGroups, ncol = 1)								# Initiate eta^2 vector
for (i in 1:ageGroups)	{													# Calculate prop. variance expl.
	for (j in 1:timePeriods) {												
		num2[j] <- (mrMatrix[i,j] - exp(Mhat2[i,j]))^2
		den2[j] <- (mrMatrix[i,j] - exp(a_x[i]))^2
	}
etaSqrd2[i] = 1 - (sum(num2)/sum(den2))
}
## Confidence intervals procedure ##
total1 = matrix(nrow = timePeriods - 1, ncol = 1)							# Initiate temp vectors
total2 = matrix(nrow = timePeriods - 1, ncol = 1)	
for (i in 1:timePeriods-1) {												# Calculate see
	total1[i] <- (k_t1[i+1] - k_t1[i] - dhat1)^2
	total2[i] <- (k_t2[i+1] - k_t2[i] - dhat2)^2
}
see1 = sqrt(1/(timePeriods-2)*sum(total1))									# Define see
see2 = sqrt(1/(timePeriods-2)*sum(total2))
upperk_k1 <- matrix(nrow = ageGroups, ncol = 1)								# Initiate interval vectors
lowerk_k1 <- matrix(nrow = ageGroups, ncol = 1)
upperk_k2 <- matrix(nrow = ageGroups, ncol = 1)
lowerk_k2 <- matrix(nrow = ageGroups, ncol = 1)
lower2 <- matrix(nrow = ageGroups, ncol = forcHorizon)						# Initiate confidence matrix
upper2 <- matrix(nrow = ageGroups, ncol = forcHorizon)

Mtilde2 <- exp(Mtilde2)

for (i in 1:ageGroups) {	
	for (j in 1: forcHorizon) {													# Create confidence intervals
	upper2[i,j] <- Mtilde2[i,j]*exp(-1.96*(b_x1[i]*see1 + b_x2[i]*see2))
	lower2[i,j] <- Mtilde2[i,j]*exp(1.96*(b_x1[i]*see1 + b_x2[i]*see2))
	}
}