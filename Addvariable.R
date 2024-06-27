rm(list= ls())

# ------provide directory of the RData--------
setwd("/Users/grace/Desktop/Projectone/RealData/DATA1")

source('PMSE.R')
source('Funs.R')

library(fda)
library(dbplyr)
library(MASS)
library(splines)
library(rgl)
library(plotly)
library(psych)

windspeed = read.csv("WindSpeed.csv",header = T,row.names = 1)

CAweather = CanadianWeather
place     <- CAweather$place

Precip = CAweather$monthlyPrecip
arvida_indexx <- which(dimnames(Precip)[[2]] == "Arvida")
Precip = Precip[, -arvida_indexx]

Temp  <- CAweather$monthlyTemp
arvida_index <- which(dimnames(Temp)[[2]] == "Arvida")
Temp <- Temp[, -arvida_index]

Month = seq(1, 12, length = 12)

par(mar = c(5,5,3,1))
plot(Month, Temp[,1], type = 'l',
     xlab= "Month" , ylab="Temperature(deg C)",
     ylim = c(-31,26), cex.lab = 2, cex.axis = 1.5,
     cex.main = 2, main = "(b)")
for (i in 2: ncol(Temp)) {
  lines(Month, Temp[,i] , col = i) 
}

par(mar = c(5,5,3,1))
plot(Month, Precip[,1], type = 'l',
     xlab= "Month" , ylab="Precipitation(mm)",
     ylim = c(0,13),cex.lab = 2, cex.axis = 1.5 , 
     cex.main = 2, main = "(c)")
for (i in 2: ncol(Precip)) {
  lines(Month, Precip[,i] , col = i) 
}

par(mar = c(4.5,5,3,1))
plot(Month, windspeed[,1], type = 'l',
     xlab= "Month" , ylab="Max Wind Speed(km/h)",
     ylim = c(10 , 130), cex.lab = 2, cex.axis = 1.5, 
     cex.main = 2, main = "(a)")
for (i in 2: ncol(windspeed)) {
  lines(Month, windspeed[,i] , col = i) 
}

#--smooth three variables-----------
# -----Precp-----

monthbasis12 <- create.bspline.basis(rangeval=c(0, 12), nbasis= 12) #Creating a B-spline basis

lambda_values <- seq(1e-5, 10, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  Precip_smooth <- matrix(0, nrow = nrow(Precip), ncol = ncol(Precip))
  
  for(i in 1: 34) {
    curveData <- Precip[,i]
    fdParObj <-  fdPar(monthbasis12, Lfdobj=2, lambda = lambda_values[j]) # Adjust lambda for smoothing parameter
    smoothObj <- smooth.basis(Month, curveData, fdParObj)
    Precip_smooth[,i] <- eval.fd(Month, smoothObj$fd)
  }
  
  # Calculate the error for this lambda value
  errors[j] <- mean((Precip_smooth - as.matrix(Precip))^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))

dense  <- seq(0, 12, length.out = 200)
Precip_smooth   <- matrix(0, nrow = nrow(Precip), ncol = ncol(Precip))
Precip_smoothed <- matrix(0, nrow =  200, ncol = 34)

for(i in 1: 34) {
  curveData <- Precip[,i]
  fdParObj <-  fdPar(monthbasis12, Lfdobj=2, lambda = best_lambda) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(Month, curveData, fdParObj)
  Precip_smooth[,i] <- eval.fd(Month, smoothObj$fd)
  Precip_smoothed[,i] <- eval.fd(dense, smoothObj$fd)
}

par(mar = c(5,5,3,1))
plot(dense, Precip_smoothed[,1], type = 'l',
     xlab= "Month" , ylab="Precipitation(mm)",
     ylim = c(0,13),cex.lab = 2, cex.axis = 1.5 , 
     cex.main = 2, main = "(c)")
for (i in 2: ncol(Precip)) {
  lines(dense, Precip_smoothed[,i] , col = i) 
}

# -----Temp-----
lambda_values <- seq(1e-2, 10, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  Temp_smooth <- matrix(0, nrow = nrow(Temp), ncol = ncol(Temp))
  
   for(i in 1: 34) {
     curveData <- Temp[,i]
     fdParObj <-  fdPar(monthbasis12, Lfdobj=2, lambda = lambda_values[j]) # Adjust lambda for smoothing parameter
     smoothObj <- smooth.basis(Month, curveData, fdParObj)
     Temp_smooth[,i] <- eval.fd(Month, smoothObj$fd)
   }
   
  # Calculate the error for this lambda value
  errors[j] <- mean((Temp_smooth - as.matrix( Temp))^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))

dense  <- seq(0, 12, length.out = 200)
Temp_smooth   <- matrix(0, nrow = nrow(Temp), ncol = ncol(Temp))
Temp_smoothed <- matrix(0, nrow =  200, ncol = 34)

for(i in 1: 34) {
  curveData <- Temp[,i]
  fdParObj <-  fdPar(monthbasis12, Lfdobj=2, lambda = best_lambda) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(Month, curveData, fdParObj)
  Temp_smooth[,i] <- eval.fd(Month, smoothObj$fd)
  Temp_smoothed[,i] <- eval.fd(dense, smoothObj$fd)
}

par(mar = c(5,5,3,1))
plot(dense, Temp_smoothed[,1], type = 'l',
     xlab= "Month" , ylab="Temperature(deg C)",
     ylim = c(-31,26), cex.lab = 2, cex.axis = 1.5,
     cex.main = 2, main = "(b)")
for (i in 2: 34) {
  lines(dense, Temp_smoothed[,i] , col = i) 
}

# -----wind speed-----

lambda_values <- seq(1e-2, 10, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  windspeed_smooth <- matrix(0, nrow = nrow(windspeed), ncol = ncol(windspeed))
  
  for(i in 1: 34) {
    curveData <- windspeed[,i]
    fdParObj <-  fdPar(monthbasis12, Lfdobj=2, lambda = lambda_values[j]) # Adjust lambda for smoothing parameter
    smoothObj <- smooth.basis(Month, curveData, fdParObj)
    windspeed_smooth[,i] <- eval.fd(Month, smoothObj$fd)
  }
  
  # Calculate the error for this lambda value
  errors[j] <- mean((windspeed_smooth - as.matrix( windspeed))^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))

windspeed_smooth <- matrix(0, nrow = nrow(windspeed), ncol = ncol(windspeed))
dense  <- seq(0, 12, length.out = 200)
windspeed_smoothed <- matrix(0, nrow = 200, ncol = 34)

for(i in 1: 34) {
  curveData <- windspeed[,i]
  fdParObj <-  fdPar(monthbasis12, Lfdobj=2, lambda = best_lambda) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(Month, curveData, fdParObj)
  windspeed_smooth[,i] <- eval.fd(Month, smoothObj$fd)
  windspeed_smoothed[,i] <- eval.fd(dense, smoothObj$fd)
}

par(mar = c(4.5,5,3,1))
plot(dense, windspeed_smoothed[,1], type = 'l',
     xlab= "Month" , ylab="Max Wind Speed(km/h)",
     cex.lab = 2, cex.axis = 1.5, ylim = c(10,130),
     cex.main = 2, main = "(a)")
for (i in 2: ncol(windspeed_smooth)) {
  lines(dense, windspeed_smoothed[,i] , col = i) 
}

# write.csv(windspeed_smoothed,"wind.csv")
# write.csv(Temp_smoothed,"Temp.csv")
# write.csv(Precip_smoothed,"Precip.csv")
Precip_smoothed = read.csv('Precip.csv',header = TRUE)
windspeed_smoothed = read.csv('wind.csv',header = TRUE)
Temp_smoothed = read.csv('Temp.csv',header = TRUE)

# ------Define the number of individuals and time points T------
N      <- 34
ps     <- 12
Times  <- 12
Month = seq(1, 12, length = 12)
tobs   <- Month
time_rangeval  <- c(0, Times)
s      <- Month
dense  <- seq(0, 12, length.out = 200)
#--------------Y(t)-- X(t) ---A(s)---------

Y_mat = Precip_smoothed[,-1]
X_mat = windspeed_smoothed[,-1]
A_mat = Temp_smoothed[,-1]

xmean <- apply(X_mat ,2,mean)
ymean <- apply(Y_mat ,2,mean)
Amean <- apply(A_mat ,2,mean)

# centering data
# xmean <- apply(X_mat ,1,mean)
# ymean <- apply(Y_mat ,1,mean)
# Amean <- apply(A_mat ,1,mean)

Xmat <- X_mat - xmean
Ymat <- Y_mat - ymean
#Amat <- A_mat - Amean
Amat = A_mat
set.seed(123) # Setting a seed for reproducibility
random_columns <- sample(1:34, 34)

# Select the first 28 columns for Xmat and Ymat based on the shuffled indices
Xmat_train <- t(Xmat[, random_columns[1:28]])
Ymat_train <- t(Ymat[, random_columns[1:28]])
Amat_train <- t(Amat[, random_columns[1:28]])

# Select the remaining 6 columns for Xmat and Ymat based on the shuffled indices
Xmat_test <- t(Xmat[, random_columns[29:34]])
Ymat_test <- t(Ymat[, random_columns[29:34]])
Amat_test <- t(Amat[, random_columns[29:34]])

time_rangeval <- c(0,12)
tobs <- dense

# Define t basis function
t_basis <- basisphi.gen(d = 3, M = 4, domain = time_rangeval)
L          <- t_basis$nbasis # d+M

# evaluate basis matrix of X(t)
t_basisMat  <- compute.phi( d = 3, M = 4, domain = time_rangeval , t= tobs)

# Define temperate s basis function
As_basis  <- psi.generator (d = 3, M = 4, domain = c(min(Amat), max(Amat)))

# evaluate basis matrix of psi at A(s)

As_psi_evalist <- psi.eval.fun(n = 28, psi = As_basis, As = Amat_train)
nbasis_As   <- As_basis$nbasis

# penalty matirx
theta = theta.basis.fun(d=3,M=4,time_rangeval= c(0,12))
R2     <- eval.penalty(theta, Lfdobj =int2Lfd(2), tobs)


# basis function of theta(t)
basismat <- theta.spline.generator(d = 3, M = 4 ,time_rangeval= c(0,12), tobs)

# solve b0 begin value in iterative
#b_0       <- b.fit(Xmat_train, Ymat_train, basismat, tobs)
b_0       <- b.fit(Xmat_train, Ymat_train, basismat, tobs, k= 0.1, R = R2)
# solve Beta
Betafit_0 <- Beta.function(b = b_0 ,basismat, tobs)

Betafit_0 <- matrix(diag(Betafit_0 ), nrow =length(tobs) , ncol =length(tobs))

# Y*(t) at 1st iterative
Ystar_1   <- Ymat_train - Xmat_train %*% t(Betafit_0)

W_1    <- compute.u(Ystar_1 , d = 3, M = 4,  domain = time_rangeval)

zeta   <- computeZeta(n = 28, nbasis_As, As_psi_evalist, s=dense)

# ----gamma^(1) at 1 iterative-------------
y <- seq(min(Amat_train),max(Amat_train), length.out = 12)

R1  <- eval.penalty(As_basis, Lfdobj = int2Lfd(2),  y) # L*L
R3  <- eval.penalty(t_basis,  Lfdobj = int2Lfd(2), tobs) # K*K
Psiy  = eval.basis(y, As_basis)

gamma_1 = gamma.fun(zeta, W_1 , t_basis,Month, nbasis_As, As_basis, y,0.1,0, R1,R3)



# h^(1) at 1st iterative

hpredorder1 <- t_basisMat  %*% gamma_1 %*% t(Psiy)

plott <- plot_ly(x = ~ tobs, y = ~ y, z = ~ hpredorder1, type = "surface")
layout(plott, scene =
         list(xaxis = list(title = "t-axis"),
              yaxis = list(title = "s-axis"),
              zaxis = list(title = "h"),
              title = "H Surface Plot"))

integralH_1  <- zeta %*% t(gamma_1) %*% t(t_basisMat)

# Y^~(t)
Ytilde_2 = Ymat_train - integralH_1

# b^(1) at 2nd iterative
b_1  = b.fit(Xmat_train , Ytilde_2, basismat, tobs, k= 0.1, R = R2)

# beta^(1) at 2nd iterative
Betafit_1 = Beta.function(b = b_1 ,basismat, tobs)

plot(Betafit_1) 

Betafit_1  <- matrix(diag(Betafit_1), nrow =length(tobs) , ncol =length(tobs))

# Y^*(2) at 2nd iterative
Ystar_2  = Ymat_train - Xmat_train %*% t(Betafit_1)

W_2  = compute.u(Ystar_2 , d = 3, M = 4,  domain = time_rangeval)

# gamma^(2) at 2st iterative
gamma_2 = gamma.fun(zeta, W_2 , t_basis, tobs = Month, nbasis_As)

# h^(2) at 2nd iterative
hpredorder2 <- t_basisMat  %*% gamma_2 %*% t(Psiy)

# ---Set a convergence criterion and maximum number of iterations
convergence_criterion = 1e-3
max_iterations = 300 #136 Converged in 136 iterations

iteration = 1
converged = FALSE

while (!converged && iteration <= max_iterations) {
  # Compute b_x
  b_x = b.fit(Xmat_train , Ymat_train, basismat, tobs)
  if (iteration > 1) {
    b_x = b.fit(Xmat_train, Ytilde_x, basismat, tobs)
  }
  
  # Compute Betafit
  Betafit = Beta.function(b = b_x, basismat, tobs)
  Betafit_x = matrix(diag(Betafit), nrow = length(tobs), ncol = length(tobs))
  
  # Compute Ystar_x
  Ystar_x = Ymat_train - Xmat_train %*% t(Betafit_x)
  
  # Compute W_x
  W_x = compute.u(Ystar_x, d = 3, M = 4, domain = time_rangeval)
  
  # Compute gamma_x
  gamma_x = gamma.fun(zeta, W_x, t_basis, tobs, nbasis_As)
  
  # Compute horder_x
  horder_x = t_basisMat %*% gamma_x %*% t(Psiy)
  
  # Compute integralH_pred and update Ytilde_x
  
  integralH_pred_x = zeta %*% t(gamma_x) %*% t(t_basisMat)
  
  Ytilde_x = Ymat_train - integralH_pred_x
  
  # Check for convergence (based on the change in Ytilde_x or Ystar_x)
  if (iteration > 1 && max(abs(Ystar_x - prev_Ystar_x)) < convergence_criterion && 
      max(abs(Ytilde_x - prev_Ytilde_x)) < convergence_criterion) {
    converged = TRUE
  }
  
  # Store current Ystar_x and Ytilde_x for next iteration's comparison
  prev_Ystar_x = Ystar_x
  prev_Ytilde_x = Ytilde_x
  
  # Increment iteration counter
  iteration = iteration + 1
}

# Output the results or status
if (converged) {
  cat("Converged in", iteration, "iterations\n")
} else {
  cat("Maximum iterations reached without convergence\n")
}

# Final results are in b_x, W_x, gamma_x, Ystar_x, Ytilde_x, horder_x, Betafit_x
gamma_hat = gamma_x

plot(tobs, diag(Betafit_x),type = 'l',
     xlab = "Month",
     ylab = expression(beta))

h_pred <- t_basisMat  %*% gamma_hat %*% t(Psiy)
open3d()
persp3d(tobs , y , h_pred,
        theta = 30, phi = 20,
        xlab = 'Time', ylab = 'Temperature', zlab = 'h',
        axes = TRUE, col = 'lightblue',
        border = NA,  facets = FALSE, zlim = c(-2,2))

integralH  <- zeta %*% t(gamma_hat) %*% t(t_basisMat)
XtBeta_train <- Xmat_train %*% Betafit_x

XtBeta <- Xmat_test %*% Betafit_x 

temp_test <- seq(min(Amat_test),max(Amat_test), length.out = 12)
Psi_temp_test = eval.basis(temp_test, As_basis)
As_psi_test <- psi.eval.fun(n = 6, psi = As_basis, As = Amat_test)

zeta_test   <- computeZeta(n = 6, nbasis_As, As_psi_test, s)

integralH_train  <- zeta %*% t(gamma_hat) %*% t(t_basisMat)
integralH_test  <- zeta_test %*% t(gamma_hat) %*% t(t_basisMat)

Yhat = XtBeta + integralH_test 

PMSE = PMSE(6, Ymat_test, Yhat, tobs)

#contour(tobs , y , z = h_pred, xlab = "Time", ylab = "temperature")

library(plotly)
plottt = plot_ly(x = ~ Month, y = ~ y, z = ~ h_pred , type = "surface", 
            colorscale = 'grey')
layout(plottt, scene = list(xaxis = list(title = "Time"),
                            yaxis = list(title = "Temperature"),
                            zaxis = list(title = "h", range = c(-3, 2))))


# ---- smooth the beta-----------

bsplineBasis <- create.bspline.basis(c(0,12), 7, 4)

lambda_values <- seq(1e-4, 1, length.out = 10)

# Initialize a vector to store errors for each lambda
errors <- numeric(length(lambda_values))

for (j in 1:length(lambda_values)) {
  smoothedCurves <- matrix(0, nrow = 1, ncol = 1)
  
    curveData <- Betafit
    fdParObj <- fdPar(bsplineBasis, Lfdobj=2, lambda=lambda_values[j])
    smoothObj <- smooth.basis(tobs, curveData, fdParObj)
    smoothedCurves<- eval.fd(tobs, smoothObj$fd)
    errors <- mean((smoothedCurves - Betafit)^2)
}

# Find the lambda with the smallest error
best_lambda <- lambda_values[which.min(errors)]
cat("The best lambda value is:", best_lambda, "with an MSE of:", min(errors))

dense_points <- seq(1, 12, length.out = 50)

smoothed_beta <- matrix(0, nrow = 1, ncol = 1)
dense_smoothed_beta <- matrix(0, nrow = 1, ncol = 100)

  curveData <- Betafit
  fdParObj <-  fdPar(bsplineBasis, Lfdobj=2, lambda=1e-4) # Adjust lambda for smoothing parameter
  smoothObj <- smooth.basis(tobs, curveData, fdParObj)
  smoothed_beta <- eval.fd(tobs, smoothObj$fd)
  dense_smoothed_beta<- eval.fd(dense_points, smoothObj$fd)

plot(dense_points, dense_smoothed_beta,
       type = 'l',
       ylab=expression(beta), 
       xlab="Month")
  
plot(tobs, smoothed_beta,type = 'l',
     ylab=expression(beta), 
     xlab="Month")

# ---------------Choosing Smoothing Parameter lamada-----------------------
choosing.lambda.gcv <-function(theta.basis, Betafit, loglam , t){
  
  nlam    <- length(loglam)
  gcvsave <- rep(0, nlam)
  dfsave  <- rep(0, nlam)
  
  for (i in  1 : nlam ) {
    cat(paste("lambda",loglam[i],"\n"))
    
    lambdai <- 10^loglam[i]
    thetafdPari <- fdPar(theta.basis, lambda = lambdai)
    betafdi  <- smooth.basis(t,as.vector(Betafit),  thetafdPari)$fd
    dfi      <- smooth.basis(t, as.vector(Betafit),  thetafdPari)$df
    gcvi     <- smooth.basis(t,as.vector(Betafit),  thetafdPari)$gcv
    gcvsave[i] <- sum(gcvi)
    
    dfsave[i] <- dfi
  }
  min_gcv <- min(gcvsave)
  # Find the index of the minimum GCV value
  min_gcv_index <- which.min(gcvsave)
  # Find the lambda value corresponding to the minimum GCV
  min_lambda <- 10^loglam[min_gcv_index]
  return(list(min_lambda,min_gcv))
}

theta.basis <- theta.basis.fun(d = 3, M = 7 ,time_rangeval= c(0,12))

gcvlambda  <- choosing.lambda.gcv (theta.basis , Betafit , loglam = seq(-1, 2, by = 0.05), tobs)
lambda.min <- gcvlambda[[1]]
print(lambda.min)
#--------------------with roughness penalty parameter---------------


R2   <- fdPar( theta.basis, 2, lambda.min)
#R2   <- fdPar( theta.basis, 2 , 10)
Beta.smooth.gcv <- smooth.basis(tobs ,as.vector(Betafit) , R2)

# Rmat <- Beta.smooth.gcv$penmat

Beta.GCV.fd  <- Beta.smooth.gcv$fd
Beta.gcv     <- predict(Beta.smooth.gcv,tobs)
plot(tobs, Beta.gcv, type = 'l',
     cex.lab = 1.5, cex.axis=1.5,
     xlab = "Month", ylab = expression(beta))




#---add smooth parameter-------
## Computer Jphiphi
Jphiphi  <- eval.penalty(t_basis, Lfdobj =int2Lfd(0), tobs)
Jpsipsi  <- eval.penalty(As_basis, Lfdobj = int2Lfd(0), y)

# Create the penalty matrices
R1  <- eval.penalty(As_basis,Lfdobj = int2Lfd(2),  y) # L*L
R2  <- eval.penalty(t_basis, Lfdobj =int2Lfd(2), tobs) # K*K

PenMat_y <- kronecker(Jphiphi,R1)
PenMat_t <- kronecker(R2,Jpsipsi)

Ystar = Ymat_train - XtBeta_train
Y.Phi.Mat = compute.u(Ystar , d = 3, M = 7,  domain = time_rangeval)


lambda_t_values <- exp(seq(-10, 0,length.out = 11)) #-20 1 30
lambda_t_values
lambda_y_values <- exp(seq(-10, 0,length.out = 11)) 
fold_indices <- split(1:28, cut(1:28, breaks = 7, labels = FALSE))

# Initialize values to store results
min_cv      <- Inf
min_lambda_y <- NULL
min_lambda_t <- NULL
min_gamma    <- NULL

# Loop through lambda_t_values and lambda_y_values
for (lambda_t_idx in 1:length(lambda_t_values)) {
  for (lambda_y_idx in 1:length(lambda_y_values)) {
    
    lambda_t <- lambda_t_values[lambda_t_idx]
    lambda_y <- lambda_y_values[lambda_y_idx]
    
    cv_scores <- numeric(length(lambda_y_values))
    
    for (fold in 1: 7) {
      test_indices  = fold_indices[[fold]]
      train_indices = setdiff(1:28, test_indices)
      
      Ystar_train    = Ystar[train_indices,] # N_train x T 
      Ystar_test     = Ystar[test_indices,]  # N_test x T 
      W_trian  = Y.Phi.Mat[train_indices,]
      W_test   = Y.Phi.Mat[test_indices,]
      zeta_train = zeta[train_indices,]
      zeta_test  = zeta[test_indices,]
      
      ZetaW_Train  <- t(zeta_train) %*% W_trian
      
      vecZetaW_Train <- as.vector(ZetaW_Train)
      # Compute gamma
      inprod_zeta_train  <- t(zeta_train)%*% zeta_train
      JphiZeta_train     <- kronecker(Jphiphi, inprod_zeta_train)
      
      gamma_mat <- ginv(JphiZeta_train + lambda_y * PenMat_y + lambda_t * PenMat_t) %*%vecZetaW_Train
      
      gamma     <- matrix(gamma_mat, nrow = nrow(Jphiphi) , ncol= ncol(Jpsipsi))  # Reshape gamma back to a matrix
      
      YPred <- zeta_test %*% t(gamma) %*% t(t_basisMat)
      
      MSE <- int.Yerror(4, Ystar_test, YPred,s)
      cv_scores[fold] <- MSE
      
    }
    avg_cv <- mean(cv_scores)
    print(paste("t: ", lambda_t, "y: ", lambda_y," Avg CV: ", avg_cv))
    # Update best results if necessary
    if (avg_cv < min_cv) {
      min_cv       <- avg_cv
      min_lambda_y <- lambda_y
      min_lambda_t <- lambda_t
      min_gamma    <- gamma
    }
  }
}

# --------------min GCV, min MSE  lambda_y, lambda_t and gamma-------

print(paste("Min GCV: ", min_cv, " Lambda_t: ", min_lambda_t, " Lambda_y: ", min_lambda_y))
print(lambda_t_values)
print(lambda_y_values)
# -----------

intH_smooth <- zeta %*% t(min_gamma) %*% t(t_basisMat)

for (i in 1:28) {
  plot(tobs ,intH_smooth[i,] ,type = 'l',xlab = 't',ylab = 'h', main='h', 
       col = 'blue' , lwd = 1, ylim = c(-1,1))
  lines(tobs, Ystar[i,], type = 'l',col = 'red', lwd = 1)
}
