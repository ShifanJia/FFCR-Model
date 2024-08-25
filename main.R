# functions for
# *** FFCR method
# *** close form
# *** FFLR (Ramsay and Silverman (2005))





# main function
library(fda)
library(ggplot2)
library(MASS)
library(zoo)
library(splines)
library(rgl)
library(plotly)
library(psych)
library(pracma)
library(grDevices)
library(Matrix)
library(glmnet)
library(npreg)

theta.basis.fun <- function(d, M, time_rangeval){
  norder  = d + 1
  nknots  = M + 1 # knots number
  knots  = seq(time_rangeval[1], time_rangeval[2], length.out = nknots)# b-spline knot points 
  # basis function number
  nbasis  = nknots + norder - 2 
  theta.basis = create.bspline.basis(time_rangeval, nbasis, norder)
  return(theta.basis)
}


theta.spline.generator = function(d, M, time_rangeval,tobs ) { 
  norder  = d + 1
  nknots  = M + 1 # knots number
  knots  = seq(time_rangeval[1], time_rangeval[2], length.out = nknots)# b-spline knot points 
  # basis function number
  nbasis  = nknots + norder - 2 
  theta.basis = create.bspline.basis(time_rangeval, nbasis, norder)
  # theta(t) matrix value
  basismat = eval.basis(tobs,theta.basis) # 1000 x 53
  return(basismat)
}

#-------------------------- X(t)--------------------------------------------

data.generator.bsplines = function(n, nknots, norder, p , domain) {
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2 
  basis    = create.bspline.basis(knots, nbasis, norder)
  tobs = seq(domain[1],domain[2],length.out = p)
  basismat = eval.basis(tobs, basis) 
  x=array(NA,c(n,p))
  for(i in 1:n){
    x[i,] = rnorm(nbasis, 0, 1)%*%t(basismat)
  }
  return(x)
}


# function with beta

b.fit <- function(X, Y, Basis_mat, t, k, R) {
  bhat <- matrix(0, nrow = ncol(Basis_mat), ncol =length(t))
  for (i in 1:length(t)) {
    Xtitle <- X[, i, drop = FALSE] %*% Basis_mat[i,, drop = FALSE]
    H <- ginv(t(Xtitle) %*% Xtitle + k*R) %*% t(Xtitle)
    bhat[, i] <- H %*% Y[, i, drop = FALSE]
  }
  
  return(bhat)
}

Beta.function = function(b, basismat, t)
{
  beta.fit <- numeric(length(t))
  for (i in 1:length(t)) {
    beta_fit <- basismat[i,,drop = FALSE]%*%b[,i,drop = FALSE]
    beta.fit[i] <- beta_fit
  }
  return(beta.fit)
}

# ----function with gamma-------
basisphi.gen = function(d , M , domain ){
  norder     = d + 1
  nknotsphi  = M + 1 # knots number
  knots_phi  = seq(domain[1], domain[2], length.out = nknotsphi)
  nbasis_phi = nknotsphi + norder - 2
  basis = create.bspline.basis(domain, nbasis_phi, norder)
  return(basis)
}

# evaluate basis function phi matrix

compute.phi  = function( d , M , domain , t){
  norder     = d + 1
  nknotsphi  = M + 1 # knots number
  knots_phi  = seq(domain[1], domain[2], length.out = nknotsphi)
  nbasis_phi = nknotsphi + norder - 2  # basis function number
  K          = nbasis_phi
  basisphi   = create.bspline.basis(domain, nbasis_phi, norder)
  PhiMat     = eval.basis(t , basisphi) # 
  return(PhiMat)
}


# -----------Y.Phi.Mat ---------------------

compute.u  = function(Y, d, M, domain){
  norder   = d + 1
  nknots   = M + 1
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots, nbasis , norder)
  p        = dim(Y)[2]
  tobs     = seq(domain[1], domain[2], length.out = p)
  basismat = eval.basis(tobs, basis) # p x M+d
  
  cef      = c(1,rep(c(4,2),(p-1)/2))
  cef[p]   = 1
  h        = 1/(  p- 1)
  
  u   = h/3* Y %*% diag(cef) %*% basismat
  return(u)
}

# ------------------ basis  psi--------------------------------

psi.generator = function(d , M , domain){
  norder    = d +1 
  nknotspsi = M + 1
  knots_psi = seq(domain[1], domain[2], length.out = nknotspsi)
  nbasis_psi = nknotspsi + norder - 2 #
  L      =  nbasis_psi
  psi.fun <- create.bspline.basis(domain , nbasis_psi, norder)
  return(psi.fun)
}

psi.eval.fun <- function(n, psi, As){ 
  psievalist <- list()
  
  for (i in 1 : n) {
    # Evaluate the basis function psi on the i-th row of the As matrix
    Psieval <- eval.basis(As[i, ], psi)
    psievalist[[i]] <- Psieval
  }
  return(psievalist)
}

computeZeta = function(n, nbasis_psi, psievalist, s){
  zeta = matrix(0, nrow = n, ncol = nbasis_psi)
  for (i in 1 : n) {
    for (j in 1 : nbasis_psi) {
      Psieval <- psievalist[[i]]
      delta   <- 1/length(s) 
      nA0     <- length(s)  
      coeffs  <- c(1, rep(c(4, 2), (nA0 - 1) / 2))
      coeffs[nA0]  <- 1
      quadwts      <- coeffs * delta / 3
      zeta[i , j]  <- quadwts %*% Psieval[ , j]
    }
  }
  return(zeta)
}


gamma.fun <- function(zeta, Y.Phi.Mat, basisphi, tobs, nbasis_psi, psi,y_order,lambda_y,lambda_t, R1,R3){
  ZetaJYPhi    <- t(zeta) %*% Y.Phi.Mat   # Compute the in-product of zeta with Y.Phi.Mat
  vec.ZetaYPhi <- as.vector(ZetaJYPhi)    # vec in col
  inprod_zeta  <- t(zeta)%*% zeta
  # Evaluate the penalty on the basis functions. This penalty could be related to 
  # the roughness or other characteristics of the function.
  Jphiphi      <- eval.penalty(basisphi, Lfdobj =int2Lfd(0), tobs) 
  Jpsipsi  <- eval.penalty(psi, Lfdobj = int2Lfd(0),  y_order)
  PenMat_y <- kronecker(Jphiphi,R1)
  PenMat_t <- kronecker(R3,Jpsipsi)
  B            <- kronecker(Jphiphi,inprod_zeta)
  Vec.gamma    <- ginv(B + lambda_y * PenMat_y + lambda_t * PenMat_t) %*% vec.ZetaYPhi
  gamma_Matrix <- matrix(Vec.gamma, nrow = basisphi$nbasis , ncol= nbasis_psi)
  return(gamma_Matrix)
}

Yfit.fun = function(t, X, beta){
  Ymat <- matrix(0, nrow = nrow(X), ncol =length(t))
  for (i in 1:length(t)) {
    Y_hat     <- X[,i,drop = FALSE]%*%beta[i]
    Ymat[,i]  <- Y_hat
  }
  return(Ymat)
}

Y2bMap = function(X , Basis_mat, tobs)
{
  Hmat <- list()
  for (i in 1:length(tobs)){
    Xtitle    <- X[,i, drop = FALSE]%*%Basis_mat[i,,drop = FALSE] #  X~(t) = X(t)theta^T(t)
    H.mat      <- ginv(t(Xtitle)%*%Xtitle)%*%t(Xtitle)
    Hmat[[i]] <-  H.mat  
  }
  Hmat<- do.call(cbind,Hmat)
  return(Hmat) 
}


YMSE = function(n, YTrue , Ypred, s){
  MSE = matrix(0, nrow = 1, ncol = n)
  for (i in 1 : n ) {
    Yerror <- YTrue - Ypred
    delta   <- 1/length(s) 
    nA0     <- length(s)  
    coeffs  <- c(1, rep(c(4, 2), (nA0 - 1) / 2))
    coeffs[nA0]  <- 1
    quadwts      <- coeffs * delta / 3
    MSE[i]  <- quadwts %*% Yerror[i,]^2
  }
  pmse = mean(MSE)
  return(pmse)
}

int.Yerror = function(test_indices, Ystar_test , YPred, s){
  ye = matrix(0, nrow = 1, ncol = test_indices)
  for (i in 1 : test_indices ) {
    Yerror <- Ystar_test - YPred
    delta   <- 1/length(s) 
    nA0     <- length(s)  
    coeffs  <- c(1, rep(c(4, 2), (nA0 - 1) / 2))
    coeffs[nA0]  <- 1
    quadwts      <- coeffs * delta / 3
    ye[i]  <- quadwts %*% (Yerror[i,])^2
  }
  pmse = mean(ye)
  return(pmse)
}

SSEY = function( Ypred, Ytilde, n , nbasis){ 
  SSE = numeric(length(dim(Y)[1]))
  for (i in 1 : n) {
    sse    <- t(Y[i,] - Ytilde[i,])%*%(Y[i,] - Ytilde[i,])
    SSE[i] <- sse
  }
  MSE = sum(SSE)/ (n - nbasis)
  return(MSE)
}

PMSE = function(n, YTrue , Ypred, s){
  MSE = matrix(0, nrow = 1, ncol = n)
  for (i in 1 : n ) {
    Yerror <- YTrue - Ypred
    delta   <- 1/length(s) 
    nA0     <- length(s)  
    coeffs  <- c(1, rep(c(4, 2), (nA0 - 1) / 2))
    coeffs[nA0]  <- 1
    quadwts      <- coeffs * delta / 3
    MSE[i]  <- quadwts %*% Yerror[i,]^2
  }
  pmse = mean(MSE)
  return(pmse)
}


ISE.fun <- function(m, n, h_hat, h ){
  ISE <- 0
  delta_t <- 1 / (m - 1)
  delta_s <- 1 / (n - 1)
  
  for (i in 1: m) {
    for (j in 1: n) {
      ISE <- ISE + get_coeff(i) * get_coeff(j) * (h_hat[i,j] - h[i,j])^2
    }
  }
  ISE <- (delta_t * delta_s) / 9 * ISE
  return(ISE)
}

ISE_beta <- function(tobs,beta_true,beta_smooth_value ){
  
  quadpts <- tobs
  delta  <- 1/length(tobs)
  nquadpts <- length(quadpts)
  quadwts <- as.vector(c(1,rep(c(4,2),(nquadpts-2)/2),4,1),mode="any")
  quadwts <- c(1,rep(c(4,2),(nquadpts-1)/2))
  quadwts[nquadpts] <- 1
  quadwts <- quadwts*delta/3
  
  ISE <- sum(quadwts*(beta_true - beta_smooth_value)^2)
  return(ISE)
}







#---- GIVE phi psi R1 R2-------

basisphi <- basisphi.gen(d , M , domain = time_rangeval)
PhiMat   <- compute.phi ( d , M , domain = time_rangeval , t= tobs)
psi      <- psi.generator (d  , M , domain = c( min(As), max(As)))
Psiyorder  = eval.basis(y, psi)
psievalist <- psi.eval.fun(n =n, psi = psi, As) # N dim:P*7
nbasis_psi <- psi$nbasis

s    <- seq(0, 1, length.out = p)
zeta <- computeZeta(n = n, nbasis_psi, psievalist, s) # n*7

R1  <- eval.penalty(psi, Lfdobj = int2Lfd(2),  y) # L*L
R3  <- eval.penalty(basisphi, Lfdobj =int2Lfd(2), tobs) # K*K

R2       <- eval.penalty(theta.basis, Lfdobj =int2Lfd(2), tobs)

#---Step 1. Begin with an initial estimate b0-----------

b_coef0  <- b.fit(XTrue , YTrue, basismat, tobs, k= 0.1, R = R2)
Betafit  <- Beta.function(b = b_coef0 ,basismat, tobs)

Ystar1 = Yfit.fun(tobs,  XTrue , Betafit) 
nbasis <- dim(basismat)[2]
#MSE    <- t(Y - Ytilde)%*%(Y - Ytilde)/(n - nbasis)
#MSE = SSEY(Y, Ytilde, n , nbasis )

# for (i in 1:100) {
#   plot(Ytilde[i,], Y[i,], xlab = 'Y observe',ylab = 'Y Fit',
#        main='Fitted vs. Observed')
#   abline(a = 0,  b = 1)
# }

# -------step 2 solve gamma-----------------------------------------------

W_1  <- compute.u(Ystar1 , d = 3, M = 20,  domain = time_rangeval)

gamma_1 = gamma.fun(zeta, W_1 , basisphi, tobs, nbasis_psi, psi, y, 0.01, 0,01, R1,R3)

hpredorder1 <- PhiMat  %*% gamma_1 %*% t(Psiyorder)

# -------step 3 update beta-----------------------------------------------

Ytilde_2 <- YTrue - zeta %*% t(gamma_1) %*% t(PhiMat)
b_coef1  <- b.fit(XTrue , Ytilde_2 , basismat, tobs, k= 0.1, R = R2)
Betafit1  <- Beta.function(b = b_coef1 ,basismat, tobs)

# --------plot smooth beta and no smooth beta------------------
plot(tobs, Betafit, type='l', )
lines(tobs, Beta_true, col='red')

#---------next iertation----------------------------------
Ystar2 = Yfit.fun(tobs,  XTrue , Betafit1) 

W_2  <- compute.u(Ystar2 , d = 3, M = 20,  domain = time_rangeval)

gamma_2 = gamma.fun(zeta, W_2 , basisphi, tobs, nbasis_psi, psi, y, 0.01, 0.01, R1, R3)

hpredorder2 <- PhiMat  %*% gamma_2 %*% t(Psiyorder)

open3d()
persp3d(t , y_order, hpredorder1, 
        col = 'lightgreen',
        xlab = 't', 
        ylab = 'z',
        zlab = 'h(t,z)',
        axes = TRUE,
        border = NA, 
        facets = FALSE,
        xlim = c(0, 2*pi/10),
        zlim  = c(min(hsurface2)+1, max(hsurface2)))



fold_indices <- split(1:n, cut(1:n, breaks = 10, labels = FALSE))
results <- list()  # Initialize a list to store results
average_rmse <- numeric()  # Vector to store the average RMSE of each combination
params_labels <- character()  # Vector to store the corresponding parameter labels

#---------------given b0--------------
b_coef0   <- b.fit(XTrue , YTrue, basismat, tobs, k= 0.1, R = R2)
Betafit_0 <- Beta.function(b = b_coef0 ,basismat, tobs)
Betafit_0 <- matrix(diag(Betafit_0 ), nrow =length(tobs) , ncol =length(tobs))
Ystar1    = Yfit.fun(tobs,  XTrue , Betafit_0) 


fold_indices <- split(1:n, cut(1:n, breaks = 10, labels = FALSE))
candidate_values <- c(0.01, 0.1, 0.5, 1)
results <- list()  # Initialize a list to store results
average_rmse <- numeric()  # Vector to store the average RMSE of each combination
params_labels <- character()  # Vector to store the corresponding parameter labels

for (fold in 1: 10) {
  test_indices  <- fold_indices[[fold]]
  train_indices <- setdiff(1:n, test_indices)
  
  XTrain = XTrue[train_indices,] # N_train x T
  YTrain = YTrue[train_indices,] # N_train x T
  ZetaTrain = zeta[train_indices,]
  XTest = XTrue[test_indices,] # N_test x T
  YTest = YTrue[test_indices,] # N_test x T
  ZetaTest = zeta[test_indices,]
  
  # Iterate over each combination of parameter values
    for (lambda_y in candidate_values) {
      for (lambda_t in candidate_values) {
        Ystar1  = YTrain - XTrain %*% t(Betafit_0)
        
        W_1     <- compute.u(Ystar1 , d = 3, M = 20,  domain = time_rangeval)
        
        gamma_1 = gamma.fun(ZetaTrain, W_1 , basisphi, tobs, nbasis_psi, psi, y, lambda_y,lambda_t, R1,R3)
        
        #---------get cv score-------
        integralH_pred  <- ZetaTest %*% t(gamma_1) %*% t(PhiMat)
        Ypred           <- XTest %*% BetaTruemat + integralH_pred
        RMSE            <- YMSE(n=20, YTest, Ypred, s)

        # Store RMSE with parameter combination information
        params <- paste( "lambda_y =", lambda_y, "lambda_t =", lambda_t, sep=",")
        if (!is.list(results[[params]])) {
          results[[params]] <- numeric(10)  # Initialize numeric vector if not already initialized
        }
        results[[params]][fold] <- RMSE
      }
    }
  }

# Calculate average RMSE for each parameter combination
for (param in names(results)) {
  average_rmse <- c(average_rmse, mean(results[[param]], na.rm = TRUE))
  params_labels <- c(params_labels, param)
}

# Identify the minimum RMSE and corresponding parameters
min_index <- which.min(average_rmse)
min_rmse <- average_rmse[min_index]
best_params <- params_labels[min_index]

cat("The minimum average RMSE is", min_rmse, "with parameters", best_params, "\n")

# Calculate average RMSE for each parameter combination and find the best
for (param in names(results)) {
  avg_rmse = mean(results[[param]], na.rm = TRUE)
  average_rmse <- c(average_rmse, avg_rmse)
  params_labels <- c(params_labels, param)
  
  if (!exists("min_rmse") || avg_rmse < min_rmse) {
    min_rmse <- avg_rmse
    best_params <- param
    best_gamma  <- gamma_1  # Update the best gamma
  }
}

# --------------iterantion update b---------
Ytilde_2 <- YTrue - zeta %*% t(best_gamma) %*% t(PhiMat)

fold_indices <- split(1:n, cut(1:n, breaks = 10, labels = FALSE))
results <- list()  # Initialize a list to store results
average_rmse <- numeric()  # Vector to store the average RMSE of each combination
params_labels <- character()  # Vector to store the corresponding parameter labels

for (fold in 1: 10) {
  test_indices  <- fold_indices[[fold]]
  train_indices <- setdiff(1:n, test_indices)
  
  XTrain    = XTrue[train_indices,] # N_train x T
  YTrain    = Ytilde_2[train_indices,] # N_train x T
  XTest     = XTrue[test_indices,] # N_test x T
  YTest     = Ytilde_2[test_indices,] # N_test x T
  # Iterate over each combination of parameter values
  for (k in candidate_values) {
        b_s = b.fit(XTrain , YTrain, basismat, tobs, k, R = R2)
        
        Betafit_s = Beta.function(b = b_s ,basismat, tobs)
        
        Betafit_s   <- matrix(diag(Betafit_s ), nrow =length(tobs) , ncol =length(tobs))
        #---------get cv score-------
        Ypred           <- XTest %*% BetaTruemat
        RMSE            <- YMSE(n = 20, YTest, Ypred, s)

        # Store RMSE with parameter combination information
        params <- paste("k=", k, sep=",")
        if (!is.list(results[[params]])) {
          results[[params]] <- numeric(10)  # Initialize numeric vector if not already initialized
        }
        results[[params]][fold] <- RMSE
      }
    }

# Calculate average RMSE for each parameter combination
for (param in names(results)) {
  average_rmse <- c(average_rmse, mean(results[[param]], na.rm = TRUE))
  params_labels <- c(params_labels, param)
}

# Identify the minimum RMSE and corresponding parameters
min_index <- which.min(average_rmse)
min_rmse <- average_rmse[min_index]
best_params <- params_labels[min_index]

cat("The minimum average RMSE is", min_rmse, "with parameters", best_params, "\n")



# ------------Initialization--------------
# --------------Full iterantion ---------

max_iterations <- 50  # Limit on number of iterations to prevent infinite loops
tolerance <- 1e-3     # Convergence criterion for changes in Beta
converged <- FALSE
iteration <- 0
candidate_values <- c(1, 10, 100)  # Candidate values for tuning parameters

# Initial computation of b and beta
b_0 <- b.fit(XTrue, YTrue, basismat, tobs, k = 1, R = R2)
Betafit <- Beta.function(b = b_coef, basismat, tobs)
Betafit_0 <- matrix(diag(Betafit), nrow = length(tobs), ncol = length(tobs))

# Main loop for iterative update
while (!converged && iteration < max_iterations) {
  iteration <- iteration + 1
  previous_Betafit <- Betafit_0  # Save the previous Beta for convergence check
  
  # Cross-validation for gamma selection
  best_gamma_value <- NULL
  min_rmse_gamma   <- Inf
  
  for (lambda_y in candidate_values) {
    for (lambda_t in candidate_values) {
      fold_rmse <- c()
      for (fold in 1:10) {
        train_indices <- setdiff(1:n, fold_indices[[fold]])
        test_indices <- fold_indices[[fold]]
        # Training data preparation
        ZetaTrain <- zeta[train_indices,]
        ZetaTest  <- zeta[test_indices,]
        Ystar_Train    <- YTrue[train_indices,] - XTrue[train_indices,] %*% previous_Betafit
        Ystar_Test     <- YTrue[test_indices,]  - XTrue[test_indices,] %*% previous_Betafit
        
        W         <- compute.u(YTrain, d = 3, M = 20, domain = time_rangeval)
        # Gamma computation
        gamma_temp <- gamma.fun(ZetaTrain, W, basisphi, tobs, nbasis_psi, psi, y, lambda_y, lambda_t, R1, R3)
        # Validation data preparation and RMSE calculation
       
        integralH_pred <- ZetaTest %*% t(gamma_temp) %*% t(PhiMat)
        Ypred      <- XTest %*% Betafit_0 + integralH_pred
        RMSE       <- sqrt(mean((YTest - Ypred)^2))
        fold_rmse  <- c(fold_rmse, RMSE)
      }
      
      avg_rmse <- mean(fold_rmse)
      if (avg_rmse < min_rmse_gamma) {
        min_rmse_gamma   <- avg_rmse
        best_gamma_value <- list(lambda_y = lambda_y, lambda_t = lambda_t, gamma = gamma_temp)
      }
    }
  }
  
  # Update gamma based on best parameters found in CV
  gamma <- best_gamma_value$gamma
  
  # Cross-validation for beta selection
  best_k_value <- NULL
  min_rmse_beta <- Inf
  
  for (k in candidate_values) {
    fold_rmse <- c()
    for (fold in 1:10) {
      train_indices <- setdiff(1:n, fold_indices[[fold]])
      test_indices <- fold_indices[[fold]]
      YTide_Train  <- YTrue[train_indices,] - ZetaTrain %*% t(gamma) %*% t(PhiMat)
      YTide_Test   <- YTrue[test_indices,] - ZetaTest %*% t(gamma) %*% t(PhiMat)
      
      # Training and prediction under current k
      
      b_s <- b.fit(XTrain, YTide_Train , basismat, tobs, k, R = R2)
      Betafit_s <- Beta.function(b = b_s, basismat, tobs)
      Betafit_s <- matrix(diag(Betafit_s), nrow = length(tobs), ncol = length(tobs))
      
      Ypred <- XTest %*% Betafit_s
      
      # RMSE calculation
      RMSE      <- ISE_beta(tobs, Beta_true, diag(Betafit_s))
      fold_rmse <- c(fold_rmse, RMSE)
    }
    
    avg_rmse <- mean(fold_rmse)
    if (avg_rmse < min_rmse_beta) {
      min_rmse_beta <- avg_rmse
      best_k_value  <- k
    }
  }
  
  # Update Beta based on the best k found in CV
  b_coef <-  b.fit(XTrue, YTrue - zeta %*% t(gamma) %*% t(PhiMat), basismat, tobs, best_k_value, R = R2)
  Betafit <- Beta.function(b = b_coef, basismat, tobs)
  Betafit <- matrix(diag(Betafit), nrow = length(tobs), ncol = length(tobs))
  
  # Convergence check
  if (max(abs(Betafit - previous_Betafit)) < tolerance) {
    converged <- TRUE
  }
}

# Output final parameters and check convergence status
cat("Iteration completed: ", iteration, "\n")
cat("Convergence Status: ", ifelse(converged, "Converged", "Not Converged"), "\n")

