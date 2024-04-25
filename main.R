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

