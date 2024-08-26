
library(fda)
library(MASS)
library(plotly)
library(glmnet)
library(dplyr)

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
    Psieval <- eval.basis(As[i, ], psi)
    psievalist[[i]] <- Psieval
  }
  return(psievalist)
}

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
  #return(list(U=u))
  return(u)
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

gamma.fun <- function(zeta,Y.Phi.Mat,basisphi,tobs,nbasis_psi){
  ZetaJYPhi    <- t(zeta) %*% Y.Phi.Mat
  vec.ZetaYPhi <- as.vector(ZetaJYPhi) # vec in col
  inprod_zeta  <- t(zeta)%*% zeta
  Jphiphi      <- eval.penalty(basisphi, Lfdobj =int2Lfd(0), tobs)
  B            <- kronecker(Jphiphi,inprod_zeta)
  Vec.gamma    <- ginv(B) %*% vec.ZetaYPhi
  gamma_Matrix <- matrix(Vec.gamma, nrow = basisphi$nbasis , ncol= nbasis_psi)
  return(gamma_Matrix)
}

# Prediton Mean Square Error 
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



speed_Men   = read.csv( "MenSpeed.csv", header = TRUE, row.name = 1)
#speed_Women = read.csv("WomenSpeed.csv",header = TRUE, row.name = 1)
stroke_Men  = read.csv("MenStroke.csv",header = TRUE, row.name = 1)
#stroke_Wonmen = read.csv("WomenStroke.csv",header = TRUE, row.name = 1)

speed_Men = as.matrix(speed_Men)
stroke_Men = as.matrix(stroke_Men)


rangeval <- c(0, 2000)
races_N      <- dim(speed_Men)[1]
meter <- seq(50, 2000, length.out = 500)
# generate phi basis function of m
dis_basisphi <- basisphi.gen(d = 3, M = 197, domain = rangeval)
dis_PhiMat   <- compute.phi ( d = 3, M = 197, domain = rangeval , t = meter)

print(paste0("dimension of Phi is: ", dim(dis_PhiMat)))

# generate psi basis function of stroke(s)
stroke_psi <- psi.generator (d = 3, M = 197, domain = c( min(stroke_Men), max(stroke_Men)))

# Eval stroke
stroke_psievalist <- psi.eval.fun(n = races_N , psi = stroke_psi , stroke_Men) # 41 x basis

nbasis_psi <- stroke_psi$nbasis
speed.Phi.Mat  <- compute.u(speed_Men , d = 3, M = 197,  domain = rangeval)

s    <- meter
zeta <- computeZeta(n = races_N, stroke_psi$nbasis , stroke_psievalist, s) # n*7

# solving gamma
Gamma = gamma.fun(zeta, speed.Phi.Mat, dis_basisphi, meter, stroke_psi$nbasis)

write.csv(Gamma, file = "Gamma.csv")

#  --------------------- H integral ---------------------------

int_H_pred  <- zeta %*% t(Gamma) %*% t(dis_PhiMat) # 2000 x 41
write.csv(int_H_pred, file = "int_H_pred.csv")

# ------------------ Plot H 3D predict surface----------------------

y <- seq(min(stroke_Men),max(stroke_Men), length.out = 500)
stroke_psi_y      <-  eval.basis(y, stroke_psi)

H   <- dis_PhiMat %*% Gamma %*% t(stroke_psi_y)

write.csv(H, file = "H.csv")

plott <- plot_ly(x = ~ meter, y = ~ y, z = ~ H, type = "surface")
layout_plot = layout(plott, scene = list(xaxis = list(title = "distance"),
                           yaxis = list(title = "Y"),
                           zaxis = list(title = "H"), title = "H"))


HMSE = PMSE(races_N , speed_Men, int_H_pred, s)
print(HMSE)

