# INPUTs
#   n: The observation times for Y_i X_i and A_i(s), i= 1....n
#   T: Time (0 ,1 , by = 0.01)
#   p: length of T
#   Y: matrix of n x p 
#   X: matrix of n x p
#   Beta: matrix of J x J, one covariate
#   Gamma: matrix of K x L
#   J: integer, the number of basis function of Beta(t)
#   d: integer, the degree of B-splines (cubic) using 4
#   K: integer, the number of basis function of h(t,A(s)) for t
#   L: integer, the number of basis function of h(t,A(s)) for A(s)
#   lambda: real number, tuning parameter for roughness penalty

# ------Define the number of individuals and time points T------
n  <- 200 
p  <- 500
T  <- 1
tobs <- seq(0, T, length.out = p)
time_rangeval  <- c(0, T)


#-------------------------- simulate X(t)--------------------------------------------

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

XTrue <- data.generator.bsplines(n = n, nknots= 40, norder = 4, p = p,domain = time_rangeval)

XTrue <- XTrue -  apply(XTrue,1,mean)

# ------------------------ simulate A(s)-----------------------------

A.fun = function(s) { 
  c1 <- runif(1, 0, 0.05)
  c2 <- runif(1, 0, 0.5)
  c3 <- runif(1, 0, 0.5)
  c4 <- runif(1, 0, 0.5)
  c5 <- runif(1, 0, 0.5)
  # return(c1*2 + c2*sin(pi * s/4) + c3* cos(pi*s/4)
  #        + c4* sin(pi*s/2)  + c5* cos(pi*s/2))
  return(c1 + cos(c2)*sin(pi * s/2) + sin(c3)* cos(pi*s/2)
         + cos(c4)* sin(pi*s)  + sin(c5)* cos(pi*s))
  
}

s          <- seq(0,T,length.out = p)

simulate.As <- function(num_samples, num_points, s) {
  data <- matrix(0, nrow = num_samples, ncol = num_points)
  for (i in 1:num_samples) {
    data[i, ] <- A.fun(s)
  }
  return(data)
} # N*P

As <- simulate.As(num_samples = n, num_points = p, s) # simulate A(s) in (0,1)

# ----------------------Simulate Beta(t)--------------------

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


# Here is two simulation setting I and II

#----------------------scenario I: beta2(t) ---------------- 
Beta_true <- cos(2*pi*tobs)

#----------------------scenario II: beta1(t) ---------------- 
Beta_true <- sin(2*pi*tobs)

# use FDA package to change beta into functional object
Beta.fd <- fd(Beta_true)

# Plot Beta(t)
plot(tobs, Beta_true, xlab = "t", ylab = expression(beta), type = 'l',col = 'red')
#plot.fd(Beta.fd, xlab = "t", ylab = expression(beta))

# Beta(t) Matrix

BetaTruemat   <- matrix(diag(Beta_true), nrow =length(tobs) , ncol =length(tobs))

# -----------------------simulate h(t, y) ------------------------------

# simulation I 
h.fun = function(t, y){
  h = cos(10*t - 5 - (5*(y - 0.5)))
  return(h)
}

# simulation II 

h.fun <- function(t,y){
  z =  t*exp(y)
  return(z)
}

# --------------plot H surface---------------
t <- tobs  
#y <- As[1,]
y <- seq(min(As),max(As),length.out = p)
y_order   <- sort(y)
hsurface  <- outer(t, y , h.fun)
hsurface2 <- outer(t, y_order, h.fun)
hsurface2 <- hsurface2 - apply(hsurface2, 1, mean) # centering data

# plot <- plot_ly(x = ~ t, y = ~ y_order, z = ~ hsurface2, type = "surface")
# x_range <- c(0, 2*pi/10)
# layout(plot, scene = list(xaxis = list(title = "t", range = x_range), 
#                           yaxis = list(title = "z"), 
#                           zaxis = list(title = "h(t,z)")))

plot <- plot_ly(x = ~ t, y = ~ y_order, z = ~ hsurface2, type = "surface")
layout(plot, scene = list(xaxis = list(title = "t"), 
                          yaxis = list(title = "z"), 
                          zaxis = list(title = "h(t,z)")))

# ---------------------- 3D Surface Plot ---------------------

# png("H2.png", width = 1600, height = 1200)
# plot3D::persp3D(t, y_order, hsurface2, theta = -40, phi = 20, col = NULL,
#        xlab = 't',ylab = 'A(s)', zlab = 'h(t,A(s))', ticktype = "detailed",
#        colkey = list(side = 1, length = 0.5, dist = -0.1),
#        plot= TRUE)
# dev.off()


open3d()
persp3d(t , y_order, hsurface2,theta = 30, phi = 20,col = 'lightblue',
        xlab = 't', ylab = 'z', zlab = 'h(t,z)',
        axes = TRUE, 
        cex.lab = 1.5,  
        zlim  = c(-3, 5),
        border = NA, facets = FALSE)



open3d()
persp3d(t , y_order, hsurface2, 
        col = 'lightgreen',
        xlab = 't', 
        ylab = 'z',
        zlab = 'h(t,z)',
        axes = TRUE,
        border = NA, 
        facets = FALSE,
        xlim = c(0, 2*pi/10),
        zlim  = c(min(hsurface2), max(hsurface2)))


# png 1 xlim = c(0,0.8) simulation 1 

# -------------------Fix a specific y to plot h -----------
for (i in 1:100) {
  h_fixed_y        <- hsurface2[, i]
  plot(tobs, h_fixed_y, type = "l", xlab = "t", ylab = "h")
}
# ------------------Fix a specific t  plot h ---------------

for (j in 200 : 256) {
  h_t        <- hsurface2[j,]
  plot(y_order, h_t, type = "l", xlab = "A1(s)", ylab = "h",
       main = paste("h vs y for fixed_t"))
}

# -----------------h min and h max--------------------
print(paste("h in ", min(hsurface2), "to", max(hsurface2)))

#---------------compute the integral of h surface---------------
hintergral.fun <- function(n, p, tobs, As){
  U  <- matrix(0, nrow = n, ncol = p) # n * p
  hh <- matrix(0, nrow = n, ncol = p )
  for (i in 1:length(tobs)) {
    for (j in 1: n) {
      hh[j,]   = h.fun(t = tobs[i], y = As[j,])
      cef      = c(1,rep(c(4,2),(p-1)/2))
      cef[p]   = 1
      h        = 1/(  p- 1)
      U[j,i]   = h/3* t(cef) %*% hh[j,]
    }
  }
  return(U)
}

integralHTrue <- hintergral.fun(n = n, p, tobs, As) # N*P


# ---------------Signal-to-Noise ratio-----------------------------------------
y0   = XTrue %*% BetaTruemat + integralHTrue 

#y0- apply(y0,1,mean)

#eps0 = apply(y0- apply(y0,1,mean), 1 , sd)

eps0 = apply(y0, 1 , sd)
SNR1  = 5
SNR2  = 10


#--------------------------noise-----------------------------------------------

#noise1 = matrix(rnorm(n * length(tobs), mean = 0, sd = eps0/sqrt(SNR1)), nrow = n)
#noise = matrix(rnorm(n * length(tobs), mean = 0, sd = y0_sd/sqrt(SNR1)), nrow = n)

noise1 <- matrix(nrow = nrow(y0), ncol = length(tobs))

# Generate data for each row using rnorm with sd[i]
for (i in 1:nrow(y0)) {
  noise1[i, ] <- rnorm(length(tobs), 0, eps0[i]/sqrt(SNR1))
}

noise2 <- matrix(nrow = nrow(y0), ncol = length(tobs))
for (i in 1:nrow(y0)) {
  noise2[i, ] <- rnorm(length(tobs), 0, eps0[i]/sqrt(SNR2))
}


# --Y(t) = X(t)Beta(t) + \int h(t,A(s))ds + error(t) 
# ------------------ Y(t)----------------------------------------------

YTrue = XTrue %*% BetaTruemat + integralHTrue + noise1 # N*P under snr=5



basisphi   <- basisphi.gen(d = 3, M = 47, domain = time_rangeval)
PhiMat   <- compute.phi (d = 3,  M = 47, domain = time_rangeval , t= tobs)
psi      <- psi.generator(d = 3, M = 47, domain = c( min(As), max(As)))
Psiyorder  = eval.basis(y,psi)
psievalist <- psi.eval.fun(n =n, psi = psi, As) # N dim:P*7
nbasis_psi <- psi$nbasis

zeta <- computeZeta(n = n, nbasis_psi, psievalist, s) # n*7

R1  <- eval.penalty(psi, Lfdobj = int2Lfd(2),  y) # L*L
R3  <- eval.penalty(basisphi, Lfdobj =int2Lfd(2), tobs) # K*K

# -------------given b0----------------
theta.basis <- theta.basis.fun(d = 3, M = 47 ,c(0, 1))
basismat    <- theta.spline.generator(d = 3, M = 47 ,time_rangeval= c(0,1), tobs = tobs)
R2          <- eval.penalty(theta.basis, Lfdobj =int2Lfd(2), tobs)

b_0       = b.fit(XTrue , YTrue, basismat, tobs, k = 1, R2)
Betafit_0 = Beta.function(b = b_0 , basismat, tobs)

Betafit_0  <- matrix(diag(Betafit_0 ), nrow =length(tobs) , ncol =length(tobs))

# -------------s = 1 given b solve r----------

Ystar  = YTrue - XTrue%*% t(Betafit_0)
fold_indices <- split(1:n, cut(1:n, breaks = 10, labels = FALSE))
candidate_values <- seq(1e-10,10000, length=10)
results <- list()  
average_rmse <- numeric()  
params_list <- list()
params_labels <- character()

for (fold in 1: 10) {
  test_indices  <- fold_indices[[fold]]
  train_indices <- setdiff(1:n, test_indices)
  
  XTrain    = XTrue[train_indices,] # N_train x T
  YTrain    = YTrue[train_indices,] # N_train x T
  ZetaTrain = zeta[train_indices,]
  XTest     = XTrue[test_indices,] # N_test x T
  YTest     = YTrue[test_indices,] # N_test x T
  ZetaTest  = zeta[test_indices,]

    for (lambda_y in candidate_values) {
    for (lambda_t in candidate_values) {
      Ystar_1  = YTrain - XTrain %*% t(Betafit_0)
      
      W_1  <- compute.u(Ystar_1 , d = 3, M = 47,  domain = c(0,1))
      
      gamma_1 = gamma.fun(ZetaTrain, W_1 , basisphi, tobs, nbasis_psi, psi, y, lambda_y,lambda_t, R1,R3)
      
      #---------get cv score-------
      integralH_pred  <- ZetaTest %*% t(gamma_1) %*% t(PhiMat)
      Ypred           <- XTest %*% BetaTruemat + integralH_pred
      RMSE            <- YMSE(n = 20, YTest, Ypred, s)
      
      # Store RMSE with parameter
      params <- c(lambda_y, lambda_t)  # Store parameters as a numeric vector
      
      params_name <- paste(params, collapse=",")  # Create a string name for list indexing
      if (!is.list(results[[params_name]])) {
        results[[params_name]] <- numeric(10)
      }
      results[[params_name]][fold] <- RMSE
      params_list[[params_name]] <- params  # Store parameters in a list for easy retrieval
    }
    }
}

for (param in names(results)) {
  average_rmse <- c(average_rmse, mean(results[[param]], na.rm = TRUE))
  params_labels <- c(params_labels, param)
}

min_index <- which.min(average_rmse)
min_rmse <- average_rmse[min_index]
best_params <- params_list[[params_labels[min_index]]]

cat("The minimum average RMSE is", min_rmse, "with parameters lambda_y =", best_params[1], "lambda_t =", best_params[2], "\n")
gamma_1 = gamma.fun(ZetaTrain, W_1, basisphi, tobs, nbasis_psi, psi, y, best_params[1], best_params[2], R1, R3)


#---------- Next iteration given r solve b--------
Ytilde_1 <- YTrue - zeta %*% t(gamma_1) %*% t(PhiMat)
tau <- seq(1e-20, 10000, length = 20)

results <- list()  # Initialize a list to store results
average_rmse <- numeric()  # Vector to store the average RMSE of each combination
params_list <- list()
params_labels <- character()

for (fold in 1: 10) {
  test_indices  <- fold_indices[[fold]]
  train_indices <- setdiff(1:n, test_indices)
  
  XTrain    = XTrue[train_indices,]    # N_train x T
  YTrain    = Ytilde_1[train_indices,] # N_train x T
  XTest     = XTrue[test_indices,]     # N_test x T
  YTest     = Ytilde_1[test_indices,]  # N_test x T
  # Iterate over each combination of parameter values
  for (k in tau) {
    b_s = b.fit(XTrain , YTrain, basismat, tobs, k, R = R2)
    
    Betafit_s <- Beta.function(b = b_s ,basismat, tobs)
    
    Betafit_s   <- matrix(diag(Betafit_s ), nrow = length(tobs) , ncol =length(tobs))
    
    #---------get cv score-------
    
    IMSE = PMSE(n = 20, YTest , XTest%*%Betafit_s, s)
    params_name <- paste("k", k, sep=",")  # Used for list naming
    if (!is.list(results[[params_name]])) {
      results[[params_name]] <- numeric(10)
    }
    results[[params_name]][fold] <- IMSE
    params_list[[params_name]] <- k  # Directly store the numeric value of k
  }
}


for (param_name in names(results)) {
  average_rmse <- c(average_rmse, mean(results[[param_name]], na.rm = TRUE))
  params_labels <- c(params_labels, param_name)  # Carry over names for referencing
}

min_index <- which.min(average_rmse)
min_rmse <- average_rmse[min_index]
best_params <- params_list[[params_labels[min_index]]]  # Directly retrieve the best k value

cat("The minimum average RMSE is", min_rmse, "with parameter k =", best_params, "\n")

b_1       = b.fit(XTrue , Ytilde_1, basismat, tobs, k = best_params , R = R2)
Betafit_1 = Beta.function(b = b_1 ,basismat, tobs)
Betafit_1  <- matrix(diag(Betafit_1 ), nrow =length(tobs) , ncol =length(tobs))


#-- next s=2 given b solve r ----------
candidate_values <- seq(1e-10, 1000, length=10)

results <- list()  # Initialize a list to store results
average_rmse <- numeric()  # Vector to store the average RMSE of each combination
params_list <- list()
params_labels <- character()

for (fold in 1: 10) {
  test_indices  <- fold_indices[[fold]]
  train_indices <- setdiff(1:n, test_indices)
  
  XTrain    = XTrue[train_indices,] # N_train x T
  YTrain    = YTrue[train_indices,] # N_train x T
  ZetaTrain = zeta[train_indices,]
  XTest     = XTrue[test_indices,] # N_test x T
  YTest     = YTrue[test_indices,] # N_test x T
  ZetaTest  = zeta[test_indices,]
  # Iterate over each combination of parameter values
  for (lambda_y in candidate_values) {
    for (lambda_t in candidate_values) {
      Ystar_2  = YTrain - XTrain %*% t(Betafit_1)
      
      W_2  <- compute.u(Ystar_2 , d = 3, M = 47,  domain = c(0,1))
      
      gamma_2 = gamma.fun(ZetaTrain, W_2 , basisphi, tobs, nbasis_psi, psi, y, lambda_y,lambda_t, R1,R3)
      
      #---------get cv score-------
      integralH_pred  <- ZetaTest %*% t(gamma_2) %*% t(PhiMat)
      Ypred           <- XTest %*% BetaTruemat + integralH_pred
      RMSE            <- YMSE(n=20, YTest, Ypred, s)
      
      # Store RMSE with parameter
      params <- c(lambda_y, lambda_t)  # Store parameters as a numeric vector
      
      params_name <- paste(params, collapse=",")  # Create a string name for list indexing
      if (!is.list(results[[params_name]])) {
        results[[params_name]] <- numeric(10)
      }
      results[[params_name]][fold] <- RMSE
      params_list[[params_name]] <- params  # Store parameters in a list for easy retrieval
    }
  }
}

for (param in names(results)) {
  average_rmse <- c(average_rmse, mean(results[[param]], na.rm = TRUE))
  params_labels <- c(params_labels, param)
}

min_index <- which.min(average_rmse)
min_rmse <- average_rmse[min_index]
best_params <- params_list[[params_labels[min_index]]]  # Directly retrieve the numeric vector

cat("The minimum average RMSE is", min_rmse, "with parameters lambda_y =", best_params[1], "lambda_t =", best_params[2], "\n")
gamma_2 = gamma.fun(ZetaTrain, W_2, basisphi, tobs, nbasis_psi, psi, y, best_params[1], best_params[2], R1, R3)

hpredorderx     <- PhiMat  %*% gamma_2 %*% t(Psiyorder)

ISEh <- ISE.fun(m = dim(hpredorderx)[1], n = dim(hpredorderx)[1], h_hat = hpredorderx, h = hsurface2 )
print(ISEh)
ISEh

ISEbeta = ISE_beta(tobs, Beta_true , diag(Betafit_1))
print(ISEbeta)
ISEbeta

integralH_pred  <- ZetaTest %*% t(gamma_2) %*% t(PhiMat)
Ypred           <- XTest %*% BetaTruemat + integralH_pred
RMSE = YMSE(n = 20, YTest, Ypred, s)
print(RMSE)
RMSE
# Data frame creation for saving
results_df <- data.frame(
  Metric = c("ISEh", "ISEbeta", "RMSE"),
  Value = c(ISEh, ISEbeta, RMSE)
)

results_df
# Save to CSV file with dynamic name
write.csv(results_df, file_name, row.names = FALSE)



