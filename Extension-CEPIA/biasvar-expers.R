library(mvCircular)
library(ggplot2)
library(cowplot)
library(pracma)
library(prallel)

# EXPER RESULT STRG
expset_A <- list()
expset_A_reg <- list()
expset_B <- list()
expset_B_reg <- list()
expset_C <- list()
expset_C_reg <- list()

### Function that executes the experiment
mseExp <- function(kappa, lambda, 
                   phi = NULL, #diag(kappa) - lambda,
                   H = NULL, #matrix(1, ncol = length(kappa), nrow = length(kappa)),
                   nrep = 100, 
                   nmin = 10, 
                   nmax = 5000, 
                   nfact = 2, 
                   cl = NULL){

  # Generate N values
  nlist <- nmin * nfact ^ (0:floor(log(nmax,nfact) - log(nmin,nfact)))
  # Add nmax if it is not in the list
  if (nlist[length(nlist)] != nmax) nlist <- c(nlist, nmax)

  ##
  # Create Distribution
  ##
  x <- mvVonMises( rep(0,length(kappa)), kappa = kappa, lambda)
  
  ##
  # Generate samples
  ##
  x.samples <- getSamples(x, nrep*nmax)
  
  # Given parameter vector
  paramVector <- c(kappa, lambda[upper.tri(lambda)])

  # Execute experiment
  ret <- parallel::parLapplyLB(cl,nlist, function(n, nrep, data, real, phi, H){
  #ret <- lapply(nlist, function(n, nrep, data, real){
    
    # Compute Penalized estimator
    vals <- sapply(1:nrep, function(x, data, size){
      fit <- fit_mvvonmises(data[sample(nrow(data),size = size),], phi = phi, H = H)
      return(c(fit$kappa, fit$lambda[upper.tri(fit$lambda)]))
    },data, n)
    
    # Transpose vals
    vals <- t(vals)
    
    # Expected Bias value
    bias <- base::colMeans(vals) - real
    
    # Exp. Var. matrix
    variance <- circular::var.default(vals)
    
    mse <- sum(diag(variance)) + sqrt(sum(bias ^ 2))
    
    # Return mean and variance
    return(list(n = n, bias = bias, var = variance, mse = mse ))
  },nrep = nrep,data = x.samples, real = paramVector, phi = phi, H = H )
  

  results <- as.data.frame(t( sapply(ret, function(x){
    return( c( n = x$n, bias = sqrt(sum(x$bias ^ 2)), var = sum(diag(x$var)), mse = x$mse  ))
  })))
  
  colnames(results) <- c("n","bias","var", "mse")
  
  return(list(
    dist = x,
    results = ret,
    results_table = results
  ))
}





################################
####### EXPERIMENT SETUP A #####
# 
# Eigenvalues range is small but with different centers (10,1,0.1)
# Lambda matrix is 0
#
################################

# CONFIG
nrep = 100
nmax = 5000
nmin = 10
nfact = 2
p = 3
ncores = 6

cl <- parallel::makeForkCluster(nnodes = ncores)

# Kappa with close-to-0 inf. norm
kappa_0 = abs(rnorm(p,0,0.1)) + 1E-6

# Kappa with close-to-1 inf. norm
kappa_1 = rnorm(p,0,0.1) + 1

# Kappa with greater than 1 inf. norm
kappa_g = rnorm(p,0,1) + 10

lambda <- matrix(0, ncol = p, nrow = p)

## A.1
expset_A[[1]] <- mseExp(kappa_0, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)
expset_A_reg[[1]] <- mseExp(kappa_0, lambda,
                            phi = diag(kappa_0) - lambda,
                            H = matrix(1, ncol = p, nrow = p),
                        nrep, nmin, nmax,
                        nfact, cl)

## A.2
expset_A[[2]] <- mseExp(kappa_1, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)
expset_A_reg[[2]] <- mseExp(kappa_1, lambda,
                            phi = diag(kappa_1) - lambda,
                            H = matrix(1, ncol = p, nrow = p),
                        nrep, nmin, nmax,
                        nfact, cl)

## A.3
expset_A[[3]] <- mseExp(kappa_g, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)
expset_A_reg[[3]] <- mseExp(kappa_g, lambda,
                            phi = diag(kappa_g) - lambda,
                            H = matrix(1, ncol = p, nrow = p),
                        nrep, nmin, nmax,
                        nfact, cl)

parallel::stopCluster(cl)

################################

#########
## Plots
#########

bias_plot <- ggplot( expset_A[[1]]$results_table, aes(n) ) + 
   geom_line( data =  expset_A[[1]]$results_table, mapping = aes(y = bias, color = "Zero"), size = 1.2) +
   geom_line( data =  expset_A[[2]]$results_table, mapping = aes(y = bias, color = "One"), size = 1.2) +
   geom_line( data =  expset_A[[3]]$results_table, mapping = aes(y = bias, color = "Large"), size = 1.2) +
  geom_line( data =  expset_A_reg[[1]]$results_table, mapping = aes(y = bias, color = "Zero_R"), size = 1.2) +
  geom_line( data =  expset_A_reg[[2]]$results_table, mapping = aes(y = bias, color = "One_R"), size = 1.2) +
  geom_line( data =  expset_A_reg[[3]]$results_table, mapping = aes(y = bias, color = "Large_R"), size = 1.2) +
   scale_x_log10() + 
   theme_cowplot() +
   guides(color = guide_legend(title = "Experimental set"))
 
var_plot <- ggplot( expset_A[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_A[[1]]$results_table, mapping = aes(y = var, color = "Zero"), size = 1.2) +
  geom_line( data =  expset_A[[2]]$results_table, mapping = aes(y = var, color = "One"), size = 1.2) +
  #geom_line( data =  expset_A[[3]]$results_table, mapping = aes(y = var, color = "Large"), size = 1.2) +
  geom_line( data =  expset_A_reg[[1]]$results_table, mapping = aes(y = var, color = "Zero_R"), size = 1.2) +
  geom_line( data =  expset_A_reg[[2]]$results_table, mapping = aes(y = var, color = "One_R"), size = 1.2) +
  geom_line( data =  expset_A_reg[[3]]$results_table, mapping = aes(y = var, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))


mse_plot <- ggplot( expset_A[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_A[[1]]$results_table, mapping = aes(y = mse, color = "Zero"), size = 1.2) +
  geom_line( data =  expset_A[[2]]$results_table, mapping = aes(y = mse, color = "One"), size = 1.2) +
  geom_line( data =  expset_A[[3]]$results_table, mapping = aes(y = mse, color = "Large"), size = 1.2) +
  geom_line( data =  expset_A_reg[[1]]$results_table, mapping = aes(y = mse, color = "Zero_R"), size = 1.2) +
  geom_line( data =  expset_A_reg[[2]]$results_table, mapping = aes(y = mse, color = "One_R"), size = 1.2) +
  geom_line( data =  expset_A_reg[[3]]$results_table, mapping = aes(y = mse, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))


#############################################
####### EXPERIMENT SETUP B (Prev. C)    ##### 
#
# Eigenvalues range is small but with different centers (10,1,0.1)
# Lambda matrix is not 0
#
#############################################

# CONFIG
nrep = 100
nmax = 5000
nmin = 10
nfact = 2
p = 3
ncores = 6


cl <- parallel::makeForkCluster(nnodes = ncores)

# Random orthogonal matrix
ortho <- pracma::rortho(p)

# Kappa with close-to-0 inf. norm
diag_0 = abs(rnorm(p,0,0.1)) + 1E-6

# Kappa with close-to-1 inf. norm
diag_1 = rnorm(p,0,0.1) + 1

# Kappa with greater than 1 inf. norm
diag_g = rnorm(p,0,1) + 10

## B.1
# Generate P matrix
P <- ortho %*% diag(diag_0) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_B[[1]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

expset_B_reg[[1]] <- mseExp(kappa, lambda,
                            phi = P,
                            H = matrix(1, ncol = p, nrow = p),
                            nrep, nmin, nmax,
                            nfact, cl)

## B.2
P <- ortho %*% diag(diag_1) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
lambda <- (lambda + t(lambda)) /2

expset_B[[2]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)
expset_B_reg[[2]] <- mseExp(kappa, lambda,
                            phi = P,
                            H = matrix(1, ncol = p, nrow = p),
                            nrep, nmin, nmax,
                            nfact, cl)

## B.3
P <- ortho %*% diag(diag_g) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
lambda <- (lambda + t(lambda)) /2

expset_B[[3]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)
expset_B_reg[[3]] <- mseExp(kappa, lambda,
                            phi = P,
                            H = matrix(1, ncol = p, nrow = p),
                            nrep, nmin, nmax,
                            nfact, cl)

parallel::stopCluster(cl)

################################

#########
## Plots
#########

bias_plot <- ggplot( expset_B[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_B[[1]]$results_table, mapping = aes(y = bias, color = "Zero"), size = 1.2) +
  geom_line( data =  expset_B[[2]]$results_table, mapping = aes(y = bias, color = "One"), size = 1.2) +
  #geom_line( data =  expset_B[[3]]$results_table, mapping = aes(y = bias, color = "Large"), size = 1.2) +
  geom_line( data =  expset_B_reg[[1]]$results_table, mapping = aes(y = bias, color = "Zero_R"), size = 1.2) +
  geom_line( data =  expset_B_reg[[2]]$results_table, mapping = aes(y = bias, color = "One_R"), size = 1.2) +
  geom_line( data =  expset_B_reg[[3]]$results_table, mapping = aes(y = bias, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

var_plot <- ggplot( expset_B[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_B[[1]]$results_table, mapping = aes(y = var, color = "Zero"), size = 1.2) +
  geom_line( data =  expset_B[[2]]$results_table, mapping = aes(y = var, color = "One"), size = 1.2) +
  #geom_line( data =  expset_B[[3]]$results_table, mapping = aes(y = var, color = "Large"), size = 1.2) +
  geom_line( data =  expset_B_reg[[1]]$results_table, mapping = aes(y = var, color = "Zero_R"), size = 1.2) +
  geom_line( data =  expset_B_reg[[2]]$results_table, mapping = aes(y = var, color = "One_R"), size = 1.2) +
  geom_line( data =  expset_B_reg[[3]]$results_table, mapping = aes(y = var, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))


mse_plot <- ggplot( expset_B[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_B[[1]]$results_table, mapping = aes(y = mse, color = "Zero"), size = 1.2) +
  geom_line( data =  expset_B[[2]]$results_table, mapping = aes(y = mse, color = "One"), size = 1.2) +
  geom_line( data =  expset_B[[3]]$results_table, mapping = aes(y = mse, color = "Large"), size = 1.2) +
  geom_line( data =  expset_B_reg[[1]]$results_table, mapping = aes(y = mse, color = "Zero_R"), size = 1.2) +
  geom_line( data =  expset_B_reg[[2]]$results_table, mapping = aes(y = mse, color = "One_R"), size = 1.2) +
  geom_line( data =  expset_B_reg[[3]]$results_table, mapping = aes(y = mse, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

###########################


#############################################
####### EXPERIMENT SETUP C              #####
#
# eigenvalues range varies from [10,10] to [0.1,10] with underlying uniform distribution
# lambda is not 0
#
#############################################

# CONFIG
nrep = 1000
nmax = 1000
nmin = 10
nfact = 2
p = 3
ncores = 6


cl <- parallel::makeForkCluster(nnodes = ncores)

diagonal_expc <- function(p, emax, emin){
  return( c(emin, runif(p - 2)*(emax - emin) + emin ,emax))
}

# Random orthogonal matrix
ortho <- pracma::rortho(p)

# Cases
diag_1 <- diagonal_expc(p, 10, 10)
diag_2 <- diagonal_expc(p, 10, 5)
diag_3 <- diagonal_expc(p, 10, 1)
diag_4 <- diagonal_expc(p, 10, 0.1)


## C.1
# Generate P matrix
P <- ortho %*% diag(diag_1) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C[[1]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

# expset_C_reg[[1]] <- mseExp(kappa, lambda,
#                             phi = P,
#                             H = matrix(1, ncol = p, nrow = p),
#                             nrep, nmin, nmax,
#                             nfact, cl)

## C.2
# Generate P matrix
P <- ortho %*% diag(diag_2) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C[[2]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

# expset_C_reg[[2]] <- mseExp(kappa, lambda,
#                             phi = P,
#                             H = matrix(1, ncol = p, nrow = p),
#                             nrep, nmin, nmax,
#                             nfact, cl)

## C.3
# Generate P matrix
P <- ortho %*% diag(diag_3) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C[[3]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

# expset_C_reg[[3]] <- mseExp(kappa, lambda,
#                             phi = P,
#                             H = matrix(1, ncol = p, nrow = p),
#                             nrep, nmin, nmax,
#                             nfact, cl)

## C.4
# Generate P matrix
P <- ortho %*% diag(diag_4) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C[[4]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

# expset_C_reg[[4]] <- mseExp(kappa, lambda,
#                             phi = P,
#                             H = matrix(1, ncol = p, nrow = p),
#                             nrep, nmin, nmax,
#                             nfact, cl)



parallel::stopCluster(cl)

################################

#########
## Plots
#########

bias_plot <- ggplot( expset_C[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_C[[1]]$results_table, mapping = aes(y = bias, color = "[10,10]"), size = 1.2) +
  geom_line( data =  expset_C[[2]]$results_table, mapping = aes(y = bias, color = "[5,10]"), size = 1.2) +
  geom_line( data =  expset_C[[3]]$results_table, mapping = aes(y = bias, color = "[1,10]"), size = 1.2) +
  geom_line( data =  expset_C[[4]]$results_table, mapping = aes(y = bias, color = "[0.1,10]"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[1]]$results_table, mapping = aes(y = bias, color = "Zero_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[2]]$results_table, mapping = aes(y = bias, color = "One_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[3]]$results_table, mapping = aes(y = bias, color = "Large_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[4]]$results_table, mapping = aes(y = bias, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

var_plot <- ggplot( expset_C[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_C[[1]]$results_table, mapping = aes(y = var, color = "[10,10]"), size = 1.2) +
  geom_line( data =  expset_C[[2]]$results_table, mapping = aes(y = var, color = "[5,10]"), size = 1.2) +
  geom_line( data =  expset_C[[3]]$results_table, mapping = aes(y = var, color = "[1,10]"), size = 1.2) +
  geom_line( data =  expset_C[[4]]$results_table, mapping = aes(y = var, color = "[0.1,10]"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[1]]$results_table, mapping = aes(y = var, color = "Zero_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[2]]$results_table, mapping = aes(y = var, color = "One_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[3]]$results_table, mapping = aes(y = var, color = "Large_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[4]]$results_table, mapping = aes(y = var, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))


mse_plot <- ggplot( expset_C[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_C[[1]]$results_table, mapping = aes(y = mse, color = "[10,10]"), size = 1.2) +
  geom_line( data =  expset_C[[2]]$results_table, mapping = aes(y = mse, color = "[5,10]"), size = 1.2) +
  geom_line( data =  expset_C[[3]]$results_table, mapping = aes(y = mse, color = "[1,10]"), size = 1.2) +
  geom_line( data =  expset_C[[4]]$results_table, mapping = aes(y = mse, color = "[0.1,10]"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[1]]$results_table, mapping = aes(y = mse, color = "Zero_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[2]]$results_table, mapping = aes(y = mse, color = "One_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[3]]$results_table, mapping = aes(y = mse, color = "Large_R"), size = 1.2) +
  #geom_line( data =  expset_C_reg[[4]]$results_table, mapping = aes(y = mse, color = "Large_R"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

###########################

#############################################
####### EXPERIMENT SETUP C - BIS        #####
#
# Compare three opts (dist with:)
#  - Same max, same sum
#  - Dist max, same sum
#
# We need three dists
# A and B with same Max and same sum
# C with different max but same sum
#
# lambda is not 0
#
#############################################


expset_C_bis <- list()

# CONFIG
nrep = 1000
nmax = 1000
nmin = 10
nfact = 2
p = 4 # First is not valid with p=3 -> no choice
ncores = 6

cl <- parallel::makeForkCluster(nnodes = ncores)

# Random orthogonal matrix
ortho <- pracma::rortho(p)

# MAX IS 10 SUM IS 20
# Cases
diag_1 <- c(10, 3, 3, 4) # SUM 20 MAX 10
diag_2 <- c(10, 4.5, 4.5, 1) # SUM 20 MAX 10
diag_3 <- c(5, 5, 5, 5) # SUM 20 MAX 15
diag_4 <- c(19, 0.4, 0.3, 0.3) # SUM 20 MAX 15
diag_5 <- c(10, 0.5, 0.5, 9) # SUM 20 MAX 10

## C.1
# Generate P matrix
P <- ortho %*% diag(diag_1) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C_bis[[1]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

## C.2
# Generate P matrix
P <- ortho %*% diag(diag_2) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C_bis[[2]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

## C.3
# Generate P matrix
P <- ortho %*% diag(diag_3) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C_bis[[3]] <- mseExp(kappa, lambda, NULL, NULL,
                        nrep, nmin, nmax,
                        nfact, cl)

## C.4
# Generate P matrix
P <- ortho %*% diag(diag_4) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C_bis[[4]] <- mseExp(kappa, lambda, NULL, NULL,
                            nrep, nmin, nmax,
                            nfact, cl)

## C.5
# Generate P matrix
P <- ortho %*% diag(diag_5) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_C_bis[[5]] <- mseExp(kappa, lambda, NULL, NULL,
                            nrep, nmin, nmax,
                            nfact, cl)


#########
## Plots
#########

bias_plot <- ggplot( expset_C_bis[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_C_bis[[1]]$results_table, mapping = aes(y = bias, color = "A(10,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[2]]$results_table, mapping = aes(y = bias, color = "B(10,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[3]]$results_table, mapping = aes(y = bias, color = "C(05,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[4]]$results_table, mapping = aes(y = bias, color = "D(19,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[5]]$results_table, mapping = aes(y = bias, color = "E(10,20)"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

var_plot <- ggplot( expset_C_bis[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_C_bis[[1]]$results_table, mapping = aes(y = var, color = "A(10,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[2]]$results_table, mapping = aes(y = var, color = "B(10,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[3]]$results_table, mapping = aes(y = var, color = "C(05,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[4]]$results_table, mapping = aes(y = var, color = "D(19,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[5]]$results_table, mapping = aes(y = var, color = "E(10,20)"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))


mse_plot <- ggplot( expset_C_bis[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_C_bis[[1]]$results_table, mapping = aes(y = mse, color = "A(10,20)"), size = 1.2) +
  geom_line( data =  expset_C_bis[[2]]$results_table, mapping = aes(y = mse, color = "B(10,20)"), size = 1.2) +
  #geom_line( data =  expset_C_bis[[3]]$results_table, mapping = aes(y = mse, color = "C(05,20)"), size = 1.2) +
  #geom_line( data =  expset_C_bis[[4]]$results_table, mapping = aes(y = mse, color = "D(19,20)"), size = 1.2) +
  #geom_line( data =  expset_C_bis[[5]]$results_table, mapping = aes(y = mse, color = "E(10,20)"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

parallel::stopCluster(cl)

#############################################
####### EXPERIMENT SETUP D        #####
#
# Compare three params steups with same eigenvals 2-norm
#
#
# lambda is not 0
#
#############################################
expset_D <- list()

# CONFIG
nrep = 1000
nmax = 1000
nmin = 10
nfact = 2
p = 3
ncores = 6
norm2 = 10

cl <- parallel::makeForkCluster(nnodes = ncores)

# Random orthogonal matrix
ortho <- pracma::rortho(p)

# Cases (3 different vectors with same 2-nrom)
theta <- runif(1) * 2 * pi
phi <- runif(1) * pi
diag_1 <- abs(c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))*norm2)
theta <- runif(1) * 2 * pi
phi <- runif(1) * pi
diag_2 <- abs(c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))*norm2)
theta <- runif(1) * 2 * pi
phi <- runif(1) * pi
diag_3 <- abs(c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))*norm2)

## D.1
# Generate P matrix
P <- ortho %*% diag(diag_1) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_D[[1]] <- mseExp(kappa, lambda, NULL, NULL,
                            nrep, nmin, nmax,
                            nfact, cl)

## D.2
# Generate P matrix
P <- ortho %*% diag(diag_2) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_D[[2]] <- mseExp(kappa, lambda, NULL, NULL,
                            nrep, nmin, nmax,
                            nfact, cl)

## D.3
# Generate P matrix
P <- ortho %*% diag(diag_3) %*% t(ortho)
kappa <- diag(P)
lambda <- diag(kappa) - P
# Sometimes lambda is not symmetric (due to approx. errors)
lambda <- (lambda + t(lambda)) /2

expset_D[[3]] <- mseExp(kappa, lambda, NULL, NULL,
                            nrep, nmin, nmax,
                            nfact, cl)


#########
## Plots
#########

bias_plot <- ggplot( expset_D[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_D[[1]]$results_table, mapping = aes(y = bias, color = "A"), size = 1.2) +
  geom_line( data =  expset_D[[2]]$results_table, mapping = aes(y = bias, color = "B"), size = 1.2) +
  geom_line( data =  expset_D[[3]]$results_table, mapping = aes(y = bias, color = "C"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

var_plot <- ggplot( expset_D[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_D[[1]]$results_table, mapping = aes(y = var, color = "A"), size = 1.2) +
  geom_line( data =  expset_D[[2]]$results_table, mapping = aes(y = var, color = "B"), size = 1.2) +
  geom_line( data =  expset_D[[3]]$results_table, mapping = aes(y = var, color = "C"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))


mse_plot <- ggplot( expset_D[[1]]$results_table, aes(n) ) + 
  geom_line( data =  expset_D[[1]]$results_table, mapping = aes(y = mse, color = "A"), size = 1.2) +
  geom_line( data =  expset_D[[2]]$results_table, mapping = aes(y = mse, color = "B"), size = 1.2) +
  geom_line( data =  expset_D[[3]]$results_table, mapping = aes(y = mse, color = "C"), size = 1.2) +
  scale_x_log10() + 
  theme_cowplot() +
  guides(color = guide_legend(title = "Experimental set"))

parallel::stopCluster(cl)
