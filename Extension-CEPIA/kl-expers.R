library(mvCircular)
library(ggplot2)
library(cowplot)
library(prallel)

# CONFIG
nrep = 100
nmax = 5000
nmin = 10
nfact = 2
k = 2
p = 10

cl <- parallel::makeCluster(6, type = "FORK")

##
# P=3
##
x <- mvVonMises( rep(0,p), kappa = rep(1,p), matrix(0, ncol = p, nrow = p))
x <- normalize(x)
y <- x
z.mu <- rep(0,p)
z.mu[ sample.int(p,floor(p/3)) ] <- pi/2
z.kappa = rep(1,p)
z.kappa[ sample.int(p,floor(p/3)) ] <- 3
z.kappa[ sample.int(p,floor(p/3)) ] <- 0.5
z <- mvVonMises( z.mu, kappa = z.kappa, matrix(0, ncol = p, nrow = p))
z <- normalize(z)

##
# Generate samples
##
x.samples <- getSamples(x, nrep*nmax)
y.samples <- getSamples(y, nrep*nmax)
z.samples <- getSamples(z, nrep*nmax)

# Generate N values
nlist <- nmin * nfact ^ (0:floor(log(nmax,nfact) - log(nmin,nfact)))

# Add nmax if it is not on the list
if (nlist[length(nlist)] != nmax) nlist <- c(nlist, nmax)

# EXP 1: X -Y
ret <- parallel::parLapplyLB(cl,nlist, function(n, nrep,A,B,k){
  
 # Emp. KL divergence
 vals <- sapply(1:nrep, function(x, data_1, data_2, size){
   return(empKL.circular(data_1[sample(nrow(data_1),size = n),],data_2[sample(nrow(data_2),size=n),])$kl)
 },A,B,n)
 
 # Return mean and variance
 return(list(n = n, mean = mean(vals), var = var(vals) ))
},nrep = nrep,x.samples,y.samples, k = k)

ret <- as.data.frame(t(sapply(ret,as.numeric)))
colnames(ret) <- c("n","mean","var")

ret$lwr <- ret$mean - ret$var
ret$upr <- ret$mean + ret$var

# Plot
pl <- ggplot( ret, aes(n,mean) ) + 
  geom_line( colour = "darkred") +
  geom_ribbon(aes(ymin = lwr,ymax = upr), alpha = 0.3, fill = "red")+ 
  scale_x_log10() + 
  scale_y_log10() + theme_cowplot()

save_plot("kl_p10_k2_div0.png",pl)

realkl <- kl(x,z, N = 1E5)
ret2 <-  parallel::parLapplyLB(cl,nlist, function(n, nrep,A,B,k){
  
  # Emp. KL divergence
  vals <- sapply(1:nrep, function(x, data_1, data_2, size){
    return(empKL.circular(data_1[sample(nrow(data_1),size = n),],data_2[sample(nrow(data_2),size = n),])$kl)
  },A,B,n)
  
  # Return mean and variance
  return(list(n = n, mean = mean(vals), var = var(vals) ))
},nrep = nrep,x.samples,z.samples, k = k)

ret2 <- as.data.frame(t(sapply(ret2,as.numeric)))
colnames(ret2) <- c("n","mean","var")

ret2$lwr <- ret2$mean - ret2$var
ret2$upr <- ret2$mean + ret2$var
ret2$real <- realkl

# Plot
pl <- ggplot( ret2, aes(n,mean) ) + 
  geom_line( colour = "darkred") +
  geom_ribbon(aes(ymin = lwr,ymax = upr), alpha = 0.3, fill = "red")+ 
  scale_x_log10() + 
  geom_hline(yintercept = realkl) +
  theme_cowplot()

save_plot("kl_p10_k2_divnot0.png",pl)

parallel::stopCluster(cl)
