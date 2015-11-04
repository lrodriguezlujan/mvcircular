
#' Von mises multivariate sampler based on the Gibbs sampler with multiple chains, thinning and burnin options
#'
#' Samples from a multivariate von mises distribution
#'
#'@export
#'@useDynLib mvCircular
#'
#' @importFrom bit64 as.integer64
#'
#'@examples
#'rmvvonmises_rs(100,rep(0,3),rep(1,3),matrix(rep(0,9),nrow=3,ncol=3))
rmvvonmises_rs <- function(n,
                           mu,kappa,lambda,
                           burnin=0.5, # Burn-in factor
                           burnin.tol=0.1,
                           burnin.chains=chains,
                           thinning=0.01, # Thinning parameter
                           chains=5, # Number of chaines
                           alpha=rep(1/length(mu),length(mu)), # Selection probabilities
                           alphaPrec = 4, # Alpha lookup table precision
                           seed=runif(64,0,.Machine$integer.max),
                           finalSize=T,
                           verbose=F,
                           type = c("auto","gibbs","rejection")){
  
  ## Check type
  type <- match.arg(type)
  
  ## Check values
  p <- length(mu);
  
  # Check params
  if(p==0)
    stop("Empty parameters provided")
  else if(p!=length(kappa) || p!=nrow(lambda))
    stop("Parameter length do not match")
  else if(nrow(lambda)!=ncol(lambda))
    stop("Lambda is not a square matrix")
  else if( any(t(lambda) != lambda) )
    warning("Lambda is not symmetric")
  
  # Check for rejection (auto with low dim or selected)
  if( (type=="auto" && p<=25) || type=="rejection"){
    
    # Check that lambda is positive semidefinite
    aux <- -lambda
    diag(aux) <- kappa
    eigVals <- eigen(aux,symmetric = T,only.values = T)$values
    lmin <- min(eigVals)
    
    ## Compute approx acceptance ratio
    if(lmin>0)
      accRatio = (1/(2^p)) * (sqrt((lmin^p) / prod(eigVals)))
    else
      accRatio=0
    
    if( lmin>0 && (accRatio>1E-4) ){
      #cat(sprintf("START Sampling REJECTION, maxKappa: %f , sum(ABS)Lambda: %f\n",max(kappa),sum(abs(lambda))))
      ret <- .C("__R_mvvm_sampler_safe_rejection",
                n=as.integer(n),
                p=as.integer(p),
                mu = as.double(mu),
                kappa = as.double(kappa),
                lambda = as.double(lambda),
                lmin = lmin,
                seeds = bit64::as.integer64(seed),
                ret = double(n*p),
                errno = integer(1),
                PACKAGE="mvCircular")
      #cat(sprintf("END Sampling REJECTION, maxKappa: %f , sum(ABS)Lambda: %f\n",max(kappa),sum(abs(lambda))))
      return(circular(matrix(ret$ret,ncol=p,byrow=T)))
    }
    #else{
    #  # If type was rejection create same sample using gibbs
    #  if(type=="rejection")
    #    warning("Matrix is not definite positive. Using Gibbs sampler instead")
    #}
  }
  # Use gibbs sampler
  
  # Compute alpha array
  aux <- round(alpha,digits = alphaPrec)*(10^alphaPrec)
  alphaArray <- rep( seq(1,p)-1,times=aux)
  
  # Compute thinning value
  if(thinning<=0)
    thinningVal <- 1
  else if(thinning>=1)
    thinningVal <- thinning
  else
    # Añadimos la suma de lambda y kappa para evitar mucha autocorrelacion hasta 10
    thinningVal <- ceiling(log(thinning/p)/log(1-min(alpha)))*min(5,ceiling(1+log(1+sum(abs(lambda))+sum(kappa))))
  
  # Compute initial seeds as a list of seeds
  theta <- split( runif(p*chains,min = 0,max=(2*pi) ), rep(1:chains, each = p) )
  
  # Perform fixed burn-in if required
  if(burnin > 0 ){
    
    # First case: Burn in is a ratio of total
    if( burnin <1 )
      # if FinalSize is N we have to return a total of n samples per chain
      # so real "size" is N =  n / (1-burnin) -> total iterations to discard: n * burnin / (1-burnin)
      if(finalSize)
        burnin_iters <- (n * burnin) / (1-burnin)
      else
      {
        burnin_iters <- n * burnin
        n <- n*(1-burnin)
      }
    
    # Second case: Just a number is given == Number of iters to be discarded
    else
      burnin_iters <- burnin
    
    # Now call C procedure to actually do it (one call per chain)
    theta_n <- lapply(theta,function(x){
      # Execute burn-in (with thinning)
      #cat(sprintf("START Sampling BURNIN GIBBS, maxKappa: %f , sum(ABS)Lambda: %f\n",max(kappa),sum(abs(lambda))))
      ret <- .C("__R_mvvm_sampler_safe_rsgibbs",
                n=as.integer(burnin_iters),
                p=as.integer(p),
                mu = as.double(mu),
                kappa = as.double(kappa),
                lambda = as.double(lambda),
                prec = as.integer(length(alphaArray)),
                alphaTable = as.integer(alphaArray),
                thinning = as.integer(burnin_iters + 1),
                theta = as.double(x),
                seeds = bit64::as.integer64(seed),
                ret = double(p),
                errno = integer(1),
                PACKAGE="mvCircular")
      #cat(sprintf("END Sampling BURNIN GIBBS, maxKappa: %f , sum(ABS)Lambda: %f\n",max(kappa),sum(abs(lambda))))
      seed <- as.integer64(ret$seeds)
      return(ret$theta)
    })
    theta <- theta_n
  }
  # Auto burning until convergence
  else if(burnin<0){
    
      # Check number of chains
     if(burnin.chains < chains)
     {
       warning("Burnin chains should be greater or equal to general chains")
       burnin.chains <- chains
     }
     
     # Create burnin seeds
     burnin.theta <- split( runif(p*burnin.chains,min = 0,max=(2*pi) ), rep(1:burnin.chains, each = p) )
     convVal <- burnin.tol
     blocksize <- -burnin
     
     #Run burnin till convergece
     while(convVal>=burnin.tol){
  
       # TODO: This can be done much more efficient (reusing samples) but since
       # sampling is "fast"-ish is ok at the moment
       
       # Create samples
       burnin.samples <- lapply(burnin.theta,function(x){
         ret <- .C("__R_mvvm_sampler_safe_rsgibbs",
                   n=as.integer(blocksize*thinningVal),
                   p=as.integer(p),
                   mu = as.double(mu),
                   kappa = as.double(kappa),
                   lambda = as.double(lambda),
                   prec = as.integer(length(alphaArray)),
                   alphaTable = as.integer(alphaArray),
                   thinning = as.integer(thinningVal),
                   theta = as.double(x),
                   seeds = as.integer64(seed),
                   ret = double(p*blocksize),
                   errno = integer(1),
                   PACKAGE="mvCircular")
         seed <- as.integer64(ret$seeds)
         return(circular(matrix(ret$ret,ncol=p,byrow=T)))
      })

     # Create MCMC 
     gelman <- gelman.diag(mcmc.list(lapply(burnin.samples,mcmc)),autoburnin = F)
     convVal <- abs(gelman$mpsrf-1)
     
     # Increment block size
     blocksize <- blocksize - burnin
    }
    if(verbose)
      write(sprintf("Final Burnin size: %d\n",blocksize+burnin), stdout())

    # Burning completed. Copy n seeds
    # OJO! Que no está cambiando esto!!! burnin.theta es el original!
    theta <- burnin.theta[1:length(theta)]
  }
  
  # Create samples
  samples <- lapply(theta,function(x){
    #cat(sprintf("START Sampling GIBBS, maxKappa: %f , sum(ABS)Lambda: %f\n",max(kappa),sum(abs(lambda))))
    ret <- .C("__R_mvvm_sampler_safe_rsgibbs",
              n=as.integer(n*thinningVal),
              p=as.integer(p),
              mu = as.double(mu),
              kappa = as.double(kappa),
              lambda = as.double(lambda),
              prec = as.integer(length(alphaArray)),
              alphaTable = as.integer(alphaArray),
              thinning = as.integer(thinningVal),
              theta = as.double(x),
              seeds = as.integer64(seed),
              ret = double(p*n),
              errno = integer(1),
              PACKAGE="mvCircular")
    #cat(sprintf("END Sampling GIBBS, maxKappa: %f , sum(ABS)Lambda: %f\n",max(kappa),sum(abs(lambda))))
    seed <- bit64::as.integer64(ret$seeds)
    return(matrix(ret$ret,ncol=p,byrow=T))
  })
  
  # Return value
  return(samples)
}