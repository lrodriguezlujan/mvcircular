
#' Von mises multivariate sampler based on the Gibbs sampler
#'
#' Samples from a multivariate von mises distribution
#'
#'@export
#'@useDynLib mvCircular
#'
#'@examples 
#'rmvvonmises(10,rep(0,3),rep(1,3),matrix(0,nrow=3,ncol=3))
rmvvonmises <- function(n,mu,kappa,lambda,k_init=500,k_inter=100,seed=runif(64,0,.Machine$integer.max),type=c("auto","gibbs","rejection")){

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
  else if(!is.symmetric.matrix(lambda))
    warning("Lambda is not symmetric")
  
  # Call C procedure
  if(type=="auto" || type=="rejection"){
    aux <- -lambda
    diag(aux) <- kappa
    lmin <- min(eigen(aux,symmetric = T,only.values = T)$values)
    if( lmin<=0 || (p>10 && type!="rejection") ){
      if(type=="rejection")
        warning("Matrix is not definite positive. Using Gibbs sampler instead")
      ret <- .C("__R_mvvm_sampler_safe_gibbs",
                n=as.integer(n),
                p=as.integer(p),
                mu = as.double(mu),
                kappa = as.double(kappa),
                lambda = as.double(lambda),
                initk = as.integer(k_init),
                k = as.integer(k_inter),
                seeds = as.integer64(seed),
                ret = double(n*p),
                errno = integer(1),
                PACKAGE="mvCircular")  
    }
    else{
      ret <- .C("__R_mvvm_sampler_safe_rejection",
                n=as.integer(n),
                p=as.integer(p),
                mu = as.double(mu),
                kappa = as.double(kappa),
                lambda = as.double(lambda),
                lmin = lmin,
                seeds = as.integer64(seed),
                ret = double(n*p),
                errno = integer(1),
                PACKAGE="mvCircular")
    }
  }
  else{
    ret <- .C("__R_mvvm_sampler_safe_gibbs",
              n=as.integer(n),
              p=as.integer(p),
              mu = as.double(mu),
              kappa = as.double(kappa),
              lambda = as.double(lambda),
              initk = as.integer(k_init),
              k = as.integer(k_inter),
              seeds = as.integer64(seed),
              ret = double(n*p),
              errno = integer(1),
              PACKAGE="mvCircular")  
  }
    
  
  # Return value
  return(matrix(ret$ret%%(2*pi),ncol=p,nrow=n,byrow=T))
}
