#' Fits a multivariate von Mises  distribution
#'
#' Given a set of multivariate circular samples, compute ml mu, kappa and lambda parameters
#' based on the pseudolikelihood function. Descent is made using lbfgs-b algorithm.
#' 
#' @param samples Matrix or df of circular samples
#' @param mu Initial value for mean vector. It should be Circular.
#' @param kappa Initial value for concentration vector
#' @param lambda Initial value for lambda matrix
#' @param phi Prior parameter matrix. P matrix is Diag(kappa) - LAMBDA
#' @param H confidence matrixs
#' @param verbose Verbose degree parameter. -1 for silent mode.
#' @param prec Precision parameter of the LBFGS-B method
#' @param tolerance Tolerance parameter of the LBFGS-B mehtod
#' @param mprec Number of steps that the pseudoinverse calculation remembers
#' @param retry Number of retries (Convergence can fail sometimes)
#' @param type Wheter to use optimized expresion or regular
#'
#' @return List with mu, kappa and lambda parameters, data samples and loss function value. If method don't reach convergence loss will be NA
#' 
#' @useDynLib mvCircular
#'
#' @examples 
#' fit_mvvonmises(rmvVonMises(100,rep(0,3),rep(1,3),matrix(0,nrow=3,ncol=3)))
#' 
#' @export
fit_mvvonmises <- function(samples,
                           mu = rep(0,ncol(samples)),
                           kappa = rep(1,ncol(samples)),
                           lambda = matrix(0,ncol = ncol(samples), nrow = ncol(samples)),
                           phi = NULL,
                           H = NULL,
                           verbose = -1,
                           prec = 1E5,
                           tolerance = 1E-7,
                           mprec = 30,
                           retry = 10){
  

  
  # VALIDATE INPUTS
  
  if (!is.mvCircular(samples))
    stop("Samples should be multivariate circular")
  
  # Get sample size
  p <- ncol(samples)
  n <- nrow(samples)

  # Check params
  if (length(mu) != p || length(kappa) != p || nrow(lambda) != p || ncol(lambda) != p)
    stop("Initial parameters length do not match")
  
  if (any(H < 0) )
    stop("Confidence should be positive")
  
  # Ensure that phi is a matrix and that it is symmetric
  if (!is.null(phi)) {
    phi <- as.matrix(phi)
    phi <- (phi + t(phi))/2
    
    if (is.null(H))
      H <- matrix(1, nrow = p, ncol = p )
    
    phi <- t(phi)
    withPenalization <- 1
  }
  else
    withPenalization <- 0
  
  # Ensure that H is upper triangular
  if (!is.null(H)) {
    H[lower.tri(H)] <- 0
    H <- t(H) # C matrices are row-ordered
  }
    
  if (!is.numeric(verbose))
    stop("Verbose should be numeric")
  
  if (!is.numeric(prec) || prec < 0 )
    stop("Preccision should be a strictly positive number")
  
  if (!is.numeric(tolerance) || tolerance < 0 )
    stop("Tolerance should be a strictly positive number ")
  
  if (!is.numeric(mprec) || floor(mprec) != mprec || mprec < 0 )
    stop("Matrix inversion steps should be a strictly positive integer")
  
  if (!is.numeric(retry) || floor(retry) != retry || retry < 0 )
    stop("Number of retrials should be a strictly positive integer")
  
  #### END VALIDATION ### 
  
  # Number of lambda variables
  lambda.nvars <- ((p*p) - p)/2
  kappa.nvars <- p
  
  # Set bounds 
  lower <- rep(0,kappa.nvars + lambda.nvars )
  upper <- double(kappa.nvars + lambda.nvars)
  bounded <- rep(0,kappa.nvars + lambda.nvars )
  
  # Only kappa is lower bounded
  bounded[1:kappa.nvars] < -1 
  
  # Call procedure
  i <- 0;
  fail <- T
  while (fail && i < retry) {
    
    # Reset failure flag
    fail <- F
    
    # Optimized procedure
    ret <- .C("__R_mvvonmises_lbfgs_fit",
              p = as.integer(p),
              mu = as.double(mu),
              kappa = as.double(kappa),
              lambda = as.double(lambda),
              n = as.integer(n),
              samples = as.double(t(samples)),
              penalized = as.integer(withPenalization),
              phi = as.double(phi),
              H = as.double(H),
              verbose = as.integer(verbose),
              prec = as.double(prec),
              tol = as.double(tolerance),
              mprec = as.integer(mprec),
              lower = as.double(lower),
              upper = as.double(upper),
              bounded = as.integer(bounded),
              loss = double(1),
              PACKAGE = "mvCircular")
    
      # Check for failure
      if ( is.nan(ret$loss) ) {
        fail <- T
        warning(sprintf(" Failure at iteration %d",i))
      }
      i <- i + 1
  }
  
  # Last exec was successfull
  if ( !fail )
    # Return parameters
    return(list(samples = samples,
              mu = ret$mu,
              kappa = ret$kappa,
              lambda = matrix(ret$lambda,ncol = p,nrow = p),
              loss = ret$loss))
  else
    # Return call
    return(list( samples = samples,
                 mu = mu,
                 kappa = kappa,
                 lambda = lambda,
                 loss = NA ))
}
  
  
  
  