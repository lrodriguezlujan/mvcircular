#'
#' #@param par parameter vector (length p + p(p-1)/2 ) in kappa - upper.tri(lambda) order
#'
#'@nord
mvVonMises.logplik.pen <- function(par, phi = NULL, H = NULL, alpha = NULL, data.s, data.c ){

  # data.s and data.c are assumed to be already centered
  # coerce to matrices
  data.s <- as.matrix(data.s)
  data.c <- as.matrix(data.c)
  
  # data.s is a Nxp matrix (so it is data.c)
  p <- ncol(data.s)
  N <- nrow(data.s)
  
  # Rebuild kappa and lambda
  kappa <- par[1:p]
  lambda <- matrix(0, ncol = p, nrow = p)
  lambda[upper.tri(lambda)] <- par[(p + 1):length(par)]
  lambda <- lambda + t(lambda)
  
  # Now compute the (normalized) log pseudo likelihood

  # Fix term
  fixed <- -p * log(2*pi) 
  
  # compute sin lambda product
  ro <- Matrix::Matrix(data.s) %*% Matrix::Matrix(lambda)
  
  # For each instance
  pll <- sapply(1:N, function(i){
    
    # Conditional kappa value
    kappa.conditional <- sqrt(kappa ^ 2 + ro[i,] ^ 2)
    
    #Expon scaled bessel function (exp(-x)I_0) to avoid overflow
    # Then, we need to add +x to log() to have log(I_0)
    # "normalization" side of the log pseudo likelihood
    return( -(sum(log(besselI(kappa.conditional, 0 , expon.scaled = T)) + kappa.conditional) +
    # the "cosine" (linear) part (cos(theta) * k_j)
    as.numeric(crossprod( data.c[i,], kappa)) +
    # thee "sine" (quadratic) part
    as.numeric(crossprod( data.s[i,], ro[i,])))/N )
  })
  
  if ( is.null(phi) || is.null(alpha) || is.null(H) ) {
    penalization <- 0
  }
  else {
    # Force H to be upper triangular
    H[ lower.tri(H) ] <- 0
    # F norm of ( P - Phi) H
    # Alpha is the diminishing factor / function
    if ( is.numeric(alpha) ) 
      penalization <- (N ^ (-abs(alpha)) ) * norm(((diag(kappa) - lambda) - phi) * H, "F")
    else if (is.function(alpha))
      penalization <- alpha(N) * norm(((diag(kappa) - lambda) - phi) * H, "F")
    else 
      penalization <- 0
  }
  
  return(fixed + sum(pll) - penalization)
}

#'
#'
#'@nord
mvVonMises.logplik.pen.gr <- function(par, phi = NULL, H = NULL, alpha = NULL, data.s, data.c ){
  
  # data.s is a Nxp matrix (so it is data.c)
  p <- ncol(data.s)
  N <- nrow(data.s)
  
  # coerce to matrices
  data.s <- as.matrix(data.s)
  data.c <- as.matrix(data.c)
  
  # Rebuild kappa and lambda
  kappa <- par[1:p]
  lambda <- matrix(0, ncol = p, nrow = p)
  lambda[upper.tri(lambda)] <- par[(p + 1):length(par)]
  lambda <- lambda + t(lambda)
  
  # compute sin lambda product
  ro <- Matrix::Matrix(data.s) %*% Matrix::Matrix(lambda)
  
  # For each instance
  pll <- sapply(1:N, function(i){
    
    # Conditional kappa value
    kappa.conditional <- sqrt(kappa ^ 2 + ro[i,] ^ 2)
    
    # Auxiliar: A_1/k_cond
    a1.auxiliar <- (besselI(kappa.conditional, 1, expon.scaled = T)/besselI(kappa.conditional, 0, expon.scaled = T))/ kappa.conditional 
    
    # Kappa gradient
    kappa.gr <-  data.c[i,] - a1.auxiliar * kappa
    
    # Lambda gradient
    lambda.gr <- tcrossprod( data.s[i, ] - a1.auxiliar *  ro[i, ], data.s[i,])
    
    # lambda_jk = lambda_kj
    lambda.gr <- lambda.gr + t(lambda.gr)
    
    # Return kappa + upper tri lambda
    return( c(kappa.gr, lambda.gr[ upper.tri(lambda.gr) ])/N )
  })
  
  # Penalization
  if ( is.null(phi) || is.null(alpha) || is.null(H) ) {
    return( as.numeric(rowSums(pll)) ) # No penalization
  }
  else {
    
    # Force H to be upper triangular
    H[ lower.tri(H) ] <- 0
    # F norm of ( P - Phi) H
    # Alpha is the diminishing factor / function
    if ( is.numeric(alpha) ) 
      penalization <- (N ^ (-abs(alpha)) ) * (1 / norm(((diag(kappa) - lambda) - phi) * H, "F"))
    else if (is.function(alpha))
      penalization <- alpha(N) * (1 / norm(((diag(kappa) - lambda) - phi) * H, "F"))
    else 
      return( as.numeric(rowSums(pll)) ) # No penalization
    
    # For each element 
    penalization <- penalization * (1 * (H ^ 2) * (((diag(kappa) - lambda)) - phi ))
    penalization <- c(diag(penalization), penalization[upper.tri(penalization)])
    
    return( as.numeric(rowSums(pll)) - penalization)
  }
  
}

#test <- optim(rep(0.5,6), mvCircular:::mvVonMises.logplik.pen , mvCircular:::mvVonMises.logplik.pen.gr,
      #data.s = sin(x.samples), data.c = cos(x.samples), method = "L-BFGS-B", lower = c(0,0,0,-Inf, -Inf, -Inf) )
