#############################################################
#                                                           #
#   rayleigh.qsample.test                                  #
#                                                           #
#############################################################
#' @export
rayleigh.qsample.test <- function(..., alpha = 0.05) {
  y <- list(...)
  # Unify data in x. each samples as a list
  x <- list()
  for (i in 1:length(y)) {
    if (is.data.frame(y[[i]])) {
      x <- c(x, as.list(y[[i]]))
    } else if (is.matrix(y[[i]])) {
      for (j in 1:ncol(y[[i]])) {
        x <- c(x, list(y[[i]][,j]))
      }
    } else if (is.list(y[[i]])) {
      x <- c(x, y[[i]])
    } else {
      x <- c(x, list(y[[i]]))
    }
  }
  if (length(x) < 2)
    stop("There must be at least two samples")
  for (i in 1:length(x)) {
    x[[i]] <- conversion.circular(x[[i]], units = "radians", zero = 0, 
                                  rotation = "counter", modulo = "2pi")
    attr(x[[i]], "circularp") <- attr(x[[i]], "class") <- NULL                           
  } 
  
  # Handling missing values
  x <- lapply(x, na.omit)                                                                 
  
  result <- RayqsampleTestRad(x)                                                          
  result$call <- match.call() 
  result$alpha <- alpha
  class(result) <- "rayleigh.qsample.test"
  return(result)
}

RayqsampleTestRad <- function(x) {                                                          
  # Get each sample length and number of samples                                          
  n <- vapply(x, length, FUN.VALUE = numeric(1))                                          
  q <- length(x)
  
  # Add sample number to every component                                                  
  x <- lapply(1:q, function(i) return( cbind(x[[i]],i)))                                  
  
  # All data together + ranks
  xx <- Reduce(rbind, x)                                                                  
  rank <- order(xx[,1])
  
  # Angles
  beta  <-  (rank * 2 * pi)/sum(n)
  
  #CI and si
  ci <- aggregate( beta, by = list(sample = xx[,2]), function(x) sum(cos(x)) ^ 2 )
  si <- aggregate( beta, by = list(sample = xx[,2]), function(x) sum(sin(x)) ^ 2 )
  ri <- ci[,2] + si[,2]
  df <- 2*(q - 1)
  W <- 2*sum(ri/n)
  
  result <- list()
  result$statistic <- W
  result$df <- df
  result$p.value <- (1 - pchisq(W, df))
  return(result)
}

#' @export
print.rayleigh.qsample.test <- function(x, digits=4, ...) {
  rbar <- x$statistic
  p.value <- x$p.value
  alpha <- x$alpha
  mu <- x$mu
  cat("\n", "      Rayleigh q-sample Test \n")
  cat("Test Statistic: ", round(rbar, digits = digits), "\n")
  cat("P-value: ", round(p.value, digits = digits), "\n")
  if( p.value < alpha ) cat(sprintf(" Null hypothesis F1 = F2 = ... Fq -- REJECTED -- (significance %.4f) \n\n", alpha ))
  else cat(sprintf(" Null hypothesis F1 = F2 = ... Fq -- NOT REJECTED -- (significance %.4f) \n\n", alpha ))
  invisible(x)
}


