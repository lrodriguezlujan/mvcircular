#' Own bessi0
#'
#' TODO
#'
#'@export
#'@useDynLib mvCircular
#'
custom_bessi0 <- function(x,expon.scaled=F){
  
  if(expon.scaled)
    ex <- 1
  else
    ex <- 0
  
  # Call procedure
  ret <- .C("__R_bessi0",
            x = as.double(x),
            expon = as.integer(ex),
            value = double(1),
            PACKAGE="mvCircular")
  
  # Get parameters
  return(ret$value)
}


#' Own bessi1
#'
#' TODO
#'
#'@export
#'@useDynLib mvCircular
#'
custom_bessi1 <- function(x,expon.scaled=F){
  
  if(expon.scaled)
    ex <- 1
  else
    ex <- 0
  
  # Call procedure
  ret <- .C("__R_bessi1",
            x = as.double(x),
            expon = as.integer(ex),
            value = double(1),
            PACKAGE="mvCircular")
  
  # Get parameters
  return(ret$value)
}