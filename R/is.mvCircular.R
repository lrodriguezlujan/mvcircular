#' 
#' @importFrom circular is.circular
#'
is.mvCircular <- function(x){
  
  # Admit either a circular df or a numeric matrix
  if ( is.data.frame(x) ) return( all(vapply(x, circular::is.circular, FUN.VALUE = logical(1) )) )
  else if (is.matrix(x) ) return( is.numeric(x) )
  else return(FALSE)
}

#'
#' @importFrom circular as.circular
#'
as.mvCircular <- function(x, control.circular = list(), ...){
  
  if ( is.data.frame(x) && all(vapply(x, circular::is.circular, FUN.VALUE = logical(1) )) )
    return(x)
  # Double to circular 
  else if (is.data.frame(x) && all(vapply(x, circular::is.numeric, FUN.VALUE = logical(1) )) )
    return( as.data.frame(lapply( x, circular::as.circular, control.circular = control.circular )) )
  else if (is.matrix(x) && is.numeric(x) ) 
    return( as.data.frame(apply( x, 2, circular::as.circular, control.circular = control.circular )) )
  else
    stop("Cannot convert to circular")
}