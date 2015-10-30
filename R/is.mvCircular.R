#' Multivariate circular data test
#' 
#' Checks if an object is multivariate circular data
#' 
#' @param x The object to test
#' 
#' @return True/False
#' 
#' @importFrom circular is.circular
#'
#' @export
is.mvCircular <- function(x){
  
  # Admit either a circular df or a numeric matrix
  if ( is.data.frame(x) ) return( all(vapply(x, circular::is.circular, FUN.VALUE = logical(1) )) )
  else if (is.matrix(x) ) return( is.numeric(x) )
  else return(FALSE)
}

#' Multivariate circular data conversor
#' 
#' Transforms a given object into a standard multivariate circular dataframe
#' 
#' @param x The object to convert
#' @param control.ciruclar parameter for \link{circular::as.circular}
#' 
#' @return Data frame with circular columns
#' 
#' @importFrom circular as.circular
#' 
#' @export
as.mvCircular <- function(x, control.circular = list(modulo = "2pi", zero = 0, template = "none",
                                                     type = "angles", units = "radians", rotation = "counter"), ...){
  
  if ( is.data.frame(x) && all(vapply(x, circular::is.circular, FUN.VALUE = logical(1) )) )
    return(x)
  # Double to circular 
  else if (is.data.frame(x) && all(vapply(x, is.numeric, FUN.VALUE = logical(1) )) )
    return( as.data.frame(lapply( x, circular::as.circular, control.circular = control.circular )) )
  else if (is.matrix(x) && is.numeric(x) ) 
    return( as.data.frame( lapply( lapply(as.data.frame(x), as.numeric) , 
                                   circular::as.circular, control.circular = control.circular )) )
  else
    stop("Cannot convert to multivariate circular")
}