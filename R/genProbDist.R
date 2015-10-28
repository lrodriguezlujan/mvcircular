#' genProbDist
#' 
#' Generic probability distribution class
#' 
#' Defines a generic probability distribution with common methods such as sample, fit, f, F, getParams, setParams
#' that should be defined by its instantiations. It is like a Java Abstract class (it cannot be istantiated) and will
#' throw an error in that case. 
#' 
#' To be completed with more functions / attributes / etc.
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords probability distribution
#' 
#' @section Methods:
#' 
#' \describe{
#'   \item{\code{fit(x,...)}}{This method uses data matrix or dataframe \code{x} to compute the parameters of the distribution}
#'   \item{\code{sample(n,...)}}{This method generates \code{n} samples of the distribution with current parameters}
#'   \item{\code{plot(...)}}{This method plots a representative graph of the prob. distribution, e.g. the density function for univariate}
#'   \item{\code{getParams()}}{This method returns the list of parameters}
#'   \item{\code{setParams(...)}}{This method sets the parameter given the arguments}
#'   \item{\code{f(...)}}{This method returns the density in the given points}
#'   \item{\code{F(...)}}{This method sets the cdf in the given points}
#'   \item{\code{print()}}{This method prints a string that describes de distribution and its parameters}
#' }
#' 
#' @section Attributes:
#' 
#' \describe{
#'   \item{\code{name}}{<PRIVATE> Distribution name. To be set during instantiation}
#'   \item{\code{parameters}}{<PRIVATE> List of parameters}
#'   \item{\code{fitted}}{<PRIVATE> Fitted from data flag}
#'   \item{\code{fitData}}{<PRIVATE> Data used to fit parameters}
#' }
#' 
genProbDist <- R6Class("probDist",
    public = list(
      # Initializer
      initialize = function(){ stop("probDist is an abstract class. It cannot be instantiated")},
      # Fit from data
      fit = function(){ warning("Not implemented")},
      # Sampling
      sample = function(n){ warning("Not implemented")},      
      # Density value
      f = function(x){ warning("Not implemented")},
      # Cummulative density function
      F = function(x){ warning("Not implemented")},
      # Plot distribution (TBD in subclases)
      plot = function(){ warning("Not implemented")},
      # Returns parameter list
      getParams = function(){
        return(private$parameters)
      },
      # Changes parameters
      setParams = function(...){
        args <- list(...)
        # Lets keep it simple, ... is the param list 
        private$parameters <- args
        private$fitted <- F
        private$fitData <- NULL
        return(self)
      },
      # Print distribution and parameters
      print = function(...){
        str <- sprintf("<%s> (%s) - Parameters:\n","probDist",private$name)
        for(i in 1:length(private$parameters)){
          str <- paste0(str,attributes(private$parameters)$names[[i]],":",toString(private$parameters[[i]]),"\n", collapse = " ")
        }
        
        cat(str)
        invisible(self) 
      }
      
    ),
    private = list(
      parameters = list(),
      name = NULL,
      fitted = F,
      fitData = NULL
    )
)