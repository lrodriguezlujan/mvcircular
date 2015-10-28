#' MvVmDist
#' 
#' Multivariate von Mises distribution class definition
#' 
#' The definition of the R6 class that corresponds to the multivariate von Mises distribution.
#' The idea behind this class is to be extended to other classes in the future in order to
#' reuse external methods via public methods, such as sample, fit, f, F, etc. although specific
#' parameters, fit or sampling methods may differ.
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords vonMises multivariate distribution
#' @inheritParams genProbDist
#' @seealso genProbDist
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
#' }
#' 
MvVmDist <- R6Class("MvVmDist",
            inherit = genProbDist,
            public = list(
              # Initializer
              initialize = function(mu,kappa,lambda){
                            
                # Set distribution name
                private$name <- "Multivariate von Mises"
                
                # If params are missing
                if( (!missing(mu)) && (!missing(kappa)) && (!missing(lambda)) ){
                
                  # Check parameters length (mu == kappa == nrow lambda == ncol lambda)
                  if(length(mu)!=length(kappa)){stop("Parameters length do not match")}
                  else if(length(mu)!=nrow(lambda)){stop("Parameters length do not match")}
                  else if(nrow(lambda) != ncol(lambda)){stop("Lambda is not square")}
                  
                  # Check mu is circular
                  if(!is.circular(mu)){stop("Mu must be circular (see ?as.circular())")}
                  else{
                    # Circular conversion data
                    mu <- conversion.circular(mu,modulo="2pi")
                  }
                  
                  # Check kappa is positive
                  if(any(kappa<=0)) stop("Kappa must be positive")
                  
                  # Fill: is.positive(-Lambda diag(kappa...))
                  private$parameters <- list( mu=mu, kappa=kappa, lambda = lambda)
                }
                else
                {
                  private$parameters <- list( mu=NULL, kappa=NULL, lambda = NULL)
                }
                self
              },
              # Fit from data
              fit = function(samples,...){
                #cat(sprintf("START Vm FIT. nsamples: %d\n",nrow(samples)))
                ret_fit <- fit_mvvonmises(samples,...)
                # FLag as fitted
                private$fitted <- T
                private$fitData <- samples
                # Fill parameters
                private$parameters <- list(mu=ret_fit$mu, kappa=ret_fit$kappa, 
                                           lambda=matrix(ret_fit$lambda,nrow=ncol(samples),ncol=ncol(samples),byrow=T))
                # Return loss value
                #cat(sprintf("END Vm FIT. nsamples: %d\n",nrow(samples)))
                return(ret_fit$loss)
              },
              # Sampling
              sample = function(n,mixChains=length(private$parameters$mu),...){
                
                #cat(sprintf("START Sampling, maxKappa: %f , sum(ABS)Lambda: %f\n",max(private$parameters$kappa),sum(abs(private$parameters$lambda))))
                ## Create p mixChains
                res <- rmvvonmises_rs(n,private$parameters$mu,private$parameters$kappa,private$parameters$lambda,chains = mixChains,...)
                #cat(sprintf("END Sampling, maxKappa: %f , sum(ABS)Lambda: %f\n",max(private$parameters$kappa),sum(abs(private$parameters$lambda))))
                
                ## Join them
                if(is.list(res)){
                  if(length(res)==1)
                    return(circular(res[[1]]))
                  else
                    return(circular(Reduce(rbind,res)[sample.int(n*mixChains,size=n,replace = F),]))
                }
                else
                  return (circular(res,modulo="2pi"))
                #return(conversion.circular(rmvvonmises_rs(n,private$parameters$mu,private$parameters$kappa,private$parameters$lambda,...,chains = 1)[[1]],modulo="2pi"))
              },
              # Density value
              f = function(x){ 
                warning("not normalized")
                aux <- (x - private$parameters$mu)
                return ( exp(as.numeric( crossprod(private$parameters$kappa,cos(aux)) + 
                           0.5* ( sin(aux) %*% private$parameters$lambda %*% sin(aux) ) ))) },
              # Cummulative density function
              F = function(x){ warning("Not implemented")},
              # Plot distribution (TBD in subclases)
              plot = function(data=NULL,n=1000,...){ 
                if(is.null(data)){
                  if(is.null(private$fitData))
                    plot_circularVm_results(self,self$sample(n),...)
                  else
                    plot_circularVm_results(self,private$fitData,...)
                }
                else{
                  plot_circularVm_results(self,data,....)
                }
                return(self)
              }
            )
)
            
                    