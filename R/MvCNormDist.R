#' MvCircNormDist
#' 
#' Multivariate circular Normal distribution
#' 
#' The definition of the R6 class that corresponds to the multivariate normal distribution wrapped.
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
#' @seealso genProbDist tmvtnorm
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
MvCNormDist <- R6Class("MvCNormDist",
                    inherit = genProbDist,
                    public = list(
                      # Initializer
                      initialize = function(mu,sigma){
                        
                        # Set distribution name
                        private$name <- "Multivariate circular Normal"
                        
                        # If params are missing
                        if( (!missing(mu)) && (!missing(sigma))){
                          
                          # Check parameters length (mu == kappa == nrow lambda == ncol lambda)
                          if(length(mu)!=nrow(sigma)){stop("Parameters length do not match")}
                          else if(nrow(sigma) != ncol(sigma)){stop("Sigma is not square")}
                          
                          # Check mu is circular
                          if(!is.circular(mu)){stop("Mu must be circular (see ?as.circular())")}
                          else{
                            # Circular conversion data
                            mu <- conversion.circular(mu,modulo="2pi")
                          }
                          
                          # Fill: is.positive(-Lambda diag(kappa...))
                          private$parameters <- list( mu=mu, sigma = sigma, lower = mu - pi, upper = mu + pi)
                        }
                        else
                        {
                          private$parameters <- list(mu=NULL, sigma=NULL, lower =NULL , upper = NULL)
                        }
                        self
                      },
                      f = function(x, ...){
                        return(dtmvnorm(x,
                                        mean =  as.numeric(private$parameters$mu),
                                        sigma = private$parameters$sigma,
                                        lower = as.numeric(private$parameters$lower),
                                        upper = as.numeric(private$parameters$upper),
                                        ...))
                      },
                      F = function(x, ... ){
                        return(ptmvnorm(as.numeric(private$parameters$lower),x ,
                                        mean =  as.numeric(private$parameters$mu),
                                        sigma = private$parameters$sigma,
                                        lower = as.numeric(private$parameters$lower),
                                        upper = as.numeric(private$parameters$upper),
                                        ...))
                      },
                      # Fit from data
                      fit = function(samples,...){
                        
                        # Number of variables
                        p <- ncol(samples)
                        n <- nrow(samples)
                        
                        # Samples should be circular
                        if(!is.circular(samples))
                          stop("Provided samples are not circular")
                        else
                          mean <- apply(samples,2,function(x){mean.circular(circular(x))})
                        
                        # Center sample
                        centeredSample<-circular(sweep(samples,2,mean),modulo="2pi")
                        centeredSample[centeredSample > pi] <- centeredSample[centeredSample > pi]- 2*pi
                        
                        # Fit (TODO: use other method than the default BFGS which can be too expensive in memory for high dims)
                        # FIXME: mle.tmvnorm IS REALLY SLOW (like 10x10 matrix takes forever)
                        #mle <- mle.tmvnorm(matrix(as.double(centeredSample),ncol=p,nrow=n),cholesky = T,
                        #                   lower=rep(-pi,p),
                        #                   upper=rep(pi,p))

                        # Get sigma coef and build matrix
                        #sigma <- matrix(NA,ncol=p,nrow=p);
                        #sigma_coefs <- mle@fullcoef[-p:-1]
                        #sigma[lower.tri(sigma,diag=T)] <- sigma_coefs
                        #sigma[upper.tri(sigma,diag=F)] <- sigma[lower.tri(sigma,diag=F)]
                        
                        if(n<=p){
                          sigma <- matrix(cov.shrink(centeredSample),nrow=p,ncol=p)
                        }
                        else
                        {
                          sigma <- stats::var(centeredSample)
                        }
                        
                        # Fill parameters
                        private$parameters <- list(mu=mean, sigma = sigma, 
                                                   lower = mean - pi, upper = mean + pi);
                        # Return loss value
                        #return(logLik(mle))
                        return(0)
                      },
                      # Sampling
                      sample = function(n,burn.in.samples=1E5,thinning=1000,start.value=private$parameters$mu,...){
                          if(length(private$parameters$mu)<100)
                            return(circular(rtmvnorm(n,mean=as.numeric(private$parameters$mu),sigma=private$parameters$sigma,
                                        lower=as.numeric(private$parameters$lower),
                                        upper=as.numeric(private$parameters$upper),algorithm="rejection",...),modulo="2pi"))
                          else
                            return(circular(rtmvnorm(n,mean=as.numeric(private$parameters$mu),sigma=private$parameters$sigma,
                                                     lower=as.numeric(private$parameters$lower),
                                                     upper=as.numeric(private$parameters$upper),algorithm="gibbs",
                                                     burn.in.samples=burn.in.samples,thinning=thinning,start.value=start.value
                                                     ,...),modulo="2pi"))
                      },
                      # Plot distribution (TBD in subclases)
                      plot = function(data=NULL,n=1000,...){ 
                        if(is.null(data)){
                          if(is.null(private$fitData))
                            plot_circularNorm(self,self$sample(n),...)
                          else
                            plot_circularNorm(self,private$fitData,...)
                        }
                        else{
                          plot_circularNorm(self,data,...)
                        }
                        return(self)
                      }
                    )
)
