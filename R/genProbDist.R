PROBDIST_CLASS <- "mvCircularProbDist"
#' genProbDist Methods
#' 

#'
#'
#' @export
getSamples <- function(obj, n, ...){
  UseMethod("getSamples", obj)
}

#'
#'
#' @export
fval <- function(obj, x, ...){
  UseMethod("fval", obj)
}

#'
#'
#' @export
Fval <- function(obj, x, ...){
  UseMethod("Fval", obj)
}

#'
#'
#'@export
normalize <- function(obj, ...){
  UseMethod("normalize",obj)
}

#'
#'
#'@export
circMarginal <- function(obj, x, i, ...){
  UseMethod("circMarginal",obj)
}

#'
#'
#'@export
circCor <- function(obj, i, j, ...){
  UseMethod("circCor",obj)
}

#'
#'
#'@export
circMarginalMean <- function(obj, i, ...){
  UseMethod("circMarginalMean",obj)
}

#'
#'
#'@export
circMarginalConcentration <- function(obj, i, ...){
  UseMethod("circMarginalConcentration",obj)
}

#'
#'
#'@export
contourPlot <- function(obj, ...){
  UseMethod("contourPlot",obj)
}


#' Plot Bivariate contours
#'
#'@examples
#'  obj <- mvVonMises(c(3*pi/2, pi/2), rep(1.5,2), matrix( c(0,2,2,0) ,ncol=2,nrow=2) )
#'  obj <- normalize(obj)
#'  contourPlot(obj)
#'  
#'@export
contourPlot.mvCircularProbDist <- function(obj,  npoints = 1E4, color = T, k = 4, ...) {
  
  if ( obj$dim != 2 ) stop("Contour plot only works for 2-d distributions")
  
  # Generate lattice
  dims <- obj$dim 
  # Compute number point per dim
  samplesPerDim <- floor(npoints ^ (1/dims))
  
  # Since 2pi == 0 we create  n+1 samples per dim and remove the last one (2pi) to avoid double-evaluating the same point
  dimvals <- seq(from = 0, to = 2*pi, length.out = (samplesPerDim + 1))[-(samplesPerDim + 1)] 
  points  <- expand.grid( split(rep(dimvals, each = dims),1:dims) )
  
  z <- matrix( fval(obj, points, k), ncol = samplesPerDim, nrow = samplesPerDim, byrow = F )
  if (!color)
    contour(x = dimvals, y = dimvals, z = z, ...)
  else
    filled.contour(x = dimvals, y = dimvals, z = z, ...)
  
}


#'Plot multivariate circular distribution
#'
#' @param obj
#' @param data
#' @param ...
#' @param colors
#' @param cex.params
#' @param rose.bins
#' @param cex.sigma
#' @param pch
#' 
#' @importFrom MASS ginv
#' @importFrom circular rose.diag
#' 
#'@export
plot.mvCircularProbDist <- function(obj, data = NULL, n = 100 ,
                             ...,
                             color = "red",
                             cex.params = 2,
                             rose.bins = 25,
                             cex.sigma = 2, 
                             pch = 18){
  
  # If data is null and n != 0 we create n samples 
  if ( is.null(data) && n > 0 )
    data <- getSamples(obj, n)
  
  # Number of columns
  p <- obj$dim
  
  # p x p+1 plot , 
  par(mfrow = c(p,p + 1),
      oma = c(2.5,0,1.2,0))
  
  options(warn = -1)
  
  # Plot time!
  for (i in 1:p) {
    for (j in 1:(p + 1)) {
      
      #Plot rose diag. in the diagonal
      if (j == 1) {
        par(mar = c(0, 0, 0, 0))
        plot.new()
        m = circMarginalMean(obj,i)
        k = circMarginalConcentration(obj,i)
        txt <- sprintf("θ[%d]\n\nμ = %.3f\nκ = %.3f",i,m,k)
        text(0.5, 0.5, txt, cex =  cex.params)
      }
      else if (i == (j - 1)) {
        
        par(mar = c(0, 0, 0, 0))
        den_x <- seq( circMarginalMean(obj,i) - pi,  circMarginalMean(obj,i) + pi, 0.01)
        den_y <- circMarginal(obj, den_x, i, ... )
        # Scale it
        den_y <- den_y / (2.25*max(den_y) )
        
        # Create circular density from data
        circDensity <- list(data = data[,i], x = circular(den_x),
                            y = den_y, bw = 1, N = 0, call = nrow(data),
                            data.name = "", has.na = F)
        class(circDensity) <- "density.circular"
        
        # Plot circular density
        plot(circDensity,shrink = 1.25, col = color, main = "", xlab = "", ylab = "")
        circular::rose.diag(as.circular(data[,i],control.circular = list(0)), bins = rose.bins, 
                            col = color, main = "", add = T)
        box()
        
      }
      else {
        if (i < (j - 1)) {
          par(mar = c(0, 0, 0, 0))
          plot.new()
          r = circCor(obj, i, j - 1 )
          txt <- format(c(r, 0.123456789), digits = 3)[1]
          text(0.5, 0.5, txt, cex = cex.sigma)
        }
        else{
          par(mar = c(1, 1, 1, 1))
          plot(as.numeric(data[,j - 1])*(180/pi), as.numeric(data[, i])*(180/pi),
               type = "p", xlim = c(0, 360), ylim = c(0, 360), pch = pch, col = color, lwd = 0.1)
        }
      }
    }
  }
  options(warn = 0)
}