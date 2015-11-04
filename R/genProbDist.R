#' @name mvCircularProbDist
#' @rdname mvCircularProbDist
#' @title Multivariate circular distribution
#' 
#' Generic defnition for a multivariate ciruclar distribution
#' 
#' @author Luis Rodriguez Lujan 
#' 
#' @keywords multivariate circular
#' 
#' @export
NULL

# Class name
PROBDIST_CLASS <- "mvCircularProbDist"

#' @param dim  Number of variables
#' @param \dots (\code{mvCircularProbDist}) Named list with aditional parameters to include in the object
#' 
#' @export
mvCircularProbDist <- function(dim, ... ){
  
  ret <- list(
    dim = dim
  )
  ret <- c(ret, list(...))
  
  class(ret) <- PROBDIST_CLASS
  
  return(ret)
}

#' @return \code{getSamples} returns a circular dataframe with n samples from the given circular distribution
#' 
#' @param obj A multivariate circular distribution object
#' @param n The number of samples to return
#' 
#' @rdname mvCircularProbDist
#' 
#' @export
getSamples <- function(obj, n, ...){
  UseMethod("getSamples", obj)
}

#' @param x A circular point (p-vector)
#' 
#' @return \code{fval} evaluates the density function in the point \code{x}
#'
#' @rdname mvCircularProbDist
#'
#' @export
fval <- function(obj, x, ...){
  UseMethod("fval", obj)
}

#'
#' @return \code{Fval} returns the cumulative distribution function in the point \code{x}
#' 
#' @rdname mvCircularProbDist
#' 
#' @export
Fval <- function(obj, x, ...){
  UseMethod("Fval", obj)
}

#'
#' @return \code{normalize} returns the object with the normalization term set
#' 
#' @rdname mvCircularProbDist
#' 
#' @export
normalize <- function(obj, ...){
  UseMethod("normalize",obj)
}

# Default behaviour: No normalization
normalize.default <- function(obj){
  return(obj)
}


#' @param i Component index
#' 
#' @return \code{circMarginal} returns the \code{i}-th compontent marginal density function evaluated at \code{x}
#'
#' @rdname mvCircularProbDist
#'
#' @export
circMarginal <- function(obj, x, i, ...){
  UseMethod("circMarginal",obj)
}

#' @rdname mvCircularProbDist
#' 
#' @param j Component index
#'
#' @return \code{circCor} returns the dependency between \code{i}-th and \code{j}-th compontents
#' @export
circCor <- function(obj, i, j, ...){
  UseMethod("circCor",obj)
}

#' @rdname mvCircularProbDist
#' 
#' @return \code{circMarginalMean} returns the marginal mean of the \code{i}-th compontent
#' 
#' @export
circMarginalMean <- function(obj, i, ...){
  UseMethod("circMarginalMean",obj)
}

#'
#' 
#' @return \code{circMarginalConcentration} returns the marginal concentration of the \code{i}-th compontent
#' 
#' @rdname mvCircularProbDist
#' 
#' @export
circMarginalConcentration <- function(obj, i, ...){
  UseMethod("circMarginalConcentration",obj)
}
 
#' \code{contourPlot} creates a 2D plot where the level based on on the density function value. 
#' Only available for bivariate distributions
#' 
#' @rdname mvCircularProbDist
#'
#' @export
contourPlot <- function(obj, ...){
  UseMethod("contourPlot",obj)
}


#' 
#' \code{torusPlot} creates a 3D interactive plot using \link{threejs} library. The surface of a torus is coloured
#' based on the density function value. Only available for bivariate distributions
#'
#' @rdname mvCircularProbDist
#'
#' @export
torusPlot <- function(obj, ...){
  UseMethod("torusPlot",obj)
}

#' @rdname mvCircularProbDist
#' 
#' @param npoints Number of points to evaluate
#' @param \dots (\code{torusPlot}) Additional params for \code{feval}
#' 
#' @examples
#'  obj <- mvVonMises(c(3*pi/2, pi/2), rep(1.5,2), matrix( c(0,2,2,0) ,ncol=2,nrow=2) )
#'  obj <- normalize(obj)
#'  torusPlot(obj)
#'  
#'  @export
torusPlot.mvCircularProbDist <- function(obj,  npoints = 1E4,  ...) {
  
  if ( obj$dim != 2 ) stop("Contour plot only works for 2-d distributions")
  
  # Generate lattice
  dims <- obj$dim 
  # Compute number point per dim
  samplesPerDim <- floor(npoints ^ (1/dims))
  
  # Since 2pi == 0 we create  n+1 samples per dim and remove the last one (2pi) to avoid double-evaluating the same point
  dimvals <- sort(seq(from = 0, to = 2*pi, length.out = (samplesPerDim + 1))[-(samplesPerDim + 1)] )
  points  <- expand.grid( split(rep(dimvals, each = dims),1:dims) )
  z <- fval(obj, points, ...)
  
  # Create colorscale
  colorscale <- threejs::colorscaleThreeJS(max = max(z), min = min(z), legend.on = T, legend.labels.on = T, 
                                  legend.labels.title = 'f(x)')
  
  # Lighting
  light <- list( threejs::ambientLight("#555555"), threejs::directionalLight("#FFFFFF",intensity = 0.5,from = c(0,0,1)) )
  
  # Helpers
  helpers <- list( threejs::axisHelper(), threejs::gridHelper(rotation = c(0,pi/2,pi/2), size = 100) )
  
  # Color function
  color.extra <- list( npoints = samplesPerDim, r = 7.5, tuber = 2.5, z = z )
  colorFunction <- htmlwidgets::JS("function(v,g,i,e,cs){ 
  // Compute angles
  var theta = Math.atan2(v.y,v.x);
  theta += (theta < 0 ? 2 * Math.PI : 0)
  var phi = Math.asin(v.z/e.tuber);
  phi = ( (Math.sqrt((v.y*v.y) + (v.x*v.x)) < e.r) ? Math.PI-phi : phi );
  phi += (phi < 0 ? 2 * Math.PI : 0)
  // Get poistion
  return( cs.getColor( e.z[ e.npoints * Math.floor( theta * e.npoints / (2*Math.PI) ) + Math.floor( phi * e.npoints / (2*Math.PI) )] ) );
  }")
  
  # Object (torus)
  obj <- threejs::objThreeJS( 
    threejs::torusGeometry( inner = 5, outer = 10, color.extra = color.extra, vertex.color.fun = colorFunction),
    threejs::lambertMaterial(vertexColors = htmlwidgets::JS("THREE.VertexColors" ), side = htmlwidgets::JS("THREE.DoubleSide") ))
  
  # Plot
  threejs::composition( objects = obj , lights = light, helpers = helpers, colorscale = colorscale)
}

#' @rdname mvCircularProbDist
#' 
#'  @param \dots (contourPlot) Additional parameters for \code{\link{contour}} or \code{\link{filled.contour}}
#'
#' @examples
#'  obj <- mvVonMises(c(3*pi/2, pi/2), rep(1.5,2), matrix( c(0,2,2,0) ,ncol=2,nrow=2) )
#'  obj <- normalize(obj)
#'  contourPlot(obj)
#'  
#'  @export
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

#'
#' @rdname mvCircularProbDist
#'
#' @param data Data to plot in the lower triagonal segment. Usually data used for fitting.
#' @param \dots (plot) Additional arguments for \link{circMarginal}
#' @param color Plot main color 
#' @param cex.params Parameters font size
#' @param rose.bins Number of bins in the rose plot
#' @param cex.sigma Upper triangular dep. font size
#' @param pch Point character to use in the 2D plots
#' 
#' @importFrom circular rose.diag
#' 
#' @examples
#'  obj <- mvVonMises(c(3*pi/2, pi/2, 0, 0), c(1.5,2,1,5), matrix( c(0,2,-1,0, 2,0,1,0, -1,1,0,1, 0,0,1,0) ,ncol=4,nrow=4) )
#'  plot(obj)
#' 
#' @export
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