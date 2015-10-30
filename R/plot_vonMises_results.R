#' Plot multivariate von Mises distribution
#' 
#' @param mu
#' @param kappa
#' @param lambda
#' @param data
#' @param color
#' @param cex.params
#' @param rose.bins
#' @param cex.sigma
#' @param pch
#' 
#' @importFrom circular density.circular as.circular rose.diag
#'  
#'@export
circularMvVm.plot <- function(mu,kappa,lambda,data,color="red",
                              cex.params = 2, rose.bins = 25, cex.lambda = 2, pch = 18){
  
  # Number of columns
  p <- length(mu)
  
  # pxp+1 plot , 
  par(mfrow = c(p,p + 1),
      oma = c(2.5,0,1.2,0))
  
  # Generate samples (for kernel dens estimation)
  samples_mvvm <- rmvVonMises(1E4, mu ,kappa ,lambda);
  
  options(warn = -1)
  
  # Plot time!
  for (i in 1:p) { #rows
    for (j in 1:(p + 1)) { # Cols
      
      # Plot parameters in the first column
      if (j == 1) {
        par(mar = c(0, 0, 0, 0))
        plot.new()
        m = mu[i]
        k = kappa[i]
        txt <- sprintf("θ[%d]\n\nμ = %.3f\nκ = %.3f",i,m,k)
        text(0.5, 0.5, txt, cex =  cex.params)
      }
      
      #Plot rose diag. in the diagonal
      else if ( i == (j - 1) ) {
        par(mar = c(0, 0, 0, 0))
        x2 <- circular::density.circular(samples_mvvm[,i],kernel = "vonmises",bw = kappa[i])
        plot(x2, shrink = 1.2,col = color,main = "",xlab = "",ylab = "",tcl = c(0,pi))
        circular::rose.diag(circular::as.circular(data[,i],shrink = 0.9,control.circular = list(0)), bins = rose.bins,
                            col = color, main = "", axes = F, add = T)
        box()
      }

      else if (i < (j - 1)) {
        # Plot lambda values in the upper triang.
        par(mar = c(0, 0, 0, 0))
        plot.new()
        r = lambda[i,j - 1]
        txt <- format(c(r, 0.123456789), digits = 3)[1]
        text(0.5, 0.5, txt, cex = cex.lambda )
      }
      else{
        # scatterplot in the lower triang
        par(mar = c(1, 1, 1, 1))
        plot( as.numeric(data[,i])*(180/pi), as.numeric(data[,j - 1]) * (180 / pi), type = "p",
             xlim = c(0, 360),ylim = c(0, 360),pch = pch,col = color,lwd = 0.1)
      }
    }
  }
  options(warn = 0)
}
