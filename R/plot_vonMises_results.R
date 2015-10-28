#'plot multivariate von mises results
#'@export
plot_circularVm_results <- function(mvDist,data,color="red",cex.params=2,rose.bins=25,cex.lambda=2, pch=18){
  
  #Get dist params
  params <- mvDist$getParams()
  mu <- params$mu
  kappa <- params$kappa
  lambda <- params$lambda
  
  # Number of columns
  p<-length(mu)
  
  # pxp plot , 
  par(mfrow=c(p,p+1),
      oma=c(2.5,0,1.2,0))
  
  # Generate samples
  samples_mvvm <- rmvvonmises(1E4, mu ,kappa ,lambda);
  
  options(warn=-1)
  # Plot time!
  for(i in 1:p){
    for(j in 1:(p+1)){
      
      #Plot rose diag. in the diagonal
      if(j==1){
        par(mar=c(0, 0, 0, 0))
        plot.new()
        m = mu[i]
        k = kappa[i]
        #txt_m <- format(c(m, 0.123456789), digits=3)[1]
        #txt_k <- format(c(k, 0.123456789), digits=3)[1]
        txt <- sprintf("θ[%d]\n\nμ = %.3f\nκ = %.3f",i,m,k)
        text(0.5, 0.5, txt, cex =  cex.params)
      }
      else if(i==(j-1)){
        par(mar=c(0, 0, 0, 0))
        x2<-density.circular(samples_mvvm[,i],kernel="vonmises",bw=kappa[i]) #obtenemos las densidades de esos valores
        plot(x2,shrink=1.2,col=color,main="",xlab="",ylab="",tcl = c(0,pi))
        rose.diag(as.circular(data[,i],shrink=0.9,control.circular=list(0)), bins=rose.bins, col=color,main="",axes=F, add=T)
        box()
      }
      else{
        if(i<(j-1)){
          par(mar=c(0, 0, 0, 0))
          plot.new()
          r = lambda[i,j-1]
          txt <- format(c(r, 0.123456789), digits=3)[1]
          text(0.5, 0.5, txt, cex = cex.lambda )
        }
        else{
          par(mar=c(1, 1, 1, 1))
          plot(data[,i]*(180/pi),data[,j-1]*(180/pi),type="p",xlim=c(0, 360),ylim=c(0, 360),pch=pch,col=color,lwd=0.1)
        }
      }
    }
  }
  options(warn=-0)
}
