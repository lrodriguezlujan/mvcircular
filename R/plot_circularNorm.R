#'plot multivariate circular normal results
#'@export
plot_circularNorm <- function(mvDist,data,color="red",cex.params=2,rose.bins=25,cex.sigma=2, pch=18){
  
  #Get dist params
  params <- mvDist$getParams()
  mu <- as.numeric(params$mu)
  sigma_inv <- ginv(params$sigma)
  
  # Number of columns
  p<-length(mu)
  
  # pxp plot , 
  par(mfrow=c(p,p+1),
      oma=c(2.5,0,1.2,0))
  
  options(warn=-1)
  # Plot time!
  for(i in 1:p){
    for(j in 1:(p+1)){
      
      #Plot rose diag. in the diagonal
      if(j==1){
        par(mar=c(0, 0, 0, 0))
        plot.new()
        m = mu[i]
        k = sigma_inv[i,i]
        #txt_m <- format(c(m, 0.123456789), digits=3)[1]
        #txt_k <- format(c(k, 0.123456789), digits=3)[1]
        txt <- sprintf("θ[%d]\n\nμ = %.3f\nκ = %.3f",i,m,k)
        text(0.5, 0.5, txt, cex =  cex.params)
      }
      else if(i==(j-1)){
        par(mar=c(0, 0, 0, 0))
        den_x <- seq(as.numeric(params$lower[i]),as.numeric(params$upper[i]),0.01)
        den_y<-dtmvnorm.marginal(den_x,mean=mu,sigma=params$sigma,n=i,
                              lower=as.numeric(params$lower),upper=as.numeric(params$upper))
        circDensity <- list(data=data,x=circular(den_x),y=den_y,bw=1,N=0,call=nrow(data),data.name="",has.na=F)
        class(circDensity) <- "density.circular"
        plot(circDensity,shrink=1.25,col=color,main="",xlab="",ylab="")
        rose.diag(as.circular(data[,i],control.circular=list(0)), bins=rose.bins, col=color,main="", add=T)
        box()
      }
      else{
        if(i<(j-1)){
          par(mar=c(0, 0, 0, 0))
          plot.new()
          r = sigma_inv[i,j-1]
          txt <- format(c(r, 0.123456789), digits=3)[1]
          text(0.5, 0.5, txt, cex = cex.sigma)
        }
        else{
          par(mar=c(1, 1, 1, 1))
          plot(as.numeric(data[,i])*(180/pi),as.numeric(data[,j-1])*(180/pi),type="p",xlim=c(0, 360),ylim=c(0, 360),pch=pch,col=color,lwd=0.1)
        }
      }
    }
  }
  options(warn=-0)
}