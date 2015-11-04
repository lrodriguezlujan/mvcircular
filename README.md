# mvCircular

Multivariate circular  R package 

This is an R package that includes some basic functions to use multivariate circular distributions. The current version implements four S3 clases:
- `mvVonMises` for the multivariate von Mises distribution
- `mvWrappedNormal` for the multivariate wrapped-normal distribution
- `mvWrappedCauchy` for the multivariate wrapped-Cauchy distribution
- `mvCircTruncNormal` for the multivariate normal distribution truncated to [μ - π ,μ + π )

**mvCircular** also implements functions to compute the geometric median for mv. circular data, a Kullback-Leibler divergence approximation, and three types of plots `plot`, `contourPlot` and `toursPlot`. Torus plot requires some functionality part of a [custom development](https://github.com/lrodriguezlujan/rthreejs/tree/add_composition) based on [rthreejs](http://bwlewis.github.io/rthreejs/).

# Installation

From github:

``` R 

# Install rthreejs branch (optional)
install_github("rthreejs", username = "lrodriguezlujan", ref = "add_composition")

# Install mvCircular
install_github("mvCircular", username = "lrodriguezlujan", ref = "master")

```

# Usage

``` R

# Load lib
library(mvCircular)

# Bivariate plots
vmbiv <- mvVonMises(c(3*pi/2, pi/2), rep(1.5,2), matrix( c(0,2,2,0) ,ncol=2,nrow=2) )
vmbiv <- normalize(vmbiv)
# Contour
contourPlot(vmbiv)
# Torus (needs add_composition branch)
torusPlot(vmbiv)

# Create samples from a distribution
samples <- rmvVonMises(1E4, rep(0,4), rep(1,4), matrix(0,ncol=4,nrow=4) )

# Median
geomedian.circular(samples)

# Fit
vmfit <- mvVonMises.fit(samples)

# Empiric KL divergence
vmfit.kl <- empKL.circular(samples, vmfit, m = 1E4)
vmfit.kl$kl

# Compute normalization term
vmfit <- normalize(vmfit)

# Plot
plot(vmfit, data = vmfit$fitted.data[1:100,])

```

# License

L-BFGS-B is released under the BSD 3-clause license

This package uses [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C):
- L-BFGS-B is released under the BSD 3-clause license (J.L. Morales and J. Nocedal. L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization (2011), ACM Transactions on Mathematical Software, Vol 38, Num. 1.)
- [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C) is released under the BSD 3-clause license

# Authors

This R package is written by Luis Rodríguez, [luis.rodriguezl@upm.es](mailto:luis.rodriguezl@upm.es)


