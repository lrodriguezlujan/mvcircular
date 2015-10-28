/*******************
 * Gibbs sampler for the multivariate von Mises distribution
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: mvvonmises_sampler.c 19 2015-05-11 15:05:35Z lrodriguez $
 **/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "mvvonmises_sampler.h"
#include "prng.h"
#include "vonmises_sampler.h"
#include "mvvonmises_marginalCond.h"


/*******************
 *
 * REAL METHODS
 *
 */

/**
 * Gibbs sampler for multivariate von Mises distribution
 *
 * @param [in] n Number of samples to generate
 * @param [in] p Number of variables
 * @param [in] mu pdim Real vector. Mu parameter of the mv-vm dist
 * @param [in] kappa pdim Real vector. Kappa parameter of the mv-vm dist
 * @param [in] lambda pxp Real matrix on row-leading order. Lambda parameter of the mv-vm dist
 * @param [in] initk initial number of iterations before accept
 * @param [in] k Number of iterations between samples.
 *
 * @returns nxp Real matrix on row-leading order with the samples
 */
double* mvvm_sampler_gibbs(int n, int p, double* mu, double* kappa, double* lambda,int initk,int k){

    // Sample and return
    double *theta; 
    double *ret;

    // Conditional mu and kappa
    double mu_c,kappa_c;

    //iterator
    int i,j;

    // Initialize theta (unique sample)
    theta = malloc(sizeof(double)*p);
    for(i=0;i<p;i++) theta[i]=1;

    // Initialize return matrix
    ret = malloc(sizeof(double) * n * p );
    
    // initial rejection
     for(i=0;i<=initk;i++){
        // Generate samples
        for(j=0;j<p;j++){
            condKappaMu(p,j,theta,mu,kappa,lambda,&mu_c,&kappa_c);
            theta[j] = rvm_s(mu_c,kappa_c);
        }
     }

    // Start iterating and filling ret (NOTE: initial rejection is initk + k)
    for(i=0;i<=(k*n);i++){
        // Generate one sample
        for(j=0;j<p;j++){
            condKappaMu(p,j,theta,mu,kappa,lambda,&mu_c,&kappa_c);
            theta[j] = rvm_s(mu_c,kappa_c);
        }
        if(i!=0 && i%k == 0 )
            memcpy(ret+((i/k)-1)*p,theta,sizeof(double)*p);
    }
    free(theta);
    return ret;
}

/* Thread safe version */
int mvvm_sampler_safe_gibbs(int n, int p, double* mu, double* kappa, double* lambda,int initk,int k,uint64_t *seeds,double *ret){

    // Sample and return
    double *theta; 

    // Conditional mu and kappa
    double mu_c,kappa_c;

    //iterator
    int i,j;

    // Initialize theta (unique sample)
    theta = malloc(sizeof(double)*p);
    for(i=0;i<p;i++) theta[i]=1;

    // initial rejection
     for(i=0;i<=initk;i++){
        // Generate samples
        for(j=0;j<p;j++){
            condKappaMu(p,j,theta,mu,kappa,lambda,&mu_c,&kappa_c);
            theta[j] = rvm_s_r(mu_c,kappa_c,seeds);
        }
     }

    // Start iterating and filling ret (NOTE: initial rejection is initk + k)
    for(i=0;i<=(k*n);i++){
        // Generate one sample
        for(j=0;j<p;j++){
            condKappaMu(p,j,theta,mu,kappa,lambda,&mu_c,&kappa_c);
            theta[j] = rvm_s_r(mu_c,kappa_c,seeds);
        }
        if(i!=0 && i%k == 0 )
            memcpy(ret+((i/k)-1)*p,theta,sizeof(double)*p);
    }
    free(theta);
}


/**
 * Rejection sampler for multivariate von Mises distribution.
 *
 * Rejection-based sampler as proposed by Mardia 2014. Only for
 * moderate and low values of p.
 *
 * @param [in] n Number of samples to generate
 * @param [in] p Number of variables
 * @param [in] mu pdim Real vector. Mu parameter of the mv-vm dist
 * @param [in] kappa pdim Real vector. Kappa parameter of the mv-vm dist
 * @param [in] lambda pxp Real matrix on row-leading order. Lambda parameter of the mv-vm dist
 *
 * @returns nxp Real matrix on row-leading order with the samples
 */
inline int mvvm_sampler_safe_rejection(int n, int p, double* mu, double* kappa, double* lambda,double lmin, uint64_t *seeds,double *ret){

    // iterators
    uint_fast16_t i,j,k;
    // Auxiliar values
    double *auxv;
    double v;

    auxv = malloc(sizeof(double) * p);

    i=0; 
    while (i<n) {

        // Generate p samples from VM(0,lmin/4) (step 1)
        rvm_r(ret+(i*p),p,0,lmin/4,seeds);
        
        for(j=0;j<p;j++) 
            if(xorshift4096_r(seeds) > (UINT64_MAX/2)) ret[i*p+j]= ret[i*p+j]/2  + M_PI ;
            else ret[i*p+j]/=2;

        // Compute terms.
        v = 0;

        for(j=0;j<p;j++)
            // compute sin(t_i) in aux
            auxv[j] = sin(ret[i*p +j]);

        for(j=0;j<p;j++){
            // k_i (cos(t_i -1 )
            v+= kappa[j]*(cos(ret[i*p + j]) -1);

            // add lmin to lambda diagonal
            lambda[j*(p+1)] = lmin;

            // Add sin term
            for(k=0;k<p;k++) // Lambda is symmetric, so instead of div. by 2 we can just compute "one triangle"
                v+=0.5*auxv[j]*auxv[k]*lambda[j*p + k];

            // Restore lambda
            lambda[j * (p+1)]=0;
        }
        if( ((xorshift4096_r(seeds)*1.0)/(UINT64_MAX)) <= exp(v) ){
            // Accept current (add mu)
            for(j=0;j<p;j++)
                ret[i*p +j] += mu[j];
            i++;
        }
    }
    free(auxv);

    return 0;
}


/**
 *
 * Random scan gibbs sampler - Implementation
 *
 */

/**
 * Random scan gibbs sampler (Base)
 *
 * Random scan gibbs sampler - Non adaptative. Generates k independent chains, choosing starting points at random, and selects
 * samples with probability alpha (prec digit precission) Marginals are von Mises @see vonmises_sampler. Thread safe version
 * 
 * No burn-in is performed so this method can be reused for more sophisticated gibbs sampling
 *
 * @param [in] n Number of samples to generate per chain
 * @param [in] p Number of variables (dimensions)
 * @param [in] mu Mean direction of the VM distribution
 * @param [in] Kappa Concentration vector of de vM distibution
 * @param [in] Lambda Correlation matrix of the vM distribution
 * @param [in] prec Decimal precision of the lookup array. Array length (10^digits)
 * @param [in] alphaTable Lookup array (length prec=10^digits) to select the next component to update
 * @param [in] thinning Value for thinnin (i.e Number of samples to skip). Set to 1 if thinning is not desired. Set to n+1 to get the last value in theta. (burn-in)
 * @param [in/out] theta Initial value for theta. Size p
 * @param [in/out] seeds Seeds vector to be used
 * @param [out] ret nxp Matrix with generated samples
 *
 * @return Error code. 0 OK.
 */
int mvvm_sampler_safe_rsgibbs(int n, int p, double* mu, double* kappa, double* lambda, int prec, int *alphaTable, int thinning, double *theta, uint64_t *seeds,double *ret){

    // Conditional mu and kappa
    double mu_c,kappa_c;

    //iterators
    uint_fast16_t i,j;

    // Start iterating and filling ret
    for(i=1;i<=n;i++){

            // Select component to update
            j = alphaTable[xorshift4096_r(seeds)%prec]; 

            // Generate one sample
            condKappaMu(p,j,theta,mu,kappa,lambda,&mu_c,&kappa_c);
            theta[j] = rvm_s_r(mu_c,kappa_c,seeds);

            // Thinning
            if( i%thinning == 0 )
                memcpy(ret+((i/thinning)-1)*p,theta,sizeof(double)*p);
    }
    return(0);
}

/***********
 *
 * R-C Wrappers
 *
 */

/**
 * Gibbs sampler for multivariate von Mises distribution (R-Wrapper)
 *
 * @param [in] n Number of samples to generate
 * @param [in] p Number of variables
 * @param [in] mu pdim Real vector. Mu parameter of the mv-vm dist
 * @param [in] kappa pdim Real vector. Kappa parameter of the mv-vm dist
 * @param [in] lambda pxp Real matrix on row-leading order. Lambda parameter of the mv-vm dist
 * @param [in] initk initial number of iterations before accept
 * @param [in] k Number of iterations between samples.
 *
 * @seeAlso mvvm_sampler_safe_gibbs
 */
void __R_mvvm_sampler_safe_gibbs(int* n, int* p, double* mu, double* kappa, double* lambda,int* initk,int* k,uint64_t *seeds, double *ret, int * errno){
   *errno= mvvm_sampler_safe_gibbs(*n,*p,mu,kappa,lambda,*initk,*k,seeds,ret);
}

/**
 * Rejection sampler for multivariate von Mises distribution. (R-wrapper)
 *
 * Rejection-based sampler as proposed by Mardia 2014. Only for
 * moderate and low values of p.
 *
 * @param [in] n Number of samples to generate
 * @param [in] p Number of variables
 * @param [in] mu pdim Real vector. Mu parameter of the mv-vm dist
 * @param [in] kappa pdim Real vector. Kappa parameter of the mv-vm dist
 * @param [in] lambda pxp Real matrix on row-leading order. Lambda parameter of the mv-vm dist
 *
 * @seeAlso mvvm_sampler_safe_rejection
 */
int __R_mvvm_sampler_safe_rejection(int *n, int *p, double* mu, double* kappa, double* lambda,double *lmin, uint64_t *seeds,double *ret,int *errno){
    *errno=mvvm_sampler_safe_rejection(*n, *p, mu, kappa, lambda, *lmin, seeds, ret);
}

/**
 * Random scan gibbs sampler (Base)
 *
 * Random scan gibbs sampler - Non adaptative. Generates k independent chains, choosing starting points at random, and selects
 * samples with probability alpha (prec digit precission) Marginals are von Mises @see vonmises_sampler. Thread safe version
 * 
 * No burn-in is performed so this method can be reused for more sophisticated gibbs sampling
 *
 * @param [in] n Number of samples to generate per chain
 * @param [in] p Number of variables (dimensions)
 * @param [in] mu Mean direction of the VM distribution
 * @param [in] Kappa Concentration vector of de vM distibution
 * @param [in] Lambda Correlation matrix of the vM distribution
 * @param [in] prec Decimal precision of the lookup array. Array length (10^digits)
 * @param [in] alphaTable Lookup array (length prec=10^digits) to select the next component to update
 * @param [in] thinning Value for thinnin (i.e Number of samples to skip). Set to 1 if thinning is not desired. Set to n+1 to get the last value in theta. (burn-in)
 * @param [in/out] theta Initial value for theta. Size p
 * @param [in/out] seeds Seeds vector to be used
 * @param [out] ret nxp Matrix with generated samples
 *
 * @return Error code. 0 OK.
 */
void __R_mvvm_sampler_safe_rsgibbs(int *n, int *p, double* mu, double* kappa, double* lambda, int *prec, int *alphaTable, int *thinning, double *theta, uint64_t *seeds,double *ret,int *errno){
    *errno = mvvm_sampler_safe_rsgibbs(*n,*p,mu,kappa,lambda,*prec,alphaTable,*thinning,theta,seeds,ret);
}

