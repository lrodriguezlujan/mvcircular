/*******************
 * Gibbs sampler for the multivariate von Mises distribution
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: mvvonmises_sampler.h 19 2015-05-11 15:05:35Z lrodriguez $
 **/

#ifndef __MVVONMISES_SAMPLER_H_
#define __MVVONMISES_SAMPLER_H_

#include <stdint.h>

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
double* mvvm_sampler_gibbs(int n, int p, double* mu, double* kappa, double* lambda,int initk,int k);
int mvvm_sampler_safe_gibbs(int n, int p, double* mu, double* kappa, double* lambda,int initk,int k,uint64_t *seeds, double *ret); // Thread safe version

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
int mvvm_sampler_safe_rejection(int n, int p, double* mu, double* kappa, double* lambda,double lmin, uint64_t *seeds,double *ret);

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

int mvvm_sampler_safe_rsgibbs(int n, int p, double* mu, double* kappa, double* lambda, int prec, int *alphaTable, int thinning, double *theta, uint64_t *seeds,double *ret);

/*************************
 *
 * R - WRAPPERS
 *
 *************************/

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
void __R_mvvm_sampler_safe_gibbs(int* n, int* p, double* mu, double* kappa, double* lambda,int* initk,int* k,uint64_t *seeds, double *ret, int * errno);

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
int __R_mvvm_sampler_safe_rejection(int *n, int *p, double* mu, double* kappa, double* lambda,double *lmin, uint64_t *seeds,double *ret,int *errno);


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
void __R_mvvm_sampler_safe_rsgibbs(int *n, int *p, double* mu, double* kappa, double* lambda, int *prec, int *alphaTable, int *thinning, double *theta,uint64_t *seeds,double *ret,int *errno);



#endif
 
