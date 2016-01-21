/*******************
 * Module to fit a mVM distribution from data samples
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id:$
 **/

#ifndef _MVVONMISES_FIT_H__
#define _MVVONMISES_FIT_H__

/**
 * Maximum full pseudolikelihood method to estimate parameters.
 *
 * Minimization done by gradient L-BFGS with wolfe
 *
 * @param [in] p Number of variables
 * @param [in,out] mu pdim Real vector. Mu parameter of the mv-vm dist
 * @param [in,out] kappa pdim Real vector. Kappa parameter of the mv-vm dist
 * @param [in,out] lambda pxp Real matrix on row-leading order. Lambda parameter of the mv-vm dist
 * @param [in] n Number of samples
 * @param [in] samples
 * @param [in] phi Prior matrix
 * @param [in] H Confidence matrix
 * @param [in] verbose
 * @param [in] prec
 * @param [in] tol 
 * @param [in] mprec
 * @param [in] lower
 * @param [in] upper
 * @param [in] bounded
 *
 * @returns Natural logarithm of the pseudolikelihood of the fitted distribution (aprox)
 */
double mvvonmises_lbfgs_fit(int p, double* mu, double* kappa, double* lambda, int n, double *samples, double* phi, double *H,
                            int verbose, double prec, double tol, int mprec, double *lower, double *upper, int *bounded );

#endif
