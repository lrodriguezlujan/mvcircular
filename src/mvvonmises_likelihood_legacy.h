/*******************
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: mvvonmises.c 13 2015-04-17 08:17:47Z lrodriguez $
 **/

#ifndef _MVVM_LIKELIHOOD_LEGACY_H_
#define _MVVM_LIKELIHOOD_LEGACY_H_

/*************
 * Likelihood
 */

/* Transforms theta to Sin(theta-mu)
 *
 * Computes Sin(theta-mu) / Cos(theta-mu) for each variable in the dataset.
 *
 * @param [in] n
 * @param [in] p
 * @param [in] theta
 * @param [in] mu
 * @param [out] S
 * @param [out] C
 *
 */
void mv_theta_cos_sinTransform(int n, int p, double*theta, double* mu, long double *S, long double *C);

double vm_loglikelihood_lambda_partial_naive(int n, int p, int R,int S, double *mu, double *kappa, double *lambda, double *theta);

double vm_loglikelihood_kappa_partial_naive(int n, int p, int R, double *mu, double *kappa, double *lambda, double *theta);

double mv_vonmises_logpseudolikelihood_naive(int n,int p, double* mu, double* kappa, double* lambda, double *theta);

#endif
