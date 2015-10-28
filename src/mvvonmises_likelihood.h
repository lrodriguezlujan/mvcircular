/*******************
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: mvvonmises.c 13 2015-04-17 08:17:47Z lrodriguez $
 **/

#ifndef _MVVM_LIKELIHOOD_H_
#define _MVVM_LIKELIHOOD_H_

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


/***
 * Multivaraite von mises density. Not normalized
 *
 * @param [in] p
 * @param [in] theta
 * @param [in] mu
 * @param [in] kappa
 * @param [in] Lambda
 *
 * @return Density val
 * 
 */
double mvvm_density_unregularized(int p, double *theta, double *mu, double *kappa, double *lambda);

/*****
 * Optimized version. Computes simplified log likelihood and its derivatives
 *
 * Removes constant terms. Uses matrix multiplications. Reduce number of trigonometric opers.
 * Reuses terms. Marginal mu is not longer needed, derivatives are simplified as well. 
 * 
 * Instead of PL we will name it "loss function"
 *
 * Lambda is taken as a full row-leading pxp matrix with 0's in the diagonal. But d_lambda is only
 * computed for the upper triangle.
 *
 * If d_kappa or d_lambda are NULL, derivatives are not computed.
 *
 * To even speed up the process a little bit more, we make the following assumptions:
 *
 *  Mu is 0 for every j
 *  (input) Theta_i_j = Sin( (real) Theta_i_j - (real)mu_j )
 *
 *
 * @param [in] n
 * @param [in] p
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta
 * @param [out] ro
 * @param [out] d_kappa
 * @param [out] d_lambda
 *
 * @return loss function value
 */
long double mv_vonmises_lossFunction(int n,int p, double* kappa, double* lambda, long double *S, long double *C, long double *ro, double *d_kappa, double *d_lambda);

double vm_loglikelihood_lambda_deriv_old(int n, int p, int R,int S, double *mu, double *kappa, double *lambda, double *theta);

double vm_loglikelihood_kappa_deriv_old(int n, int p, int R, double *mu, double *kappa, double *lambda, double *theta);

double mv_vonmises_logpseudolikelihood_nonopt(int n,int p, double* mu, double* kappa, double* lambda, double *theta);

#endif
