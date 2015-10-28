/*******************
 * Sampler for the univariate von Mises distribution
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: vonmises_sampler.h 19 2015-05-11 15:05:35Z lrodriguez $
 **/

#ifndef __VONMISES_SAMPLER_H_
#define __VONMISES_SAMPLER_H_


#include <stdint.h>



/**
 * Von mises density function. Informative
 *
 * @param [in] theta
 * @param [in] mu
 * @param [in] kappa
 *  
 * @return density
 */
double vmDensity(double theta, double mu, double kappa);

/**
 * Rejection-based sampler for the von Mises distribution. Retruns several samples per call
 *
 *  @param [out] x Real ndim vector to return the samples
 *  @param [in] n Number of samples to generate
 *  @param [in] mu Mu parameter of the vm dist
 *  @param [in] kappa Kappa parameter of the vm dist
 */
void rvm(double *x, int n, double mu, double kappa);
void rvm_r(double *x, int n, double mu, double kappa, uint64_t *seeds);

/******
 * Rejection-based sampler for the von Mises distribution. Returns one sample per call
 *
 *  @param [in] mu Mu parameter of the vm dist
 *  @param [in] kappa Kappa parameter of the vm dist*
 *
 *  @return one sample
 * *****/
double rvm_s(double mu, double kappa);
double rvm_s_r(double mu, double kappa, uint64_t * seeds);

#endif
