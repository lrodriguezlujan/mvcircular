/*******************
 * Module to fit a mVM distribution from data samples
 * with L1 Penalization
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id:$
 **/

// We should start using pragma include once...
#ifndef __MVVONMISES_FIT_L1_H__
#define __MVVONMISES_FIT_L1_H__


#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "lbfgsb.h"
#include "mvvonmises_fit.h"
#include "mvvonmises_likelihood.h"
#include "circularMvStats.h"

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
 * @param [in] Penparam
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
double mvvonmises_lbfgs_fit_l1(int p, double* mu, double* kappa, double* lambda, int n, double *samples, double penparam,
                            int verbose, double prec, double tol, int mprec, double *lower, double *upper, int *bounded );



/****
 *
 */
void __R_mvvonmises_lbfgs_fit_l1(int* p, double* mu, double* kappa, double* lambda, int* n, double *samples, double* penparam,
                            int* verbose, double* prec, double* tol, int* mprec, double *lower, double *upper, int *bounded , double *loss);

#endif
