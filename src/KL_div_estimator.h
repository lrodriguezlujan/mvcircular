/*******************
 * KL divergence estimator
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/


#ifndef _KLDIV_ESTIMATOR_H__
#define _KLDIV_ESTIMATOR_H__

/*
 * Given two set of samples, computes approx. KL divergence
 *
 * computes divergence estimation by using "euclidean" distance bw points as described in perez-cruz(2008) (4)
 *
 * @param [in] p Number of variables
 * @param [in] n Number of samples in A
 * @param [in] A nxp Matrix in row leading order
 * @param [in] m Number of samples in B
 * @param [in] B mxp Matrix in row leading order
 *
 * @return D_k(P||Q)
 */
double KL_estimator_circulardist(        int p,
                                         int k,
                                         int n,
                                         double* A,
                                         int m,
                                         double* B);

/**
 *
 * R - Wrapper
 *
 */
void   __R_KL_estimator_circulardist(    int *p,
                                         int *k,
                                         int *n,
                                         double* A,
                                         int *m,
                                         double* B,
                                         double *kl);

#endif
