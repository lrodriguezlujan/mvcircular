/*******************
 * Circular multivariate distance
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#ifndef _R__circular_geomedian
#define _R__circular_geomedian


/**
 * Angular multivariate geometric median
 *
 * Computed with Weisfelz algorithm -> Iterative schema
 * Possible improvements: faster evaluation via preprocessing.
 * See: Base et al 2003, chadrasekaran, tamir 1990
 *
 * @param [in] n Number of instances
 * @param [in] p Number of variables
 * @param [in] angles Angular data in row-order basis
 * @param [in] eps First stop condition (improvement threshold)
 * @param [in] maxiter Maximum number of iterations
 * @param [out] maxiter_flag If true (!=0) max number of iterations reached
 *
 * @returns p-dimensional vector with computed median
*/
void circularGeometricMedian(int n, int p, double *angles, double eps, char maxiter, int* maxiter_flag, double *median);

void _R__circularGeometricMedian(int *n, int *p, double *angles, double *eps, int *maxiter, int* maxiter_flag, double *median);

#endif
