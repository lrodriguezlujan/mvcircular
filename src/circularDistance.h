/*******************
 * Circular multivariate distance
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/


#ifndef _CIRCULARMVDIST_H__
#define _CIRCULARMVDIST_H__


/**
 * Compute circular distance between two angular multivariate points
 *
 * @param [in] p Number of variables
 * @param [in] a First point
 * @param [in] b Second point
 *
 * @returns Angular distance (from 0 to sqrt(d)*pi)
 */
double circularMvDistance(int p, double* a, double *b );

#endif
