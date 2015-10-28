/*******************
 * Circular multivariate statistics
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/


#ifndef _CIRCULARMVSTATS_H__
#define _CIRCULARMVSTATS_H__

/**
 * Computes the circular mean from p-variate n samples of data
 * 
 * @param [in] n number of samples
 * @param [in] p number of variables
 * @param [in] data nxp matrix of samples
 * @param [out] mean p variate mean vector
 *
 * @returns vector of length p with the marginal circular means
 */

void multiCircularMean(int n,int p, double* data, double *mean);

#endif
