/*******************
 * Linear statistics
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#ifndef _LINEARSTATS_H__
#define _LINEARSTATS_H__

/**
 * Computes the sample linear mean length
 * 
 * @param [in] n number of samples
 * @param [in] p Number of variables (set to 1 if data is a vector)
 * @param [in] j Component (set to 0 if data is a vector)
 * @param [in] data nxp matrix of samples
 * @param [out] mean linear mean
 *
 */
void linearMean(int n,int p, int j, double* data, double* mean);

/**
 * Computes the sample linear variancea
 * 
 * @param [in] n number of samples
 * @param [in] p Number of variables (set to 1 if data is a vector)
 * @param [in] j Component (set to 0 if data is a vector)
 * @param [in] data nxp matrix of samples
 * @param [out] var Linear variance (sigma)
 *
 */
void linearVariance(int n,int p, int j, double* data, double* var);

#endif
