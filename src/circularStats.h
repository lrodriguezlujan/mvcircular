/*******************
 * Circular statistics
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/


#ifndef _CIRCULARSTATS_H__
#define _CIRCULARSTATS_H__

/**
 * Computes the circular mean length (R) from samples
 * 
 * @param [in] n number of samples
 * @param [in] p Number of variables (set to 1 if data is a vector)
 * @param [in] j Component (set to 0 if data is a vector)
 * @param [in] data nxp matrix of samples
 * @param [out] R mean resultant length
 *
 */
void circularMeanLength(int n,int p, int j, double* data, double* R);



/**
 * Computes the circular variance (V) from samples
 * 
 * @param [in] n number of samples
 * @param [in] p Number of variables (set to 1 if data is a vector)
 * @param [in] j Component (set to 0 if data is a vector)
 * @param [in] data nxp matrix of samples
 * @param [out] var Circular variance
 *
 */
void circularVariance(int n,int p, int j, double* data, double* var);


#endif
