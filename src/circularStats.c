/*******************
 * Circular statistics
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <math.h>
#include "circularStats.h"

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
void circularMeanLength(int n,int p, int j, double* data, double* R){

    // Accumulators
    double c,s;

    // Iterators
    int i;

    // Iterate over data
    c=s=0.0;
    for(i=0;i<n;i++){
        c+=cos(data[i*p+j]);
        s+=sin(data[i*p+j]);
    }
    // R^2 = (1/n Sum(C) )^2 + (1/n Sum(S))^2) 
    *R = sqrt(c*c + s*s)/n;
}


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
void circularVariance(int n,int p, int j, double* data, double* var){
    //Compute mean resultant length
    circularMeanLength(n,p,j,data,var);
    // V = 1-R
    *var = 1 - *var;
}
