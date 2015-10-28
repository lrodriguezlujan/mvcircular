/*******************
 * Linear statistics
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <math.h>
#include "linearStats.h"

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
void linearMean(int n,int p, int j, double* data, double* mean){

    // Accumulators
    double temp;

    // Iterators
    int i;

    // Iterate over data
    temp=0.0;
    for(i=0;i<n;i++){
        temp+=data[i*p+j];
    }
    // R^2 = (1/n Sum(C) )^2 + (1/n Sum(S))^2) 
    *mean = temp/n;
}


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
void linearVariance(int n,int p, int j, double* data, double* var){

    // Accumulator
    double temp;

    // mean
    double mu;

    //Compute mean
    circularMeanLength(n,p,j,data,&mu);

    // Iterators
    int i;

    // Iterate over data
    temp=0.0;
    for(i=0;i<n;i++){
        temp+=pow(data[i*p+j]-mu,2);
    }

    // unbiased
    *var = temp / (n-1);
}
