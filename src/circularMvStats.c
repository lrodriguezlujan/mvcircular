/*******************
 * Circular multivariate statistics
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <math.h>
#include "circularMvStats.h"

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
void multiCircularMean(int n,int p, double* data, double *mean){

    // Accumulators
    double c,s;
    // Iterators
    int i,j;

    // Iterate over data
    for(j=0;j<p;j++){
        // Initialize
        c=s=0.0;
        for(i=0;i<n;i++){
            c+=cos(data[i*p+j]);
            s+=sin(data[i*p+j]);
        }
        // Mean angle is arctan of mean sin mean cos
        mean[j]=atan2(s/n,c/n);
        //Adjust to [0,2pi]
        if(mean[j]<0) mean[j]+=2*M_PI;
    }
}
