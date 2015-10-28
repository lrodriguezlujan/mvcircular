/*******************
 * Circular multivariate distance
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "circularDistance.h"

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
void circularGeometricMedian(int n, int p, double *angles, double eps, int maxiter, char* maxiter_flag, double *median){

    // Auxiliar memory
    double *last = calloc(sizeof(double),p);

    // Sample iterator
    int i,j;
    // trial iterator
    int iter;

    // Auxiliar point to median distance
    double auxDist;
    // Auxiliar aggregate
    double aggr;

    // EPS control
    double improv = DBL_MAX ;

    // Weisfelz iterative algorithm
    iter=0;
    while (improv > eps && maxiter > iter){
        // iter + 1
        iter++;
        // reset Aggregate
        aggr=0;
        // Copy current to last
        memcpy(last,median,sizeof(double)*p);
        // Reset current
        memset(median,0,sizeof(double)*p);
        // Sample by sample
        for(i=0;i<n;i++){

            // Candidate to sample dist
            auxDist = circularMvDistance(p,last,angles+(i*p));
            if(auxDist!=0){
                aggr+=1/auxDist;
                // Update current
                for(j=0;j<p;j++) median[j] += angles[(i*p)+j]/auxDist;
            }
        }

        // aggr div
        for(j=0;j<p;j++) {
            median[j] /= aggr;
            median[j]=fmod(median[j],2*M_PI);
        }

        // Compute eps
        improv = circularMvDistance(p,last,median);
    }

    // Free last and return current
    free(last);
    if(maxiter_flag!=NULL) *maxiter_flag = maxiter>iter;
    return(median);
}

void _R__circularGeometricMedian(int *n, int *p, double *angles, double *eps, int *maxiter, int* maxiter_flag, double *median){
    circularGeometricMedian(*n,*p, angles, *eps, *maxiter, *maxiter_flag, median);
}
