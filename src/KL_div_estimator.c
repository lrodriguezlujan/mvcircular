/*******************
 * KL divergence estimator
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <math.h>
#include <float.h>
#include <stdint.h>
#include "circularDistance.h"
#include "KL_div_estimator.h"

/* Computes distance to the 2-nd nearest neighbor within the set
 * based on indexed. (excluded himself=
 *
 */
double k2nn_dist_inSet(int p, int idx, int n, double *set){
    
    uint_fast16_t i,j;

    double lower=-1;
    double dist=0;
    double v;

    //For the first k points we search for the maximum (i.e the 2-th )
    for(i=0,j=0; j<2; i++){
        if(i!=idx){
            j++ ; 
            v = circularMvDistance(p,set + idx*p , set + i*p);
            if(v > dist){
                lower = dist;
                dist = v;
            }
            else if (v>lower)
                lower = v;
        }
    }

    // Then we just update dist and lower
    for(;i<n;i++){
        if(i!=idx){
            v = circularMvDistance(p,set+idx*p,set+i*p);
            if(v>lower){
                if( (v < dist) )
                    dist = v;
            }
            else{
                dist=lower;
                lower=v;
            }
        }
    }

    // Then return dist
    return dist;
}


/* Computes distance to the 2-nd nearest neighbor in the set
 *
 *
 */
double k2nn_dist_point(int p, double *a , int n, double *set){

    uint_fast16_t i,j;

    double lower=-1;
    double dist=0;
    double v;

    //For the first k points we search for the maximum (i.e the k-th )
    lower = circularMvDistance(p,a,set); // Take distance to the first point as lower
    v = circularMvDistance(p,a,set+p);
    if(v<lower){
        dist = lower;
        lower = v;
    }
    else{
        dist = v;
    }
    for(i=0,j=0; j<2; i++){
        v = circularMvDistance(p, a , set + i*p);
        if(v>0){
            j++ ; 
            if(v > dist){
                lower = dist;
                dist = v;
            }
            else if (v>lower)
                lower = v;
        }
    }

    // Then we just update dist and lower
    for(i=2;i<n;i++){
            v = circularMvDistance(p, a, set + i*p);
            if(v>0){
                if(v>lower){
                    if( (v < dist) )
                        dist = v;
                }
                else{
                    dist=lower;
                    lower=v;
                }
            }
    }

    // Then return dist
    return dist;
}

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
                                         double* B){

    // Iterator (over A set)
    int i;

    // Sum
    double acum;

    // Dobule r and s quantities
    double r_k,s_k;

    // Check that we have enough samples to run the test
    if( (n < k+1) || (m<k) )
        return DBL_MAX;

    // At the moment only take k=2 (eff. purposses)
    if(k!=2) return 0;

    acum=0;
    for(i=0;i<n;i++){
        // Get k-nearest point for rk
        r_k = k2nn_dist_inSet(p, i, n, A);

        // Get k-nearest point for sk
        s_k = k2nn_dist_point(p, A+(i*p) , m, B);

        // Compute sum
        acum+=log(s_k/r_k);
    }
    return ((1.0*p)/n) * acum + log(m/(n-1.0));
}

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
                                         double *kl){
    *kl=KL_estimator_circulardist(*p,*k,*n,A,*m,B);
}

