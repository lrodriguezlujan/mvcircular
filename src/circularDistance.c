/**
 * Compute circular distance between two angular multivariate points
 *
 * @param [in] p Number of variables
 * @param [in] a First point
 * @param [in] b Second point
 *
 * @returns Angular distance (from 0 to sqrt(d)*pi)
 */

#include <math.h>
#include <stdint.h>

double circularMvDistance(int p, double* a, double *b ){

    // Dim iterator
    uint_fast16_t i;

    // distance
    double dist = 0;

    // Move b point to compute minimum distance. This avoids to cross "edges"
    for(i=0;i<p;i++){
        // If not in the same half
        if(fabs(a[i]-b[i]) > M_PI){        
                
                // If a is in the right segment -> move b right
                if(a[i] > M_PI)  dist += (a[i] - b[i] -2*M_PI) *  (a[i] - b[i] -2*M_PI);
                // Else move b left
                else dist += (a[i] - b[i] + 2*M_PI) *  (a[i] - b[i] + 2*M_PI);
        }
        else
            dist += (a[i] - b[i]) *  (a[i] - b[i]);
    }
    
    // sqroot 
    return(sqrt(dist));
}

