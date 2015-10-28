/*******************
 * Sampler for the univariate von Mises distribution
 * Best and Fisher 1979
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: vonmises_sampler.c 19 2015-05-11 15:05:35Z lrodriguez $
 **/

#include <math.h>
#include "vonmises_sampler.h"
#include "bessel.h"
#include "prng.h"


/**
 * Von mises density function. Informative
 *
 * @param [in] theta
 * @param [in] mu
 * @param [in] kappa
 *  
 * @return density
 */
double vmDensity(double theta, double mu, double kappa){
    return exp(kappa * cos(theta-mu)) / (2 * M_PI * bessi0_spd(kappa));
}

/**
 * Rejection-based sampler for the von Mises distribution. Retruns several samples per call
 *
 *  @param [out] x Real ndim vector to return the samples
 *  @param [in] n Number of samples to generate
 *  @param [in] mu Mu parameter of the vm dist
 *  @param [in] kappa Kappa parameter of the vm dist
 */
void rvm(double *x, int n, double mu, double kappa){
 
    int i;
    double U2;
    register double r, f, c;

    /* SET r */
    r = 1 + sqrt(1+4*kappa*kappa);
    r = (r - sqrt(2*r))/(2*kappa);
    // At this point, r is ro as in (Best and Fisher 1979)
    r = (1 + r*r)/(2*r); // r = (1+\rho^2)/(2*\rho)

    i = 0;
    /** Fill N Samples */
    do{
        /* Set F and C */
        // Step 1
        // f = (1 + r cos(\pi u_1)) / (r+cos(\pi u_1)
        f = cos( M_PI * ((xorshift4096()*1.0)/UINT64_MAX)  );
        f = (1. + r * f)/(r + f); 

        //  c = k * (r - f)
        c = kappa * (r - f);

        // Step 2 & 3
        U2 = (xorshift4096()*1.0)/UINT64_MAX;

        //Rejection
        if( (c * (2 - c) - U2 > 0) || ( (log(c/U2) + 1 - c ) >= 0.)) {
            if ( (xorshift4096()*1.0)/UINT64_MAX > 0.50) x[i] = acos(f) + mu;
            else x[i] = -acos(f) + mu;
            i++;
        }

    } while(i < n);
}

/**
 * Rejection-based sampler for the von Mises distribution. Retruns several samples per call
 *
 *  @param [out] x Real ndim vector to return the samples
 *  @param [in] n Number of samples to generate
 *  @param [in] mu Mu parameter of the vm dist
 *  @param [in] kappa Kappa parameter of the vm dist
 *  @param [in] seeds Random seeds
 */
void rvm_r(double *x, int n, double mu, double kappa, uint64_t *seeds){
 
    int i;
    double U2;
    register double r, f, c;

    /* SET r */
    r = 1 + sqrt(1+4*kappa*kappa);
    r = (r - sqrt(2*r))/(2*kappa);
    // At this point, r is ro as in (Best and Fisher 1979)
    r = (1 + r*r)/(2*r); // r = (1+\rho^2)/(2*\rho)

    i = 0;
    /** Fill N Samples */
    do{
        /* Set F and C */
        // Step 1
        // f = (1 + r cos(\pi u_1)) / (r+cos(\pi u_1)
        f = cos( M_PI * ((xorshift4096()*1.0)/UINT64_MAX)  );
        f = (1. + r * f)/(r + f); 

        //  c = k * (r - f)
        c = kappa * (r - f);

        // Step 2 & 3
        U2 = (xorshift4096()*1.0)/UINT64_MAX;

        //Rejection
        if( (c * (2 - c) - U2 > 0) || ( (log(c/U2) + 1 - c ) >= 0.)) {
            if ( (xorshift4096()*1.0)/UINT64_MAX > 0.50) x[i] = acos(f) + mu;
            else x[i] = -acos(f) + mu;
            i++;
        }

    } while(i < n);
}

/******
 * Rejection-based sampler for the von Mises distribution. Returns one sample per call
 *
 *  @param [in] mu Mu parameter of the vm dist
 *  @param [in] kappa Kappa parameter of the vm dist*
 *
 *  @return one sample
 * *****/
double rvm_s(double mu, double kappa){

    double U2;
    double r,f,c;

    r = 1 + sqrt(1+4*kappa*kappa);
    r = (r - sqrt(2*r))/(2*kappa);
    r = (1 + r*r)/(2*r);

    do{
        /* Set F and C */
        f = cos( M_PI * ((xorshift4096()*1.0)/UINT64_MAX)  );
        f = (1. + r * f)/(r + f);
        c = kappa * (r - f);

        U2 = (xorshift4096()*1.0)/UINT64_MAX;
        if( (c * (2 - c) - U2 > 0) || ( (log(c/U2) + 1 - c ) >= 0.)) {
            if ( (xorshift4096()*1.0)/UINT64_MAX > 0.50) return (acos(f) + mu);
            else return (-acos(f) + mu);
        }
    } while(1);
}

double rvm_s_r(double mu, double kappa, uint64_t* seeds){

    double U2;
    double r,f,c;

    r = 1 + sqrt(1+4*kappa*kappa);
    r = (r - sqrt(2*r))/(2*kappa);
    r = (1 + r*r)/(2*r);

    do{
        /* Set F and C */
        f = cos( M_PI * ((xorshift4096_r(seeds)*1.0)/UINT64_MAX)  );
        f = (1. + r * f)/(r + f);
        c = kappa * (r - f);

        U2 = (xorshift4096_r(seeds)*1.0)/UINT64_MAX;
        if( (c * (2 - c) - U2 > 0) || ( (log(c/U2) + 1 - c ) >= 0.)) {
            if ( (xorshift4096_r(seeds)*1.0)/UINT64_MAX > 0.50) return (acos(f) + mu);
            else return (-acos(f) + mu);
        }

    } while(1);
}
