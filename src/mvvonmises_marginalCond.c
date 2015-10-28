/*******************
 * Set of functions related to a multivariate von Mises distribution
 * Includes pseudolikelihood calculations and functions to calculate 
 * the parameters of the conditional distributions.
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: mvvonmises.c 19 2015-05-11 15:05:35Z lrodriguez $
 **/

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>

/**
 * Computes univariate conditional mu coefficient of mvmdist for the variable j based on the ith observation
 *
 * @param [in] size Number of variable
 * @param [in] j  Index of the variable
 * @param [in] instance Sample of data
 * @param [in] mu p-vector with mu parameters of the mvm dist.
 * @param [in] kappa p-vector with kappa parameters of the mvm dist.
 * @param [in] lambda pxp matrix with lambda parameters
 *
 * @returns univariate conditional mu for j-th variable
 */

double condMu(int size,int j,double* instance,
                    double* mu, double* kappa,double* lambda){

    // Accumulator
    double sum = 0.0;

    // Iterator
    int l;

    // Sum_j!=l lambda_jl * sin(theta_li - mu_l)
    for(l=0;l<size;l++){
        // lambda[i,i] == 0
        sum += lambda[l*size+j]*sin(instance[l] - mu[l]);
    }

    // mu_j + atan(sum/kappa_j)
    return mu[j] + atan(sum/kappa[j]);
}

/**
 * Computes univariate conditional kappa coefficient of mvmdist for the variable
 * j based on the i-th observation (instance)
 *
 * @param [in] size Number of variable
 * @param [in] j  Index of the variable
 * @param [in] instance Sample of data
 * @param [in] mu p-vector with mu parameters of the mvm dist.
 * @param [in] kappa p-vector with kappa parameters of the mvm dist.
 * @param [in] lambda pxp matrix with lambda parameters

 */

double condKappa(int size,int j, double* instance,
                    double* mu, double* kappa, double* lambda){

    // Accumulator
    double sum = 0.0;

    // Iterator
    int l;

    // sum_l!=j lambda_jl*sin(theta_li - mu_l)
    for(l=0;l<size;l++){
        // lambda[i,i] == 0
        sum+= lambda[l*size+j]* sin( instance[l] - mu[l] );
    }

    // return (kappa_j ^2 + sum^2)^(1/2)
    return sqrt( ( kappa[j] * kappa[j] ) + ( sum * sum ) );
}

/**
 * Computes both mu and kappa marginal conditional coeficients of mv von Mises given the
 * ovservation instance.
 *
 * @param [in] size Number of variable
 * @param [in] j  Index of the variable
 * @param [in] instance Sample of data
 * @param [in] mu p-vector with mu parameters of the mvm dist.
 * @param [in] kappa p-vector with kappa parameters of the mvm dist.
 * @param [in] lambda pxp matrix with lambda parameters
 * @param [out] muCond Conditional mu value
 * @param [out] kappaCond Conditional Kappa value
 */
void condKappaMu(int size, int j, double* instance, double *mu, double *kappa,
                        double *lambda, double *muCond, double *kappaCond){

    // Accumulator
    double sum = 0.0;

    // Iterator
    int l;

    // sum_l!=j lambda_jl*sin(theta_li - mu_l)
    for(l=0;l<size;l++){
        // lambda[i,i] == 0
        sum+= lambda[l*size+j]* sin( instance[l] - mu[l] );
    }
    
    // Values
    *kappaCond = sqrt( ( kappa[j] * kappa[j] ) + ( sum * sum ) );
    *muCond = mu[j] + atan(sum/kappa[j]);
}
