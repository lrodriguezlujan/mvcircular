/*******************
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: mvvonmises.c 13 2015-04-17 08:17:47Z lrodriguez $
 **/

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include "bessel.h"
#include "mvvonmises_likelihood.h"
#include "matrixOps.h"

/*************
 * Likelihood
 */

/* Transforms theta to Sin(theta-mu)
 *
 * Computes Sin(theta-mu) / Cos(theta-mu) for each variable in the dataset.
 *
 * @param [in] n
 * @param [in] p
 * @param [in] theta
 * @param [in] mu
 * @param [out] S
 * @param [out] C
 *
 */
void mv_theta_cos_sinTransform(int n, int p, double *theta, double* mu, long double *S, long double *C){
    
    int i,j;

    for(i=0;i<n;i++){
        for(j=0;j<p;j++){
            S[i*p + j] = sinl(theta[i*p + j] - mu[j]);
            C[i*p + j] = cosl(theta[i*p + j] - mu[j]);
        }
    }
}


/***
 * Multivaraite von mises density. Not normalized
 *
 * @param [in] p
 * @param [in] theta
 * @param [in] mu
 * @param [in] kappa
 * @param [in] Lambda
 *
 * @return Density val
 * 
 */
double mvvm_density_unregularized(int p, double *theta, double *mu, double *kappa, double *lambda){



    double aux=0;
    long double *s = malloc(sizeof(long double) * p);
    long double *v = malloc(sizeof(long double) * p);
    int i;

    for(i=0;i<p;i++) s[i] = sin(theta[i]-mu[i]);

    // Bilinear form
    generalMatrixMatrix(1,p,p,p,'n','n',s,lambda,v);
    generalMatrixMatrix(1,p,p,1,'n','n',v,s,&aux);
    aux /= 2;

    // Cos prod
    for(i=0;i<p;i++)
        aux += kappa[i] * cos(theta[i]-mu[i]);

    return (exp(aux));
}

/*****
 * Optimized version. Computes simplified log likelihood and its derivatives
 *
 * Removes constant terms. Uses matrix multiplications. Reduce number of trigonometric opers.
 * Reuses terms. Marginal mu is not longer needed, derivatives are simplified as well. 
 * 
 * Instead of PL we will name it "loss function"
 *
 * Lambda is taken as a full row-leading pxp matrix with 0's in the diagonal. But d_lambda is only
 * computed for the upper triangle.
 *
 * If d_kappa or d_lambda are NULL, derivatives are not computed.
 *
 * To even speed up the process a little bit more, we make the following assumptions:
 *
 *  Mu is 0 for every j
 *  (input) Theta_i_j = Sin( (real) Theta_i_j - (real)mu_j )
 *
 *
 * @param [in] n
 * @param [in] p
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta
 * @param [out] ro
 * @param [out] d_kappa
 * @param [out] d_lambda
 *
 * @return loss function value
 */
long double mv_vonmises_lossFunction(int n,int p, double* kappa, double* lambda, long double *S, long double *C, long double *ro, double *d_kappa, double *d_lambda){

    // iterators
    int i,j,s;

    // Auxiliar terms
    long double fx;    // Loss function value
    long double k_i_j; // Marginal concentration por i,j
    long double bes_kij; // Bessel I0(k_ij)
    long double bes1_kij; // I1(k_ij)
    long double a_kij; // A_0(k_j_i)/k_j_i

    // Compute ro as ro = Theta * Lambda (Lambda is symmetric)
    generalMatrixMatrix(n,p,p,p,'n','n',S,lambda,ro);
    
    // Initialize
    fx=0;

    
    if(d_lambda != NULL && d_kappa != NULL){

        // Initialize d_lambda d_kappa
        memset(d_lambda,0,sizeof(double)*p*p);
        memset(d_kappa,0,sizeof(double)*p);

        for(j=0;j<p;j++){

            for(i=0;i<n;i++){

                // Compute k_i_j and a_kij
                k_i_j = sqrtl( (kappa[j] * kappa[j]) + (ro[i*p + j] * ro[i*p + j]) );

                // Bessel kij (Expon scaled)
                // We get exp(-k_i_j) * I0(k_i_j)
                bes_kij = bessi0_expon(k_i_j);
                
                // A0/k_i_j
                a_kij = ((bessi1_expon(k_i_j)) / (bes_kij * k_i_j)); // Actually this is A0(k_j_i)/ k_j_i but getting rid of one div.

                // Loss function part
                //
                // logl(exp(-kij)*I0) = -kij + logl(I0) -> We add k_ij to get logl(I0)
                fx += (k_i_j + logl(bes_kij) ) - ( C[i*p + j] * kappa[j] ) - ( S[i*p + j] * ro[i*p + j] );

                // Sum the part that coresponds to d_kappa
                d_kappa[j] += a_kij * kappa[j] - C[i*p + j];
               
                // Update lambda partial
                for(s=0;s<j;s++)
                    d_lambda[s*p + j] += S[i*p +s] * (a_kij * ro[i*p + j] - S[i*p + j] ); 

                for(s=j+1;s<p;s++)
                    d_lambda[j*p + s] += S[i*p +s] * (a_kij * ro[i*p + j] - S[i*p + j] ); 

            }  
            // Copy to full matrix
            for(s=j+1;s<p;s++){
                   d_lambda[s*p + j] = d_lambda[j*p + s];
            }
        }
    }
    else{
        for(j=0;j<p;j++){
            for(i=0;i<n;i++){

                // Compute k_i_j and a_kij
                k_i_j = sqrtl( (kappa[j] * kappa[j]) + (ro[i*p + j] * ro[i*p + j]) );

                // Bessel kij (Expon scaled)
                bes_kij = bessi0_expon(k_i_j);

                // Loss function part
                fx += (k_i_j + logl(bes_kij)) - (C[i*p + j] * kappa[j]) - (S[i*p + j] * ro[i*p + j]);
            }
        }
    }
    
    return(fx);
}
