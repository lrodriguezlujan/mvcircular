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
#include "mvvonmises_likelihood_legacy.h"
#include "mvvonmises_likelihood.h"
#include "matrixOps.h"

/*************
 * Likelihood
 */


/*****************************************************
 *
 * Old version
 *
 */

/**
 * Computes j-nth marginal mean given a sample and distribution parameters
 *
 * @param [in] p
 * @param [in] j
 * @param [in] mu
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta_i
 *
 * @return marginal mu for jth variable
 */
double marginal_j_mean(int p, int j, double* mu, double* kappa, double* lambda, double *theta_i){

    // iterator
    int l;

    // Sum term
    double acum =0.;

    for(l=0;l<p;l++){
        if(l!=j){
            acum+=lambda[j*p+l]*sin(theta_i[l]-mu[l]);
        }
    }

    return mu[j] + atan(acum/kappa[j]);
}

/**
 * Computes j-nth marginal concentration given a sample and distribution parameters
 *
 * @param [in] p
 * @param [in] j
 * @param [in] mu
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta_i
 *
 * @return marginal concentration for jth variable
 */
double marginal_j_concentration(int p, int j, double* mu, double* kappa, double* lambda, double *theta_i){

    // iterator
    int l;

    // Sum term
    double acum =0.;

    for(l=0;l<p;l++){
        if(l!=j){
            acum+=lambda[j*p+l]*sin(theta_i[l]-mu[l]);
        }
    }
    return sqrt(kappa[j]*kappa[j] + acum*acum);
}

/*************
 * Likelihood
 */

/**
 * Computes log pseudolikelihood for multivariate von mises
 *
 * @param [in] n
 * @param [in] p
 * @param [in] mu
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta
 *
 * @return negative log pseudolikelihood
 */
double mv_vonmises_logpseudolikelihood_naive(int n,int p, double* mu, double* kappa, double* lambda, double *theta){

    // iterator
    int j,i;

    // Sum term
    double acum =0.;

    // Auxiliar kappa
    double k_rest,m_rest;

    // Double sum
    for(j=0;j<p;j++){
        for(i=0;i<n;i++){
            k_rest = marginal_j_concentration(p,j,mu,kappa,lambda,theta+i*p);
            m_rest = marginal_j_mean(p,j,mu,kappa,lambda,theta+i*p);
            acum+= k_rest* cos(theta[i*p + j] - m_rest) - log(bessi0(k_rest));
        }
    }

    //return n*p*log(2*M_PI) - acum ;
    return -acum ;
}

double vm_loglikelihood_kappa_partial_naive(int n, int p, int R, double *mu, double *kappa, double *lambda, double *theta){

    // auxiliar marginals
    double k_rest, m_rest;

    // sum
    double acum=0.,acum_scnd;

    //iterator
    int i,l;

    for(i=0;i<n;i++){
        // compute marginals
        k_rest = marginal_j_concentration(p,R,mu,kappa,lambda,theta+i*p);
        m_rest = marginal_j_mean(p,R,mu,kappa,lambda,theta+i*p);

        //first term
        acum+=(cos(theta[i*p+R]-m_rest)-(bessi1_spd(k_rest)/bessi0_spd(k_rest)))/k_rest;

        // second acum
        acum_scnd=0.;
        for(l=0;l<p;l++){
            if(l!=R){
                acum_scnd+=lambda[R*p + l]*sin(theta[i*p+l]-mu[l]);
            }
        }
        acum+= (sin(theta[i*p+R]-m_rest)*acum_scnd)/k_rest;
    }

    return -acum*kappa[R];
}


/**
 * Partial derivative of the marginal mean wrt Lambda_RS 
 *
 * @param [in] p
 * @param [in] j
 * @param [in] R
 * @param [in] S
 * @param [in] mu
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta_i
 *
 * @return partial deriv
 */

double vm_marginalMu_lambda_partial_naive(int p,int j, int R, int S, double *mu, double *kappa, double *lambda, double *theta_i){

    // iterator
    int l;

    // accumulator
    double acum=0.;

    if( j != R ) return 0;
    else{
        for(l=0;l<p;l++){
            if(l!=j){
                acum+= lambda[j*p+l]*sin(theta_i[l]-mu[l]);
            }
        }
        acum/=kappa[j];
        acum = (1+acum*acum)*kappa[j];
        return sin(theta_i[S]-mu[S])/acum;
    }
}

/**
 * Partial derivative of the marginal concentration wrt Lambda_RS 
 *
 * @param [in] p
 * @param [in] j
 * @param [in] R
 * @param [in] S
 * @param [in] mu
 * @param [in] kappa
 * @param [in] lambda
 * @param [in] theta_i
 *
 * @return partial deriv
 */
double vm_marginalKapa_lambda_partial_naive(int p,int j, int R, int S, double *mu, double *kappa, double *lambda, double *theta_i){

    // iterator
    int l;

    // accumulator
    double acum=0.;

    if( j != R ) return 0;
    else{
        for(l=0;l<p;l++){
            if(l!=j){
                acum+= lambda[R*p+l]*sin(theta_i[l]-mu[l]);
            }
        }
        return acum*sin(theta_i[S]-mu[S])/marginal_j_concentration(p,R,mu,kappa,lambda,theta_i);

    }
}

double vm_loglikelihood_lambda_partial_naive(int n, int p, int R,int S, double *mu, double *kappa, double *lambda, double *theta){

    //iterators
    int j,i;

    // Accum
    double acum=0.;

    // auxiliar
    double m_rest,k_rest;

    for(j=0;j<p;j++){
        if(R==j){
            for(i=0;i<n;i++){
                // compute marginals
                k_rest = marginal_j_concentration(p,j,mu,kappa,lambda,theta+i*p);
                m_rest = marginal_j_mean(p,j,mu,kappa,lambda,theta+i*p);

                // First term
                acum+=(cos(theta[i*p+j] - m_rest)-(bessi1_spd(k_rest)/bessi0_spd(k_rest)) )* vm_marginalKapa_lambda_partial_naive(p,j,R,S,mu,kappa,lambda,theta+i*p);

                // Second term
                acum+=k_rest*sin(theta[i*p+j]-m_rest)*vm_marginalMu_lambda_partial_naive(p,j,R,S,mu,kappa,lambda,theta+i*p);
            }
        }
    }
    return -acum;
}
