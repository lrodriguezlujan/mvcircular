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
double mv_vonmises_logpseudolikelihood_nonopt(int n,int p, double* mu, double* kappa, double* lambda, double *theta){

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

double vm_loglikelihood_kappa_deriv_old(int n, int p, int R, double *mu, double *kappa, double *lambda, double *theta){

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

double vm_marginalMu_lambda_deriv_old(int p,int j, int R, int S, double *mu, double *kappa, double *lambda, double *theta_i){

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
double vm_marginalKapa_lambda_deriv(int p,int j, int R, int S, double *mu, double *kappa, double *lambda, double *theta_i){

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

double vm_loglikelihood_lambda_deriv_old(int n, int p, int R,int S, double *mu, double *kappa, double *lambda, double *theta){

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
                acum+=(cos(theta[i*p+j] - m_rest)-(bessi1_spd(k_rest)/bessi0_spd(k_rest)) )* vm_marginalKapa_lambda_deriv(p,j,R,S,mu,kappa,lambda,theta+i*p);

                // Second term
                acum+=k_rest*sin(theta[i*p+j]-m_rest)*vm_marginalMu_lambda_deriv_old(p,j,R,S,mu,kappa,lambda,theta+i*p);
            }
        }
    }
    return -acum;
}
