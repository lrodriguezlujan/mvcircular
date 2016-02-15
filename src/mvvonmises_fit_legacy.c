/*******************
 * Module to fit a mVM distribution from data samples
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id:$
 **/
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "lbfgsb.h"
#include "mvvonmises_fit_legacy.h"
#include "mvvonmises_fit.h"
#include "mvvonmises_likelihood_legacy.h"
#include "circularMvStats.h"


double mvvonmises_lbfgs_fit_naive(int p, double* mu, double* kappa, double* lambda, int n, double *samples,int verbose, double prec, double tol, int mprec, double *lower, double *upper, int *bounded ){

    uint_fast16_t i,j,k;

    /* Von Mises computation parameters*/
    long double *S,*C,*ro;
    double *d_kappa,*d_lambda;

    S = malloc(sizeof(long double)*n*p*3);
    C = S + (n*p);
    ro = C + (n*p);

    /* Set up instance */
    multiCircularMean(n,p,samples,mu);

    // Do the theta transformation
    // mv_theta_cos_sinTransform(n,p,samples,mu,S,C);

    /* Derivatives */
    d_kappa = malloc(sizeof(double)*(p+(p*p)));
    d_lambda = d_kappa + p;


        /* LBFGS variables */
    /* integer and logical types come from lbfgsb.h */

    /* LBFGS static parameters - To be included as inputs! */
    integer iprint = verbose; // No output
    double  factr = prec; // Moderate prec
    double  pgtol = tol; // Gradient tolerance
    integer m = mprec;       // Number of corrections

    /* Fixd workspaces */
    integer taskValue;
    integer *task = &taskValue;
    integer csaveValue;
    integer *csave = &csaveValue;
    integer isave[44];
    double  dsave[29];
    logical lsave[4];

    /* Dynamic parameters (Given, but we might have to cast)*/
    integer nvar = ((p*p)-p)/2 + p;  // Number of variables
    double  f = DBL_MAX;   // Eval value
    double *g = calloc(sizeof(double),nvar*4); // Gradient value
    double *l = g+nvar; // lower bounds
    double *u = l+nvar; // upper bounds
    double *x = u+nvar; // Point value
    integer *nbd = calloc(sizeof(integer),nvar);

    /* Dynamic Workspaces*/
    double  *wa  = calloc(sizeof(double),( (2*m + 5)*nvar + 11 * m * m + 8 * m));
    integer *iwa = calloc(sizeof(integer) , 3 * nvar );

    /* Copy parameters (this way so casting is done)*/
    for(i=0;i<nvar; i++){
        l[i]   = lower[i] ;
        u[i]   = upper[i];
        nbd[i] = bounded[i];
    }


    /* Copy initial values to X*/
    // Kappa values
    for(i=0;i<p;i++){
        x[i] = kappa[i];
    }

    // Lambda values
    for(i=0, k=0; i < (p-1) ; i++) {
        for( j=(i+1) ; j < p ; j++, k++){
            x[p+k] = lambda[ i*p + j];
        }
    } 
   
    /* Set task to START*/
    *task = (integer)START;

    do{
        setulb(&nvar,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave);

        // If F and G comp. is requiered
        if( IS_FG(*task) ){

            // Copy kappa
            memcpy(kappa,x,sizeof(double)*p);

            // Copy lambda
            matrixUpperToFull(p,x+p,lambda);

            // Call loss function with current x
            f = mv_vonmises_logpseudolikelihood_naive(n,p,mu,kappa,lambda,samples);
            
            // kappa
            for(i=0;i<p;i++) d_kappa[i] = vm_loglikelihood_kappa_partial_naive(n,p,i,mu,kappa,lambda,samples);

            // Approx DF for lambda
            for(i=0;i<(p-1);i++){
                    for(j=i+1;j<p;j++){
                        d_lambda[i*p + j] = vm_loglikelihood_lambda_partial_naive(n,p,i,j,mu,kappa,lambda,samples);
                        d_lambda[j*p + i] = d_lambda[i*p + j];
                    }
            }

            // Kappa partials
            memcpy(g,d_kappa,sizeof(double)*p);

            // Lambda partials (to vector)
            for(i=0, k=0; i < (p-1) ; i++) {
                for( j=(i+1) ; j < p ; j++, k++){
                    g[p+k] = d_lambda[ i*p + j];
                }
            }
        }

        // Additional stopping criteria to avoid infinite looping

        else if( *task == NEW_X ){

            /* Control number of iterations */
            if(isave[33] >= 100 ){
                *task = STOP_ITER;
            }
            /* Terminate if |proj g| / (1+|f|) < 1E-10 */
            else if ( dsave[12] <= (fabs(f) + 1) * 1E-10 ){
                *task = STOP_GRAD;
            }
        }

    }while( IS_FG(*task) || *task==NEW_X);


    /* Copy back x */
    for(i=0;i<p;i++) kappa[i]=x[i];
    matrixUpperToFull(p,x+p,lambda);

    /** FreE*/

    free(S); free(d_kappa); free(g);
    free(wa); free(iwa); free(nbd);

    if( IS_ERROR(*task) || IS_WARNING(*task) )
        return NAN;
    else
        return (double)f;
}

void __R_mvvonmises_lbfgs_fit_naive(int* p, double* mu, double* kappa, double* lambda, int* n, double *samples, int* verbose, double* prec, double* tol, int* mprec, double *lower, double *upper, int *bounded , double *loss){
    *loss = mvvonmises_lbfgs_fit_naive(*p, mu, kappa, lambda, *n, samples, *verbose, *prec, *tol, *mprec, lower, upper, bounded );
}
