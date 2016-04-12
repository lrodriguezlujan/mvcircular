/*******************
 * Module to fit a mVM distribution from data samples
 * with L1 Penalization
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
#include "mvvonmises_fit_l1.h"
#include "mvvonmises_likelihood.h"
#include "circularMvStats.h"

/**
 * Maximum full pseudolikelihood method to estimate parameters.
 *
 * Minimization done by gradient L-BFGS with wolfe
 *
 * @param [in] p Number of variables
 * @param [in,out] mu pdim Real vector. Mu parameter of the mv-vm dist
 * @param [in,out] kappa pdim Real vector. Kappa parameter of the mv-vm dist
 * @param [in,out] lambda pxp Real matrix on row-leading order. Lambda parameter of the mv-vm dist
 * @param [in] n Number of samples
 * @param [in] samples
 * @param [in] Penparam
 * @param [in] verbose
 * @param [in] prec
 * @param [in] tol 
 * @param [in] mprec
 * @param [in] lower
 * @param [in] upper
 * @param [in] bounded
 *
 * @returns Natural logarithm of the pseudolikelihood of the fitted distribution (aprox)
 */
double mvvonmises_lbfgs_fit_l1(int p, double* mu, double* kappa, double* lambda, int n, double *samples, double penparam,
                            int verbose, double prec, double tol, int mprec, double *lower, double *upper, int *bounded ){

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
    mv_theta_cos_sinTransform(n,p,samples,mu,S,C);

    /* Derivatives */
    d_kappa = malloc(sizeof(double)*(p+(p*p)));
    d_lambda = d_kappa + p;


    /* LBFGS variables */
    /* integer and logical types come from lbfgsb.h */
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
            f = mv_vonmises_lossFunction(n,p,kappa,lambda,S,C,ro,d_kappa,d_lambda);

            // Kappa partials
            memcpy(g,d_kappa,sizeof(double)*p);

            // Apply penalization
            if( (penparam > 0) ){

                double l1pen = 0;

                // Compute l1
                 for (i=0; i < (p-1); i++) {
                    for (j=i+1; j < p; j++){
                        l1pen += abs( lambda[i * p + j]);
                    }
                 }
            
                 // We also need to include kappa
                 for (i=0; i<p; i++){
                        l1pen += kappa[ i ]++ ;
                 }
                 
                 // Add norm to F
                 f += penparam * l1pen;

                 if(penparam > 0){
                    // Modify lambda partials
                    for (i=0; i < (p-1); i++) {
                        for (j=i+1; j < p; j++){
                            // RECALL: P is actually diag(kappa) - lambda
                            d_lambda[i*p + j] += penparam * (lambda[i*p + j] > 0)?1:-1; 
                            d_lambda[j*p + i] = d_lambda[i*p + j];
                        }
                    }
    
                    // Modify kappa partials
                    for (i=0; i<p; i++){
                            d_kappa[i] += penparam ;
                     }
                }
            }

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



/****
 *
 */
void __R_mvvonmises_lbfgs_fit_l1(int* p, double* mu, double* kappa, double* lambda, int* n, double *samples, double* penparam,
                            int* verbose, double* prec, double* tol, int* mprec, double *lower, double *upper, int *bounded , double *loss){
    if((*penparam) > 0)
        *loss = mvvonmises_lbfgs_fit_l1(*p, mu, kappa, lambda, *n, samples, *penparam, *verbose, *prec, *tol, *mprec, lower, upper, bounded );
    else
        *loss = mvvonmises_lbfgs_fit(*p, mu, kappa, lambda, *n, samples, NULL, NULL, *verbose, *prec, *tol, *mprec, lower, upper, bounded );
}
