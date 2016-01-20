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
#include "mvvonmises_fit.h"
#include "mvvonmises_likelihood.h"
#include "circularMvStats.h"

/** Auxiliar diagonal to full matrix */
void matrixDiagToFull(int p, double* diags, double*full){

    int row,col;
    int col_shift,l;
    int i,j; 

    for(i=0,col_shift=1,l=(p-1),j=0;j<(p*(p-1))/2;j++,i++){

            if(l==i){
                l--;
                i=0;
                col_shift++;
            }

            row = i;
            col = i+col_shift;
            full[row*p+col]=diags[j];
            full[col*p+row]=diags[j];
    }
}

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
 * @param [in] phi Prior matrix
 * @param [in] H Confidence matrix
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
double mvvonmises_lbfgs_fit(int p, double* mu, double* kappa, double* lambda, int n, double *samples, double* phi, double *H,
                            int verbose, double prec, double tol, int mprec, double *lower, double *upper, int *bounded ){

    uint_fast16_t i,j,l_int,col_shift;

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

#ifdef DEBUG
    double *d_kappa_fd,*d_lambda_fd;

    d_kappa_fd = malloc(sizeof(double)*(p+(p*p)));
    d_lambda_fd = d_kappa_fd + p;
#endif


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
    for(i=0,col_shift=1,l_int=(p-1),j=0;j<(p*(p-1))/2;j++,i++){
        if(l_int==i){
            l_int--;i=0;col_shift++;
        }
        /* i is row , i+colshift is col */
        x[j+p]=lambda[ i*p + (i+col_shift) ];
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
            matrixDiagToFull(p,x+p,lambda);

            // Call loss function with current x
            f = mv_vonmises_lossFunction(n,p,kappa,lambda,S,C,ro,d_kappa,d_lambda);
/*******************************************************/
#ifdef DEBUG
            // Approx DF for kappa
            for(i=0;i<p;i++){
                kappa[i]+=1E-9;
                d_kappa_fd[i] = (mv_vonmises_lossFunction(n,p,kappa,lambda,S,C,ro,NULL,NULL)-f)/1E-9  ;
                kappa[i]-=1E-9;
                printf("KAPPA(%d):  calc %f \t  fd approx: %f \t DIFF: %f \n",i,d_kappa[i],d_kappa_fd[i],fabs(d_kappa[i]-d_kappa_fd[i]));
            }

            // Approx DF for lambda
            for(i=0;i<(p-1);i++){
                    for(j=i+1;j<p;j++){
                        lambda[i*p + j] += 1E-9 ; 
                        lambda[j*p + i] += 1E-9 ; 
                        d_lambda_fd[i*p + j] = (mv_vonmises_lossFunction(n,p,kappa,lambda,S,C,ro,NULL,NULL)-f)/1E-9  ;
                        d_lambda_fd[j*p + i] = d_lambda_fd[i*p + j];
                        printf("Lambda(%d,%d):  calc %f \t  fd approx: %f \t DIFF: %f \n",i,j,d_lambda[i*p + j],d_lambda_fd[i*p +j],fabs(d_lambda[i*p +j ]-d_lambda_fd[i*p + j]));
                        lambda[i*p + j] -= 1E-9 ; 
                        lambda[j*p + i] -= 1E-9 ; 
                    }
            }

#endif
/**************************************************************/

            // Kappa partials
            memcpy(g,d_kappa,sizeof(double)*p);

            // Apply penalization
            if( (phi!=NULL) && (H != NULL) ){

                double fnorm = 0;
                double fnorm_term;

                // Compute f norm
                 for (i=0; i < (p-1); i++) {
                    for (j=i+1; j < p; j++){
                        fnorm_term = (lambda[ i*p + j] - phi[ i*p + j]) * H[ i*p + j];
                        fnrom += fnorm_term * fnorm_term;
                    }
                 }
            
                 // We also need to include kappa
                 for (i=0; i<p; i++){
                        fnorm_term = (kappa[ i ] - phi[ i*p + i]) * H[ i*p + i];
                        fnrom += fnorm_term * fnorm_term;
                 }
                 // Add norm to F
                 fnorm = sqrtl(fnorm);
                 f += fnorm;

                // Modify lambda partials
                for (i=0; i < (p-1); i++) {
                    for (j=i+1; j < p; j++){
                        d_lambda[i*p + j]+= H[i*p + j] * H[i*p + j] * ( lambda[i*p + j] - phi[ i*p + j] ) / fnorm ;
                        d_lambda[j*p + i]=d_lambda[i*p + j];
                    }
                }

                // Modify kappa partials
                for (i=0; i<p; i++){
                        d_kappa[i]+= H[i*p + i] * H[i*p + i] * ( kappa[i] - phi[ i*p + i] ) / fnorm ;
                 }
            }

            // Lambda partials (to vector)
            for(i=0,col_shift=1,l_int=(p-1),j=0;j<(p*(p-1))/2;j++,i++){
                if(l_int==i){
                    l_int--;i=0;col_shift++;
                }
                /* i is row , i+colshift is col */
                g[j+p]=d_lambda[ i*p + (i+col_shift) ];
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
    matrixDiagToFull(p,x+p,lambda);

    /** FreE*/

#ifdef DEBUG
            free(d_kappa_fd);
#endif
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
void __R_mvvonmises_lbfgs_fit(int* p, double* mu, double* kappa, double* lambda, int* n, double *samples, double* penMatrix,
                            int* verbose, double* prec, double* tol, int* mprec, double *lower, double *upper, int *bounded , double *loss){
    *loss = mvvonmises_lbfgs_fit(*p, mu, kappa, lambda, *n, samples, penMatrix, *verbose, *prec, *tol, *mprec, lower, upper, bounded );
}


double mvvonmises_lbfgs_fit_OLD(int p, double* mu, double* kappa, double* lambda, int n, double *samples, double* penMatrix,
                            int verbose, double prec, double tol, int mprec, double *lower, double *upper, int *bounded ){

    uint_fast16_t i,j,l_int,col_shift;

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
    for(i=0,col_shift=1,l_int=(p-1),j=0;j<(p*(p-1))/2;j++,i++){
        if(l_int==i){
            l_int--;i=0;col_shift++;
        }
        /* i is row , i+colshift is col */
        x[j+p]=lambda[ i*p + (i+col_shift) ];
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
            matrixDiagToFull(p,x+p,lambda);

            // Call loss function with current x
            f = mv_vonmises_logpseudolikelihood_nonopt(n,p,mu,kappa,lambda,samples);
            
            // kappa
            for(i=0;i<p;i++) d_kappa[i] = vm_loglikelihood_kappa_deriv_old(n,p,i,mu,kappa,lambda,samples);

            // Approx DF for lambda
            for(i=0;i<(p-1);i++){
                    for(j=i+1;j<p;j++){
                        d_lambda[i*p + j] = vm_loglikelihood_lambda_deriv_old(n,p,i,j,mu,kappa,lambda,samples);
                        d_lambda[j*p + i] = d_lambda[i*p + j];
                    }
            }

            // Kappa partials
            memcpy(g,d_kappa,sizeof(double)*p);

            // Apply penalization
            if(penMatrix!=NULL){
                for(i=0;i<(p-1);i++){
                    for(j=i+1;j<p;j++){
                        f += penMatrix[i*p+j]*fabs(lambda[i*p+j]);
                        d_lambda[i*p + j]+=penMatrix[i*p + j]*(lambda[i*p + j]<0?-1:(lambda[i*p +j]>0?1:0));
                        d_lambda[j*p + i]=d_lambda[i*p + j];
                    }
                }
            }

            // Lambda partials
            for(i=0,col_shift=1,l_int=(p-1),j=0;j<(p*(p-1))/2;j++,i++){
                if(l_int==i){
                    l_int--;i=0;col_shift++;
                }
                /* i is row , i+colshift is col */
                g[j+p]=d_lambda[ i*p + (i+col_shift) ];
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
    matrixDiagToFull(p,x+p,lambda);

    /** FreE*/

    free(S); free(d_kappa); free(g);
    free(wa); free(iwa); free(nbd);

    if( IS_ERROR(*task) || IS_WARNING(*task) )
        return NAN;
    else
        return (double)f;
}

void __R_mvvonmises_lbfgs_fit_OLD(int* p, double* mu, double* kappa, double* lambda, int* n, double *samples, double* penMatrix,
                            int* verbose, double* prec, double* tol, int* mprec, double *lower, double *upper, int *bounded , double *loss){
    *loss = mvvonmises_lbfgs_fit_OLD(*p, mu, kappa, lambda, *n, samples, penMatrix, *verbose, *prec, *tol, *mprec, lower, upper, bounded );
}
