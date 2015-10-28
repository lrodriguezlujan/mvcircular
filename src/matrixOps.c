/*******************
 * Matrix and vector auxiliar operations
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <stdint.h>
#include <string.h>
#include "matrixOps.h"


/**
 * Computes c = v * op(A)
 *
 * A  is  n x p matrix
 * v  is  1xn vector   || 1 x p vector if transA is T or t
 * c  is  1 x p vector || 1 x n vector if transA is T or t
 */
void leftMatrixVector(int n, int p, char transA, double *A, double *v, double *C){

    // Iterators
    uint_fast16_t i,j;

    memset(C,0,sizeof(double)*p);

    // Change loop order for maximizing cache usage
    if(transA=='T' || transA=='t'){
         for(i=0;i<n;i++)
            for(j=0;j<p;j++)
                C[j]+=v[i]*A[i*p+j];
    }
    else{
        for(j=0;j<p;j++)
            for(i=0;i<n;i++)
                C[i]+=v[j]*A[i*p+j];
    }
}

/**
 * Computes c = op(A) * v
 *
 * A  is  n x p matrix
 * v  is  p x 1 vector || n x 1 if transA is T or t
 * c  is  n x 1 vector || p x 1 if transA is T or t
 */
void rightMatrixVector(int n, int p, char transA, double *A, double *v, double *C){

    // Iterators
    uint_fast16_t i,j;

    memset(C,0,sizeof(double)*n);

    // Change loop order for maximizing cache usage
    if(transA=='T' || transA=='t'){
        for(i=0;i<n;i++)
            for(j=0;j<p;j++)
                C[j]+=v[i]*A[i*p+j];
    }
    else
    {
         for(j=0;j<p;j++)
            for(i=0;i<n;i++)
                C[i]+=v[j]*A[i*p+j];
    }
}

/**
 * Computes C = op(A) * op(B)
 *
 * A is m x n matrix
 * B is p x k matrix
 * C depends on op(A) and op(B)
 */
int generalMatrixMatrix(int m, int n, int p, int k, char transA, char transB, long double *A, double *B, long double *C){

    // iterators
    uint_fast16_t i,j,l;

    char opA = 0 , opB = 0 ;

    // Check parameters in each case
    opA = (transA=='T') || (transA=='t') ;
    opB = (transB=='T') || (transB=='t') ;

    // First case: C = A * B
    if( !opA && !opB ){
        /* n should be equal to p */
        if(n!=p) return -1;
        // C is m x k matrix
        else
        {
            memset(C,0,sizeof(long double)*m*k);
            for(j=0;j<k;j++){
                for(l=0;l<n;l++){
                    for(i=0;i<m;i++){
                        C[i*k + j] += B[l*k + j] *  A[i*n + l] ; 
                    }
                }
            }
        }
    }
    // Second case C = A^t * B
    else if (opA && !opB){
        // m shoud be equal to p */
        if(m!=p) return -1;
        // C is n x k matrix
        else{
            memset(C,0,sizeof(double)*n*k);
            for(j=0;j<k;j++){
                for(l=0;l<m;l++){
                    for(i=0;i<n;i++){
                        C[i*k + j] += B[l*k + j] *  A[i*m + l] ; 
                    }
                }
            }
        }
    }
    // Third case  C = A * B^t  
    else if( !opA && opB ){
        // n should be equal to k
        if(n!=k) return -1;
        // C matrix is m x p 
        else{
            memset(C,0,sizeof(double)*m*p);
            for(j=0;j<p;j++){
                for(l=0;l<n;l++){
                    for(i=0;i<m;i++){
                        C[i*p + j] += B[l*p + j] *  A[i*n + l] ; 
                    }
                }
            }
        }
    }
    // Last case C = A^t * B^t
    // m should be equal to k
    else if(m!=k) return -1;
    // C matrix is n x p 
    else{
        memset(C,0,sizeof(double)*n*p);
        for(j=0;j<p;j++){
            for(l=0;l<m;l++){
                for(i=0;i<n;i++){
                    C[i*p + j] += B[l*p + j] *  A[i*m + l] ; 
                }
            }
        }
    }
    return 0;
}

/******************************** END ***********************************/
