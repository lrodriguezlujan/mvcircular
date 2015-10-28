/*******************
 * Matrix and vector auxiliar operations
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#ifndef _MATRIXOPS_H__
#define _MATRIXOPS_H__
#include <stdint.h>


/**
 * Computes c = v * op(A)
 *
 * A  is  n x p matrix
 * v  is  1xn vector   || 1 x p vector if transA is T or t
 * c  is  1 x p vector || 1 x n vector if transA is T or t
 */
void leftMatrixVector(int n, int p, char transA, double *A, double *v, double *C);


/**
 * Computes c = op(A) * v
 *
 * A  is  n x p matrix
 * v  is  p x 1 vector || n x 1 if transA is T or t
 * c  is  n x 1 vector || p x 1 if transA is T or t
 */
void rightMatrixVector(int n, int p, char TransA, double *A, double *v, double *C);
/**
 * Computes C = op(A) * op(B)
 *
 * A is m x n matrix
 * B is p x k matrix
 * C depends on op(A) and op(B)
 */
int generalMatrixMatrix(int m, int n, int p, int k, char transA, char transB, long double *A, double *B, long double *C);


/******************************** END ***********************************/
#endif
