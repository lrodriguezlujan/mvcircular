/*******************
 * Bessel functions
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/


#ifndef _BESSEL_H__
#define _BESSEL_H__

/* Use different approaches depending on the value of x */
#define BESSEL_THRESHOLD 1E31
#define bessi0(x) (x)<(BESSEL_THRESHOLD)?(bessi0_spd(x)):(bessi0_salahat(x))
#define bessi0_expon(x) (x)<(BESSEL_THRESHOLD)?(bessi0_spd_exponScaled(x)):(bessi0_salahat_exponScaled(x))
#define bessi1(x) (x)<(BESSEL_THRESHOLD)?(bessi1_spd(x)):(bessi1_salahat(x))
#define bessi1_expon(x) (x)<(BESSEL_THRESHOLD)?(bessi1_spd_exponScaled(x)):(bessi1_salahat_exponScaled(x))


/***************************************************
 * Bessel function of order 0
 **************************************************/

/**
 * Computes bessel I function of zero order applying the formula
 *
 *  @param [in] x Bessel function argument
 *  @param [in] tol Tolerance
 *
 *  @return value
 */
long double bessi0_exact(long double x, long double tol);

/**
 * Computes bessel I0 function approximately based on precomputed coefficients
 *
 * @param [in] x Bessel function argument
 *
 * @return approx(I_0) 
 * @seealso bessi0_exact
 */
long double bessi0_spd(long double x);

/**
 * Exponential scaled version of the bessel I0 function
 *
 * Approximated,  computed using only the first 5 coefficients.
 *
 * @param [in] x Bessel function argument
 *
 * @return exp(-x)I_0(x)
 * @seealso bessi0_spd
 */
long double bessi0_spd_exponScaled(long double x);

/**
 * Computes bessel I function of zero order
 *
 * Computes bessel modified function of first kind, zero order
 * applying the approximation given by salahat et al (2013)
 * Provides good absolute relative error in all ranges. To be used
 * for very high values of x
 * 
 * @param [in] x Bessel function argument (>0)
 *
 * @return I_0(x)
 */
long double bessi0_salahat(double x);
/**
 * Computes exponential scaled bessel I function of zero order
 *
 * Computes bessel modified function of first kind, zero order
 * applying the approximation given by salahat et al (2013)
 * Provides good absolute relative error in all ranges. To be used
 * for very high values of x
 * 
 * @param [in] x Bessel function argument (>0)
 *
 * @return exp(-x)I_0(x)
 * @seealso bessi0_salahat
 */
long double bessi0_salahat_exponScaled(double x);


/************************************************************
 *
 * Modified bessel function order 1
 *
 ************************************************************/


/**
 * Computes bessel I function of first order
 *
 * Computes bessel modified function of first kind, first order
 * applying the approximation given by salahat et al (2013)
 * Provides good absolute relative error in all ranges)
 * 
 * @param [in] x Bessel function argument (>0)
 *
 * @return I_1(x)
 */
long double bessi1_salahat(double x);

/**
 * Computes exponential scaled bessel I function of first order
 *
 * Computes bessel modified function of first kind, first order
 * applying the approximation given by salahat et al (2013)
 * Provides good absolute relative error in all ranges. To be used
 * for very high values of x
 * 
 * @param [in] x Bessel function argument (>0)
 *
 * @return exp(-x)I_1(x)
 * @seealso bessi0_salahat
 */
long double bessi1_salahat_exponScaled(double x);

/**
 * Computes bessel I1 function approximately based on precomputed coefficients
 *
 * @param [in] x Bessel function argument
 *
 * @return approx(I_1) 
 * @seealso bessi1_exact
 */

long double bessi1_spd(long double x);

/**
 * Exponential scaled version of the bessel I1 function
 *
 * Approximated,  computed using only the first 5 coefficients.
 *
 * @param [in] x Bessel function argument
 *
 * @return exp(-x)I_1(x)
 * @seealso bessi1_spd
 */

long double bessi1_spd_exponScaled(long double x);


/******************************** END ***********************************/
#endif
