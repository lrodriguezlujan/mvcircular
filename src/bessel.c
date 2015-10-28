/*******************
 * Bessel functions
 *
 * $Date: 2015-04-20 17:09:57 +0200 (lun 20 de abr de 2015) $
 * $Revision: 14 $
 * $Author: lrodriguez $
 * $HeadUrl:$
 * $Id: stats.h 14 2015-04-20 15:09:57Z lrodriguez $
 **/

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "bessel.h"


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
long double bessi0_exact(long double x, long double tol){

    // Return value and current
    long double acum, term;

    // Iterator
    register uint_fast64_t i;

    uintmax_t factorial=1;

    // Initialization
    acum = 1 ; term = 1 ; i = 1;

    // While last computed term is greater than given tolerance
    while( term > tol ){
        factorial*=i; // Update factorial

        // Compute term
        term = powl(x/2.0,2*i) / (factorial*factorial);

        if(!isfinite(term)) break;
        acum+=term;

        i++;
    }
    return acum;
}

/**
 * Computes bessel I0 function approximately based on precomputed coefficients
 *
 * @param [in] x Bessel function argument
 *
 * @return approx(I_0) 
 * @seealso bessi0_exact
 */
long double bessi0_spd(long double x ){
   long double ax;
   long double y;

   if ((ax=fabsl(x)) < 3.75) {
      y=x/3.75,y=y*y;
      return (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))));
   } else {
      y=3.75/ax;
      return ((expl(ax)/sqrtl(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2)))))))));
   }
}

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
long double bessi0_spd_exponScaled(long double x){

    long double y;
    if (x < 3.75) {
      y=x/3.75,y=y*y;
      return expl(-x)*(1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))));
   } else {
      y=3.75/x;
      return ((0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))))/sqrtl(x));
   }

}

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
long double bessi0_salahat(double x){

    if(x<11.5) 
        // First case: alfa = 0.1682 , 0.1472 , 0.4450 , 0.2382 ;
        //             beta = 0.7536 , 0.9739 , -0.715 , 0.2343
        return 0.1682 * expl(x * 0.7536) + 0.1472 * expl(x * 0.9739) + 0.4450 * expl(x * -0.715) + 0.2382 * expl(x * 0.2343);
    else if (x<20)
         // Second case: alfa =  ;
        //             beta = 
        return 0.2667 * expl(x * 0.4710) + 0.4916 * expl(x * -163.4) + 0.1110 * expl(x * 0.9852) + 0.1304 * expl(x * 0.8554);
    else if (x<37.25)
         // third case: alfa =  ;
        //             beta = 
        return 0.1121 * expl(x * 0.9807) + 0.1055 * expl(x * 0.8672) + -0.00018 * expl(x * 1.0795) + 0.00326 * expl(x * 1.0385);
    else
        // Last one
        return 2.41E-9 * expl(x * 1.144) + 0.06745 * expl(x * 0.995) + 0.05471 * expl(x * 0.5686) + 0.07869 * expl(x * 0.946);

}

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
long double bessi0_salahat_exponScaled(double x){

    if(x<11.5) 
        // First case: alfa = 0.1682 , 0.1472 , 0.4450 , 0.2382 ;
        //             beta = 0.7536 , 0.9739 , -0.715 , 0.2343
        return 0.1682 * expl(x * (0.7536 - 1) ) + 0.1472 * expl(x * (0.9739 - 1) ) + 0.4450 * expl(x * (-0.715 - 1) ) + 0.2382 * expl( x * (0.2343 - 1) );
    else if (x<20)
         // Second case: alfa =  ;
        //             beta = 
        return 0.2667 * expl(x * (0.4710 - 1) ) + 0.4916 * expl(x * (-163.4 - 1) ) + 0.1110 * expl(x * (0.9852 - 1) ) + 0.1304 * expl( x * (0.8554 - 1) );
    else if (x<37.25)
         // third case: alfa =  ;
        //             beta = 
        return 0.1121 * expl(x * (0.9807 - 1) ) + 0.1055 * expl(x * (0.8672 - 1) ) + -0.00018 * expl(x * (1.0795 - 1 )) + 0.00326 * expl(x * (1.0385 - 1));
    else
        // Last one
        return 2.41E-9 * expl(x * (1.144 - 1) ) + 0.06745 * expl(x * (0.995 - 1) ) + 0.05471 * expl(x * (0.5686 - 1) ) + 0.07869 * expl( x * (0.946 - 1) );
}



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
long double bessi1_salahat(double x){

    if(x<11.5) 
        // First case: alfa = 0.1682 , 0.1472 , 0.4450 , 0.2382 ;
        //             beta = 0.7536 , 0.9739 , -0.715 , 0.2343
        return 0.1682 * 0.7536 * expl(x * 0.7536) + 0.1472 * 0.9739 * expl(x * 0.9739) + 0.4450 * -0.715 * expl(x * -0.715) + 0.2382 * 0.2343 * expl(x * 0.2343);
    else if (x<20)
         // Second case: alfa =  ;
        //             beta = 
        return 0.2667 * 0.4710 * expl(x * 0.4710) + 0.4916 * -163.4 * expl(x * -163.4) + 0.1110 * 0.9852 * expl(x * 0.9852) + 0.1304 * 0.8554 * expl(x * 0.8554);
    else if (x<37.25)
         // third case: alfa =  ;
        //             beta = 
        return 0.1121 * 0.9807 * expl(x * 0.9807) + 0.1055 * 0.8672 * expl(x * 0.8672) + -0.00018 * 1.0795 * expl(x * 1.0795) + 0.00326 * 1.0385 * expl(x * 1.0385);
    else
        // Last one
        return 2.41E-9 * 1.144 * expl(x * 1.144) + 0.06745 * 0.995 * expl(x * 0.995) + 0.05471 * 0.5686 * expl(x * 0.5686) + 0.07869 * 0.946 * expl(x * 0.946);
}


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
long double bessi1_salahat_exponScaled(double x){

    if(x<11.5) 
        // First case: alfa = 0.1682 , 0.1472 , 0.4450 , 0.2382 ;
        //             beta = 0.7536 , 0.9739 , -0.715 , 0.2343
        return 0.1682 * (0.7536-1) * expl(x * (0.7536 - 1) ) + 0.1472 * (0.9739-1) * expl(x * (0.9739 - 1) ) + 0.4450 * (-0.715-1) * expl(x * (-0.715 - 1) ) + 0.2382 * (0.2343-1) * expl( x * (0.2343 - 1) );
    else if (x<20)
         // Second case: alfa =  ;
        //             beta = 
        return 0.2667 * (0.4710-1) * expl(x * (0.4710 - 1) ) + 0.4916 * (-163.4-1) * expl(x * (-163.4 - 1) ) + 0.1110 * (0.9852-1) * expl(x * (0.9852 - 1) ) + 0.1304 * (0.8554-1) * expl( x * (0.8554 - 1) );
    else if (x<37.25)
         // third case: alfa =  ;
        //             beta = 
        return 0.1121 * (0.9807-1) * expl(x * (0.9807 - 1) ) + 0.1055 * (0.8672-1) * expl(x * (0.8672 - 1) ) + -0.00018 * (1.0795-1) * expl(x * (1.0795 - 1 )) + 0.00326 * (1.0385-1) * expl(x * (1.0385 - 1));
    else
        // Last one
        return 2.41E-9 * (1.144-1) * expl(x * (1.144 - 1) ) + 0.06745 * (0.995-1) * expl(x * (0.995 - 1) ) + 0.05471 * (0.5686-1) * expl(x * (0.5686 - 1) ) + 0.07869 * (0.946-1) * expl( x * (0.946 - 1) );
}


/**
 * Computes bessel I1 function approximately based on precomputed coefficients
 *
 * @param [in] x Bessel function argument
 *
 * @return approx(I_1) 
 * @seealso bessi1_exact
 */

long double bessi1_spd(long double x)
{
   long double ax,ans;
   long double y;


   if ((ax=fabsl(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (expl(ax)/sqrtl(ax));
   }
   return x < 0.0 ? -ans : ans;
}

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

long double bessi1_spd_exponScaled(long double x)
{
   long double ans;
   long double y;


   if ( x < 3.75) {
      y=x/3.75,y=y*y;
      ans=expl(-x)*x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/x;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans /= sqrtl(x);
   }
   return  ans;
}


/******************************** END ***********************************/


/* R Wrappers*/

void __R_bessi0(double *x,int *expon,double *val){
    if(*expon)
        *val = bessi0_expon(*x);
    else
        *val = bessi0(*x);
}

void __R_bessi1(double *x,int *expon,double *val){
    if(*expon)
        *val = bessi1_expon(*x);
    else
        *val = bessi1(*x);
}
