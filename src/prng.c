/*******************
 * Implementation of Pseudo random number generators
 * 
 * Implemented Methods: Xorshift*1024, Xorshift*4096
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: prng.c 10 2015-02-24 13:37:13Z lrodriguez $
 **/

#include <string.h>
#include <stdint.h>
#include "prng.h"

/** Xorshift* Seeds circular array */
/* Random numbers from random.org (Default seeds) */
uint64_t xs_seeds[64]={396211266,318936545,239782012,327524388,262111400,138645947,129773377,388466163,112746286,213523003,275358881,282949203,178572974,\
        90135547,421770271,303986047,309203197,255247428,268867007,318765536,304687789,405614684,269760794,236489007,411620509,166766312,\
            361210002,283644967,119641924,69485002,377870531,113980432,364456658,287448779,141061818,198138578,116842912,307991489,201984000,\
            350364459,92809236,164348060,148959980,326501783,132148331,24694032,34476500,274413237,164035902,11016039,329087019,194669370,\
            190407088,121086105,26299711,367248076,92649296,79189179,256366023,204156770,2324851,391323303,102893852,11160450};


/*********************************************/

/**
 * xorshift1024* method as in Vigna[2014a] a:31 b:11 c:30
 */
uint64_t xorshift1024_r(uint64_t *s){

    static int p=0; /** Current status pointer (circular array) **/
    register uint64_t s0 = s[p];
    register uint64_t s1 = s[ p = (p+1) & 15U ] ;/* Circular pointer, &&15 <-> mod 16 */
    
    s1 ^= s1 << 31;
    s1 ^= s1 >> 11;
    s0 ^= s0 >> 30;
    
    return ( s[p] = s0 ^ s1 ) * XORSHIFT_MULT1024;
}

uint64_t xorshift1024(){
    return xorshift1024_r(xs_seeds);
}

/**
 * xorshift4096* method as in Vigna[2014a] a:25 b:3 c:49
 **/
uint64_t xorshift4096_r(uint64_t *s){

    static int p=0; /** Current status pointer (circular array) **/
    register uint64_t s0 = s[p];
    register uint64_t s1 = s[ p = ((p+1) & 63U) ] ;/* Circular pointer, &&15 <-> mod 64 */
    
    s1 ^= s1 << 25;
    s1 ^= s1 >> 3;
    s0 ^= s0 >> 49;
    
    return ( s[p] = s0 ^ s1 ) * XORSHIFT_MULT4096;
}

uint64_t xorshift4096(){
    return xorshift4096_r(xs_seeds);
}

/**
 * Sets first n seed values (up to 64)
 **/
void xorshift4096_seed(uint64_t *v){
    memcpy(xs_seeds,v,sizeof(uint64_t)*64);
}

void xorshift1024_seed(uint64_t *v){
    memcpy(xs_seeds,v,sizeof(uint64_t)*16);
}

