/*******************
 * Implementation of Pseudo random number generators
 * 
 * Implemented Methods: Xorshift*1024, Xorshift*4096
 *
 * $Date:$
 * $Revision:$
 * $Author:$
 * $HeadUrl:$
 * $Id: prng.h 6 2015-02-17 15:19:00Z lrodriguez $
 **/
#ifndef __PRNG_H_
#define __PRNG_H_

#include <stdint.h>

/**
 * XORSHIFT * CONSTS
 */
#define XORSHIFT_MULT1024 1181783497276652981LL
#define XORSHIFT_MULT4096 8372773778140471301LL

/**
 *  Definition:  PRNG Function pointers
 **/

 /* Non reentrant prng */ 
 typedef uint64_t (*pfunc_prng) ();

 /* Reentrant prng */
 typedef uint64_t (*pfunc_prng_r) (void *status);

 /* Seed setter for non reentrant */
 typedef void (*pfunc_prng_set_seed) (void *status);

/*********************************************/

/**
 * xorshift1024* method as in Vigna[2014a] a:31 b:11 c:30
 */
uint64_t xorshift1024(void);
uint64_t xorshift1024_r(uint64_t *status);

/**
 * xorshift4096* method as in Vigna[2014a] a:25 b:3 c:49
 **/
uint64_t xorshift4096(void);
uint64_t xorshift4096_r(uint64_t *status);

/**
 * Sets first n seed values (64 or 16)
 **/
void xorshift4096_seed(uint64_t *v);
void xorshift1024_seed(uint64_t *v);

#endif
