/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project
 * 
 * File: BigInt.cpp
 *
 * Synopsis: a wrapper class of mpz in GMP
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#include <CORE/BigInt.h>

CORE_BEGIN_NAMESPACE

// get exponent of power k
void BigInt::getKaryExpo(BigInt& m, int& e, unsigned long k) const {
  BigInt q, r;
  m = *this; e = 0;
  mpz_tdiv_qr_ui(q.mpz(), r.mpz(), m.mpz(), k);
  while (r.sign() == 0) {
    e ++; m = q; 
    mpz_tdiv_qr_ui(q.mpz(), r.mpz(), m.mpz(), k);
  }
} 

// return a gmp_randstate_t structure
gmp_randstate_t* getRandstate() {
  static gmp_randstate_t rstate;
  static bool initialized = false;
  if (!initialized) {
    gmp_randinit(rstate, GMP_RAND_ALG_DEFAULT, 32L);
    initialized = true;
  }
  return &rstate;
}

CORE_END_NAMESPACE

