/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/); you may
 * redistribute it under the terms of the Q Public License version 1.0.
 * See the file LICENSE.QPL distributed with CORE.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: CoreAux.cpp
 * Synopsis:
 *       Auxiliary routines such as ceiling of log_2, etc. 
 *       they are not specific to any Core classes.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Source$
 * $Revision$ $Date$
 ***************************************************************************/

#include "CORE/CoreAux.h"
#include <gmp.h>

CORE_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////
//  More useful functions to implement:
//
//  To convert digits into bits:
//      given X, compute X * log_2(10)
//  To convert bits into digits:
//      given X, compute X * log_10(2)
//
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// flrLg(x)
//      returns floor log base 2 of abs(x)
// CONVENTION lg(0) = -1	(Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
int flrLg(long x) {
  if (x == LONG_MIN) {
    // special treatment as -LONG_MIN would be not representable as "long"
    return LONG_BIT - 1;
  } else {
    //  1 <= |x| <= LONG_MAX
    if (x < 0)
      x = - x;

    int lg = -1;
    while (x > 0) {
      lg++;
      x >>= 1;
    }
    return lg;
  }
}

////////////////////////////////////////////////////////////
// floor log base 2 of unsigned long x
// CONVENTION lg(0) = -1	(Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
int flrLg(unsigned long x) {
  int lg = -1;
  while (x > 0) {
    lg++;
    x >>= 1;
  }
  return lg;
}

////////////////////////////////////////////////////////////
// ceiling log base 2 of abs(x)
// CONVENTION lg(0) = -1	(Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
int clLg(long x) {
  if (x == LONG_MIN)
    return LONG_BIT - 1;
  if (x < 0)
    x = -x; 		// use absolute value
  if (x > (LONG_MAX >> 1)) 	// the leading data bit is 1
    return (LONG_BIT - 1);	// exclude the sign bit
  if (x >= 2)
    return flrLg((unsigned long)((x << 1) - 1));
  // SINCE ceilLog_2(x) = floorLog_2(2x-1) for x>=2
  if (x == 1)
    return 0;
  return -1;			// x must be 0 here
}

////////////////////////////////////////////////////////////
// ceiling log base 2 of unsigned long x
// CONVENTION lg(0) = -1
////////////////////////////////////////////////////////////
int clLg(unsigned long x) {
  if (x > (ULONG_MAX >> 1))	// the leading bit is 1
    return LONG_BIT;
  if (x >= 2)
    return flrLg((x << 1) - 1);
  // SINCE ceilLog_2(x) = floorLog_2(2x-1) for x>=2.
  if (x == 1)
    return 0;
  return -1;	// x must be equal to 0
}

/// gcd for machine type long
/** This is needed when we construct Polynomials with int or long coefficients */
long gcd(long m, long n) {
  if (m == 0)
    return core_abs(n);
  if (n == 0)
    return core_abs(m);
  long p = core_abs(m);
  long q = core_abs(n);
  if (p<q)
    core_swap(p, q);
  while (q>0) {
    long r = p % q;
    p = q;
    q = r;
  }
  return p;
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

// char* core_itoa(int n, char* buffer)
//      returns a pointer to the null-terminated string in buffer
// NOTES:
// 0. Buffer size should be 17 bytes (resp., 33 bytes, 65 bytes) on 16-bit
//      (resp., 32-bit, 64-bit) machines.  Formula: 1+sizeof(int)*8 bytes.
// 1. itoa(...) is available on some stdlib.h, but it is
//      not defined by ANSI-C and so not all compilers support it.
// 2. Our use of sprintf(...) to do the job is known to
//      be inefficient, but this is hardly critical for our usage.
// 3. A more general program should take a 3rd argument (the radix of
//      output number).  We assume radix 10.
char * core_itoa(int n, char* buffer) {
	sprintf(buffer, "%d", n);
	return buffer;
}

/// implements the "integer mantissa" function
//      (See CORE_PATH/progs/ieee/frexp.cpp for details)
double IntMantissa(double d) {
	int e;
	return ldexp(frexp(d, &e), 53);
}

/// implements the "integer exponent" function
//      (See CORE_PATH/progs/ieee/frexp.cpp for details)
int IntExponent(double d) {
	int e;
	frexp(d, &e);
	return e-53;
}

// core_error is the method to write Core Library warning or error messages
// 	Warnings and Errors are written to a file called CORE_DIAGFILE
// 	But Errors are also written on std:cerr (via the std::perror()).
//
// file where core_error(...) writes its output
char* CORE_DIAGFILE = "Core_Diagnostics"; // global file name

// error function for general use in Core Library
void core_error(std::string msg, std::string file, int lineno, bool err) {
  std::ofstream outFile(CORE_DIAGFILE, std::ios::app);  // open to append
  if (!outFile) {
     // perror("CORE ERROR: cannot open Core Diagnostics file");
     std::cerr << "CORE ERROR: can't open Core Diagnostics file"<<std::endl;
     exit(1); //Note: do not call abort()
  }
  outFile << "CORE " << (err? "ERROR" : "WARNING")
     << " (at " << file << ": " << lineno << "): "
     << msg << std::endl;
  outFile.close();
  if (err) {
     char buf[65];
     // perror((std::string("CORE ERROR") + " (file " + file + ", line "
     //        + core_itoa(lineno,buf) +"):" + msg + "\n").c_str());
     std::cerr << (std::string("CORE ERROR") + " (file " + file + ", line "
             + core_itoa(lineno,buf) +"):" + msg + "\n").c_str();
     exit(1); //Note: do not call abort()
  }
}


CORE_END_NAMESPACE
