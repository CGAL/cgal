/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: BigRat.h
 * Synopsis:
 *                 a wrapper class for mpq from GMP
 *
 * Written by
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: https://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

#ifndef _CORE_BIGRAT_H_
#define _CORE_BIGRAT_H_

#include <CGAL/CORE/BigInt.h>

namespace CORE {
#ifdef CGAL_CORE_USE_GMP_BACKEND
  typedef boost::multiprecision::mpq_rational BigRat;
#else
  typedef  boost::multiprecision::cpp_rational BigRat;
#endif


  inline BigInt numerator(const BigRat& q)
  {
    return boost::multiprecision::numerator(q);
  }

  inline BigInt denominator(const BigRat& q)
  {
    return boost::multiprecision::denominator(q);
  }

  // Chee (3/19/2004):
//   The following definitions of div_exact(x,y) and gcd(x,y)
//   ensures that in Polynomial<NT>
/// divisible(x,y) = "x | y"
  inline BigRat div_exact(const BigRat& x, const BigRat& y) {
    BigRat z = x / y;
    return z;
  }

  inline BigRat gcd(const BigRat& x , const BigRat& y)
  {
    //        return BigRat(1);  // Remark: we may want replace this by
                           // the definition of gcd of a quotient field
                           // of a UFD [Yap's book, Chap.3]
  //Here is one possible definition: gcd of x and y is just the
  //gcd of the numerators of x and y divided by the gcd of the
  //denominators of x and y.
  BigInt n = gcd(numerator(x), numerator(y));
  BigInt d = gcd(denominator(x), denominator(y));
  return BigRat(n,d);
  }

  /// BigIntValue
inline BigInt BigIntValue(const BigRat& br)
{
  BigInt r, rem;
  divide_qr(numerator(br), denominator(br), r, rem);
  return r;
}

} // namespace CORE


#endif // _CORE_BIGRAT_H_
