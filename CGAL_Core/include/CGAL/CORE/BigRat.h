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
 * WWW URL: http://cs.nyu.edu/exact/
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


  /// BigIntValue
inline BigInt BigIntValue(const BigRat& br)
{
  BigInt r, rem;
  divide_qr(numerator(br), denominator(br), r, rem);
  return r;
}

} // namespace CORE


#endif // _CORE_BIGRAT_H_
