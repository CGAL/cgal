// Copyright (c) 2002  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Susan Hert, Sylvain Pion
 
#ifndef CGAL_TO_RATIONAL_H
#define CGAL_TO_RATIONAL_H

#include <CGAL/Number_type_traits.h>

CGAL_BEGIN_NAMESPACE

template <class Rational>
Rational
to_rational(double x)
{ 
    typename Rational_traits<Rational>::RT num = 0;
    typename Rational_traits<Rational>::RT den = 1;

    if (x != 0.0)
    { bool neg = (x < 0);
      if (neg) x = -x;

      const unsigned shift = 15;   // a safe shift per step
      const int shift_pow = 32768; // = 2^shift
      const double width = 32768;  // = 2^shift
      const int maxiter = 20;      // ought not be necessary, but just in
                                   // case, max 300 bits of precision
      int expt;
      double mantissa = frexp(x, &expt);
      long exponent = expt;
      double intpart;
      int k = 0;

      while (mantissa != 0.0 && k++ < maxiter)

      { mantissa *= width; // shift double mantissa
        mantissa = CGAL_CLIB_STD::modf(mantissa, &intpart);
        num *= shift_pow;
        num += (int)intpart;
        exponent -= shift;
      }
      int expsign = (exponent>0 ? +1 : (exponent<0 ? -1 : 0));
      exponent *= expsign;
      typename Rational_traits<Rational>::RT twopot(2);
      typename Rational_traits<Rational>::RT exppot(1);
      while (exponent!=0) {
        if (exponent & 1)
          exppot *= twopot;
        exponent >>= 1;
        twopot *= twopot;
      }

      if (expsign > 0)
        num *= exppot;
      else if (expsign < 0)
        den *= exppot;
      if (neg)
        num = -num;
    }
    return Rational(num,den);
}

CGAL_END_NAMESPACE

#endif // CGAL_TO_RATIONAL_H
