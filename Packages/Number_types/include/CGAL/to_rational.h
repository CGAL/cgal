// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/to_rational.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Susan Hert, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
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
      const unsigned int shift_pow = 32768; // = 2^shift
      const double width = 32768;  // = 2^shift
      const int maxiter = 20;      // ought not be necessary, but just in case,
                                   // max 300 bits of precision
      int expt;
      double mantissa = frexp(x, &expt);
      long exponent = expt;
      double intpart;
      int k = 0;
      
      while (mantissa != 0.0 && k++ < maxiter)

      { mantissa *= width; // shift double mantissa
        mantissa = CGAL_CLIB_STD::modf(mantissa, &intpart);
        num *= shift_pow;
        num += (long)intpart;
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
