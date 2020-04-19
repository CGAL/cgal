// Copyright (c) 2005,2007
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion, Michael Hemmer

#ifndef CGAL_LONG_DOUBLE_H
#define CGAL_LONG_DOUBLE_H

#include <utility>
#include <cmath>
#ifdef CGAL_CFG_IEEE_754_BUG
#  include <CGAL/IEEE_754_unions.h>
#endif

// #include <CGAL/FPU.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/assertions.h>
#include <CGAL/functional.h>

namespace CGAL {

// Is_valid moved to top, since used by is_finite
#ifdef CGAL_CFG_IEEE_754_BUG

#define CGAL_EXPONENT_DOUBLE_MASK   0x7ff00000
#define CGAL_MANTISSA_DOUBLE_MASK   0x000fffff

inline
bool
is_finite_by_mask_long_double(unsigned int h)
{
  unsigned int e = h & CGAL_EXPONENT_DOUBLE_MASK;
  return ( ( e ^ CGAL_EXPONENT_DOUBLE_MASK ) != 0 );
}

inline
bool
is_nan_by_mask_long_double(unsigned int h, unsigned int l)
{
  if ( is_finite_by_mask_long_double(h) )
      return false;
  return (( h & CGAL_MANTISSA_DOUBLE_MASK ) != 0) || (( l & 0xffffffff ) != 0);
}

template<>
class Is_valid< long double >
  : public CGAL::cpp98::unary_function< long double, bool > {
  public :
    bool operator()( const long double& x ) const {
      double d = x;
      IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
      return ! ( is_nan_by_mask_long_double( p->c.H, p->c.L ));
    }
};

#else

template<>
class Is_valid< long double >
  : public CGAL::cpp98::unary_function< long double, bool > {
  public :
    bool operator()( const long double& x ) const {
      return (x == x);
    }
};

#endif




template <> class Algebraic_structure_traits< long double >
  : public Algebraic_structure_traits_base< long double,
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false            Is_exact;
    typedef Tag_true             Is_numerical_sensitive;

    class Sqrt
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return std::sqrt( x );
        }
    };

    class Kth_root
      :public CGAL::cpp98::binary_function<int, Type, Type > {
      public:
        Type operator()( int k,
                                        const Type& x) const {
          CGAL_precondition_msg( k > 0,
                                    "'k' must be positive for k-th roots");
          return std::pow(x, (long double)1.0 / (long double)(k));
        };
    };

};

template <> class Real_embeddable_traits< long double >
  : public INTERN_RET::Real_embeddable_traits_base< long double , CGAL::Tag_true > {
  public:

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          // The conversion long double to double does not always follow the
          // rounding mode (e.g. on MacIntel/g++-4.0.2).  So let's do a naive
          // conversion to double, then widen the interval.
          return (Interval_nt<>((double)x)+Interval_nt<>::smallest()).pair();
#if 0
          // We hope that the long double -> double conversion
          // follows the current rounding mode.
          Protect_FPU_rounding<true> P(CGAL_FE_UPWARD);
          volatile long double mx = -x; // needed otherwise the conversion can
                                        // get factorized between d and -d...
          if (x>0)
            return std::make_pair(- static_cast<double>(CGAL_IA_FORCE_TO_DOUBLE(mx)),
                                  static_cast<double>(CGAL_IA_FORCE_TO_DOUBLE(x)));
          else
            return std::make_pair((double) CGAL_IA_FORCE_TO_DOUBLE(x),
                                  - (double) CGAL_IA_FORCE_TO_DOUBLE(mx));
#endif
        }
    };

// Is_finite depends on platform
    class Is_finite
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& x ) const {
#ifdef CGAL_CFG_IEEE_754_BUG
          Type d = x;
          IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
          return is_finite_by_mask_long_double( p->c.H );
#else
          return (x == x) && (is_valid(x-x));
#endif
        }
    };

};

} //namespace CGAL

#endif // CGAL_LONG_DOUBLE_H
