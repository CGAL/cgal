// Copyright (c) 1999,2007
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
// Author(s)     : Geert-Jan Giezeman, Michael Hemmer


#ifndef CGAL_FLOAT_H
#define CGAL_FLOAT_H

#include <CGAL/utils.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>

#include <cmath> // std::sqrt, std::pow

#ifdef CGAL_CFG_IEEE_754_BUG
#  include <CGAL/IEEE_754_unions.h>
#endif

namespace CGAL {

#ifdef CGAL_CFG_IEEE_754_BUG

#define CGAL_EXPONENT_FLOAT_MASK   0x7f800000
#define CGAL_MANTISSA_FLOAT_MASK   0x007fffff

inline
bool
is_finite_by_mask_float(unsigned int u)
{
  unsigned int e = u & CGAL_EXPONENT_FLOAT_MASK;
  return ( (e ^ CGAL_EXPONENT_FLOAT_MASK) != 0);
}

inline
bool
is_nan_by_mask_float(unsigned int u)
{
  if ( is_finite_by_mask_float(u) ) return false;
  // unsigned int m = u & CGAL_MANTISSA_FLOAT_MASK;
  return ( (u & CGAL_MANTISSA_FLOAT_MASK) != 0);
}

template<>
class Is_valid< float >
  : public CGAL::cpp98::unary_function< float, bool > {
  public :
    bool operator()( const float& x ) const {
      float f = x;
      IEEE_754_float* p = reinterpret_cast<IEEE_754_float*>(&f);
      return !is_nan_by_mask_float( p->c );
    }
};

#else

template<>
class Is_valid< float >
  : public CGAL::cpp98::unary_function< float, bool > {
  public :
    bool operator()( const float& x ) const {
      return (x == x);
    }
};

#endif

template <> class Algebraic_structure_traits< float >
  : public Algebraic_structure_traits_base< float,
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
      : public CGAL::cpp98::binary_function<int, Type, Type> {
      public:
        Type operator()( int k, const Type& x) const {
          CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
          return (Type) std::pow(double(x), 1.0 / double(k));
        };
    };

};

template <> class Real_embeddable_traits< float >
  : public INTERN_RET::Real_embeddable_traits_base< float , CGAL::Tag_true> {
public:
// Is_finite depends on platform
    class Is_finite
      : public CGAL::cpp98::unary_function< Type, bool > {
      public:
        bool operator()( const Type& x ) const {

#if defined CGAL_CFG_IEEE_754_BUG
          Type f = x;
          IEEE_754_float* p = reinterpret_cast<IEEE_754_float*>(&f);
          return is_finite_by_mask_float( p->c );
#else
          return std::isfinite(x);
#endif
        }
    };

};

} //namespace CGAL

#endif // CGAL_FLOAT_H
