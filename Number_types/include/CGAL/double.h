// Copyright (c) 1999,2007  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Geert-Jan Giezeman, Michael Hemmer

#ifndef CGAL_DOUBLE_H
#define CGAL_DOUBLE_H

#include <CGAL/utils.h>
#include <CGAL/utils_classes.h>
#include <CGAL/number_utils.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <utility>
#include <cmath>
#include <math.h> // for nextafter
#include <limits>

#ifdef CGAL_USE_SSE2_FABS
#include <CGAL/sse2.h>
#endif

#ifdef _MSC_VER
#include <cfloat>
#endif

#ifdef CGAL_CFG_IEEE_754_BUG
#  include <CGAL/IEEE_754_unions.h>
#endif

namespace CGAL {

#ifdef CGAL_CFG_IEEE_754_BUG

#define CGAL_EXPONENT_DOUBLE_MASK   0x7ff00000
#define CGAL_MANTISSA_DOUBLE_MASK   0x000fffff

inline
bool
is_finite_by_mask_double(unsigned int h)
{
  unsigned int e = h & CGAL_EXPONENT_DOUBLE_MASK;
  return ( ( e ^ CGAL_EXPONENT_DOUBLE_MASK ) != 0 );
}

inline
bool
is_nan_by_mask_double(unsigned int h, unsigned int l)
{
  if ( is_finite_by_mask_double(h) )
      return false;
  return (( h & CGAL_MANTISSA_DOUBLE_MASK ) != 0) || (( l & 0xffffffff ) != 0);
}

template<>
class Is_valid< double >
  : public std::unary_function< double, bool > {
  public :
    bool operator()( const double& x ) const{
      double d = x;
      IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
      return ! ( is_nan_by_mask_double( p->c.H, p->c.L ));
    }
};

#else

template<>
class Is_valid< double >
  : public std::unary_function< double, bool > {
  public :
    bool operator()( const double& x ) const {
#ifdef _MSC_VER
      return ! _isnan(x);
#else
      return (x == x);
#endif
    }
};

#endif

template <> class Algebraic_structure_traits< double >
  : public Algebraic_structure_traits_base< double,
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false            Is_exact;
    typedef Tag_true             Is_numerical_sensitive;

    class Sqrt
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return std::sqrt( x );
        }
    };

    class Kth_root
      : public std::binary_function<int, Type, Type> {
      public:
        Type operator()( int k,
                                        const Type& x) const {
          CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
          return std::pow(x, 1.0 / double(k));
        }
    };

};


  
#ifdef CGAL_USE_SSE2_FABS   
inline double sse2fabs(double a)
{
  static CGAL_ALIGN_16 const union{
    __int64 i[2];
    __m128d m;
  } absMask = {0x7fffffffffffffff, 0x7fffffffffffffff};
	  
  __m128d temp = _mm_set1_pd(a);
  
  temp = _mm_and_pd(temp, absMask.m);
  _mm_store_sd(&a, temp);
  return a;
}

#endif

template <> class Real_embeddable_traits< double >
  : public INTERN_RET::Real_embeddable_traits_base< double , CGAL::Tag_true> {
  public:


// GCC is faster with std::fabs().
#if defined(__GNUG__) || defined(CGAL_MSVC_USE_STD_FABS) || defined(CGAL_USE_SSE2_FABS)
    class Abs
      : public std::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
#ifdef CGAL_USE_SSE2_FABS
          return sse2fabs(x);
#else         
          return std::fabs( x );
#endif
        }
    };
#endif

// Is_finite depends on platform
    class Is_finite
      : public std::unary_function< Type, bool > {
      public :
        bool operator()( const Type& x ) const {
#ifdef CGAL_CFG_IEEE_754_BUG
          Type d = x;
          IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
          return is_finite_by_mask_double( p->c.H );
#elif defined CGAL_CFG_NUMERIC_LIMITS_BUG
          return (x == x) && (is_valid(x-x));
#else
          return (x != std::numeric_limits<Type>::infinity())
              && (-x != std::numeric_limits<Type>::infinity())
              && is_valid(x);
#endif
      }
    };
};

inline
double
nextafter(double d1, double d2)
{
#ifdef CGAL_CFG_NO_NEXTAFTER
  return _nextafter(d1, d2); // works at least for VC++-7.1
#else
  return ::nextafter(d1,d2);
#endif
}

inline
bool
is_integer(double d)
{
  return CGAL::is_finite(d) && (std::ceil(d) == d);
}

// Returns a pair of integers <num,den> such that d == num/den.
inline
std::pair<double, double>
split_numerator_denominator(double d)
{
  // Note that it could probably be optimized.
  double num = d;
  double den = 1.0;
  while (std::ceil(num) != num)
  {
    num *= 2.0;
    den *= 2.0;
  }
  CGAL_postcondition(d == num/den);
  return std::make_pair(num, den);
}

} //namespace CGAL

#endif // CGAL_DOUBLE_H
