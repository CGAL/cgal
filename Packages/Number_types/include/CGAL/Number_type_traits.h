// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Susan Hert, Michael Hoffmann
 

#ifndef CGAL_NUMBER_TYPE_TRAITS_H
#define CGAL_NUMBER_TYPE_TRAITS_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_NEW_NT_TRAITS

template<class NT> class Number_type_traits;

namespace CGALi {

  template < typename NT >
  struct Default_ring_number_type_traits
  {
    typedef Tag_false      Has_gcd;
    typedef Tag_false      Has_division;
    typedef Tag_false      Has_sqrt;
    typedef Tag_false      Has_rational_traits;

    //    typedef typename NT::Has_exact_ring_operations Has_exact_ring_operations;
    //    typedef typename NT::Has_exact_division        Has_exact_division;
    //    typedef typename NT::Has_exact_sqrt            Has_exact_sqrt;

    //    typedef typename NT::Has_simplify          Has_simplify;


    static inline bool is_zero(const NT& x) {
      return x == 0;
    }

    static inline bool is_one(const NT& x) {
      return x == 1;
    }

    static inline bool is_negative(const NT& x) {
      return x < 0;
    }

    static inline bool is_positive(const NT& x) {
      return x > 0;
    }

    static inline Sign sign(const NT& x) {
      return (x < 0) ? NEGATIVE : (0 < x) ? POSITIVE : ZERO;
    }

    static inline NT abs(const NT& x) {
      return (x < 0) ? (-x) : x;
    }

    static inline Comparison_result
    compare(const NT& x1, const NT& x2) {
      return (x1 < x2) ? SMALLER : (x2 < x1) ? LARGER : EQUAL;
    }

    static inline NT square(const NT& x) {
      return x * x;
    }

    static inline io_Operator io_tag(const NT&) {
      return io_Operator();
    }
  };


  template < typename NT >
  struct Default_euclidean_ring_number_type_traits
    : public Default_ring_number_type_traits<NT>
  {
    typedef  Tag_true   Has_gcd;

    static inline NT gcd(const NT& n1, const NT& n2) {
      CGAL_precondition( !is_zero(n2) );
      NT x = abs(n1);
      NT y = abs(n2);
      do {
	x %= y;
	if ( is_zero(x) ) { return y; }
	y %= x;
      } while ( is_positive(y) );
      return x;
    }

    static inline NT div(const NT& n1, const NT& n2) {
      return n1 / n2;
    }
  };


  template < typename NT >
  struct Default_field_number_type_traits
    : public Default_ring_number_type_traits<NT>
  {
    typedef Tag_true  Has_division;
    //    typedef typename NT::Has_rational_traits  Has_rational_traits;
  };

} // namespace CGALi

#else // CGAL_NEW_NT_TRAITS

template < class NT >
struct Number_type_traits {
  typedef typename NT::Has_gcd       Has_gcd;
  typedef typename NT::Has_division  Has_division;
  typedef typename NT::Has_sqrt      Has_sqrt;
};

#endif // CGAL_NEW_NT_TRAITS

template < class Rational >
struct Rational_traits {
  typedef typename Rational::NT RT;

  RT numerator   (const Rational & r) const { return r.numerator(); }
  RT denominator (const Rational & r) const { return r.denominator(); }
  
  Rational make_rational(const RT & n, const RT & d) const
  { return Rational(n, d); } 
};

// number type tags
struct Ring_tag {};
struct Euclidean_ring_tag {};
struct Field_tag {};
struct Sqrt_field_tag {};

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_TYPE_TRAITS_H
