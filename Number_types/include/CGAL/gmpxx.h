// Copyright (c) 2002,2003  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_GMPXX_H
#define CGAL_GMPXX_H

#include <CGAL/basic.h>

#include <gmpxx.h>
#include <utility>

#include <CGAL/mpz_class.h>
#include <CGAL/mpq_class.h>
#include <CGAL/gmpxx_coercion_traits.h>

//#include <CGAL/Root_of_traits.h>
//#include <CGAL/Root_of_2.h>

#include <CGAL/functional_base.h> // Unary_function, Binary_function

// This file gathers the necessary adaptors so that the following
// C++ number types that come with GMP can be used by CGAL :
// - mpz_class
// - mpq_class

// - mpf_class support is commented out until to_interval() is implemented.
//   It is probably not very useful with CGAL anyway.

// Note that GMP++ use the expression template mechanism, which makes things
// a little bit complicated in order to make square(x+y) work for example.
// Reading gmpxx.h shows that ::__gmp_expr<T, T> is the mp[zqf]_class proper,
// while ::__gmp_expr<T, U> is the others "expressions".

CGAL_BEGIN_NAMESPACE


/* FIX ME: THERE IS NO CONSTRUCTOR FT(x,y)
           AVALIABLE FOR THIS TYPE

namespace CGALi {

inline
Root_of_2< ::mpz_class >
make_root_of_2_gmpxx(const ::mpz_class & a,
                     const ::mpz_class & b,
                     const ::mpz_class & c,
                     bool d)
{
  return Root_of_2< ::mpz_class >(a, b, c, d);
}

inline
Root_of_2< ::mpz_class >
make_root_of_2_gmpxx(const ::mpq_class & a,
                     const ::mpq_class & b,
                     const ::mpq_class & c,
                     bool d)
{
    typedef CGAL::Rational_traits< ::mpq_class > Rational;

    Rational r;
    CGAL_assertion( r.denominator(a) > 0 );
    CGAL_assertion( r.denominator(b) > 0 );
    CGAL_assertion( r.denominator(c) > 0 );
*/
/*   const RT lcm = ( r.denominator(a) * r.denominator(b) * r.denominator(c)          )/
               ( gcd( r.denominator(a), gcd(r.denominator(b), r.denominator(c)) ) );

      RT a_ = r.numerator(a) * ( lcm / r.denominator(a) );
      RT b_ = r.numerator(b) * ( lcm / r.denominator(b) );
      RT c_ = r.numerator(c) * ( lcm / r.denominator(c) );
*/
/*    ::mpz_class a_ = r.numerator(a) * r.denominator(b) * r.denominator(c);
    ::mpz_class b_ = r.numerator(b) * r.denominator(a) * r.denominator(c);
    ::mpz_class c_ = r.numerator(c) * r.denominator(a) * r.denominator(b);

    return Root_of_2< ::mpz_class >(a, b, c, d);
}

} // CGALi

template < typename T, typename U1, typename U2, typename U3 >
inline
typename Root_of_traits< ::__gmp_expr<T, T> >::RootOf_2
make_root_of_2(const ::__gmp_expr< T, U1> & a,
               const ::__gmp_expr< T, U2> & b,
               const ::__gmp_expr< T, U3> & c,
               bool d)
{
  return CGALi::make_root_of_2_gmpxx(a, b, c, d);
}

template < typename T, typename U >
struct Root_of_traits< ::__gmp_expr<T, U> >
{
  typedef ::mpq_class               RootOf_1;
  typedef Root_of_2< ::mpz_class >  RootOf_2;
};*/


CGAL_END_NAMESPACE
// XXX : These seem necessary.
// I don't know why doing them in namespace CGAL is not enough.
// using CGAL::to_double;
// using CGAL::is_valid;

#endif // CGAL_GMPXX_H
