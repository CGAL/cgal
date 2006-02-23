// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Athanasios Kakargias

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ROOT_OF_GMPXX_H
#define CGAL_ROOT_OF_GMPXX_H

#include <CGAL/gmpxx.h>
#include <CGAL/Root_of/Root_of_traits.h>
#include <CGAL/Root_of/Root_of_2.h>

namespace CGAL {

    template <  typename U1, typename U2, typename U3 >
    inline
    Root_of_2< ::mpz_class >
    make_root_of_2(const ::__gmp_expr< __gmpz_value, U1> & a,
                   const ::__gmp_expr< __gmpz_value, U2> & b,
                   const ::__gmp_expr< __gmpz_value, U3> & c,
                   bool d)
    {
      return Root_of_2< ::mpz_class >(a, b, c, d);
    }

    template < typename U1, typename U2, typename U3 >
    inline
    Root_of_2< ::mpz_class >
    make_root_of_2
                  (const ::__gmp_expr< __gmpq_value, U1> & a,
                   const ::__gmp_expr< __gmpq_value, U2> & b,
                   const ::__gmp_expr< __gmpq_value, U3> & c,
                   bool d)
    {
      typedef CGAL::Rational_traits< ::mpq_class > Rational;

      Rational r;
      assert( r.denominator(a) > 0 );
      assert( r.denominator(b) > 0 );
      assert( r.denominator(c) > 0 );

/*   const RT lcm = ( r.denominator(a) * r.denominator(b) * r.denominator(c)          )/
               ( gcd( r.denominator(a), gcd(r.denominator(b), r.denominator(c)) ) );

      RT a_ = r.numerator(a) * ( lcm / r.denominator(a) );
      RT b_ = r.numerator(b) * ( lcm / r.denominator(b) );
      RT c_ = r.numerator(c) * ( lcm / r.denominator(c) );
*/
      ::mpz_class a_ = r.numerator(a) * r.denominator(b) * r.denominator(c);
      ::mpz_class b_ = r.numerator(b) * r.denominator(a) * r.denominator(c);
      ::mpz_class c_ = r.numerator(c) * r.denominator(a) * r.denominator(b);

       return Root_of_2< ::mpz_class >(a, b, c, d); 
    }
    
    template < typename T, typename U >
    struct Root_of_traits< ::__gmp_expr<T, U> >
    {
      typedef ::mpq_class  RootOf_1;
      typedef Root_of_2< ::mpz_class > RootOf_2;
      typedef Root_of_3< ::mpz_class > RootOf_3;
      typedef Root_of_4< ::mpz_class > RootOf_4;
    };

}

// XXX : These seem necessary.
// I don't know why doing them in namespace CGAL is not enough.
using CGAL::to_double;
using CGAL::is_valid;

#endif // CGAL_ROOT_OF_GMPXX_H
