// Copyright (c) 2003  INRIA Sophia-Antipolis (France) and
//                     Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//           Athanasios Kakargias <grad0460@di.uoa.gr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (CGAL - Effective Computational Geometry for Curves and Surfaces) 

// file : include/CGAL/Root_of/gmpxx.h

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
