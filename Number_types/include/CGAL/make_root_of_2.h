// Copyright (c) 2005,2006  INRIA Sophia-Antipolis (France)
// All rights reserved.
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
// Author(s)     : Sylvain Pion, Monique Teillaud, Athanasios Kakargias

#ifndef CGAL_MAKE_ROOT_OF_2_H
#define CGAL_MAKE_ROOT_OF_2_H

#include <iostream>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Root_of_2_fwd.h>

namespace CGAL {

namespace CGALi {

    // This version is internal and can be re-used for
    // number types which also support division and sqrt().
    template < typename NT >
    NT
    make_root_of_2_sqrt(const NT &a, const NT &b, const NT &c, bool smaller)
    {
      CGAL_assertion( a != 0 );
      NT discriminant = CGAL_NTS square(b) - a*c*4;
      CGAL_assertion( discriminant >= 0 );
      NT d = CGAL_NTS sqrt(discriminant);
      if ((smaller && a>0) || (!smaller && a<0))
        d = -d;
      return (d-b)/(a*2);
    }

    // This version is internal and can be re-used for
    // number types which also support division and sqrt().
    template < typename NT >
    NT
    make_root_of_2_sqrt(const NT &a, const NT &b, const NT &c)
    {
      CGAL_assertion( c >= 0 );
      return a + b * sqrt(c);
    }

    // This version is internal and can be re-used for
    // number types which also support division and sqrt().
    template < typename NT >
    NT
    make_root_of_2_sqrt(const NT &a, const int &b, const NT &c)
    {
      CGAL_assertion( c >= 0 );
      return a + (NT(b)) * sqrt(c);
    }

    // This version is internal and can be re-used for
    // number types which are rational.
    template < typename RT, typename FT >
    Root_of_2< RT >
    make_root_of_2_rational(const FT &a, const FT &b, const FT &c, bool smaller)
    {
      typedef Rational_traits< FT > Rational;

      Rational r;
      // CGAL_assertion( r.denominator(a) > 0 );
      // CGAL_assertion( r.denominator(b) > 0 );
      // CGAL_assertion( r.denominator(c) > 0 );

/*   const RT lcm = ( r.denominator(a) * r.denominator(b) * r.denominator(c)          )/
               ( gcd( r.denominator(a), gcd(r.denominator(b), r.denominator(c)) ) );

      RT a_ = r.numerator(a) * ( lcm / r.denominator(a) );
      RT b_ = r.numerator(b) * ( lcm / r.denominator(b) );
      RT c_ = r.numerator(c) * ( lcm / r.denominator(c) );
*/
      RT a_ = r.numerator(a) * r.denominator(b) * r.denominator(c);
      RT b_ = r.numerator(b) * r.denominator(a) * r.denominator(c);
      RT c_ = r.numerator(c) * r.denominator(a) * r.denominator(b);

      return make_root_of_2(a_,b_,c_,smaller);
    }

    

// automatic dispatcher between the 2 generic versions (using Root_of_2 or
// sqrt()), if sqrt() exists (checking Has_sqrt).
template < typename RT,
           typename Has_sqrt /*= typename Number_type_traits<RT>::Has_sqrt*/ >
struct Make_root_of_2_helper
{
  typedef Root_of_2<RT> result_type;

  result_type operator()(const RT& a, const RT& b, const RT& c, bool smaller)
  const
  {
    return Root_of_2<RT>(a, b, c, smaller);
  }

  result_type operator()(const RT& a, const RT& b, const RT& c)
  const
  {
    return Root_of_2<RT>(a, b, c);
  }

  result_type operator()(const RT& a, const int &b, const RT& c)
  const
  {
    const RT ba = RT(b);
    return make_root_of_2(a, ba, c);
  }

};

// Specialization for Has_sqrt == Tag_true
template < typename RT >
struct Make_root_of_2_helper <RT, Tag_true>
{
  typedef RT result_type;

  result_type operator()(const RT& a, const RT& b, const RT& c, bool smaller)
  const
  {
    return CGALi::make_root_of_2_sqrt(a, b, c, smaller);
  }

  result_type operator()(const RT& a, const RT& b, const RT& c)
  const
  {
    return CGALi::make_root_of_2_sqrt(a, b, c);
  }

  result_type operator()(const RT& a, const int &b, const RT& c)
  const
  {
    return CGALi::make_root_of_2_sqrt(a, b, c);
  }

};

} // namespace CGALi

// Template default version generating a Root_of_2<>.
template < typename RT >
inline
typename CGALi::Make_root_of_2_helper<RT>::result_type
make_root_of_2(const RT &a, const RT &b, const RT &c, bool smaller)
{
  CGAL_assertion( a != 0 );
  return CGALi::Make_root_of_2_helper<RT>()(a, b, c, smaller);
}

// Template default version generating a Root_of_2<>.
template < typename RT >
inline
typename CGALi::Make_root_of_2_helper<RT>::result_type
make_root_of_2(const RT &a, const RT &b, const RT &c)
{
  return CGALi::Make_root_of_2_helper<RT>()(a, b, c);
}

// Template default version generating a Root_of_2<>.
template < typename RT >
inline
typename CGALi::Make_root_of_2_helper<RT>::result_type
make_root_of_2(const RT &a, const int &b, const RT &c)
{
  return CGALi::Make_root_of_2_helper<RT>()(a, b, c);
}

} // namespace CGAL

#endif // CGAL_MAKE_ROOT_OF_2_H
