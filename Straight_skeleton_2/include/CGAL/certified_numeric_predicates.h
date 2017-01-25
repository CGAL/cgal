// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_CERTIFIED_NUMERIC_PREDICATES_H
#define CGAL_CERTIFIED_NUMERIC_PREDICATES_H

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/number_utils.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Uncertain.h>

namespace CGAL {

inline Uncertain<bool> logical_or ( Uncertain<bool> a, Uncertain<bool> b ) { return a | b ; }
inline Uncertain<bool> logical_and( Uncertain<bool> a, Uncertain<bool> b ) { return a & b ; }

inline Uncertain<bool> logical_or ( Uncertain<bool> a, Uncertain<bool> b, Uncertain<bool> c ) { return a | b | c ; }
inline Uncertain<bool> logical_and( Uncertain<bool> a, Uncertain<bool> b, Uncertain<bool> c ) { return a & b & c ; }

template <class NT>
inline Uncertain<bool> certified_is_zero(const NT& x)
{
  return CGAL_NTS is_valid(x) ? make_uncertain(CGAL_NTS is_zero(x)) : Uncertain<bool>::indeterminate() ; 
}

template <class NT>
inline Uncertain<bool> certified_is_not_zero(const NT& x)
{
  return CGAL_NTS is_valid(x) ? make_uncertain(!CGAL_NTS is_zero(x)) : Uncertain<bool>::indeterminate() ; 
}

template <class NT>
inline Uncertain<bool> certified_is_one(const NT& x)
{
  return CGAL_NTS is_valid(x) ? make_uncertain(CGAL_NTS is_one(x) ) : Uncertain<bool>::indeterminate() ;
}

template <class NT>
inline Uncertain<bool> certified_is_negative(const NT& x)
{
  return CGAL_NTS is_valid(x) ? make_uncertain(CGAL_NTS is_negative(x) ) : Uncertain<bool>::indeterminate()  ;
}

template <class NT>
inline Uncertain<bool> certified_is_positive(const NT& x)
{
  return CGAL_NTS is_valid(x) ? make_uncertain(CGAL_NTS is_positive(x) ) : Uncertain<bool>::indeterminate() ;
}

template <class NT>
inline Uncertain<Sign> certified_sign(const NT& x)
{
  return CGAL_NTS is_valid(x) ? make_uncertain(CGAL_NTS sign(x)) : Uncertain<Sign>::indeterminate() ;
}

template <class NT1, class NT2>
inline Uncertain<Comparison_result> certified_compare(const NT1& n1, const NT2& n2)
{
  return CGAL_NTS is_valid(n1) && CGAL_NTS is_valid(n2) ? make_uncertain(CGAL_NTS compare(n1,n2)) 
                                                        : Uncertain<Comparison_result>::indeterminate()  ;
}

inline Uncertain<bool> certified_is_smaller( Uncertain<Comparison_result> c )
{
  return c == SMALLER;
}

inline Uncertain<bool> certified_is_equal( Uncertain<Comparison_result> c )
{
  return c == EQUAL;
}

inline Uncertain<bool> certified_is_larger( Uncertain<Comparison_result> c )
{
  return c == LARGER;
}

inline Uncertain<bool> certified_is_smaller_or_equal( Uncertain<Comparison_result> c )
{
  return logical_or( c == SMALLER , c == EQUAL ) ;
}

inline Uncertain<bool> certified_is_larger_or_equal( Uncertain<Comparison_result> c )
{
  return logical_or( c == LARGER , c == EQUAL ) ;
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_smaller(const NT1& n1, const NT2& n2)
{
  return certified_is_smaller(certified_compare(n1,n2));
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_equal(const NT1& n1, const NT2& n2)
{
  return certified_is_equal(certified_compare(n1,n2));
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_larger(const NT1& n1, const NT2& n2)
{
  return certified_is_larger(certified_compare(n1,n2));
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_smaller_or_equal(const NT1& n1, const NT2& n2)
{
  return certified_is_smaller_or_equal(certified_compare(n1,n2)) ;
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_larger_or_equal(const NT1& n1, const NT2& n2)
{
  return certified_is_larger_or_equal(certified_compare(n1,n2)) ;
}

template <class NT>
inline Uncertain<Sign> certified_sign_of_determinant2x2( const NT& a00
                                                       , const NT& a01
                                                       , const NT& a10
                                                       , const NT& a11
                                                       )
{
  return certified_compare(a00*a11, a10*a01) ;
}

} // end namespace CGAL

#endif // CGAL_CERTIFIED_NUMERIC_PREDICATES_H

