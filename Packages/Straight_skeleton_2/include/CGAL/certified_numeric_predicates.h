// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/certified_numeric_predicates.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_CERTIFIED_NUMERIC_PREDICATES_H
#define CGAL_CERTIFIED_NUMERIC_PREDICATES_H

#include <CGAL/number_utils.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Uncertain.h>

CGAL_BEGIN_NAMESPACE

template <class NT>
inline Uncertain<bool> certified_is_finite(const NT& x)
{ 
  return make_uncertain(CGAL_NTS is_finite(x) ) ; 
}

template <class NT>
inline Uncertain<bool> certified_is_zero(const NT& x)
{ 
  return make_uncertain(CGAL_NTS is_zero(x) ) ; 
}

template <class NT>
inline Uncertain<bool> certified_is_one(const NT& x)
{
  return make_uncertain(CGAL_NTS is_one(x) ) ; 
}

template <class NT>
inline Uncertain<bool> certified_is_negative(const NT& x)
{ 
  return make_uncertain(CGAL_NTS is_negative(x) ) ; 
}

template <class NT>
inline Uncertain<bool> certified_is_positive(const NT& x)
{ 
  return make_uncertain(CGAL_NTS is_positive(x) ) ; 
}

template <class NT>
inline Uncertain<Sign> certified_sign(const NT& x)
{
  return make_uncertain(CGAL_NTS sign(x)); 
}

template <class NT1, class NT2>
inline Uncertain<Comparison_result> certified_compare(const NT1& n1, const NT2& n2)
{ 
  return make_uncertain(CGAL_NTS compare(n1,n2)); 
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_smaller(const NT1& n1, const NT2& n2)
{ 
  return certified_compare(n1,n2) == make_uncertain(SMALLER); 
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_equal(const NT1& n1, const NT2& n2)
{ 
  return certified_compare(n1,n2) == make_uncertain(EQUAL); 
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_larger(const NT1& n1, const NT2& n2)
{ 
  return certified_compare(n1,n2) == make_uncertain(LARGER); 
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_smaller_or_equal(const NT1& n1, const NT2& n2)
{ 
  Uncertain<Comparison_result> r = certified_compare(n1,n2) ;
  return r == make_uncertain(SMALLER) || r == make_uncertain(EQUAL) ;
}

template <class NT1, class NT2>
inline Uncertain<bool> certified_is_larger_or_equal(const NT1& n1, const NT2& n2)
{ 
  Uncertain<Comparison_result> r = certified_compare(n1,n2) ;
  return r == make_uncertain(LARGER) || r == make_uncertain(EQUAL) ;
}


CGAL_END_NAMESPACE

#endif // CGAL_CERTIFIED_NUMERIC_PREDICATES_H
 
