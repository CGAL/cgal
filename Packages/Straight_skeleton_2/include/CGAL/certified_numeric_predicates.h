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
inline
Uncertain<bool>
certified_is_zero(const NT& x)
{ 
  return make_uncertain( x == 0 ) ; 
}

template <class NT>
inline
Uncertain<bool>
certified_is_one(const NT& x)
{
  return make_uncertain( x == 1 ) ; 
}

template <class NT>
inline
Uncertain<bool>
certified_is_negative(const NT& x)
{ 
  return make_uncertain( x < 0 ) ; 
}

template <class NT>
inline
Uncertain<bool>
certified_is_positive(const NT& x)
{ 
  return make_uncertain( 0 < x ) ; 
}

template <class NT>
inline
Uncertain<Sign>
certified_sign(const NT& x)
{
  return make_uncertain(CGAL_NTS sign(x)); 
}

template <class NT1, class NT2>
inline
Uncertain<Comparison_result>
certified_compare(const NT1& n1, const NT2& n2)
{ 
  return make_uncertain(CGAL_NTS compare(n1,n2)); 
}


//
// Specialization for Quotient<>
//
template <class NT1, class NT2>
CGAL_MEDIUM_INLINE
Uncertain<Comparison_result>
certified_quotient_compare(const Quotient<NT1>& x, const Quotient<NT2>& y)
{ 
    // No assumptions on the sign of  den  are made

    // code assumes that SMALLER == - 1;
    CGAL_precondition( SMALLER == static_cast<Comparison_result>(-1) );

    Uncertain<Sign> xnumsign = CGAL_NTS certified_sign(x.num) ;
    Uncertain<Sign> xdensign = CGAL_NTS certified_sign(x.den) ;
    Uncertain<Sign> ynumsign = CGAL_NTS certified_sign(y.num) ;
    Uncertain<Sign> ydensign = CGAL_NTS certified_sign(y.den) ;
    
    if (  is_indeterminate(xnumsign) 
       || is_indeterminate(xdensign) 
       || is_indeterminate(ynumsign) 
       || is_indeterminate(ydensign) 
       ) 
    {
      return make_uncertain(SMALLER,LARGER);
    } 
    else
    {
      int xsign = xnumsign * xdensign ;
      int ysign = ynumsign * ydensign ;
      if (xsign == 0) return make_uncertain(static_cast<Comparison_result>(-ysign));
      if (ysign == 0) return make_uncertain(static_cast<Comparison_result>(xsign));
      // now (x != 0) && (y != 0)
      int diff = xsign - ysign;
      if (diff == 0)
      {
          int msign = xdensign * ydensign;
          NT leftop  = x.num * y.den * msign;
          NT rightop = y.num * x.den * msign;
          return certified_compare(leftop, rightop);
      }
      else
      {
          return make_uncertain((xsign < ysign) ? SMALLER : LARGER);
      }
    }  
}

template <class NT1, class NT2>
inline
Uncertain<Comparison_result>
certified_compare(const Quotient<NT1>& n1, const Quotient<NT2>& n2)
{ 
  return certified_quotient_compare(n1,n2); 
}

template <class NT1, class NT2>
inline
Uncertain<bool>
certified_is_smaller(const NT1& n1, const NT2& n2)
{ 
  return certified_compare(n1,n2) == SMALLER; 
}

template <class NT1, class NT2>
inline
Uncertain<bool>
certified_is_equal(const NT1& n1, const NT2& n2)
{ 
  return certified_compare(n1,n2) == EQUAL; 
}

template <class NT1, class NT2>
inline
Uncertain<bool>
certified_is_larger(const NT1& n1, const NT2& n2)
{ 
  return certified_compare(n1,n2) == LARGER; 
}

template <class NT1, class NT2>
inline
Uncertain<bool>
certified_is_smaller_or_equal(const NT1& n1, const NT2& n2)
{ 
  Uncertain<Comparison_result> r = certified_compare(n1,n2) ;
  return r == make_uncertain(SMALLER,EQUAL) || r == LARGER ;
}

template <class NT1, class NT2>
inline
Uncertain<bool>
certified_is_larger_or_equal(const NT1& n1, const NT2& n2)
{ 
  Uncertain<Comparison_result> r = certified_compare(n1,n2) ;
  return r == make_uncertain(EQUAL,LARGER) || r == SMALLER ;
}

CGAL_END_NAMESPACE

#endif // CGAL_CERTIFIED_NUMERIC_PREDICATES_H
 
