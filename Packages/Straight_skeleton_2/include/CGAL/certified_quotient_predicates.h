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
// file          : include/CGAL/certified_quotient_predicates.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_CERTIFIED_QUOTIENT_PREDICATES_H
#define CGAL_CERTIFIED_QUOTIENT_PREDICATES_H

#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/Quotient.h>

CGAL_BEGIN_NAMESPACE

template <class NT>
inline Uncertain<bool> certified_quotient_is_positive(const Quotient<NT>& x)
{
  return CGAL_NTS certified_sign(x.num) == CGAL_NTS certified_sign(x.den) ;
}

template <class NT>
inline Uncertain<bool> certified_quotient_is_negative(const Quotient<NT>& x)
{
  return CGAL_NTS certified_sign(x.num) != CGAL_NTS certified_sign(x.den) ;
}

template <class NT>
inline Uncertain<bool> certified_quotient_is_zero(const Quotient<NT>& x)
{ 
  return CGAL_NTS certified_is_zero(x.num) ; 
}

template <class NT>
CGAL_MEDIUM_INLINE
Uncertain<Sign> certified_quotient_sign(const Quotient<NT>& x)
{
  Uncertain<Sign> signum = CGAL_NTS certified_sign(x.num) ;
  Uncertain<Sign> sigden = CGAL_NTS certified_sign(x.den) ;
  Uncertain<bool> same = signum == sigden ;
  if ( is_indeterminate(same) )
       return make_uncertain(NEGATIVE,POSITIVE);
  else return same ? signum : make_uncertain(NEGATIVE);
}

template <class NT1, class NT2>
CGAL_MEDIUM_INLINE
Uncertain<Comparison_result> certified_quotient_compare(const Quotient<NT1>& x, const Quotient<NT2>& y)
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
        NT1 leftop  = x.num * y.den * msign;
        NT1 rightop = y.num * x.den * msign;
        return certified_compare(leftop, rightop);
    }
    else
    {
        return make_uncertain((xsign < ysign) ? SMALLER : LARGER);
    }
  }  
}

template <class NT>
inline Uncertain<bool> certified_is_zero(const Quotient<NT>& n)
{ 
  return certified_quotient_is_zero(n); 
}
template <class NT>
inline Uncertain<bool> certified_is_positive(const Quotient<NT>& n)
{ 
  return certified_quotient_is_positive(n); 
}
template <class NT>
inline Uncertain<bool> certified_is_negative(const Quotient<NT>& n)
{ 
  return certified_quotient_is_negative(n); 
}
template <class NT>
inline Uncertain<Sign> certified_sign(const Quotient<NT>& n)
{ 
  return certified_quotient_sign(n); 
}

template <class NT1, class NT2>
inline Uncertain<Comparison_result> certified_compare(const Quotient<NT1>& n1, const Quotient<NT2>& n2)
{ 
  return certified_quotient_compare(n1,n2); 
}

CGAL_END_NAMESPACE

#endif // CGAL_CERTIFIED_QUOTIENT_PREDICATES_H
 
