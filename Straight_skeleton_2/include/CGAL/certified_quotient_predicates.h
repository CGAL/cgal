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
#ifndef CGAL_CERTIFIED_QUOTIENT_PREDICATES_H 
#define CGAL_CERTIFIED_QUOTIENT_PREDICATES_H
 
#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/Quotient.h>

namespace CGAL { 

template <class NT>
inline Uncertain<bool> certified_quotient_is_positive(const Quotient<NT>& x) 
{
  Uncertain<Sign> signum = CGAL_NTS certified_sign(x.num) ;
  Uncertain<Sign> sigden = CGAL_NTS certified_sign(x.den) ;
  Uncertain<Sign> zero(ZERO);
  return ( signum != zero ) & ( signum == sigden );
}

template <class NT> 
inline Uncertain<bool> certified_quotient_is_negative(const Quotient<NT>& x)
{
  Uncertain<Sign> signum = CGAL_NTS certified_sign(x.num) ;
  Uncertain<Sign> sigden = CGAL_NTS certified_sign(x.den) ;
  Uncertain<Sign> zero(ZERO);

  return ( signum != zero ) & ( signum != sigden );
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
  // No assumptions on the sign of  den  are made

  return CGAL_NTS certified_sign(x.num) * CGAL_NTS certified_sign(x.den);
}

template <class NT1, class NT2>
CGAL_MEDIUM_INLINE
Uncertain<Comparison_result> certified_quotient_compare(const Quotient<NT1>& x, const Quotient<NT2>& y)
{
  Uncertain<Comparison_result> r = Uncertain<Comparison_result>::indeterminate();

  // No assumptions on the sign of  den  are made

  // code assumes that SMALLER == - 1;
  CGAL_precondition( SMALLER == static_cast<Comparison_result>(-1) );

  Uncertain<Sign> xnumsign = CGAL_NTS certified_sign(x.num) ;
  Uncertain<Sign> xdensign = CGAL_NTS certified_sign(x.den) ;
  Uncertain<Sign> ynumsign = CGAL_NTS certified_sign(y.num) ;
  Uncertain<Sign> ydensign = CGAL_NTS certified_sign(y.den) ;

  if (  is_certain(xnumsign)
     && is_certain(xdensign)
     && is_certain(ynumsign)
     && is_certain(ydensign)
     )
  {
    int xsign = xnumsign * xdensign ;
    int ysign = ynumsign * ydensign ;
    if (xsign == 0) return static_cast<Comparison_result>(-ysign);
    if (ysign == 0) return static_cast<Comparison_result>(xsign);
    // now (x != 0) && (y != 0)
    int diff = xsign - ysign;
    if (diff == 0)
    {
      int msign = xdensign * ydensign;
      NT1 leftop  = x.num * y.den * msign;
      NT1 rightop = y.num * x.den * msign;
      r = certified_compare(leftop, rightop);
    }
    else
    {
      r = (xsign < ysign) ? SMALLER : LARGER;
    }
  }

  return r ;
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

} // end namespace CGAL

#endif // CGAL_CERTIFIED_QUOTIENT_PREDICATES_H

