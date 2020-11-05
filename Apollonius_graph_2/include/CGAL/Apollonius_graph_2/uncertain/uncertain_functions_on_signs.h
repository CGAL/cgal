// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_UNCERTAIN_FUNCTIONS_ON_SIGNS_H
#define CGAL_UNCERTAIN_FUNCTIONS_ON_SIGNS_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/enum.h>
#include <CGAL/Uncertain.h>

namespace CGAL {

template < class RT >
Uncertain<Sign>
uncertain_sign_a_plus_b_x_sqrt_c(const RT &a, const RT &b, const RT &c)
{
  // computes the sign of quantity: a + b * sqrt(c)

  CGAL_assertion( !(CGAL::is_negative(c)) );

  Uncertain<Sign> sa = CGAL::sign(a);
  if ( is_indeterminate(sa) ) {
    return Uncertain<Sign>::indeterminate();
  }

  Uncertain<Sign> sc = CGAL::sign(c);
  if ( is_indeterminate(sc) ) {
    return Uncertain<Sign>::indeterminate();
  }
  if ( sc == ZERO )  return sa;

  Uncertain<Sign> sb = CGAL::sign(b);
  if ( is_indeterminate(sb) ) {
    return Uncertain<Sign>::indeterminate();
  }
  if ( sa == sb )  return sa;
  if ( sa == ZERO )  return sb;

  Uncertain<Comparison_result> cr =
    CGAL::compare( CGAL::square(a), c * CGAL::square(b) );
  if ( is_indeterminate(cr) ) {
    return Uncertain<Sign>::indeterminate();
  }

  return sa * cr;
}

#if 0
template < class RT >
Sign
sign_a_x_sqrt_c_plus_b_x_sqrt_d(const RT &a, const RT &b,
                                const RT &c, const RT &d)
{
  // computes the sign of quantity: a * sqrt(c) + b * sqrt(d)

  CGAL_assertion( !(CGAL::is_negative(c)) );
  CGAL_assertion( !(CGAL::is_negative(d)) );

  Sign sb = CGAL::sign(b);
  if ( CGAL::sign(d) == ZERO )  return CGAL::sign(a * c);
  if ( CGAL::sign(c) == ZERO )  return sb;

  Sign sa = CGAL::sign(a);
  if ( sa == sb )  return sa;
  if ( sa == ZERO )  return sb;

  return sa * CGAL::compare( CGAL::square(a) * c,
                             CGAL::square(b) * d );
}

template < class RT >
Sign
sign_a_plus_b_x_sqrt_e_plus_c_x_sqrt_f(const RT &a, const RT &b,
                                       const RT &c, const RT &e,
                                       const RT &f)
{
  // computes the sign of quantity: a + b * sqrt(e) + c * sqrt(f)

  CGAL_assertion( !(CGAL::is_negative(e)) );
  CGAL_assertion( !(CGAL::is_negative(f)) );

  Sign s_a_plus_b_x_sqrt_e = sign_a_plus_b_x_sqrt_c(a, b, e);
  if ( CGAL::sign(f) == ZERO )  return s_a_plus_b_x_sqrt_e;

  Sign sc = CGAL::sign(c);
  if ( s_a_plus_b_x_sqrt_e == sc )  return sc;
  if ( s_a_plus_b_x_sqrt_e == ZERO )  return sc;

  return s_a_plus_b_x_sqrt_e *
    sign_a_plus_b_x_sqrt_c(CGAL::square(a) + CGAL::square(b) * e
                           - CGAL::square(c) * f,
                           RT(2) * a * b, e);
}

template < class RT >
Sign
sign_a_plus_b_x_sqrt_e_plus_c_x_sqrt_f_plus_d_sqrt_e_x_f(const RT &a,
                                                         const RT &b,
                                                         const RT &c,
                                                         const RT &d,
                                                         const RT &e,
                                                         const RT &f)
{
  // computes the sign of quantity:
  //           a + b * sqrt(e) + c * sqrt(f) + d * sqrt(e * f)

  CGAL_assertion( !(CGAL::is_negative(e)) );
  CGAL_assertion( !(CGAL::is_negative(f)) );

  Sign s_a_plus_b_sqrt_e = sign_a_plus_b_x_sqrt_c(a, b, e);
  Sign s_c_plus_d_sqrt_e = sign_a_plus_b_x_sqrt_c(c, d, e);

  if ( s_a_plus_b_sqrt_e == s_c_plus_d_sqrt_e )
    return s_a_plus_b_sqrt_e;

  if ( s_a_plus_b_sqrt_e == ZERO )
    return s_a_plus_b_sqrt_e;

  return s_a_plus_b_sqrt_e *
    sign_a_plus_b_x_sqrt_c(CGAL::square(a) + CGAL::square(b) * e
                           - CGAL::square(c) * f
                           - CGAL::square(d) * e * f,
                           RT(2) * (a * b - c * d * f),
                           e);
}
#endif

} //namespace CGAL


#endif // CGAL_UNCERTAIN_FUNCTIONS_ON_SIGNS_H
