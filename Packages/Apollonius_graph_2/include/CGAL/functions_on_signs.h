// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/functions_on_signs.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_FUNCTIONS_ON_SIGNS_H
#define CGAL_FUNCTIONS_ON_SIGNS_H

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

inline
Sign
operator*(const Sign &s1, const Sign &s2)
{
  if ( s1 == ZERO || s2 == ZERO )  return ZERO;
  if ( s1 == s2 )  return POSITIVE;
  return NEGATIVE;
}

template < class RT >
Sign
sign_a_plus_b_x_sqrt_c(const RT &a, const RT &b, const RT &c)
{
  // computes the sign of quantity: a + b * sqrt(c)

  CGAL_assertion( !(CGAL_NTS is_negative(c)) );

  Sign sa = CGAL_NTS sign(a);
  if ( CGAL_NTS sign(c) == ZERO )  return sa;

  Sign sb = CGAL_NTS sign(b);
  if ( sa == sb )  return sa;
  if ( sa == ZERO )  return sb;

  return Sign( sa * CGAL_NTS compare( CGAL_NTS square(a),
				      c * CGAL_NTS square(b) )
	       );
}

template < class RT >
Sign
sign_a_x_sqrt_c_plus_b_x_sqrt_d(const RT &a, const RT &b,
				const RT &c, const RT &d)
{
  // computes the sign of quantity: a * sqrt(c) + b * sqrt(d)

  CGAL_assertion( !(CGAL_NTS is_negative(c)) );
  CGAL_assertion( !(CGAL_NTS is_negative(d)) );

  Sign sb = CGAL_NTS sign(b);
  if ( CGAL_NTS sign(d) == ZERO )  return CGAL_NTS sign(a * c);
  if ( CGAL_NTS sign(c) == ZERO )  return sb;

  Sign sa = CGAL_NTS sign(a);
  if ( sa == sb )  return sa;
  if ( sa == ZERO )  return sb;

  return Sign( sa * CGAL_NTS compare( CGAL_NTS square(a) * c,
				      CGAL_NTS square(b) * d )
	       );
}

template < class RT >
Sign
sign_a_plus_b_x_sqrt_e_plus_c_x_sqrt_f(const RT &a, const RT &b,
				       const RT &c, const RT &e,
				       const RT &f)
{
  // computes the sign of quantity: a + b * sqrt(e) + c * sqrt(f)
  
  CGAL_assertion( !(CGAL_NTS is_negative(e)) );
  CGAL_assertion( !(CGAL_NTS is_negative(f)) );

  Sign s_a_plus_b_x_sqrt_e = sign_a_plus_b_x_sqrt_c(a, b, e);
  if ( CGAL_NTS sign(f) == ZERO )  return s_a_plus_b_x_sqrt_e;

  Sign sc = CGAL_NTS sign(c);
  if ( s_a_plus_b_x_sqrt_e == sc )  return sc;
  if ( s_a_plus_b_x_sqrt_e == ZERO )  return sc;

  return s_a_plus_b_x_sqrt_e * 
    sign_a_plus_b_x_sqrt_c(CGAL_NTS square(a) + CGAL_NTS square(b) * e
			   - CGAL_NTS square(c) * f,
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
  
  CGAL_assertion( !(CGAL_NTS is_negative(e)) );
  CGAL_assertion( !(CGAL_NTS is_negative(f)) );

  Sign s_a_plus_b_sqrt_e = sign_a_plus_b_x_sqrt_c(a, b, e);
  Sign s_c_plus_d_sqrt_e = sign_a_plus_b_x_sqrt_c(c, d, e);

  if ( s_a_plus_b_sqrt_e == s_c_plus_d_sqrt_e )
    return s_a_plus_b_sqrt_e;

  if ( s_a_plus_b_sqrt_e == ZERO )
    return s_a_plus_b_sqrt_e;

  return s_a_plus_b_sqrt_e *
    sign_a_plus_b_x_sqrt_c(CGAL_NTS square(a) + CGAL_NTS square(b) * e
			   - CGAL_NTS square(c) * f
			   - CGAL_NTS square(d) * e * f,
			   RT(2) * (a * b - c * d * f),
			   e);
}

CGAL_END_NAMESPACE

#include <CGAL/more_functions_on_signs.h>

#endif // CGAL_FUNCTIONS_ON_SIGNS_H
