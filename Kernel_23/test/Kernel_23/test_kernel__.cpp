// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Stefan Schirra

#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Simple_cartesian.h>

#include <cassert>

#include <CGAL/_test_cls_object.h>
#include <CGAL/_test_cls_quotient.h>
#include <CGAL/_test_fct_point_conversion.h>
#include <CGAL/_test_fct_determinant.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Precise_numbers.h>
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#endif
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
#endif

int
main()
{
  std::cout << "Testing miscellaneous" << std::endl;
  typedef   CGAL::Homogeneous<Precise_integer>                 H_Cls;
  typedef   CGAL::Simple_homogeneous<Precise_integer>          SH_Cls;
  typedef   CGAL::Cartesian<Precise_rational>                  C_Cls;
  typedef   CGAL::Simple_cartesian<Precise_rational>           SC_Cls;
  _test_cls_object( C_Cls() );
  _test_cls_object( H_Cls() );
  _test_cls_object( SC_Cls() );
  _test_cls_object( SH_Cls() );
  _test_cls_quotient( int(1) );
  _test_cls_quotient( double(1.0) );
  _test_cls_quotient( long(1) );
  _test_cls_quotient( short(1) );
  _test_fct_determinant( Precise_integer(1) );
  _test_fct_point_conversion( int(1) );
  _test_fct_point_conversion( double(1.0) );
#ifdef CGAL_USE_GMP
  _test_cls_quotient( CGAL::Gmpz(1) );
  _test_fct_point_conversion( CGAL::Gmpz(1) );
#endif
#ifdef CGAL_USE_LEDA
  _test_cls_quotient( leda_integer(1) );
  _test_fct_point_conversion( leda_integer(1) );
#endif
  _test_cls_quotient( CGAL::MP_Float(1) );
  _test_fct_point_conversion( CGAL::MP_Float(1) );
  return 0;
}
