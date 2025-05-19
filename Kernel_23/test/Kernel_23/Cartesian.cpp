// Copyright (c) 2001,2002
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
// Author(s)     : Sylvain Pion


#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <cassert>

#include "CGAL/Precise_numbers.h"
#include "CGAL/_test_cls_kernel.h"
#define TEST_FILENAME "Test-Cartesian-IO.out"
#include "CGAL/_test_io.h"
#include "CGAL/_test_2.h"
#include "CGAL/_test_3.h"

#include "CGAL/_test_new_2.h"
#include "CGAL/_test_new_3.h"

#include "CGAL/_test_fct_points_implicit_sphere.h"
#include "CGAL/_test_orientation_and_bounded_side.h"
#include "CGAL/_test_fct_constructions_2.h"
#include "CGAL/_test_fct_constructions_3.h"
#include "CGAL/_test_fct_point_3.h"
#include "CGAL/_test_fct_coplanar_3.h"
#include "CGAL/_test_cls_iso_cuboid_3.h"
#include "CGAL/_test_angle.h"
#include "CGAL/_test_cls_circle_3.h"

#include "CGAL/_test_mf_plane_3_to_2d.h"

int
main()
{
  typedef   CGAL::Cartesian<CGAL::Quotient<Precise_integer> >     Cls;
  typedef   CGAL::Cartesian<double>     Clsd;


  std::cout << "Testing nested types with Cartesian<double> :" << std::endl;
  _test_kernel( Clsd() );

  std::cout << "Testing IO with Cartesian<double> :" << std::endl;
  _test_io( Clsd() );

  std::cout << "Testing 2d with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_3( Cls() );
  _test_cls_circle_3( Cls() );

  std::cout << "Testing new 2d with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  test_new_2( Cls() );
  _test_cls_new_2( Cls() );

  std::cout << "Testing new 3d with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  test_new_3( Cls() );

  std::cout << "Testing new parts with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d with Cartesian<Quotient<Precise_integer> > :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );

  std::cout << "All tests done" << std::endl;
  return 0;
}



