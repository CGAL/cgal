// Copyright (c) 2001,2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion
 

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include <cassert>

#include "CGAL/Precise_numbers.h"
#define TEST_FILENAME "Test-Simple_cartesian-IO.out"
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

#include "CGAL/_test_mf_plane_3_to_2d.h"

int
main()
{
  typedef   CGAL::Simple_cartesian<double>     Clsd;
  std::cout << "Testing IO with Simple_cartesian<double> :" << std::endl;
  _test_io( Clsd() );

  typedef   CGAL::Simple_cartesian<CGAL::Quotient<Precise_integer> >     Cls;
  std::cout << "Testing 2d with Simple_cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d with Simple_cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_3( Cls() );

  std::cout << "Testing new 2d with Simple_cartesian<Quotient<Precise_integer>>:";
  std::cout << std::endl;
  test_new_2( Cls() );
  std::cout << "Testing new 3d with Simple_cartesian<Quotient<Precise_integer>>:";
  std::cout << std::endl;
  test_new_3( Cls() );

  std::cout << "Testing new parts with Simple_cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d with Simple_cartesian<Quotient<Precise_integer>>:";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );

  return 0;
}
