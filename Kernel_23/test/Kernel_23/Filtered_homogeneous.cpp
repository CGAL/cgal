// Copyright (c) 2003  
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
 

#include <CGAL/Homogeneous.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <cassert>

#include "CGAL/Precise_numbers.h"

#define TEST_FILENAME "Test-Filtered_homogeneous-IO.out"
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
#if 1
  std::cout << " NOT TESTED AT THE MOMENT !!! " << std::endl;
#else
  typedef   CGAL::Homogeneous<double>                             Clsdb;
  // typedef   CGAL::Homogeneous<Precise_integer>     Clsb;
  typedef   CGAL::Homogeneous<CGAL::MP_Float>                     Clsb;
  typedef   CGAL::Filtered_kernel<Clsdb>                          Clsd;
  typedef   CGAL::Filtered_kernel<Clsb>                           Cls;

  std::cout <<
   // "Testing with Filtered_kernel<Homogeneous<Precise_integer>>:"
   "Testing with Filtered_kernel<Homogeneous<MP_Float>>:"
            << std::endl;
  std::cout << "Testing IO with F_k<Homogeneous<double>>:" << std::endl;
  _test_io( Clsd() );

  std::cout << "Testing 2d :";
  std::cout << std::endl;
  _test_2( Cls() );

  std::cout << "Testing 3d :";
  std::cout << std::endl;
  _test_3( Cls() );

  std::cout << "Testing new 2d :";
  std::cout << std::endl;
  test_new_2( Cls() );

  std::cout << "Testing new 3d :";
  std::cout << std::endl;
  test_new_3( Cls() );

  std::cout << "Testing new parts :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  _test_cls_iso_cuboid_3( Cls() );
  _test_angle( Cls() );

  std::cout << "Testing 3d-2d :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( Cls() );
  
#endif
  return 0;
}
