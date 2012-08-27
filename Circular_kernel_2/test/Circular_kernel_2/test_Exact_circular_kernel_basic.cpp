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
// $URL: svn+ssh://pmachado@scm.gforge.inria.fr/svn/cgal/trunk/Kernel_23/test/Kernel_23/Cartesian.cpp $
// $Id: Cartesian.cpp 43806 2008-06-26 14:26:49Z pmachado $
// 
//
// Author(s)     : Sylvain Pion
 
#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Filtered_bbox_circular_kernel_2.h>
#include <CGAL/Quotient.h>
#include <cassert>

#include "CGAL/Precise_numbers.h"
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
	typedef CGAL::Quotient<Precise_integer>                      NT1;
	typedef double                                               NT2;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Cartesian<NT2>                                 Linear_k2;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT2>          Algebraic_k2;
  typedef CGAL::Circular_kernel_2<Linear_k1,Algebraic_k1>      Clsu;
  typedef CGAL::Circular_kernel_2<Linear_k2,Algebraic_k2>      Clsdu;
  typedef CGAL::Filtered_bbox_circular_kernel_2<Clsu>          Cls;
  typedef CGAL::Filtered_bbox_circular_kernel_2<Clsdu>         Clsd;

  std::cout << "Testing IO with Cartesian<double> :" << std::endl;
  _test_io( Clsd() );

  std::cout << "Testing 2d with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_2( Cls() );


  std::cout << "Testing new 2d with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  test_new_2( Cls() );
  _test_cls_new_2( Cls() );


  std::cout << "Testing new parts with Cartesian<Quotient<Precise_integer>> :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_angle( Cls() );

  std::cout << "All tests done" << std::endl;
  return 0;
}



