// ============================================================================
//
// Copyright (c) 2001,2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : test_all.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>

#include "../Kernel/include/CGAL/Precise_numbers.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include "../Kernel/include/CGAL/_test_io.h"

#include "../Kernel/include/CGAL/_test_2.C"
#include "../Kernel/include/CGAL/_test_3.C"

#include "../Kernel/include/CGAL/_test_new_2.h"
#include "../Kernel/include/CGAL/_test_new_3.h"

#include "../Kernel/include/CGAL/_test_fct_points_implicit_sphere.h"
#include "../Kernel/include/CGAL/_test_orientation_and_bounded_side.h"
#include "../Kernel/include/CGAL/_test_fct_constructions_2.h"
#include "../Kernel/include/CGAL/_test_fct_constructions_3.h"
#include "../Kernel/include/CGAL/_test_fct_point_3.h"
#include "../Kernel/include/CGAL/_test_fct_coplanar_3.h"
#include "../Kernel/include/CGAL/_test_cls_iso_cuboid_3.h"
#include "../Kernel/include/CGAL/_test_angle.h"

#include "../Kernel/include/CGAL/_test_mf_plane_3_to_2d.h"

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
