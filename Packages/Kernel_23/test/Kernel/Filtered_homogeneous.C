// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : test/Kernel/Filtered_homogeneous.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ============================================================================
 

#include <CGAL/Homogeneous.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <cassert>

#include "CGAL/Precise_numbers.h"

#include "CGAL/_test_io.h"

#include "CGAL/_test_2.C"
#include "CGAL/_test_3.C"

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
  
  return 0;
}
