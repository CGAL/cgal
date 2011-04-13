// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : test_new_partsS_.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Precise_numbers.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_fct_points_implicit_sphere.h>
#include <CGAL/_test_orientation_and_bounded_side.h>
#include <CGAL/_test_fct_constructions_2.h>
#include <CGAL/_test_fct_constructions_3.h>
#include <CGAL/_test_fct_point_3.C>
#include <CGAL/_test_fct_coplanar_3.h>
#include <CGAL/_test_cls_iso_cuboid_3.C>
#include <CGAL/_test_angle.h>

int
main()
{
  typedef   CGAL::Simple_cartesian<CGAL::Quotient<Precise_integer> >     Cls;
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
  
  return 0;
}
