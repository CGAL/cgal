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
// source        : news.fw
// file          : test_new_partsH_.C
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>
#ifdef CGAL_USE_GMP
# include <CGAL/Gmpz.h>
typedef CGAL::Gmpz    Precise_integer;
#else
# ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
typedef leda_integer  Precise_integer;
# endif // CGAL_USE_LEDA
#endif // CGAL_USE_GMP


#include <CGAL/Homogeneous.h>
#include <CGAL/_test_fct_points_implicit_sphere.h>
#include <CGAL/_test_orientation_and_bounded_side.h>
#include <CGAL/_test_fct_constructions_2.h>
#include <CGAL/_test_fct_constructions_3.h>
#include <CGAL/_test_fct_point_3.C>
#include <CGAL/_test_fct_coplanar_3.h>
#ifndef CGAL_STRICT21
#include <CGAL/_test_cls_iso_cuboid_3.C>
#endif // CGAL_STRICT21


int
main()
{
  typedef   CGAL::Homogeneous<Precise_integer>     Cls;
  std::cout << "Testing new parts with Homogeneous<Precise_integer> :";
  std::cout << std::endl;
  _test_orientation_and_bounded_side( Cls() );
  _test_fct_points_implicit_sphere( Cls() );
  _test_fct_constructions_2( Cls() );
  _test_fct_constructions_3( Cls() );
  _test_fct_point_3( Cls() );
  _test_fct_coplanar_3( Cls() );
  #ifndef CGAL_STRICT21
  _test_cls_iso_cuboid_3( Cls() );
  #endif // CGAL_STRICT21
  
  return 0;
}
