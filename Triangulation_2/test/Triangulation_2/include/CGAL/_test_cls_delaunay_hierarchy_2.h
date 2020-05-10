// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// source        : $URL$
// file          : include/CGAL/_test_cls_delaunay_hierarchy_2.h
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : INRIA Sophia-Antipolis <Mariette Yvinec@sophia.inria.fr>
// ============================================================================

#include <CGAL/_test_cls_delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

template <class Dh>
void
_test_cls_delaunay_hierarchy_2( const Dh & )
{
  typedef Dh  Delaunay_hierarchy;
  typedef typename Delaunay_hierarchy::Point   Point;
  typedef CGAL::Creator_uniform_2<double,Point>  Creator;


  // makes i686_CYGWINNT-5.0-1.1.4-0.26-3-2_CL.EXE-12.00.8804 crash
  _test_cls_delaunay_triangulation_2( Delaunay_hierarchy() );


  std::cout << "    insertion removal of 1000 points" << std::endl;
  Delaunay_hierarchy dh;
  CGAL::Random_points_in_square_2<Point,Creator> g(1.);
  std::copy_n( g, 1000, std::back_inserter(dh));

  dh.locate(Point(0.,0.));

  while( dh.number_of_vertices() >0) {
    dh.remove(dh.finite_vertices_begin());
  }

  return;
}
