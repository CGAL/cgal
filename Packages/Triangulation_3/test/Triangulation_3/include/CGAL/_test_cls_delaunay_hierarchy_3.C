// ============================================================================
//
// Copyright (c) 1998, 2001 The CGAL Consortium
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
// file          : test/Triangulation3/include/CGAL/
//                 _test_cls_delaunay_hierarchy_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec, Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis <Mariette Yvinec@sophia.inria.fr>
// ============================================================================

#include <CGAL/_test_cls_delaunay_3.C>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>

template <class Dh>
void
_test_cls_delaunay_hierarchy_3( const Dh & )
{
  typedef Dh                                     Delaunay_hierarchy;
  typedef typename Delaunay_hierarchy::Point     Point;
  typedef CGAL::Creator_uniform_3<double,Point>  Creator;

  std::cout << "    insertion removal of 1000 points" << std::endl;
  Delaunay_hierarchy dh;
  CGAL::Random_points_in_cube_3<Point,Creator> g(1.);
  CGAL::copy_n( g, 1000, std::back_inserter(dh));

  dh.locate(Point(0.,0.,0.));

  std::cout << " REMOVE IS NOT TESTED !!!!!!!!!!!!!!!!!!!! " << std::endl;
  // while( dh.number_of_vertices() >0)
    // dh.remove(dh.finite_vertices_begin());
}
