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
// source        : 
// file          : include/CGAL/_test_cls_geom_traits.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <assert.h>
#include <CGAL/_test_cls_distance.C>

// In order to test the class properly, we need some points whose
// geometry is defined accordingly. The following must be true:
//    compare_x  p[0] < p[1] == p[2]
//    compare_y  p[3] < p[4] == p[5]
//    compare    p[6] == p[7] != p[8]
//    orientation:
//               p[9],p[10],p[11] CCW
//               p[12],p[13],p[14] CW
//               p[15],p[16],p[17] COLLINEAR

template <class Traits, class Point >

void
CGAL__test_cls_geom_traits( Point p[34], const Traits & )
{
  typedef typename Traits::Point             Point;
  typedef typename Traits::Segment           Segment;
  typedef typename Traits::Triangle          Triangle;

  // Test Point
  Point p0(p[0]);
  Point p1 = p[1];
  Point p2 = p[2];
  
  // Test Segment
  Segment s1(p0,p1);
  assert( &s1 == &s1 ); // to avoid a warning

  // Test Triangle
  Triangle t1(p0,p1,p2);
  assert( &t1 == &t1 ); // to avoid a warning
  
  // Test constructor
  Traits gt;

  // Test compare_x()
  assert( gt.compare_x(p[0],p[1]) == CGAL_SMALLER );
  assert( gt.compare_x(p[1],p[0]) == CGAL_LARGER );
  assert( gt.compare_x(p[1],p[2]) == CGAL_EQUAL );

  // Test compare_y()
  assert( gt.compare_y(p[3],p[4]) == CGAL_SMALLER );
  assert( gt.compare_y(p[4],p[3]) == CGAL_LARGER );
  assert( gt.compare_y(p[4],p[5]) == CGAL_EQUAL );
  
  // Test compare()
  assert( gt.compare(p[6],p[7]) );
  assert( !gt.compare(p[7],p[8]) );

  // Test orientation()
  assert( gt.orientation(p[9], p[10],p[11]) == CGAL_COUNTERCLOCKWISE );
  assert( gt.orientation(p[12],p[13],p[14]) == CGAL_CLOCKWISE );
  assert( gt.orientation(p[15],p[16],p[17]) == CGAL_COLLINEAR );
}

// For Delaunay geometric traits:
//    side_of_oriented_circle
//               p[18], p[19], p[20],  p[21] POSITIVE
//               p[22], p[23], p[24],  p[25] NEGATIVE
//               p[26], p[27], p[28],  p[29] ON_BDRY
//    circumcenter
//               p[30], p[31], p[32] --> p[33]

template <class Traits, class Point >
void
CGAL__test_cls_delaunay_geom_traits( Point p[34], const Traits & )
{
  // All the other requirements
  CGAL__test_cls_geom_traits(p, Traits() );

  // Constructor
  Traits gt;

  // Distance
  CGAL__test_cls_distance(p, Traits());

  // Test side_of_oriented_circle()
  assert( gt.side_of_oriented_circle(p[18],p[19],p[20],p[21]) == CGAL_ON_NEGATIVE_SIDE );
  assert( gt.side_of_oriented_circle(p[22],p[23],p[24],p[25]) == CGAL_ON_POSITIVE_SIDE );
  assert( gt.side_of_oriented_circle(p[26],p[27],p[28],p[29]) == CGAL_ON_ORIENTED_BOUNDARY );

  // Test circumcenter()
  assert( gt.compare(p[33], gt.circumcenter(p[30],p[31],p[32])) );
}
