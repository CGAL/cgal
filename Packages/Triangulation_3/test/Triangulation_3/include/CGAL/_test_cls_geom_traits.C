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
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================
#include <cassert>

template <class Traits, class Point >


void
_test_cls_geom_traits( Point p[55], const Traits & )

{
  typedef typename Traits::Point             Point;
  typedef typename Traits::Segment           Segment;
  typedef typename Traits::Triangle          Triangle;
  typedef typename Traits::Tetrahedron       Tetrahedron;

  // test Point
  Point p0(p[0]);
  Point p1 = p[1];
  Point p2 = p[2];
  Point p3 = p[3];

  // Test Segment
  Segment s1(p0,p1);
  assert( &s1 == &s1 ); // to avoid a warning

  // Test Triangle
  Triangle t1(p0,p1,p2);
  assert( &t1 == &t1 ); // to avoid a warning
  
  // Test Tetrahedron
  Tetrahedron th(p0,p1,p2,p3);
  assert(&th==&th);

  // Test constructor
  Traits gt;

  // Test compare_x()
  assert( gt.compare_x(p[0],p[1]) == CGAL::SMALLER );
  assert( gt.compare_x(p[1],p[0]) == CGAL::LARGER );
  assert( gt.compare_x(p[1],p[2]) == CGAL::EQUAL );

  // Test compare_y()
  assert( gt.compare_y(p[3],p[4]) == CGAL::SMALLER );
  assert( gt.compare_y(p[4],p[3]) == CGAL::LARGER );
  assert( gt.compare_y(p[4],p[5]) == CGAL::EQUAL );

  // Test compare_z()
  assert( gt.compare_z(p[4],p[3]) == CGAL::SMALLER );
  assert( gt.compare_z(p[3],p[4]) == CGAL::LARGER );
  assert( gt.compare_z(p[4],p[5]) == CGAL::EQUAL );

  // Test compare()
  assert( gt.equal(p[6],p[7]) );
  assert( !gt.equal(p[7],p[8]) );


  // Test orientation()

  assert( gt.orientation(p[9], p[10],p[11],p[12]) == CGAL::POSITIVE );
  assert( gt.orientation(p[13], p[14],p[15],p[16]) == CGAL::POSITIVE );
  assert( gt.orientation(p[17], p[18],p[19],p[20]) == CGAL::POSITIVE );
  assert( gt.orientation(p[21],p[22],p[23],p[24]) == CGAL::NEGATIVE);
  assert( gt.orientation(p[25],p[26],p[27],p[28]) == CGAL::NEGATIVE);
  assert( gt.orientation(p[29],p[30],p[31],p[32]) == CGAL::NEGATIVE);
  assert( gt.orientation(p[33],p[34],p[35],p[36]) == CGAL::COLLINEAR );
  assert( gt.orientation(p[37],p[38],p[39],p[40]) == CGAL::COLLINEAR );
  assert( gt.orientation(p[41],p[42],p[43],p[44]) == CGAL::COLLINEAR );
  assert( gt.orientation(p[45],p[46],p[47],p[48]) == CGAL::COLLINEAR );

  // Test orientation in plane()

  assert(gt.orientation_in_plane(p[49],p[50],p[51],p[52])==
	 CGAL::POSITIVE);
  assert(gt.orientation_in_plane(p[51],p[50],p[53],p[49])==
	 CGAL::NEGATIVE);
  assert(gt.orientation_in_plane(p[54],p[50],p[49],p[52])== CGAL::COLLINEAR);
}
