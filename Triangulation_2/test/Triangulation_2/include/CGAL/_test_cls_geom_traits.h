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
// file          : include/CGAL/_test_cls_geom_traits.h
// revision      :
// revision_date :
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>


template <class Traits, class Point>
void
_test_cls_geom_traits(Point p[34], const Traits & )
{

  typedef typename Traits::Segment_2           Segment;
  typedef typename Traits::Triangle_2          Triangle;

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
  typename Traits::Compare_x_2 compare_x_2 =  gt.compare_x_2_object();
  assert( compare_x_2(p[0],p[1]) == CGAL::SMALLER );
  assert( compare_x_2(p[1],p[0]) == CGAL::LARGER );
  assert( compare_x_2(p[1],p[2]) == CGAL::EQUAL );

  // Test compare_y()
  typename Traits::Compare_y_2 compare_y_2 =  gt.compare_y_2_object();
  assert( compare_y_2(p[3],p[4]) == CGAL::SMALLER );
  assert( compare_y_2(p[4],p[3]) == CGAL::LARGER );
  assert( compare_y_2(p[4],p[5]) == CGAL::EQUAL );

  // Test orientation()
  typename Traits::Orientation_2 orientation = gt.orientation_2_object();
  assert( orientation(p[9], p[10],p[11]) == CGAL::COUNTERCLOCKWISE );
  assert( orientation(p[12],p[13],p[14]) == CGAL::CLOCKWISE );
  assert( orientation(p[15],p[16],p[17]) == CGAL::COLLINEAR );

  // Test side_of_oriented_circle()
  typename Traits::Side_of_oriented_circle_2
    in_circle = gt.side_of_oriented_circle_2_object();
  assert( in_circle(p[18],p[19],p[20],p[21]) == CGAL::ON_NEGATIVE_SIDE );
  assert( in_circle(p[22],p[23],p[24],p[25]) == CGAL::ON_POSITIVE_SIDE );
  assert( in_circle(p[26],p[27],p[28],p[29]) == CGAL::ON_ORIENTED_BOUNDARY );

}


template <class Traits, class Point>
void
_test_cls_delaunay_geom_traits(Point p[34], const Traits & )
{
  // All the other requirements
  _test_cls_geom_traits(p, Traits() );

  // Constructor
  Traits gt;

  // Distance
  typedef typename Traits::Compare_distance_2 Compare_distance_2;
  Compare_distance_2 closer = gt.compare_distance_2_object();
  assert( closer(p[0],p[1],p[2]) == CGAL::SMALLER);
  assert( closer(p[0],p[2],p[1]) == CGAL::LARGER);
  assert( closer(p[0],p[1],p[1]) == CGAL::EQUAL);
  assert( closer(p[0],p[0],p[0]) == CGAL::EQUAL);


  // Test circumcenter()
  typename Traits::Construct_circumcenter_2
    circumcenter = gt.construct_circumcenter_2_object();
  Point c = circumcenter(p[30],p[31],p[32]);
  assert( gt.compare_x_2_object()(p[33],c) ==         CGAL::EQUAL &&
          gt.compare_y_2_object()(p[33],c) ==  CGAL::EQUAL );

  //Test bisector()
}
