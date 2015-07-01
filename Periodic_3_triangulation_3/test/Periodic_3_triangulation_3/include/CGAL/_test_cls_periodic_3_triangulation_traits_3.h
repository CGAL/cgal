// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// $Date$

// Author(s)     :  Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                  Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#if (__GNUC__>4) || (__GNUC__ == 4 && __GNUC_MINOR__ >=6)
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include <CGAL/use.h>

template <class Traits>
void _test_for_given_domain(const Traits & traits,
    typename Traits::Point_3 * p)
{
  typedef typename Traits::Point_3                   Point;
  typedef typename Traits::Vector_3                  Vector;
  CGAL_USE_TYPE(Vector);
  typedef typename Traits::Periodic_3_offset_3       Offset;
  CGAL_USE_TYPE(typename Traits::Iso_cuboid_3);
  typedef typename Traits::Segment_3                 Segment;
  typedef typename Traits::Triangle_3                Triangle;
  typedef typename Traits::Tetrahedron_3             Tetrahedron;

  CGAL_USE_TYPE(typename Traits::Comparison_result);
  CGAL_USE_TYPE(typename Traits::Orientation);
  CGAL_USE_TYPE(typename Traits::Oriented_side);
  CGAL_USE_TYPE(typename Traits::Bounded_side);

  typedef typename Traits::Compare_xyz_3             Compare_xyz_3;
  typedef typename Traits::Orientation_3             Orientation_3;
  typedef typename Traits::Side_of_oriented_sphere_3 Side_of_oriented_sphere_3;
  typedef typename Traits::Compare_distance_3        Compare_distance_3;
  typedef typename Traits::Coplanar_orientation_3    Coplanar_orientation_3;
  typedef typename Traits::Coplanar_side_of_bounded_circle_3
    Coplanar_side_of_bounded_circle_3;
  typedef typename Traits::Side_of_bounded_sphere_3  Side_of_bounded_sphere_3;

  typedef typename Traits::Construct_point_3         Construct_point_3;
  typedef typename Traits::Construct_segment_3       Construct_segment_3;
  typedef typename Traits::Construct_triangle_3      Construct_triangle_3;
  typedef typename Traits::Construct_tetrahedron_3   Construct_tetrahedron_3;
  typedef typename Traits::Construct_circumcenter_3  Construct_circumcenter_3;

  // Operations
  Compare_xyz_3 compare_xyz = traits.compare_xyz_3_object();
  Orientation_3 orientation = traits.orientation_3_object();
  Side_of_oriented_sphere_3 soos = traits.side_of_oriented_sphere_3_object();
  Compare_distance_3 compare_distance = traits.compare_distance_3_object();
  Coplanar_orientation_3 coplanar_orientation
    = traits.coplanar_orientation_3_object();
  Coplanar_side_of_bounded_circle_3 csobc
    = traits.coplanar_side_of_bounded_circle_3_object();
  Side_of_bounded_sphere_3 sobs = traits.side_of_bounded_sphere_3_object();

  Construct_point_3 construct_point = traits.construct_point_3_object();
  Construct_segment_3 construct_segment = traits.construct_segment_3_object();
  Construct_triangle_3 construct_triangle
    = traits.construct_triangle_3_object();
  Construct_tetrahedron_3 construct_tetrahedron
    = traits.construct_tetrahedron_3_object();
  Construct_circumcenter_3 construct_circumcenter
    = traits.construct_circumcenter_3_object();

  // Create offset array for tests
  Offset o[9]={Offset( 0, 0, 0), Offset( 0, 1, 1), Offset( 1, 0, 1),
	       Offset( 1, 1, 0), Offset( 1, 1, 1), Offset( 1, 0, 0),
	       Offset(-1, 0, 0), Offset( 0, 0, 1), Offset( 1, 0,-1)};

  // Test Compare_xyz_3
  assert(compare_xyz(p[0],p[1]) == CGAL::SMALLER);
  assert(compare_xyz(p[1],p[0]) == CGAL::LARGER);
  assert(compare_xyz(p[1],p[7]) == CGAL::SMALLER);
  assert(compare_xyz(p[7],p[1]) == CGAL::LARGER);
  assert(compare_xyz(p[7],p[7]) == CGAL::EQUAL);

  assert(compare_xyz(p[0],p[1],o[0],o[0]) == CGAL::SMALLER);
  assert(compare_xyz(p[1],p[0],o[0],o[0]) == CGAL::LARGER);
  assert(compare_xyz(p[1],p[7],o[0],o[0]) == CGAL::SMALLER);
  assert(compare_xyz(p[7],p[1],o[0],o[0]) == CGAL::LARGER);
  assert(compare_xyz(p[7],p[7],o[0],o[0]) == CGAL::EQUAL);

  assert(compare_xyz(p[0],p[1],o[5],o[0]) == CGAL::LARGER);
  assert(compare_xyz(p[1],p[0],o[6],o[0]) == CGAL::SMALLER);
  assert(compare_xyz(p[1],p[7],o[0],o[6]) == CGAL::LARGER);
  assert(compare_xyz(p[7],p[1],o[0],o[5]) == CGAL::SMALLER);
  assert(compare_xyz(p[7],p[7],o[5],o[5]) == CGAL::EQUAL);
  assert(compare_xyz(p[7],p[7],o[5],o[6]) == CGAL::LARGER);
  assert(compare_xyz(p[7],p[7],o[6],o[5]) == CGAL::SMALLER);
  assert(compare_xyz(p[7],p[7],o[6],o[6]) == CGAL::EQUAL);


  // Test Orientation_3
  assert(orientation(p[1],p[2],p[3],p[0]) == CGAL::NEGATIVE);
  assert(orientation(p[1],p[2],p[3],p[4]) == CGAL::POSITIVE);
  assert(orientation(p[1],p[3],p[2],p[0]) == CGAL::POSITIVE);
  assert(orientation(p[1],p[2],p[3],p[10]) == CGAL::COPLANAR);

  assert(orientation(p[1],p[2],p[3],p[0],o[0],o[0],o[0],o[0])
      == CGAL::NEGATIVE);
  assert(orientation(p[1],p[2],p[3],p[0],o[0],o[0],o[0],o[4])
      == CGAL::POSITIVE);
  assert(orientation(p[1],p[2],p[3],p[10],o[0],o[0],o[0],o[8])
      == CGAL::COPLANAR);
  assert(orientation(p[1],p[2],p[3],p[10],o[4],o[4],o[4],o[4]+o[8])
      == CGAL::COPLANAR);
  assert(orientation(p[1],p[2],p[3],p[10],o[4],o[4],o[4],o[0])
      == CGAL::NEGATIVE);

  // Test Side_of_oriented_sphere_3
  assert(soos(p[0],p[1],p[2],p[3],p[4]) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(soos(p[0],p[1],p[2],p[3],p[5]) == CGAL::ON_POSITIVE_SIDE);
  assert(soos(p[0],p[1],p[2],p[3],p[6]) == CGAL::ON_NEGATIVE_SIDE);
  assert(soos(p[0],p[2],p[1],p[3],p[5]) == CGAL::ON_NEGATIVE_SIDE);
  assert(soos(p[0],p[2],p[1],p[3],p[6]) == CGAL::ON_POSITIVE_SIDE);

  assert(soos(p[0],p[1],p[2],p[3],p[4],o[4],o[1],o[2],o[3],o[0]) 
      == CGAL::ON_ORIENTED_BOUNDARY);
  assert(soos(p[0],p[1],p[2],p[3],p[5],o[4],o[1],o[2],o[3],o[0])
      == CGAL::ON_POSITIVE_SIDE);
  assert(soos(p[0],p[1],p[2],p[3],p[6],o[4],o[1],o[2],o[3],o[0])
      == CGAL::ON_NEGATIVE_SIDE);
  assert(soos(p[0],p[2],p[1],p[3],p[5],o[4],o[2],o[1],o[3],o[0])
      == CGAL::ON_NEGATIVE_SIDE);
  assert(soos(p[0],p[2],p[1],p[3],p[6],o[4],o[2],o[1],o[3],o[0])
      == CGAL::ON_POSITIVE_SIDE);
  assert(soos(p[0],p[1],p[2],p[3],p[5],o[4],o[1],o[2],o[3],o[4])
      == CGAL::ON_POSITIVE_SIDE);
  assert(soos(p[0],p[1],p[2],p[3],p[6],o[4],o[1],o[2],o[3],o[4])
      == CGAL::ON_POSITIVE_SIDE);
  assert(soos(p[0],p[2],p[1],p[3],p[5],o[4],o[2],o[1],o[3],o[4])
      == CGAL::ON_NEGATIVE_SIDE);
  assert(soos(p[0],p[2],p[1],p[3],p[6],o[4],o[2],o[1],o[3],o[4])
      == CGAL::ON_NEGATIVE_SIDE);

  // Test Compare_distance_3
  assert(compare_distance(p[0],p[1],p[2]) == CGAL::EQUAL);
  assert(compare_distance(p[0],p[2],p[1]) == CGAL::EQUAL);
  assert(compare_distance(p[0],p[1],p[4]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[4],p[1]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[2],p[4]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[4],p[2]) == CGAL::LARGER);
  
  assert(compare_distance(p[0],p[1],p[2],o[0],o[0],o[0]) == CGAL::EQUAL);
  assert(compare_distance(p[0],p[2],p[1],o[0],o[0],o[0]) == CGAL::EQUAL);
  assert(compare_distance(p[0],p[1],p[4],o[0],o[0],o[0]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[4],p[1],o[0],o[0],o[0]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[2],p[4],o[0],o[0],o[0]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[4],p[2],o[0],o[0],o[0]) == CGAL::LARGER);

  assert(compare_distance(p[0],p[1],p[2],o[4],o[0],o[0]) == CGAL::EQUAL);
  assert(compare_distance(p[0],p[2],p[1],o[4],o[0],o[0]) == CGAL::EQUAL);
  assert(compare_distance(p[0],p[1],p[4],o[4],o[0],o[0]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[4],p[1],o[4],o[0],o[0]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[2],p[4],o[4],o[0],o[0]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[4],p[2],o[4],o[0],o[0]) == CGAL::SMALLER);

  assert(compare_distance(p[0],p[1],p[2],o[0],o[4],o[0]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[2],p[1],o[0],o[0],o[4]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[1],p[4],o[0],o[4],o[0]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[4],p[1],o[0],o[0],o[4]) == CGAL::SMALLER);
  assert(compare_distance(p[0],p[2],p[4],o[0],o[4],o[0]) == CGAL::LARGER);
  assert(compare_distance(p[0],p[4],p[2],o[0],o[0],o[4]) == CGAL::SMALLER);

  // Test Coplanar_orientation_3
  assert(coplanar_orientation(p[1],p[3],p[0]) == CGAL::POSITIVE);
  assert(coplanar_orientation(p[1],p[3],p[7]) == CGAL::NEGATIVE);
  assert(coplanar_orientation(p[3],p[1],p[0]) == CGAL::NEGATIVE);
  assert(coplanar_orientation(p[1],p[3],p[8]) == CGAL::COLLINEAR);

  assert(coplanar_orientation(p[1],p[3],p[0],o[0],o[0],o[0]) == CGAL::POSITIVE);
  assert(coplanar_orientation(p[1],p[3],p[0],o[0],o[0],o[2]) == CGAL::NEGATIVE);
  assert(coplanar_orientation(p[1],p[3],p[8],o[0],o[0],o[8])
      == CGAL::COLLINEAR);
  assert(coplanar_orientation(p[1],p[3],p[8],o[2],o[2],o[2]+o[8])
      == CGAL::COLLINEAR);
  assert(coplanar_orientation(p[1],p[3],p[8],o[2],o[2],o[0]) == CGAL::POSITIVE);

  // Test Coplanar_side_of_bounded_circle_3
  assert(csobc(p[0],p[1],p[3],p[7]) == CGAL::ON_BOUNDARY);
  assert(csobc(p[0],p[1],p[3],p[8]) == CGAL::ON_BOUNDED_SIDE);
  assert(csobc(p[0],p[1],p[3],p[9]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(csobc(p[0],p[3],p[1],p[8]) == CGAL::ON_BOUNDED_SIDE);
  assert(csobc(p[0],p[3],p[1],p[9]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(csobc(p[0],p[1],p[3],p[7],o[2],o[7],o[5],o[0]) == CGAL::ON_BOUNDARY);
  assert(csobc(p[0],p[1],p[3],p[8],o[2],o[7],o[5],o[0])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(csobc(p[0],p[1],p[3],p[9],o[2],o[7],o[5],o[0])
      == CGAL::ON_BOUNDED_SIDE);
  assert(csobc(p[0],p[3],p[1],p[8],o[2],o[5],o[7],o[0])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(csobc(p[0],p[3],p[1],p[9],o[2],o[5],o[7],o[0])
      == CGAL::ON_BOUNDED_SIDE);
  assert(csobc(p[0],p[1],p[3],p[8],o[2],o[7],o[5],o[2])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(csobc(p[0],p[1],p[3],p[9],o[2],o[7],o[5],o[2])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(csobc(p[0],p[3],p[1],p[8],o[2],o[5],o[7],o[2])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(csobc(p[0],p[3],p[1],p[9],o[2],o[5],o[7],o[2])
      == CGAL::ON_UNBOUNDED_SIDE);

  // Test Side_of_bounded_sphere_3
  assert(sobs(p[0],p[1],p[2],p[3],p[4]) == CGAL::ON_BOUNDARY);
  assert(sobs(p[0],p[1],p[2],p[3],p[5]) == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[2],p[3],p[6]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[1],p[3],p[5]) == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[1],p[3],p[6]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(sobs(p[0],p[1],p[2],p[3],p[4],o[4],o[1],o[2],o[3],o[0])
      == CGAL::ON_BOUNDARY);
  assert(sobs(p[0],p[1],p[2],p[3],p[5],o[4],o[1],o[2],o[3],o[0])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[2],p[3],p[6],o[4],o[1],o[2],o[3],o[0])
      == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[1],p[3],p[5],o[4],o[2],o[1],o[3],o[0])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[1],p[3],p[6],o[4],o[2],o[1],o[3],o[0])
      == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[2],p[3],p[5],o[4],o[1],o[2],o[3],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[2],p[3],p[6],o[4],o[1],o[2],o[3],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[1],p[3],p[5],o[4],o[2],o[1],o[3],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[1],p[3],p[6],o[4],o[2],o[1],o[3],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);

  assert(sobs(p[0],p[2],p[4],p[3]) == CGAL::ON_BOUNDARY);
  assert(sobs(p[0],p[2],p[4],p[5]) == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[4],p[6]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[4],p[5]) == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[4],p[6]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(sobs(p[0],p[2],p[4],p[3],o[4],o[2],o[0],o[3]) == CGAL::ON_BOUNDARY);
  assert(sobs(p[0],p[2],p[4],p[5],o[4],o[2],o[0],o[0])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[4],p[6],o[4],o[2],o[0],o[0])
      == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[4],p[5],o[4],o[1],o[0],o[0])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[4],p[6],o[4],o[1],o[0],o[0]) 
      == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[4],p[5],o[4],o[2],o[0],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[2],p[4],p[6],o[4],o[2],o[0],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[4],p[5],o[4],o[1],o[0],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[1],p[4],p[6],o[4],o[1],o[0],o[4])
      == CGAL::ON_UNBOUNDED_SIDE);

  assert(sobs(p[0],p[4],p[3]) == CGAL::ON_BOUNDARY);
  assert(sobs(p[0],p[4],p[5]) == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[4],p[6]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(sobs(p[0],p[4],p[3],o[4],o[0],o[3]) == CGAL::ON_BOUNDARY);
  assert(sobs(p[0],p[4],p[5],o[4],o[0],o[0]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[4],p[6],o[4],o[0],o[0]) == CGAL::ON_BOUNDED_SIDE);
  assert(sobs(p[0],p[4],p[5],o[4],o[0],o[4]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(sobs(p[0],p[4],p[6],o[4],o[0],o[4]) == CGAL::ON_UNBOUNDED_SIDE);

  // Test Construct_point_3
  Point p0 = construct_point(p[0],o[0]);
  p0 = construct_point(p[1],o[1]);
  p0 = construct_point(p[2],o[2]);
  p0 = construct_point(p[3],o[3]);
  CGAL_USE(p0);
  // Test Construct_segment_3
  Segment s01 = construct_segment(p[0],p[1]);
  s01 = construct_segment(p[0],p[2]);
  s01 = construct_segment(p[0],p[3]);
  s01 = construct_segment(p[1],p[2]);
  s01 = construct_segment(p[1],p[3]);
  s01 = construct_segment(p[2],p[3]);

  s01 = construct_segment(p[0],p[1],o[0],o[1]);
  s01 = construct_segment(p[0],p[2],o[0],o[2]);
  s01 = construct_segment(p[0],p[3],o[0],o[3]);
  s01 = construct_segment(p[1],p[2],o[1],o[2]);
  s01 = construct_segment(p[1],p[3],o[1],o[3]);
  s01 = construct_segment(p[2],p[3],o[2],o[3]);
  CGAL_USE(s01);

  // Test Construct_triangle_3
  Triangle t0 = construct_triangle(p[1],p[2],p[3]);
  t0 = construct_triangle(p[0],p[2],p[3]);
  t0 = construct_triangle(p[0],p[1],p[3]);
  t0 = construct_triangle(p[0],p[1],p[2]);

  t0 = construct_triangle(p[1],p[2],p[3],o[1],o[2],o[3]);
  t0 = construct_triangle(p[0],p[2],p[3],o[0],o[2],o[3]);
  t0 = construct_triangle(p[0],p[1],p[3],o[0],o[1],o[3]);
  t0 = construct_triangle(p[0],p[1],p[2],o[0],o[1],o[2]);
  CGAL_USE(t0);

  // Test Construct_triangle_3
  Tetrahedron t = construct_tetrahedron(p[0],p[1],p[2],p[3]);
  t = construct_tetrahedron(p[0],p[1],p[2],p[3],
      o[0],o[1],o[2],o[3]);
  CGAL_USE(t);

  // Test of Construct_circumcenter_3
  Point c = construct_circumcenter(p[0],p[1],p[2],p[3]);
  assert(compare_distance(c,p[0],p[1]) == CGAL::EQUAL);
  assert(compare_distance(c,p[1],p[2]) == CGAL::EQUAL);
  assert(compare_distance(c,p[2],p[3]) == CGAL::EQUAL);
  assert(compare_distance(c,p[3],p[4]) == CGAL::EQUAL);
  assert(compare_distance(c,p[4],p[0]) == CGAL::EQUAL);
  assert(compare_distance(c,p[0],p[5]) == CGAL::LARGER);
  assert(compare_distance(c,p[0],p[6]) == CGAL::SMALLER);

  assert(c == construct_circumcenter(p[0],p[1],p[2],p[3]));
  assert(c == construct_circumcenter(p[0],p[1],p[2],p[4]));
  assert(c == construct_circumcenter(p[0],p[1],p[3],p[4]));
  assert(c == construct_circumcenter(p[0],p[2],p[3],p[4]));
  assert(c == construct_circumcenter(p[1],p[2],p[3],p[4]));
  assert(c == construct_circumcenter(p[1],p[2],p[4],p[3]));
  assert(c == construct_circumcenter(p[1],p[3],p[2],p[4]));
  assert(c == construct_circumcenter(p[1],p[3],p[4],p[2]));
  assert(c == construct_circumcenter(p[1],p[4],p[2],p[3]));
  assert(c == construct_circumcenter(p[1],p[4],p[3],p[2]));
	  
  c = construct_circumcenter(p[0],p[1],p[2],p[3],o[4],o[1],o[2],o[3]);
  assert(compare_distance(c,p[0],p[1],o[0],o[4],o[1]) == CGAL::EQUAL);
  assert(compare_distance(c,p[1],p[2],o[0],o[1],o[2]) == CGAL::EQUAL);
  assert(compare_distance(c,p[2],p[3],o[0],o[2],o[3]) == CGAL::EQUAL);
  assert(compare_distance(c,p[3],p[4],o[0],o[3],o[0]) == CGAL::EQUAL);
  assert(compare_distance(c,p[4],p[0],o[0],o[0],o[4]) == CGAL::EQUAL);
  assert(compare_distance(c,p[0],p[5],o[0],o[4],o[0]) == CGAL::SMALLER);
  assert(compare_distance(c,p[0],p[6],o[0],o[4],o[0]) == CGAL::LARGER);
  assert(compare_distance(c,p[0],p[5],o[0],o[4],o[4]) == CGAL::SMALLER);
  assert(compare_distance(c,p[0],p[6],o[0],o[4],o[4]) == CGAL::SMALLER);

  assert(c == construct_circumcenter(p[0],p[1],p[2],p[3],o[4],o[1],o[2],o[3]));
  assert(c == construct_circumcenter(p[0],p[1],p[2],p[4],o[4],o[1],o[2],o[0]));
  assert(c == construct_circumcenter(p[0],p[1],p[3],p[4],o[4],o[1],o[3],o[0]));
  assert(c == construct_circumcenter(p[0],p[2],p[3],p[4],o[4],o[2],o[3],o[0]));
  assert(c == construct_circumcenter(p[1],p[2],p[3],p[4],o[1],o[2],o[3],o[0]));
  assert(c == construct_circumcenter(p[1],p[2],p[4],p[3],o[1],o[2],o[0],o[3]));
  assert(c == construct_circumcenter(p[1],p[3],p[2],p[4],o[1],o[3],o[2],o[0]));
  assert(c == construct_circumcenter(p[1],p[3],p[4],p[2],o[1],o[3],o[0],o[2]));
  assert(c == construct_circumcenter(p[1],p[4],p[2],p[3],o[1],o[0],o[2],o[3]));
  assert(c == construct_circumcenter(p[1],p[4],p[3],p[2],o[1],o[0],o[3],o[2]));
}

template <class K>
void _test_cls_periodic_3_triangulation_traits_3() {
  typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K> Traits;
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::Iso_cuboid_3 Iso_cuboid;

  Traits traits;

  // Test Iso_cuboid(0,0,0,1,1,1)
  traits.set_domain(Iso_cuboid(0,0,0,1,1,1));

  Point p[11]={Point(0,0,0,8),
	       Point(4,0,0,8),
	       Point(0,4,0,8),
	       Point(0,0,4,8),
	       Point(4,4,4,8),
	       Point(2,4,2,8),
	       Point(4,6,6,8),
	       Point(4,0,4,8),
	       Point(2,0,2,8),
	       Point(4,0,6,8),
	       Point(1,1,2,8)};

  _test_for_given_domain(traits,p);

  // Test Iso_cuboid(0,0,0,8,8,8)
  traits.set_domain(Iso_cuboid(0,0,0,8,8,8));

  p[0] = Point(0,0,0);
  p[1] = Point(4,0,0);
  p[2] = Point(0,4,0);
  p[3] = Point(0,0,4);
  p[4] = Point(4,4,4);
  p[5] = Point(2,4,2);
  p[6] = Point(4,6,6);
  p[7] = Point(4,0,4);
  p[8] = Point(2,0,2);
  p[9] = Point(4,0,6);
  p[10]= Point(1,1,2);

  _test_for_given_domain(traits,p);

  // Test Iso_cuboid(-0.5,-0.5,-0.5,0.5,0.5,0.5)
  traits.set_domain(Iso_cuboid(-1,-1,-1,1,1,1,2));

  p[0] = Point(-4,-4,-4,8);
  p[1] = Point( 0,-4,-4,8);
  p[2] = Point(-4, 0,-4,8);
  p[3] = Point(-4,-4, 0,8);
  p[4] = Point( 0, 0, 0,8);
  p[5] = Point(-2, 0,-2,8);
  p[6] = Point( 0, 2, 2,8);
  p[7] = Point( 0,-4, 0,8);
  p[8] = Point(-2,-4,-2,8);
  p[9] = Point( 0,-4, 2,8);
  p[10]= Point(-3,-3,-2,8);
  
  _test_for_given_domain(traits,p);

  // Test Iso_cuboid(-4,-4,-4,4,4,4)
  traits.set_domain(Iso_cuboid(-4,-4,-4,4,4,4));
  
  p[0] = Point(-4,-4,-4);
  p[1] = Point( 0,-4,-4);
  p[2] = Point(-4, 0,-4);
  p[3] = Point(-4,-4, 0);
  p[4] = Point( 0, 0, 0);
  p[5] = Point(-2, 0,-2);
  p[6] = Point( 0, 2, 2);
  p[7] = Point( 0,-4, 0);
  p[8] = Point(-2,-4,-2);
  p[9] = Point( 0,-4, 2);
  p[10]= Point(-3,-3,-2);
  
  _test_for_given_domain(traits,p);
}
