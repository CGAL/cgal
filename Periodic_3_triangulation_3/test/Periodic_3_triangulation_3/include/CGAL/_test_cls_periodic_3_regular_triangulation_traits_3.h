// Copyright (c) 1998, 2015  INRIA Sophia-Antipolis (France).
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
// Author(s)     :  Aymeric PELLE <aymeric.pelle@sophia.inria.fr>

#if (__GNUC__>4) || (__GNUC__ == 4 && __GNUC_MINOR__ >=6)
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>

#include <CGAL/use.h>

#include <cassert>
#include <iostream>

template<class Traits>
void _test_for_given_domain (const Traits& traits,
                             typename Traits::Weighted_point_3* wp)
{
  typedef typename Traits::Weighted_point_3             Weighted_point_3;
  typedef typename Traits::Point_3                      Point_3;
  typedef typename Traits::Periodic_3_offset_3          Offset;

  CGAL_USE_TYPE(typename Traits::Vector_3);
  CGAL_USE_TYPE(typename Traits::Iso_cuboid_3);
  CGAL_USE_TYPE(typename Traits::Segment_3);
  CGAL_USE_TYPE(typename Traits::Triangle_3);
  CGAL_USE_TYPE(typename Traits::Tetrahedron_3);

  CGAL_USE_TYPE(typename Traits::Comparison_result);
  CGAL_USE_TYPE(typename Traits::Orientation);
  CGAL_USE_TYPE(typename Traits::Oriented_side);
  CGAL_USE_TYPE(typename Traits::Bounded_side);

  typedef typename Traits::Compare_power_distance_3 Compare_power_distance_3;

  CGAL_USE_TYPE(typename Traits::Compare_xyz_3);
  CGAL_USE_TYPE(typename Traits::Orientation_3);
  CGAL_USE_TYPE(typename Traits::Coplanar_orientation_3);

  CGAL_USE_TYPE(typename Traits::Construct_segment_3);
  CGAL_USE_TYPE(typename Traits::Construct_triangle_3);
  CGAL_USE_TYPE(typename Traits::Construct_tetrahedron_3);

  typedef typename Traits::Construct_weighted_circumcenter_3 Construct_weighted_circumcenter_3;

  ////////// ------------ ~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Create offset array for tests
  Offset o[9] = { Offset(0, 0, 0), Offset(0, 1, 1), Offset(1, 0, 1), Offset(1, 1, 0),
                  Offset(1, 1, 1), Offset(1, 0, 0), Offset(-1, 0, 0), Offset(0, 0, 1), Offset(1, 0, -1) };

  Weighted_point_3 wp0 = wp[0];
  Weighted_point_3 wp1 = wp[1];
  Weighted_point_3 wp2 = wp[2];
  Weighted_point_3 wp3 = wp[3];
  Weighted_point_3 wp01 = wp[4];
  Weighted_point_3 wp02 = wp[5];
  Weighted_point_3 wp03 = wp[6];
  Weighted_point_3 wp04 = wp[7];
  Weighted_point_3 wp05 = wp[8];

  Point_3 q0 = wp[9].point();
  Weighted_point_3 wq0 = wp[9];
  Weighted_point_3 wq1 = wp[10];
  Weighted_point_3 wq2 = wp[11];
  Weighted_point_3 wq3 = wp[12];
  Weighted_point_3 wq4 = wp[13];
  Weighted_point_3 wq01 = wp[14];
  Weighted_point_3 wq11 = wp[15];
  Weighted_point_3 wq21 = wp[16];
  Weighted_point_3 wq31 = wp[17];
  Weighted_point_3 wq41 = wp[18];

  // test of Construct_weighted_circumcenter_3 and compare_power_distance
  std::cout << "test of Construct_weighted_circumcenter_3" << std::endl;
  std::cout << "test of Compare_power_distance_3" << std::endl;
  Construct_weighted_circumcenter_3 weighted_circumcenter = traits.construct_weighted_circumcenter_3_object();
  Compare_power_distance_3 compare_power_distance = traits.compare_power_distance_3_object();

  {
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[0],wp[1]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[1]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[1]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[2]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[5],wp[3]) ==  CGAL::EQUAL);
  }
  {
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[0],wp[1], o[0],o[0],o[0]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[1], o[0],o[0],o[0]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[1], o[0],o[0],o[0]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[2], o[0],o[0],o[0]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[5],wp[3], o[0],o[0],o[0]) ==  CGAL::EQUAL);
  }
  {
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[0],wp[1], o[1],o[1],o[1]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[1], o[1],o[1],o[1]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[1], o[1],o[1],o[1]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[2], o[1],o[1],o[1]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[5],wp[3], o[1],o[1],o[1]) ==  CGAL::EQUAL);
  }
  {
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[0],wp[1], o[1],o[1],o[1]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[1], o[1],o[1],o[1]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[1], o[1],o[1],o[1]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[2], o[1],o[1],o[1]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[5],wp[3], o[1],o[1],o[1]) ==  CGAL::EQUAL);
  }
  {
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[1], o[0],o[0],o[1]) == CGAL::SMALLER);
    assert (compare_power_distance(c,wp[6],wp[1], o[0],o[1],o[0]) == CGAL::LARGER);
  }
}

template<class K>
void _test_cls_periodic_3_regular_triangulation_traits_3 ()
{
  typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>    Traits;
  typedef typename Traits::Weighted_point_3                     Weighted_point_3;
  typedef typename Traits::Point_3                              Point_3;
  typedef typename Traits::Iso_cuboid_3                         Iso_cuboid;

  Traits traits;

  // Test Iso_cuboid(0,0,0,1,1,1)
  std::cout << "Iso_cuboid(0,0,0,4,4,4)" << std::endl;
  {
    traits.set_domain(Iso_cuboid(0, 0, 0, 4, 4, 4));
    Point_3 p0(0.,0.,0.);
    Point_3 p1(3.,0.,0.);
    Point_3 p2(0.,3.,0.);
    Point_3 p3(0.,0.,3.);

    Point_3 q0(0.,0.,0.);
    Point_3 q1(2.,0.,0.);
    Point_3 q2(0.,2.,0.);
    Point_3 q3(0.,0.,2.);
    Point_3 q4(2.,2.,2.);

    Weighted_point_3 wp[19] =
    {
      Weighted_point_3(p0,9.),
      Weighted_point_3(p1,9.),
      Weighted_point_3(p2,9.),
      Weighted_point_3(p3,9.),
      Weighted_point_3(p0,6.),
      Weighted_point_3(p0,3.),
      Weighted_point_3(p0,12.),
      Weighted_point_3(p0,18.),
      Weighted_point_3(p0,24.),

      Weighted_point_3(q0,0.),
      Weighted_point_3(q1,0.),
      Weighted_point_3(q2,0.),
      Weighted_point_3(q3,0.),
      Weighted_point_3(q4,0.),
      Weighted_point_3(q0,2.),
      Weighted_point_3(q1,2.),
      Weighted_point_3(q2,2.),
      Weighted_point_3(q3,2.),
      Weighted_point_3(q4,2.)
    };
    _test_for_given_domain(traits, wp);
  }
  std::cout << "Iso_cuboid(-2,-2,-2,2,2,2)" << std::endl;
  {
    traits.set_domain(Iso_cuboid(-2, -2, -2, 2, 2, 2));
    Point_3 p0(-2.,-2.,-2.);
    Point_3 p1(1.,-2.,-2.);
    Point_3 p2(-2.,1.,-2.);
    Point_3 p3(-2.,-2.,1.);

    Point_3 q0(-1.,-1.,-1.);
    Point_3 q1(1.,-1.,-1.);
    Point_3 q2(-1.,1.,-1.);
    Point_3 q3(-1.,-1.,1.);
    Point_3 q4(1.,1.,1.);

    Weighted_point_3 wp[19] =
    {
      Weighted_point_3(p0,9.),
      Weighted_point_3(p1,9.),
      Weighted_point_3(p2,9.),
      Weighted_point_3(p3,9.),
      Weighted_point_3(p0,6.),
      Weighted_point_3(p0,3.),
      Weighted_point_3(p0,12.),
      Weighted_point_3(p0,18.),
      Weighted_point_3(p0,24.),

      Weighted_point_3(q0,0.),
      Weighted_point_3(q1,0.),
      Weighted_point_3(q2,0.),
      Weighted_point_3(q3,0.),
      Weighted_point_3(q4,0.),
      Weighted_point_3(q0,2.),
      Weighted_point_3(q1,2.),
      Weighted_point_3(q2,2.),
      Weighted_point_3(q3,2.),
      Weighted_point_3(q4,2.)
    };
    _test_for_given_domain(traits, wp);
  }
}
