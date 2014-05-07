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
// Author(s)     :  Aymeric PELLE <aymeric.pelle@sophia.inria.fr>

#if (__GNUC__>4) || (__GNUC__ == 4 && __GNUC_MINOR__ >=6)
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include <CGAL/use.h>

#include <cassert>

template<class Traits>
void _test_for_given_domain (const Traits & traits, typename Traits::Weighted_point* wp)
{
  typedef typename Traits::Weighted_point Weighted_point;
  typedef typename Traits::Bare_point Bare_point;
  typedef typename Traits::Vector_3 Vector;
  typedef typename Traits::Periodic_3_offset_3 Offset;
  CGAL_USE_TYPE(typename Traits::Iso_cuboid_3);
  typedef typename Traits::Segment_3 Segment;
  typedef typename Traits::Triangle_3 Triangle;
  typedef typename Traits::Tetrahedron_3 Tetrahedron;

  CGAL_USE_TYPE(typename Traits::Comparison_result);
  CGAL_USE_TYPE(typename Traits::Orientation);
  CGAL_USE_TYPE(typename Traits::Oriented_side);
  CGAL_USE_TYPE(typename Traits::Bounded_side);

  typedef typename Traits::Compare_power_distance_3 Compare_power_distance_3;

  typedef typename Traits::Compare_xyz_3 Compare_xyz_3;
  typedef typename Traits::Orientation_3 Orientation_3;
  typedef typename Traits::Coplanar_orientation_3 Coplanar_orientation_3;

  typedef typename Traits::Construct_segment_3 Construct_segment_3;
  typedef typename Traits::Construct_triangle_3 Construct_triangle_3;
  typedef typename Traits::Construct_tetrahedron_3 Construct_tetrahedron_3;

  typedef typename Traits::Construct_weighted_circumcenter_3 Construct_weighted_circumcenter_3;
  typedef typename Traits::Does_simplex_intersect_dual_support_3 Does_simplex_intersect_dual_support_3;

  // Create offset array for tests
  Offset o[9] = { Offset(0, 0, 0), Offset(0, 1, 1), Offset(1, 0, 1), Offset(1, 1, 0),
                  Offset(1, 1, 1), Offset(1, 0, 0), Offset(-1, 0, 0), Offset(0, 0, 1), Offset(1, 0, -1) };

  Weighted_point wp0 = wp[0];
  Weighted_point wp1 = wp[1];
  Weighted_point wp2 = wp[2];
  Weighted_point wp3 = wp[3];
  Weighted_point wp01 = wp[4];
  Weighted_point wp02 = wp[5];
  Weighted_point wp03 = wp[6];
  Weighted_point wp04 = wp[7];
  Weighted_point wp05 = wp[8];

  Bare_point q0 = wp[9].point();
  Weighted_point wq0 = wp[9];
  Weighted_point wq1 = wp[10];
  Weighted_point wq2 = wp[11];
  Weighted_point wq3 = wp[12];
  Weighted_point wq4 = wp[13];
  Weighted_point wq01 = wp[14];
  Weighted_point wq11 = wp[15];
  Weighted_point wq21 = wp[16];
  Weighted_point wq31 = wp[17];
  Weighted_point wq41 = wp[18];

  // test of Construct_weighted_circumcenter_3 and compare_power_distance
  std::cout << "test of Construct_weighted_circumcenter_3" << std::endl;
  std::cout << "test of Compare_power_distance_3" << std::endl;
  Construct_weighted_circumcenter_3 weighted_circumcenter = traits.construct_weighted_circumcenter_3_object();
  Compare_power_distance_3 compare_power_distance = traits.compare_power_distance_3_object();

  {
    Bare_point c;
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
    Bare_point c;
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
    Bare_point c;
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
    Bare_point c;
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
    Bare_point c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[1], o[0],o[0],o[1]) == CGAL::SMALLER);
    assert (compare_power_distance(c,wp[6],wp[1], o[0],o[1],o[0]) == CGAL::LARGER);
  }

  // test of Does_simplex_intersect_dual_support_3
  std::cout << "test of Does_simplex_intersect_dual_support_3" << std::endl;
  Does_simplex_intersect_dual_support_3 does_simplex_intersect_dual_support = traits.does_simplex_intersect_dual_support_3_object();
  assert(does_simplex_intersect_dual_support(wp0, wp1, wp2, wp3) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1, wp0, wp2, wp3) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp01, wp1, wp2, wp3) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp01, wp2, wp1, wp3) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp02, wp1, wp2, wp3) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp2, wp1, wp02, wp3) == CGAL::ON_BOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0, wp1, wp2) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp1, wp0, wp2) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp01, wp1, wp2) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp01, wp2, wp1) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp03, wp1, wp2) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp2, wp1, wp03) == CGAL::ON_UNBOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0, wp1) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1, wp0) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp04, wp1) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp1, wp04) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp05, wp1) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1, wp05) == CGAL::ON_UNBOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0, wp1, wp2, wp3, o[1], o[1], o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1, wp0, wp2, wp3, o[1], o[1], o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp01, wp1, wp2, wp3, o[1], o[1], o[1], o[1]) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp01, wp2, wp1, wp3, o[1], o[1], o[1], o[1]) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp02, wp1, wp2, wp3, o[1], o[1], o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp2, wp1, wp02, wp3, o[1], o[1], o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0, wp1, wp2, o[1], o[1], o[1]) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp1, wp0, wp2, o[1], o[1], o[1]) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp01, wp1, wp2, o[1], o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp01, wp2, wp1, o[1], o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp03, wp1, wp2, o[1], o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp2, wp1, wp03, o[1], o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0, wp1, o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1, wp0, o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp04, wp1, o[1], o[1]) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp1, wp04, o[1], o[1]) == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp05, wp1, o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1, wp05, o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp01, wp2, wp1, wp3, o[1], o[0], o[1], o[0]) == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp02, wp1, wp2, wp3, o[1], o[0], o[1], o[0]) == CGAL::ON_UNBOUNDED_SIDE);

  // test In_smallest_orthogonal_sphere_3
  // test Side_of_bounded_orthogonal_sphere_3
  std::cout << "test In_smallest_orthogonal_sphere_3" << std::endl;
  std::cout << "test Side_of_bounded_orthogonal_sphere_3" << std::endl;
  typedef typename Traits::In_smallest_orthogonal_sphere_3 In_smallest_orthogonal_sphere_3;
  typedef typename Traits::Side_of_bounded_orthogonal_sphere_3 Side_of_bounded_orthogonal_sphere_3;
  In_smallest_orthogonal_sphere_3 in_smallest_orthogonal_sphere = traits.in_smallest_orthogonal_sphere_3_object();
  Side_of_bounded_orthogonal_sphere_3 side_of_bounded_orthogonal_sphere = traits.side_of_bounded_orthogonal_sphere_3_object();

   typedef typename Traits::FT  FT;
   FT ww02 = FT(2)/FT(3);
   Weighted_point wq02(q0, ww02);

  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq0) == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq01) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq2, wq1, wq01) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq11, wq21, wq0) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq21, wq11, wq0) == CGAL::POSITIVE);
  assert(side_of_bounded_orthogonal_sphere(wq21, wq11, wq0) == CGAL::ON_UNBOUNDED_SIDE);

  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, wq3) == CGAL::POSITIVE);
  assert(side_of_bounded_orthogonal_sphere(wq1, wq0, wq2, wq3) == CGAL::ON_UNBOUNDED_SIDE);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq3, wq0) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq3, wq2, wq0) == CGAL::NEGATIVE);
//  assert(in_smallest_orthogonal_sphere(wq11, wq21, wq31, wq02) == CGAL::ZERO);
//  assert(in_smallest_orthogonal_sphere(wq31, wq21, wq11, wq02) == CGAL::ZERO);

  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, wq4) == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, wq3, wq4) == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq01, wq11, wq21, wq31, wq4) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq01, wq21, wq11, wq31, wq4) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, wq41) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq3, wq2, wq41) == CGAL::NEGATIVE);
  assert(side_of_bounded_orthogonal_sphere(wq0, wq1, wq3, wq2, wq41) == CGAL::ON_BOUNDED_SIDE);
///
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq0, o[1], o[1], o[1]) == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq01, o[1], o[1], o[1]) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq2, wq1, wq01, o[1], o[1], o[1]) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq11, wq21, wq0, o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq21, wq11, wq0, o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(side_of_bounded_orthogonal_sphere(wq21, wq11, wq0, o[1], o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);

  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, o[1], o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, wq3, o[1], o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(side_of_bounded_orthogonal_sphere(wq1, wq0, wq2, wq3, o[1], o[1], o[1], o[1]) == CGAL::ON_UNBOUNDED_SIDE);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq3, wq0, o[1], o[1], o[1], o[1]) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq3, wq2, wq0, o[1], o[1], o[1], o[1]) == CGAL::NEGATIVE);

  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, wq4, o[1], o[1], o[1], o[1], o[1]) == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, wq3, wq4, o[1], o[1], o[1], o[1], o[1]) == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq01, wq11, wq21, wq31, wq4, o[1], o[1], o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq01, wq21, wq11, wq31, wq4, o[1], o[1], o[1], o[1], o[1]) == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, wq41, o[1], o[1], o[1], o[1], o[1]) == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq3, wq2, wq41, o[1], o[1], o[1], o[1], o[1]) == CGAL::NEGATIVE);
  assert(side_of_bounded_orthogonal_sphere(wq0, wq1, wq3, wq2, wq41, o[1], o[1], o[1], o[1], o[1]) == CGAL::ON_BOUNDED_SIDE);


  // test weighted_circumcenter
  // test squared_radius_smallest_orthogonal_sphere
  // test critical_squared_radius
  std::cout << "test of squared_radius_smallest_orthogonal_sphere"
      << std::endl;
   std::cout << "test of critical_squared_radius" << std::endl;
  typedef typename Traits::Compute_power_product_3 Compute_power_product_3;
  typedef typename Traits::Compute_squared_radius_smallest_orthogonal_sphere_3 Compute_squared_radius_smallest_orthogonal_sphere_3;
  typedef typename Traits::Compute_critical_squared_radius_3 Compute_critical_squared_radius_3;
  Compute_power_product_3 power_product = traits.compute_power_product_3_object();
  Compute_squared_radius_smallest_orthogonal_sphere_3 squared_radius_smallest_orthogonal_sphere = traits.compute_squared_radius_smallest_orthogonal_sphere_3_object();
  Compute_critical_squared_radius_3 critical_squared_radius = traits.compute_critical_squared_radius_3_object();

  Weighted_point wc(weighted_circumcenter(wq11, wq21, wq31, wq41), squared_radius_smallest_orthogonal_sphere(wq11, wq21, wq31, wq41));
  Weighted_point wt(Bare_point(1., 1., 1.), 0.);
  // this test requires a weighted point with a zero weight
  assert(power_product(wc, wt) == critical_squared_radius(wq11, wq21, wq31, wq41, wt));
  assert(power_product(wc, wt) == critical_squared_radius(wq11, wq21, wq31, wq41, wt, o[1], o[1], o[1], o[1], o[1]));

  wc = Weighted_point(weighted_circumcenter(wp0, wp1, wp2, wp3), squared_radius_smallest_orthogonal_sphere(wp0, wp1, wp2, wp3));
//  assert(power_product(wc, wt) == critical_squared_radius(wp0, wp1, wp2, wp3, wt));

  wc = Weighted_point(weighted_circumcenter(wp01, wp1, wp2, wp3), squared_radius_smallest_orthogonal_sphere(wp01, wp1, wp2, wp3));
  assert(power_product(wc, wt) == critical_squared_radius(wp01, wp1, wp2, wp3,wt));
  assert(power_product(wc, wt) == critical_squared_radius(wp01, wp1, wp2, wp3,wt, o[1], o[1], o[1], o[1], o[1]));
}

template<class K>
void _test_cls_periodic_3_regular_triangulation_traits_3 ()
{
  typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
  typedef CGAL::Periodic_3_Regular_triangulation_traits_3<Regular_traits> Traits;
  typedef typename Traits::Weighted_point Weighted_point;
  typedef typename Traits::Bare_point Bare_point;
  typedef typename Traits::Iso_cuboid_3 Iso_cuboid;

  Traits traits;

  // Test Iso_cuboid(0,0,0,1,1,1)
  std::cout << "Iso_cuboid(0,0,0,4,4,4)" << std::endl;
  {
    traits.set_domain(Iso_cuboid(0, 0, 0, 4, 4, 4));
    Bare_point p0(0.,0.,0.);
    Bare_point p1(3.,0.,0.);
    Bare_point p2(0.,3.,0.);
    Bare_point p3(0.,0.,3.);

    Bare_point q0(0.,0.,0.);
    Bare_point q1(2.,0.,0.);
    Bare_point q2(0.,2.,0.);
    Bare_point q3(0.,0.,2.);
    Bare_point q4(2.,2.,2.);

    Weighted_point p[19] =
    {
        Weighted_point(p0,9.),
        Weighted_point(p1,9.),
        Weighted_point(p2,9.),
        Weighted_point(p3,9.),
        Weighted_point(p0,6.),
        Weighted_point(p0,3.),
        Weighted_point(p0,12.),
        Weighted_point(p0,18.),
        Weighted_point(p0,24.),

        Weighted_point(q0,0.),
        Weighted_point(q1,0.),
        Weighted_point(q2,0.),
        Weighted_point(q3,0.),
        Weighted_point(q4,0.),
        Weighted_point(q0,2.),
        Weighted_point(q1,2.),
        Weighted_point(q2,2.),
        Weighted_point(q3,2.),
        Weighted_point(q4,2.)
    };
    _test_for_given_domain(traits, p);
  }
  std::cout << "Iso_cuboid(-2,-2,-2,2,2,2)" << std::endl;
  {
    traits.set_domain(Iso_cuboid(-2, -2, -2, 2, 2, 2));
    Bare_point p0(-2.,-2.,-2.);
    Bare_point p1(1.,-2.,-2.);
    Bare_point p2(-2.,1.,-2.);
    Bare_point p3(-2.,-2.,1.);

    Bare_point q0(-1.,-1.,-1.);
    Bare_point q1(1.,-1.,-1.);
    Bare_point q2(-1.,1.,-1.);
    Bare_point q3(-1.,-1.,1.);
    Bare_point q4(1.,1.,1.);

    Weighted_point p[19] =
    {
        Weighted_point(p0,9.),
        Weighted_point(p1,9.),
        Weighted_point(p2,9.),
        Weighted_point(p3,9.),
        Weighted_point(p0,6.),
        Weighted_point(p0,3.),
        Weighted_point(p0,12.),
        Weighted_point(p0,18.),
        Weighted_point(p0,24.),

        Weighted_point(q0,0.),
        Weighted_point(q1,0.),
        Weighted_point(q2,0.),
        Weighted_point(q3,0.),
        Weighted_point(q4,0.),
        Weighted_point(q0,2.),
        Weighted_point(q1,2.),
        Weighted_point(q2,2.),
        Weighted_point(q3,2.),
        Weighted_point(q4,2.)
    };
    _test_for_given_domain(traits, p);
  }

//  _test_for_given_domain(traits, p);
//
//  // Test Iso_cuboid(-0.5,-0.5,-0.5,0.5,0.5,0.5)
//  traits.set_domain(Iso_cuboid(-1, -1, -1, 1, 1, 1, 2));
//
//  p[0] = Point(-4, -4, -4, 8);
//  p[1] = Point(0, -4, -4, 8);
//  p[2] = Point(-4, 0, -4, 8);
//  p[3] = Point(-4, -4, 0, 8);
//  p[4] = Point(0, 0, 0, 8);
//  p[5] = Point(-2, 0, -2, 8);
//  p[6] = Point(0, 2, 2, 8);
//  p[7] = Point(0, -4, 0, 8);
//  p[8] = Point(-2, -4, -2, 8);
//  p[9] = Point(0, -4, 2, 8);
//  p[10] = Point(-3, -3, -2, 8);
//
//  _test_for_given_domain(traits, p);
}
