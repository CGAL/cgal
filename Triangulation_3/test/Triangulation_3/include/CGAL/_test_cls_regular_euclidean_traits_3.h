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

// Author(s)     :  Mariette Yvinec
//                  Mael Rouxel-Labbé

#include <CGAL/use.h>

#include <cassert>

template <class Traits>
void
_test_cls_regular_euclidean_traits_3 (const Traits &)
{
  typedef typename Traits::Point_3           Point_3;
  typedef typename Traits::Weighted_point_3  Weighted_point_3;
  typedef typename Traits::Segment_3         Segment_3;

  typedef typename Traits::Power_side_of_oriented_power_sphere_3
                                  Power_side_of_oriented_power_sphere_3;
  typedef typename Traits::Compare_power_distance_3
                                  Compare_power_distance_3;
  typedef typename Traits::Power_side_of_bounded_power_sphere_3
                                  Power_side_of_bounded_power_sphere_3;
  typedef typename Traits::Construct_weighted_circumcenter_3
                                  Construct_weighted_circumcenter_3;
  typedef typename Traits::Compute_power_product_3
                                  Compute_power_product_3;
  typedef typename Traits::
                   Compute_squared_radius_smallest_orthogonal_sphere_3
                   Compute_squared_radius_smallest_orthogonal_sphere_3;

  typedef typename Traits::
                   Compute_power_distance_to_power_sphere_3
                   Compute_power_distance_to_power_sphere_3;

  CGAL_USE_TYPE(Segment_3);
  Traits traits;
  Power_side_of_oriented_power_sphere_3 power_test =
    traits.power_side_of_oriented_power_sphere_3_object();
  Compare_power_distance_3 compare_power_distance =
    traits.compare_power_distance_3_object();
  Power_side_of_bounded_power_sphere_3 bounded_power_test =
    traits.power_side_of_bounded_power_sphere_3_object();
  Construct_weighted_circumcenter_3 weighted_circumcenter =
    traits.construct_weighted_circumcenter_3_object();
  Compute_power_product_3 power_product =
    traits.compute_power_product_3_object();
  Compute_squared_radius_smallest_orthogonal_sphere_3
    squared_radius_smallest_orthogonal_sphere =
    traits.compute_squared_radius_smallest_orthogonal_sphere_3_object();
  Compute_power_distance_to_power_sphere_3  compute_power_distance_to_power_sphere_3 =
    traits.compute_power_distance_to_power_sphere_3_object();


  // test of Does_simplex_intersect_dual_support_3
  std::cout << "test of Does_simplex_intersect_dual_support_3" << std::endl;
  Point_3 p0(0.,0.,0.);
  Point_3 p1(3.,0.,0.);
  Point_3 p2(0.,3.,0.);
  Point_3 p3(0.,0.,3.);
  Weighted_point_3 wp0(p0,9.);
  Weighted_point_3 wp1(p1,9.);
  Weighted_point_3 wp2(p2,9.);
  Weighted_point_3 wp3(p3,9.);
  Weighted_point_3 wp01(p0,6.);
  Weighted_point_3 wp02(p0,3.);
  Weighted_point_3 wp03(p0,12.);
  Weighted_point_3 wp04(p0,18.);
  Weighted_point_3 wp05(p0,24.);

  // test of Construct_weighted_circumcenter_3 and compare_power_distance
  std::cout << "test of Construct_weighted_circumcenter_3" << std::endl;
  std::cout << "test Of Compare_power_distance_3" << std::endl;
  Point_3 c ;
  c = weighted_circumcenter(wp0,wp1,wp2,wp3);

  assert (compare_power_distance(c,wp0,wp1) == CGAL::EQUAL);

  assert (compare_power_distance(c,wp01,wp1) == CGAL::LARGER);
  assert (compare_power_distance(c,wp03,wp1) == CGAL::SMALLER);

  c = weighted_circumcenter(wp01,wp1,wp2,wp3);
  assert (compare_power_distance(c,wp01,wp2) == CGAL::EQUAL);

  c = weighted_circumcenter(wp02,wp1,wp2,wp3);
  assert (compare_power_distance(c,wp02,wp3) ==  CGAL::EQUAL);

  // test Power_side_of_bounded_power_sphere_3
  std::cout << " test Power_side_of_bounded_power_sphere_3" << std::endl;
  Point_3 q0(0.,0.,0.);
  Point_3 q1(2.,0.,0.);
  Point_3 q2(0.,2.,0.);
  Point_3 q3(0.,0.,2.);
  Point_3 q4(2.,2.,2.);
  Point_3 q5(-2.,0.,0.);
  Point_3 q6(-4.,0.,0.);

  Weighted_point_3 wq0(q0,0.);
  Weighted_point_3 wq1(q1,0.);
  Weighted_point_3 wq2(q2,0.);
  Weighted_point_3 wq3(q3,0.);
  Weighted_point_3 wq4(q4,0.);
  Weighted_point_3 wq5(q5,0.);
  Weighted_point_3 wq6(q6,0.);
  Weighted_point_3 wq01(q0,2.);
  Weighted_point_3 wq11(q1,2.);
  Weighted_point_3 wq21(q2,2.);
  Weighted_point_3 wq31(q3,2.);
  Weighted_point_3 wq41(q4,2.);

  typedef typename Traits::FT  FT;
  FT ww02 = FT(2)/FT(3);
  Weighted_point_3 wq02(q0, ww02);

  assert(bounded_power_test(wq0, wq1, wq2) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq1, wq0, wq2) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq1, wq2, wq0) == CGAL::ON_BOUNDARY);
  assert(bounded_power_test(wq1, wq2, wq01) == CGAL::ON_BOUNDED_SIDE);
  assert(bounded_power_test(wq2, wq1, wq01) == CGAL::ON_BOUNDED_SIDE);
  assert(bounded_power_test(wq11, wq21, wq0) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq21, wq11, wq0) == CGAL::ON_UNBOUNDED_SIDE);

  assert(bounded_power_test(wq0, wq1, wq2, wq3) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq1, wq0, wq2, wq3) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq1, wq2, wq3, wq0) == CGAL::ON_BOUNDED_SIDE);
  assert(bounded_power_test(wq1, wq3, wq2, wq0) == CGAL::ON_BOUNDED_SIDE);
  assert(bounded_power_test(wq11, wq21, wq31, wq02) == CGAL::ON_BOUNDARY);
  assert(bounded_power_test(wq31, wq21, wq11, wq02) == CGAL::ON_BOUNDARY);

  assert(bounded_power_test(wq0, wq1, wq2, wq3, wq4) == CGAL::ON_BOUNDARY);
  assert(bounded_power_test(wq1, wq0, wq2, wq3, wq4) == CGAL::ON_BOUNDARY);
  assert(bounded_power_test(wq01, wq11, wq21, wq31, wq4) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq01, wq21, wq11, wq31, wq4) == CGAL::ON_UNBOUNDED_SIDE);
  assert(bounded_power_test(wq0, wq1, wq2, wq3, wq41) == CGAL::ON_BOUNDED_SIDE);
  assert(bounded_power_test(wq0, wq1, wq3, wq2, wq41) == CGAL::ON_BOUNDED_SIDE);

  // test weighted_circumcenter
  // test squared_radius_smallest_orthogonal_sphere
  std::cout << "test of weighted_circumcenter" << std::endl;
   std::cout << "test of squared_radius_smallest_orthogonal_sphere" << std::endl;
  Weighted_point_3 wc(weighted_circumcenter(wq11,wq21,wq31,wq41),
                      squared_radius_smallest_orthogonal_sphere(wq11,wq21,wq31,wq41));
  Weighted_point_3 wt(Point_3(1.,1.,1.), 0.);
  // this test requires a weighted point with a zero weight
  assert( power_product(wc,wt) ==
          compute_power_distance_to_power_sphere_3(wq11,wq21,wq31,wq41,wt));

  wc = Weighted_point_3(weighted_circumcenter(wp0,wp1,wp2,wp3),
                        squared_radius_smallest_orthogonal_sphere(wp0,wp1,wp2,wp3));
  assert( power_product(wc,wt) ==
          compute_power_distance_to_power_sphere_3(wp0,wp1,wp2,wp3,wt));

  wc = Weighted_point_3(weighted_circumcenter(wp01,wp1,wp2,wp3),
                        squared_radius_smallest_orthogonal_sphere(wp01,wp1,wp2,wp3));

  assert( power_product(wc,wt) ==
          compute_power_distance_to_power_sphere_3(wp01,wp1,wp2,wp3,wt));

  // test power_test
  // null weights
  assert(power_test(wq0,wq1,wq2,wq3,wq4) ==
         traits.side_of_oriented_sphere_3_object()(q0,q1,q2,q3,q4));
  assert(power_test(wq1,wq0,wq2,wq3,wq6) ==
         traits.side_of_oriented_sphere_3_object()(q1,q0,q2,q3,q6));
  assert(power_test(wq0,wq1,wq2,wq5) == CGAL::ON_NEGATIVE_SIDE &&
         traits.coplanar_side_of_bounded_circle_3_object()(q0,q1,q2,q5) ==
                                                       CGAL::ON_UNBOUNDED_SIDE);

  // wc = (1,1,1) -3
  assert(power_test(wp01,wp1,wp2,wp3,wc) == CGAL::ON_NEGATIVE_SIDE);
  Weighted_point_3 wt2(wc.point(), wc.weight()+6);
  assert(power_test(wp01,wp1,wp2,wp3,wt2) == CGAL::ON_ORIENTED_BOUNDARY);
  Weighted_point_3 wt3(wc.point(), wc.weight()+7);
  assert(power_test(wp01,wp1,wp2,wp3,wt3) == CGAL::ON_POSITIVE_SIDE);
  Weighted_point_3 wt4(Point_3(0.,2.,2.), 6.);
  assert(power_test(wp01,wp1,wp2,wp3,wt4) == CGAL::ON_ORIENTED_BOUNDARY);

  // wc = (1,2,0) -1
  wc = Weighted_point_3(weighted_circumcenter(wp01,wp1,wq21),
                      squared_radius_smallest_orthogonal_sphere(wp01,wp1,wq21));
  assert(power_test(wp01,wp1,wq21,wc) == CGAL::ON_NEGATIVE_SIDE);
  wt = Weighted_point_3(Point_3(-1.,1.,0.),6);
  assert(power_test(wp01,wp1,wq21,wt) == CGAL::ON_ORIENTED_BOUNDARY);
  wt = Weighted_point_3(Point_3(3.,4.,0.),10);
  assert(power_test(wp01,wp1,wq21,wt) == CGAL::ON_POSITIVE_SIDE);

  // wc = (0,0,3) -9
  wc = Weighted_point_3(weighted_circumcenter(wp04,wp3),
                      squared_radius_smallest_orthogonal_sphere(wp04,wp3));
  assert(power_test(wp04,wp3,wc) == CGAL::ON_NEGATIVE_SIDE);
  wt = Weighted_point_3(Point_3(0.,0.,7.),25);
  assert(power_test(wp04,wp3,wt) == CGAL::ON_ORIENTED_BOUNDARY);
  wt = Weighted_point_3(Point_3(0.,0.,2.),12);
  assert(power_test(wp04,wp3,wt) == CGAL::ON_POSITIVE_SIDE);

  assert(power_test(wq0,wq02) == CGAL::ON_POSITIVE_SIDE);
  assert(power_test(wq01,wq01) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(power_test(wq01,wq02) == CGAL::ON_NEGATIVE_SIDE);
}
