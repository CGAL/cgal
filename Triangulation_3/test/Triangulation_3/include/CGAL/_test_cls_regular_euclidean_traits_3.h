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

template <class Traits>
void
_test_cls_regular_euclidean_traits_3 (const Traits &)
{
  typedef typename Traits::Weighted_point  Weighted_point;
  typedef typename Traits::Bare_point      Bare_point;
  typedef typename Traits::Segment_3       Segment_3;
  typedef typename Traits::Power_test_3    Power_test_3;
  typedef typename Traits::Compare_power_distance_3
                                  Compare_power_distance_3;
  typedef typename Traits::In_smallest_orthogonal_sphere_3
                                  In_smallest_orthogonal_sphere_3;
  typedef typename Traits::Side_of_bounded_orthogonal_sphere_3
                                  Side_of_bounded_orthogonal_sphere_3;
  typedef typename Traits::Does_simplex_intersect_dual_support_3 
                                  Does_simplex_intersect_dual_support_3;
  typedef typename Traits::Construct_weighted_circumcenter_3
                                  Construct_weighted_circumcenter_3;
  typedef typename Traits::Compute_power_product_3
                                  Compute_power_product_3;
  typedef typename Traits::
                   Compute_squared_radius_smallest_orthogonal_sphere_3
                   Compute_squared_radius_smallest_orthogonal_sphere_3;

  typedef typename Traits::
                   Compute_critical_squared_radius_3
                   Compute_critical_squared_radius_3;        

  Traits traits;
  Power_test_3 power_test =  traits.power_test_3_object();
  Compare_power_distance_3 compare_power_distance =
    traits.compare_power_distance_3_object();
  In_smallest_orthogonal_sphere_3 in_smallest_orthogonal_sphere =
    traits.in_smallest_orthogonal_sphere_3_object();
  Side_of_bounded_orthogonal_sphere_3 side_of_bounded_orthogonal_sphere =
    traits.side_of_bounded_orthogonal_sphere_3_object();
  Does_simplex_intersect_dual_support_3 does_simplex_intersect_dual_support =
    traits.does_simplex_intersect_dual_support_3_object();
  Construct_weighted_circumcenter_3 weighted_circumcenter =
    traits.construct_weighted_circumcenter_3_object();
  Compute_power_product_3 power_product =
    traits.compute_power_product_3_object();
  Compute_squared_radius_smallest_orthogonal_sphere_3
    squared_radius_smallest_orthogonal_sphere =
    traits.compute_squared_radius_smallest_orthogonal_sphere_3_object();
  Compute_critical_squared_radius_3  critical_squared_radius =
    traits.compute_critical_squared_radius_3_object();
  

  // test of Does_simplex_intersect_dual_support_3
  std::cout << "test of Does_simplex_intersect_dual_support_3" << std::endl;
  Bare_point p0(0.,0.,0.);
  Bare_point p1(3.,0.,0.);
  Bare_point p2(0.,3.,0.);
  Bare_point p3(0.,0.,3.);

  Weighted_point wp0(p0,9.);
  Weighted_point wp1(p1,9.);
  Weighted_point wp2(p2,9.);
  Weighted_point wp3(p3,9.);
  Weighted_point wp01(p0,6.);
  Weighted_point wp02(p0,3.);
  Weighted_point wp03(p0,12.);
  Weighted_point wp04(p0,18.);
  Weighted_point wp05(p0,24.);

  assert(does_simplex_intersect_dual_support(wp0,wp1,wp2,wp3) 
	 == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1,wp0,wp2,wp3) 
	 == CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp01,wp1,wp2,wp3) 
	 ==  CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp01,wp2,wp1,wp3) 
	 ==  CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp02,wp1,wp2,wp3) 
	 ==  CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp2,wp1,wp02,wp3) 
	 ==  CGAL::ON_BOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0,wp1,wp2)
	 == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp1,wp0,wp2)
	 == CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp01,wp1,wp2)
	 ==  CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp01,wp2,wp1)
	 ==  CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp03,wp1,wp2)
	 ==  CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp2,wp1,wp03)
	 ==  CGAL::ON_UNBOUNDED_SIDE);

  assert(does_simplex_intersect_dual_support(wp0,wp1)
	 == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1,wp0)
	 == CGAL::ON_BOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp04,wp1)
	 ==  CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp1,wp04)
	 ==  CGAL::ON_BOUNDARY);
  assert(does_simplex_intersect_dual_support(wp05,wp1)
	 ==  CGAL::ON_UNBOUNDED_SIDE);
  assert(does_simplex_intersect_dual_support(wp1,wp05)
	 ==  CGAL::ON_UNBOUNDED_SIDE);
 
  // test of Construct_weighted_circumcenter_3 and compare_power_distance
   std::cout << "test of Construct_weighted_circumcenter_3" << std::endl;
   std::cout << "test Of Compare_power_distance_3" << std::endl;
   Bare_point c ;
  c = weighted_circumcenter(wp0,wp1,wp2,wp3);
  assert (compare_power_distance(c,wp0,wp1) == CGAL::EQUAL);
  assert (compare_power_distance(c,wp01,wp1) == CGAL::LARGER);
  assert (compare_power_distance(c,wp03,wp1) == CGAL::SMALLER);

  c = weighted_circumcenter(wp01,wp1,wp2,wp3);
  assert (compare_power_distance(c,wp01,wp2) == CGAL::EQUAL);
	  
  c = weighted_circumcenter(wp02,wp1,wp2,wp3);
  assert (compare_power_distance(c,wp02,wp3) ==  CGAL::EQUAL);

  // test In_smallest_orthogonal_sphere_3
  // test Side_of_bounded_orthogonal_sphere_3
  std::cout << " test In_smallest_orthogonal_sphere_3" << std::endl;
  std::cout << " test Side_of_bounded_orthogonal_sphere_3" << std::endl;
  Bare_point q0(0.,0.,0.);
  Bare_point q1(2.,0.,0.);
  Bare_point q2(0.,2.,0.);
  Bare_point q3(0.,0.,2.);
  Bare_point q4(2.,2.,2.);
  Bare_point q5(-2.,0.,0.);
  Bare_point q6(-4.,0.,0.);

  Weighted_point wq0(q0,0.);
  Weighted_point wq1(q1,0.);
  Weighted_point wq2(q2,0.);
  Weighted_point wq3(q3,0.);
  Weighted_point wq4(q4,0.);
  Weighted_point wq5(q5,0.);
  Weighted_point wq6(q6,0.);
  Weighted_point wq01(q0,2.);
  Weighted_point wq11(q1,2.);
  Weighted_point wq21(q2,2.);
  Weighted_point wq31(q3,2.);
  Weighted_point wq41(q4,2.);

  typedef typename Traits::FT  FT;
  FT ww02 = FT(2)/FT(3);
  Weighted_point wq02(q0, ww02);
  
 
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2) 
	 == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2) 
	 == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq0) 
	 == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq01)
	 == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq2, wq1, wq01)
	 == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq11, wq21, wq0)
	 == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq21, wq11, wq0)
	 == CGAL::POSITIVE);
  assert(side_of_bounded_orthogonal_sphere(wq21, wq11, wq0)
	 == CGAL::ON_UNBOUNDED_SIDE);

  
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3)
	 == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, wq3)
	 == CGAL::POSITIVE);
  assert(side_of_bounded_orthogonal_sphere(wq1, wq0, wq2, wq3)
	 == CGAL::ON_UNBOUNDED_SIDE);
  assert(in_smallest_orthogonal_sphere(wq1, wq2, wq3, wq0)
	 == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq1, wq3, wq2, wq0)
	 == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq11, wq21, wq31, wq02)
 	 == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq31, wq21, wq11, wq02)
 	 == CGAL::ZERO);
  
  
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, wq4)
	 == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq1, wq0, wq2, wq3, wq4)
	 == CGAL::ZERO);
  assert(in_smallest_orthogonal_sphere(wq01, wq11, wq21, wq31, wq4)
	 == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq01, wq21, wq11, wq31, wq4)
	 == CGAL::POSITIVE);
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq2, wq3, wq41)
	 == CGAL::NEGATIVE);
  assert(in_smallest_orthogonal_sphere(wq0, wq1, wq3, wq2, wq41)
	 == CGAL::NEGATIVE);
  assert(side_of_bounded_orthogonal_sphere(wq0, wq1, wq3, wq2, wq41)
	 == CGAL::ON_BOUNDED_SIDE);

  
  // test weighted_circumcenter
  // test squared_radius_smallest_orthogonal_sphere
  // test critical_squared_radius
  std::cout << "test of  squared_radius_smallest_orthogonal_sphere" 
	    << std::endl;
   std::cout << "test of critical_squared_radius" << std::endl;
  Weighted_point wc(
              weighted_circumcenter(wq11,wq21,wq31,wq41),
	      squared_radius_smallest_orthogonal_sphere(wq11,wq21,wq31,wq41));
  Weighted_point wt(Bare_point(1.,1.,1.), 0.);
  // this test requires a weighted point with a zero weight 
  assert( power_product(wc,wt) == 
	  critical_squared_radius(wq11,wq21,wq31,wq41,wt));

  wc = Weighted_point(
           weighted_circumcenter(wp0,wp1,wp2,wp3),
	   squared_radius_smallest_orthogonal_sphere(wp0,wp1,wp2,wp3));
  assert( power_product(wc,wt) == 
	  critical_squared_radius(wp0,wp1,wp2,wp3,wt));

  wc = Weighted_point(
           weighted_circumcenter(wp01,wp1,wp2,wp3),
	   squared_radius_smallest_orthogonal_sphere(wp01,wp1,wp2,wp3));
  assert( power_product(wc,wt) == 
	  critical_squared_radius(wp01,wp1,wp2,wp3,wt));
  
  // test power_test
  // null weights
  assert(power_test(wq0,wq1,wq2,wq3,wq4) ==
         traits.side_of_oriented_sphere_3_object()(q0,q1,q2,q3,q4));
  assert(power_test(wq1,wq0,wq2,wq3,wq6) ==
         traits.side_of_oriented_sphere_3_object()(q1,q0,q2,q3,q6));
  assert(power_test(wq0,wq1,wq2,wq5) == CGAL::ON_NEGATIVE_SIDE &&
         traits.coplanar_side_of_bounded_circle_3_object()(q0,q1,q2,q5) ==
                                                       CGAL::ON_UNBOUNDED_SIDE);
  assert(power_test(wq5,wq1,wq0) == CGAL::ON_POSITIVE_SIDE &&
         (Segment_3(wq5,wq1)).has_on(wq0));
  assert(power_test(wq5,wq1,wq5) == CGAL::ON_ORIENTED_BOUNDARY &&
         (Segment_3(wq5,wq1)).has_on(wq5));
  assert(power_test(wq5,wq1,wq6) == CGAL::ON_NEGATIVE_SIDE &&
         !(Segment_3(wq5,wq1)).has_on(wq6));
  assert(power_test(wq3,wq2) == CGAL::ON_ORIENTED_BOUNDARY);

  // wc = (1,1,1) -3
  assert(power_test(wp01,wp1,wp2,wp3,wc) == CGAL::ON_NEGATIVE_SIDE);
  Weighted_point wt2(wc.point(), wc.weight()+6);
  assert(power_test(wp01,wp1,wp2,wp3,wt2) == CGAL::ON_ORIENTED_BOUNDARY);
  Weighted_point wt3(wc.point(), wc.weight()+7);
  assert(power_test(wp01,wp1,wp2,wp3,wt3) == CGAL::ON_POSITIVE_SIDE);
  Weighted_point wt4(Bare_point(0.,2.,2.), 6.);
  assert(power_test(wp01,wp1,wp2,wp3,wt4) == CGAL::ON_ORIENTED_BOUNDARY);

  // wc = (1,2,0) -1
  wc = Weighted_point(weighted_circumcenter(wp01,wp1,wq21),
                      squared_radius_smallest_orthogonal_sphere(wp01,wp1,wq21));
  assert(power_test(wp01,wp1,wq21,wc) == CGAL::ON_NEGATIVE_SIDE);
  wt = Weighted_point(Bare_point(-1.,1.,0.),6);
  assert(power_test(wp01,wp1,wq21,wt) == CGAL::ON_ORIENTED_BOUNDARY);
  wt = Weighted_point(Bare_point(3.,4.,0.),10);
  assert(power_test(wp01,wp1,wq21,wt) == CGAL::ON_POSITIVE_SIDE);

  // wc = (0,0,3) -9
  wc = Weighted_point(weighted_circumcenter(wp04,wp3),
                      squared_radius_smallest_orthogonal_sphere(wp04,wp3));
  assert(power_test(wp04,wp3,wc) == CGAL::ON_NEGATIVE_SIDE);
  wt = Weighted_point(Bare_point(0.,0.,7.),25);
  assert(power_test(wp04,wp3,wt) == CGAL::ON_ORIENTED_BOUNDARY);
  wt = Weighted_point(Bare_point(0.,0.,2.),12);
  assert(power_test(wp04,wp3,wt) == CGAL::ON_POSITIVE_SIDE);

  assert(power_test(wq0,wq02) == CGAL::ON_POSITIVE_SIDE);
  assert(power_test(wq01,wq01) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(power_test(wq01,wq02) == CGAL::ON_NEGATIVE_SIDE);
}
