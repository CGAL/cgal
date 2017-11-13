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
//                  Mael Rouxel-Labbé

#if (__GNUC__>4) || (__GNUC__ == 4 && __GNUC_MINOR__ >=6)
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>

#include <CGAL/use.h>
#include <CGAL/number_type_config.h>

#include <cassert>
#include <iostream>

template<class Traits>
void _test_for_given_domain (const Traits& traits,
                             typename Traits::Weighted_point_3* wp)
{
  typedef typename Traits::FT                           FT;
  typedef typename Traits::Point_3                      Point_3;
  typedef typename Traits::Periodic_3_offset_3          Offset;

  CGAL_USE_TYPE(typename Traits::Weighted_point_3);
  CGAL_USE_TYPE(typename Traits::Vector_3);
  CGAL_USE_TYPE(typename Traits::Iso_cuboid_3);

  CGAL_USE_TYPE(typename Traits::Segment_3);
  CGAL_USE_TYPE(typename Traits::Triangle_3);
  CGAL_USE_TYPE(typename Traits::Tetrahedron_3);

  CGAL_USE_TYPE(typename Traits::Comparison_result);
  CGAL_USE_TYPE(typename Traits::Orientation);
  CGAL_USE_TYPE(typename Traits::Oriented_side);
  CGAL_USE_TYPE(typename Traits::Bounded_side);

  CGAL_USE_TYPE(typename Traits::Compare_xyz_3);
  CGAL_USE_TYPE(typename Traits::Orientation_3);
  CGAL_USE_TYPE(typename Traits::Coplanar_orientation_3);

  CGAL_USE_TYPE(typename Traits::Construct_segment_3);
  CGAL_USE_TYPE(typename Traits::Construct_triangle_3);
  CGAL_USE_TYPE(typename Traits::Construct_tetrahedron_3);

  typedef typename Traits::Construct_point_3 Construct_point_3;
  typedef typename Traits::Power_side_of_oriented_power_sphere_3 Power_side_of_oriented_power_sphere_3;
  typedef typename Traits::Compare_power_distance_3 Compare_power_distance_3;
  typedef typename Traits::Construct_weighted_circumcenter_3 Construct_weighted_circumcenter_3;
  typedef typename Traits::Compare_weighted_squared_radius_3 Compare_weighted_squared_radius_3;
  typedef typename Traits::Coplanar_orientation_3 Coplanar_orientation_3;

  // Create offset array for tests
  Offset o[9] = { Offset(0, 0, 0),
                  Offset(1, 0, 0), Offset(0, 1, 0), Offset(0, 0, 1),
                  Offset(1, 1, 0), Offset(1, 0, 1), Offset(0, 1, 1),
                  Offset(1, 1, 1), Offset(10, 10, 10)};

  std::cout << "test of Construct_weighted_circumcenter_3" << std::endl;
  Construct_weighted_circumcenter_3 weighted_circumcenter = traits.construct_weighted_circumcenter_3_object();

  std::cout << "test of Compare_power_distance_3" << std::endl;
  Compare_power_distance_3 compare_power_distance = traits.compare_power_distance_3_object();

  {
    // basic version
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[0],wp[1]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[0]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[0]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[4],wp[2]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3]);
    assert (compare_power_distance(c,wp[5],wp[3]) ==  CGAL::EQUAL);
  }
  {
    // uniform (zero) offsets
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3], o[0],o[0],o[0],o[0]);
    assert (compare_power_distance(c,wp[0],wp[1], o[0],o[0],o[0]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[0], o[0],o[0],o[0]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[0], o[0],o[0],o[0]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3], o[0],o[0],o[0],o[0]);
    assert (compare_power_distance(c,wp[4],wp[2], o[0],o[0],o[0]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3], o[0],o[0],o[0],o[0]);
    assert (compare_power_distance(c,wp[5],wp[3], o[0],o[0],o[0]) ==  CGAL::EQUAL);
  }
  {
    // uniform (non-zero) offsets
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3], o[1],o[1],o[1],o[1]);
    assert (compare_power_distance(c,wp[0],wp[1], o[0],o[1],o[1]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[0], o[0],o[1],o[1]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[0], o[0],o[1],o[1]) == CGAL::SMALLER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3], o[1],o[1],o[1],o[1]);
    assert (compare_power_distance(c,wp[4],wp[2], o[0],o[1],o[1]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3], o[1],o[1],o[1],o[1]);
    assert (compare_power_distance(c,wp[5],wp[3], o[0],o[1],o[1]) ==  CGAL::EQUAL);
  }
  {
    // non-uniform offsets
    Point_3 c;
    c = weighted_circumcenter(wp[0],wp[1],wp[2],wp[3], o[1],o[2],o[3],o[2]);
    assert (compare_power_distance(c,wp[0],wp[1], o[0],o[1],o[2]) == CGAL::EQUAL);
    assert (compare_power_distance(c,wp[4],wp[0], o[0],o[2],o[2]) == CGAL::LARGER);
    assert (compare_power_distance(c,wp[6],wp[0], o[0],o[2],o[2]) == CGAL::SMALLER);
    assert (compare_power_distance(c,wp[1],wp[1], o[0],o[2],o[8]) == CGAL::SMALLER);
    assert (compare_power_distance(c,wp[1],wp[4], o[5],o[2],o[6]) == CGAL::LARGER);
    c = weighted_circumcenter(wp[4],wp[1],wp[2],wp[3], o[4],o[3],o[2],o[1]);
    assert (compare_power_distance(c,wp[4],wp[2], o[0],o[4],o[2]) == CGAL::EQUAL);
    c = weighted_circumcenter(wp[5],wp[1],wp[2],wp[3], o[2],o[1],o[3],o[1]);
    assert (compare_power_distance(c,wp[5],wp[3], o[0],o[2],o[1]) ==  CGAL::EQUAL);
  }

  std::cout << "test of Power_side_of_oriented_power_sphere_3" << std::endl;
  Power_side_of_oriented_power_sphere_3 power_test = traits.power_side_of_oriented_power_sphere_3_object();
  {
    // Triangulation_3's testsuite will check the base version in detail
    assert(power_test(wp[7],wp[8],wp[9],wp[10],wp[11]) ==
           traits.side_of_oriented_sphere_3_object()(
             wp[7].point(),wp[8].point(),wp[9].point(),wp[10].point(),wp[11].point()));

    assert(power_test(wp[8],wp[7],wp[9],wp[10],wp[13]) ==
           traits.side_of_oriented_sphere_3_object()(
             wp[8].point(),wp[7].point(),wp[9].point(),wp[10].point(),wp[13].point()));

    assert(power_test(wp[7],wp[8],wp[9],wp[13]) ==
        static_cast<CGAL::Oriented_side>(traits.coplanar_side_of_bounded_circle_3_object()(
             wp[7].point(),wp[8].point(),wp[9].point(),wp[13].point())));

    // no weight, with offsets
    assert(power_test(wp[8],wp[7],wp[9],wp[10],wp[13]) ==
           power_test(wp[8],wp[7],wp[9],wp[10],wp[13], o[2],o[2],o[2],o[2],o[2]));
    assert(power_test(wp[10],wp[9],wp[8],wp[7]) ==
           power_test(wp[10],wp[9],wp[8],wp[7], o[3],o[3],o[3],o[3]));
    assert(power_test(wp[9],wp[11],wp[13]) ==
           power_test(wp[9],wp[11],wp[13], o[4],o[4],o[4]));
    assert(power_test(wp[8],wp[14]) ==
           power_test(wp[8],wp[14],o[5],o[5]));

    assert(power_test(wp[7],wp[8],wp[9],wp[10],wp[11], o[0],o[0],o[0],o[0],o[0]) ==
           traits.side_of_oriented_sphere_3_object()(
             wp[7].point(),wp[8].point(),wp[9].point(),wp[10].point(),wp[11].point()));

    assert(power_test(wp[9],wp[13],wp[12],wp[8],wp[11],o[4],o[3],o[2],o[1],o[0]) == CGAL::ON_POSITIVE_SIDE);
    assert(power_test(wp[7],wp[8],wp[9],wp[13], o[1], o[2], o[3], o[4]) == CGAL::ON_NEGATIVE_SIDE);
    assert(power_test(wp[8],wp[11],wp[8],o[1],o[2],o[1]) == CGAL::ON_ORIENTED_BOUNDARY);

    assert(power_test(wp[0],wp[6],o[1],o[1]) == CGAL::ON_POSITIVE_SIDE);
    assert(power_test(wp[0],wp[5],o[2],o[2]) == CGAL::ON_NEGATIVE_SIDE);
    assert(power_test(wp[0],wp[0],o[3],o[3]) == CGAL::ON_ORIENTED_BOUNDARY);
  }

  std::cout << "test of Compare_weighted_squared_radius_3" << std::endl;
  Compare_weighted_squared_radius_3 compare_weighted_squared_radius = traits.compare_weighted_squared_radius_3_object();
  {
    // base, no offset
    assert(compare_weighted_squared_radius(wp[0],wp[15],wp[2],wp[3], FT(0.)) == CGAL::LARGER);
    assert(compare_weighted_squared_radius(wp[0],wp[15],wp[2], FT(1.)) == CGAL::SMALLER);
    assert(compare_weighted_squared_radius(wp[0],wp[9], FT(5.))  == CGAL::SMALLER);
    assert(compare_weighted_squared_radius(wp[15], FT(-2.)) == CGAL::EQUAL);

    // with offsets
    assert(compare_weighted_squared_radius(wp[0],wp[15],wp[2],wp[3],
                                           o[2], o[2], o[2], o[2],
                                           FT(0.))
            == compare_weighted_squared_radius(wp[0],wp[15],wp[2],wp[3], FT(0.)));
    assert(compare_weighted_squared_radius(wp[0], wp[11], wp[13],
                                           o[1], o[1], o[1],
                                           FT(1.))
            == compare_weighted_squared_radius(wp[0],wp[11],wp[13], FT(1.)));
    assert(compare_weighted_squared_radius(wp[0], wp[11], o[3], o[3], FT(1.))
            == compare_weighted_squared_radius(wp[0],wp[11], FT(1.)));
    assert(compare_weighted_squared_radius(wp[4], o[1], FT(1.))
            == compare_weighted_squared_radius(wp[4], FT(1.)));
  }

  std::cout << "test of Coplanar_orientation_3" << std::endl;
  Construct_point_3 cp = traits.construct_point_3_object();
  Coplanar_orientation_3 coplanar_orientation = traits.coplanar_orientation_3_object();
  {
    assert(coplanar_orientation(cp(wp[0]), cp(wp[1]), cp(wp[2]), cp(wp[3]),
                                o[3], o[3], o[3], o[3])
             == coplanar_orientation(cp(wp[0]), cp(wp[1]), cp(wp[2]), cp(wp[3])));
    assert(coplanar_orientation(cp(wp[0], o[3]), cp(wp[1], o[3]), cp(wp[2], o[3]), cp(wp[3], o[3]))
             == coplanar_orientation(cp(wp[0]), cp(wp[1]), cp(wp[2]), cp(wp[3])));
    assert(coplanar_orientation(cp(wp[4]), cp(wp[3]), cp(wp[2]), o[4], o[4], o[4])
             == coplanar_orientation(cp(wp[4]), cp(wp[3]), cp(wp[2])));

    assert(coplanar_orientation(cp(wp[1]), cp(wp[1]), cp(wp[1]), cp(wp[1]),
                                o[0], o[7], o[1], o[8]) == CGAL::COLLINEAR);
    assert(coplanar_orientation(cp(wp[3]), cp(wp[3]), cp(wp[3]),
                                o[0], o[7], o[8]) == CGAL::COLLINEAR);

    assert(coplanar_orientation(cp(wp[1]), cp(wp[2]), cp(wp[3]), cp(wp[4]),
                                o[6], o[7], o[4], o[1]) == CGAL::POSITIVE);
  }
}

template<class K>
void _test_cls_periodic_3_regular_triangulation_traits_3_rational ()
{
  typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>    Traits;
  typedef typename Traits::Weighted_point_3                     Weighted_point_3;
  typedef typename Traits::Point_3                              Point_3;
  typedef typename Traits::Iso_cuboid_3                         Iso_cuboid;

  Traits traits;

  std::cout << "Iso_cuboid(0,0,0,4,4,4)" << std::endl;

  traits.set_domain(Iso_cuboid(0, 0, 0, 4, 4, 4));
  Point_3 p0(0.,0.,0.);
  Point_3 p1(2.,0.,0.);
  Point_3 p2(0.,2.,0.);
  Point_3 p3(0.,0.,2.);
  Point_3 p4(2.,2.,2.);
  Point_3 p5(1.,0.,0.);
  Point_3 p6(3.,0.,0.);

  Weighted_point_3 wp[19] =
  {
    Weighted_point_3(p0, 9.),
    Weighted_point_3(p1, 9.),
    Weighted_point_3(p2, 9.),
    Weighted_point_3(p3, 9.),
    Weighted_point_3(p0, 6.),
    Weighted_point_3(p0, 3.),
    Weighted_point_3(p0, 12.),

    Weighted_point_3(p0,0.), // [7]
    Weighted_point_3(p1,0.),
    Weighted_point_3(p2,0.),
    Weighted_point_3(p3,0.),
    Weighted_point_3(p4,0.),
    Weighted_point_3(p5,0.),
    Weighted_point_3(p6,0.),
    Weighted_point_3(p0,2.),
    Weighted_point_3(p1,2.),
    Weighted_point_3(p2,2.),
    Weighted_point_3(p3,2.),
    Weighted_point_3(p4,2.)
  };

  _test_for_given_domain(traits, wp);
}

template<class K>
void _test_cls_periodic_3_regular_triangulation_traits_3_irrational ()
{
  typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>    Traits;
  typedef typename Traits::Weighted_point_3                     Weighted_point_3;
  typedef typename Traits::Point_3                              Point_3;
  typedef typename Traits::Iso_cuboid_3                         Iso_cuboid;

  Traits traits;

  std::cout << "Iso_cuboid(-pi,-sqrt(3)/2, -0.1, pi/3, pi/2, sqrt(14))" << std::endl;
  traits.set_domain(Iso_cuboid(-CGAL_PI, - CGAL::sqrt(3.) / 2., -0.1,
                                CGAL_PI / 3., CGAL_PI / 2., CGAL::sqrt(14.)));
  Point_3 p0(CGAL_PI/4., std::cos(0.7), 1./3.);
  Point_3 p1(-0.75, -1.24, CGAL::sqrt(3.));
  Point_3 p2(0.1, std::sin(1./0.002478), CGAL_PI);
  Point_3 p3(0., 0., 0.);
  Point_3 p4(1./3., 20./6., 0.);
  Point_3 p5(-2., 0., 0.);
  Point_3 p6(-4., 0., 0.);

  Weighted_point_3 wp[19] =
  {
    Weighted_point_3(p0, 9.),
    Weighted_point_3(p1, 9.),
    Weighted_point_3(p2, 9.),
    Weighted_point_3(p3, 9.),
    Weighted_point_3(p0, 6.),
    Weighted_point_3(p0, 3.),
    Weighted_point_3(p0, 12.),

    Weighted_point_3(p0,0.), // 7th
    Weighted_point_3(p1,0.),
    Weighted_point_3(p2,0.),
    Weighted_point_3(p3,0.),
    Weighted_point_3(p4,0.),
    Weighted_point_3(p5,0.),
    Weighted_point_3(p6,0.),
    Weighted_point_3(p0,2.),
    Weighted_point_3(p1,2.),
    Weighted_point_3(p2,2.),
    Weighted_point_3(p3,2.),
    Weighted_point_3(p4,2.)
  };

  _test_for_given_domain(traits, wp);
}
