// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Pedro Machado, Sebastien Loriot


#ifndef CGAL__TEST_CLS_CIRCLE_3_H
#define CGAL__TEST_CLS_CIRCLE_3_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Random.h>
#include <cassert>

// Some predicates and constructions tests related to the class Circle_3 are done here

// -------------------------- CONSTRUCTIONS

template <class K>
void _test_bounding_box_construct(const K &k)
{
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Plane_3                          Plane_3;
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Circle_3                         Circle_3;
  typedef typename K::Construct_circle_3               Construct_circle_3;
  typedef CGAL::Bbox_3                                 Bbox_3;

  Construct_circle_3 theConstruct_circle_3 = k.construct_circle_3_object();

  std::cout << "Testing the bbox of Circle_3..." << std::endl;

  Bbox_3 b;
  Circle_3 c;

  c = theConstruct_circle_3(Point_3(0,0,0), 1, Plane_3(1, 0, 0, 0));
  b = c.bbox();
  assert(b.xmin() <= 0.001);
  assert(b.xmax() >= -0.001);
  assert(b.ymin() <= -0.999);
  assert(b.ymax() >= 0.999);
  assert(b.zmin() <= -0.999);
  assert(b.zmax() >= 0.999);

  c = theConstruct_circle_3(Sphere_3(Point_3(0,0,0), 1), Plane_3(1, 0, 0, -FT(1)/FT(2)));
  b = c.bbox();
  assert(b.xmin() <= 0.501);
  assert(b.xmax() >= 0.499);
  assert(b.ymin() <= (-std::sqrt(0.5)+0.001));
  assert(b.ymax() >= (std::sqrt(0.5)-0.001));
  assert(b.zmin() <= (-std::sqrt(0.5)+0.001));
  assert(b.zmax() >= (std::sqrt(0.5)-0.001));
}

template <class K>
void _test_circle_construct(const K &k) {
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Plane_3                          Plane_3;
  typedef typename K::Circle_3                         Circle_3;
  typedef typename K::Vector_3                         Vector_3;
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Equal_3                          Equal_3;
  typedef typename K::Construct_circle_3               Construct_circle_3;
  typedef typename K::Compute_squared_distance_3       Compute_squared_distance_3;

  Equal_3 theEqual_3 = k.equal_3_object();
  Construct_circle_3 theConstruct_circle_3 = k.construct_circle_3_object();
  Compute_squared_distance_3 squared_distance = k.compute_squared_distance_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  std::cout << "Testing Construct_circle_3..." << std::endl;
  for(int i=0; i<100; i++) {
    int a,b,c,d,u,v,r;
    FT x,y,z;
    do {
      a = theRandom.get_int(random_min,random_max);
      b = theRandom.get_int(random_min,random_max);
      c = theRandom.get_int(random_min,random_max);
      d = theRandom.get_int(random_min,random_max);
    } while((a == 0) && (b == 0) && (c == 0));
    u = theRandom.get_int(random_min,random_max);
    v = theRandom.get_int(random_min,random_max);
    do {
      r = theRandom.get_int(random_min,random_max);
   } while(r <= 0);
    if(a != 0) {
      x = FT(-(b*u + c*v + d))/FT(a);
      y = FT(u);
      z = FT(v);
    } else if(b != 0) {
      x = FT(u);
      y = FT(-(a*u + c*v + d))/FT(b);
      z = FT(v);
    } else {
      x = FT(u);
      y = FT(v);
      z = FT(-(a*u + b*v + d))/FT(c);
    }
    const Plane_3 plane = Plane_3(a,b,c,d);
    const Plane_3 plane2 = Plane_3(2*a,2*b,2*c,2*d);
    const FT sqr = FT(r);
    const Point_3 p = Point_3(x,y,z);
    Circle_3 circle = theConstruct_circle_3(p,sqr,plane);

    assert(circle.supporting_plane().a() == a);
    assert(circle.supporting_plane().b() == b);
    assert(circle.supporting_plane().c() == c);
    assert(circle.supporting_plane().d() == d);
    assert(circle.center().x() == x);
    assert(circle.center().y() == y);
    assert(circle.center().z() == z);
    assert(circle.squared_radius() == sqr);
    Circle_3 circle2 = theConstruct_circle_3(p,sqr,plane2);
    Circle_3 circle3 = theConstruct_circle_3(p,sqr,Vector_3(a,b,c));
    assert(theEqual_3(circle,circle2));
    assert(theEqual_3(circle,circle3));

                Plane_3 pus(circle2);
                Sphere_3 sus(circle3);
                assert(pus == circle2.supporting_plane());
                assert(sus == circle3.diametral_sphere());
  }

  Point_3 p1, p2, p3;
  p1 = Point_3(1,0,0);
  p2 = Point_3(0,1,0);
  p3 = Point_3(0,0,1);
        Circle_3 c = theConstruct_circle_3(p1, p2, p3);
        FT r1 = squared_distance(c.center(), p1);
        FT r2 = squared_distance(c.center(), p2);
        FT r3 = squared_distance(c.center(), p3);
        assert(r1 == r2);
        assert(r2 == r3);
        assert(r3 == c.squared_radius());

        p1 = Point_3(1.3,0.2,0.1);
  p2 = Point_3(0.57,1.23,3.0);
  p3 = Point_3(9,1.2,1.3);
        c = theConstruct_circle_3(p1, p2, p3);
        r1 = squared_distance(c.center(), p1);
        r2 = squared_distance(c.center(), p2);
        r3 = squared_distance(c.center(), p3);
        assert(r1 == r2);
        assert(r2 == r3);
        assert(r3 == c.squared_radius());

  // No need to test the constructors based on intersection
  // _test_intersect_construct will test it
}

template <class K>
void _test_construct_radical_plane(const K &k) {
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Plane_3                          Plane_3;
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Circle_3                         Circle_3;
  typedef typename K::Has_on_3                         Has_on_3;
  typedef typename K::Intersect_3                      Intersect_3;
  typedef typename K::Construct_sphere_3               Construct_sphere_3;
  typedef typename K::Construct_radical_plane_3        Construct_radical_plane_3;

  Intersect_3 theIntersect_3 = k.intersect_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = k.construct_sphere_3_object();
  Construct_radical_plane_3 theConstruct_radical_plane_3 = k.construct_radical_plane_3_object();
  Has_on_3 theHas_on_3 = k.has_on_3_object();

  std::cout << "Testing radical_plane(Sphere,Sphere)..." << std::endl;
  Sphere_3 s = theConstruct_sphere_3(Point_3(0,0,0),1);
  for(int vx=-3;vx<4;vx++) {
    for(int vy=-3;vy<4;vy++) {
      for(int vz=-3;vz<4;vz++) {
        for(int vr=1;vr<6;vr++) {
          const FT x = FT(vx);
          const FT y = FT(vy);
          const FT z = FT(vz);
          const FT r = FT(vr)/FT(2);
          if(x == 0 && y == 0 && z == 0) continue;
          Sphere_3 sl_1 = theConstruct_sphere_3(Point_3(x,y,z),r*r);
          int d2 = (vx*vx + vy*vy + vz*vz);
          CGAL::Object intersection_1;
          intersection_1 = theIntersect_3(s, sl_1);
          Plane_3 p = theConstruct_radical_plane_3(s, sl_1);
          Plane_3 global_p = CGAL::radical_plane(s, sl_1);
          assert(p == global_p);
          // No intersection
          if((d2 > (r+1)*(r+1)) || (d2 < (r-1)*(r-1))) {
            assert(intersection_1.is_empty());
            CGAL::Object intersection_21, intersection_22;
            intersection_21 = theIntersect_3(p, sl_1);
            intersection_22 = theIntersect_3(p, s);
            assert(intersection_21.is_empty());
            assert(intersection_22.is_empty());
          }
          // Tangent, 1 Intersection
          else if((d2 == (r+1)*(r+1)) || (d2 == (r-1)*(r-1))) {
            Point_3 interp;
            assert(assign(interp, intersection_1));
            assert(theHas_on_3(p, interp));
          }
          // 1 Intersection Circle
          else {
            Circle_3 circle1;
            assert(assign(circle1, intersection_1));
            assert(theHas_on_3(p, circle1));
          }
        }
      }
    }
  }
}

// -------------------------- PREDICATES

template <class K>
void _test_circle_equal(const K &k) {
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Plane_3                          Plane_3;
  typedef typename K::Circle_3                         Circle_3;
  typedef typename K::Equal_3                          Equal_3;
  typedef typename K::Construct_circle_3               Construct_circle_3;

  Equal_3 theEqual_3 = k.equal_3_object();
  Construct_circle_3 theConstruct_circle_3 = k.construct_circle_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  std::cout << "Testing Equal_3 for Circle_3..." << std::endl;
  for(int i=0; i<100; i++) {
   int a,b,c,d,u,v,r;
    FT x,y,z;
    do {
      a = theRandom.get_int(random_min,random_max);
      b = theRandom.get_int(random_min,random_max);
      c = theRandom.get_int(random_min,random_max);
      d = theRandom.get_int(random_min,random_max);
    } while((a == 0) && (b == 0) && (c == 0));
    u = theRandom.get_int(random_min,random_max);
    v = theRandom.get_int(random_min,random_max);
    do {
      r = theRandom.get_int(random_min,random_max);
    } while(r <= 0);
    if(a != 0) {
      x = FT(-(b*u + c*v + d))/FT(a);
      y = FT(u);
      z = FT(v);
    } else if(b != 0) {
      x = FT(u);
      y = FT(-(a*u + c*v + d))/FT(b);
      z = FT(v);
    } else {
      x = FT(u);
      y = FT(v);
      z = FT(-(a*u + b*v + d))/FT(c);
    }
    const Plane_3 plane = Plane_3(a,b,c,d);
    const Plane_3 plane2 = Plane_3(2*a,2*b,2*c,2*d);
    const FT sqr = FT(r);
    const Point_3 p = Point_3(x,y,z);
    Circle_3 circle = theConstruct_circle_3(p,sqr,plane);
    Circle_3 circle2 = theConstruct_circle_3(p,sqr,plane2);
    assert(theEqual_3(circle,circle2));
  }
}

template <class K>
void _test_has_on_predicate(const K &k) {
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Plane_3                          Plane_3;
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Circle_3                         Circle_3;
  typedef typename K::Construct_circle_3               Construct_circle_3;
  typedef typename K::Construct_sphere_3               Construct_sphere_3;
  typedef typename K::Construct_plane_3                Construct_plane_3;
  typedef typename K::Has_on_3                         Has_on_3;

  Construct_circle_3 theConstruct_circle_3 = k.construct_circle_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = k.construct_sphere_3_object();
  Construct_plane_3 theConstruct_plane_3 = k.construct_plane_3_object();
  Has_on_3 theHas_on_3 = k.has_on_3_object();

  Sphere_3 s_1 = theConstruct_sphere_3(Point_3(0,0,0),1);
  std::cout << "Testing has_on(Sphere,Point)..." << std::endl;
  Point_3 p_1_s_1 = Point_3(1,0,0);
  Point_3 p_2_s_1 = Point_3(0,1,0);
  Point_3 p_3_s_1 = Point_3(0,0,1);
  Point_3 p_4_s_1 = Point_3(1,0,1);
  assert(theHas_on_3(s_1,p_1_s_1));
  assert(theHas_on_3(s_1,p_2_s_1));
  assert(theHas_on_3(s_1,p_3_s_1));
  assert(!theHas_on_3(s_1,p_4_s_1));

  Circle_3 c_1 = theConstruct_circle_3(Point_3(0,0,0), 1, Plane_3(1,0,0,0));
  Circle_3 c_2 = theConstruct_circle_3(Point_3(0,0,0), 1, Plane_3(1,1,1,0));
  std::cout << "Testing has_on(Circle,Point)..." << std::endl;
  Point_3 p_1_c_1 = Point_3(0,0,1);
  Point_3 p_2_c_1 = Point_3(0,1,0);
  Point_3 p_3_c_1 = Point_3(1,0,0);
  assert(theHas_on_3(c_1,p_1_c_1));
  assert(theHas_on_3(c_1,p_2_c_1));
  assert(!theHas_on_3(c_1,p_3_c_1));
  assert(!c_1.has_on(p_3_c_1));

  std::cout << "Testing has_on(Sphere,Circle)..." << std::endl;
  assert(theHas_on_3(s_1,c_1));
  assert(theHas_on_3(s_1,c_2));
  Sphere_3 s_2 = theConstruct_sphere_3(Point_3(0,0,0),2);
  assert(!theHas_on_3(s_2,c_1));
  assert(!theHas_on_3(s_2,c_2));
  Sphere_3 s_3 = theConstruct_sphere_3(Point_3(5,0,0),26);
  assert(theHas_on_3(s_3,c_1));
  assert(!theHas_on_3(s_3,c_2));

  std::cout << "Testing has_on(Plane,Circle)..." << std::endl;
  Plane_3 p_2 = theConstruct_plane_3(1,1,1,0);
  Plane_3 p_3 = theConstruct_plane_3(3,3,3,0);
  Plane_3 p_4 = theConstruct_plane_3(1,0,0,0);
  assert(theHas_on_3(p_2,c_2));
  assert(theHas_on_3(p_3,c_2));
  assert(!theHas_on_3(p_4,c_2));
  assert(!theHas_on_3(p_2,c_1));
  assert(!theHas_on_3(p_3,c_1));
  assert(theHas_on_3(p_4,c_1));
}

template <class K>
void _test_bounded_side(const K &k) {
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Construct_sphere_3               Construct_sphere_3;
  typedef typename K::Bounded_side_3                   Bounded_side_3;
  Construct_sphere_3 theConstruct_sphere_3 = k.construct_sphere_3_object();
  Bounded_side_3 theBounded_side_3 = k.bounded_side_3_object();

  std::cout << "Testing bounded_side(Sphere, Point)..." << std::endl;
  Sphere_3 s = theConstruct_sphere_3(Point_3(0,0,0),25);
  for(int x=-5; x<6; x++) {
    for(int y=-5; y<6; y++) {
      for(int z=-5; z<6; z++) {
        Point_3 p = Point_3(x,y,z);
        CGAL::Bounded_side b = theBounded_side_3(s,p);
        if((x*x + y*y + z*z) < 25) {
          assert(b == CGAL::ON_BOUNDED_SIDE);
        } else if((x*x + y*y + z*z) > 25) {
          assert(b == CGAL::ON_UNBOUNDED_SIDE);
        } else assert(b == CGAL::ON_BOUNDARY);
      }
    }
  }

  // we dont need to test bounded_side(Circle, Circular_arc_point) because
  // bounded_side(Circle, Circular_arc_point) = bounded_side(Sphere, Circular_arc_point) +
  //         has_on_3(supporting_plane, circular_arc_point) which has already been tested
  std::cout << "Testing bounded_side(Circle, Point)..." << std::endl;

  // Those predicates do not need to be tested because it is an instance of
  // bounded_side(Sphere, Circular_arc_point) and bounded_side(Circle, Circular_arc_point)
  // which have already been tested
  std::cout << "Testing has_on_bounded_side(Sphere, Point)..." << std::endl;
  std::cout << "Testing has_on_bounded_side(Circle, Point)..." << std::endl;
  std::cout << "Testing has_on_unbounded_side(Sphere, Point)..." << std::endl;
  std::cout << "Testing has_on_unbounded_side(Circle, Point)..." << std::endl;
}

// -------------------------- COMPUTATIONS

template <class K>
void _test_compute_on_circle_3(const K &k)
{
  typedef typename K::FT                                            FT;
  typedef typename K::Point_3                                       Point_3;
  typedef typename K::Plane_3                                       Plane_3;
  typedef typename K::Circle_3                                      Circle_3;
  typedef typename K::Construct_circle_3                            Construct_circle_3;
  typedef typename K::Compute_area_divided_by_pi_3                  Compute_area_divided_by_pi_3;
  typedef typename K::Compute_squared_length_divided_by_pi_square_3 Compute_squared_length_divided_by_pi_square_3;
  typedef typename K::Compute_approximate_area_3                    Compute_approximate_area_3;
  typedef typename K::Compute_approximate_squared_length_3          Compute_approximate_squared_length_3;

  std::cout << "TESTING COMPUTATIONS" << std::endl;

  Construct_circle_3 theConstruct_circle_3 = k.construct_circle_3_object();
  Compute_area_divided_by_pi_3 theCompute_area_divided_by_pi_3 = k.compute_area_divided_by_pi_3_object();
  Compute_squared_length_divided_by_pi_square_3 theCompute_squared_length_divided_by_pi_square_3 =
    k.compute_squared_length_divided_by_pi_square_3_object();
  Compute_approximate_area_3 theCompute_approximate_area_3 = k.compute_approximate_area_3_object();
  Compute_approximate_squared_length_3 theCompute_approximate_squared_length_3 =
    k.compute_approximate_squared_length_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 5;
  int random_min = -5;

  std::cout << "Testing Approximate_area of a Circle_3" << std::endl;
  std::cout << "Testing Compute_area_divided_by_pi  of a Circle_3" << std::endl;
  std::cout << "Testing Approximate_squared_length of a Circle_3" << std::endl;
  std::cout << "Testing Compute_squared_length_divided_by_pi_square of a Circle_3" << std::endl;

  for(int i=0; i<400; i++) {
    Circle_3 circle[2];
    for(int j=0; j<2; j++) {
      int a,b,c,d,u,v,r;
      FT x,y,z;
      do {
        a = theRandom.get_int(random_min,random_max);
        b = theRandom.get_int(random_min,random_max);
        c = theRandom.get_int(random_min,random_max);
        d = theRandom.get_int(random_min,random_max);
      } while((a == 0) && (b == 0) && (c == 0));
     u = theRandom.get_int(random_min,random_max);
      v = theRandom.get_int(random_min,random_max);
      do {
        r = theRandom.get_int(random_min,random_max);
      } while(r <= 0);
      if(a != 0) {
        x = FT(-(b*u + c*v + d))/FT(a);
        y = FT(u);
        z = FT(v);
     } else if(b != 0) {
        x = FT(u);
        y = FT(-(a*u + c*v + d))/FT(b);
        z = FT(v);
      } else {
        x = FT(u);
        y = FT(v);
        z = FT(-(a*u + b*v + d))/FT(c);
      }
      const Plane_3 plane = Plane_3(a,b,c,d);
      const FT sqr = FT(r*r);
      const Point_3 p = Point_3(x,y,z);
      circle[j] = theConstruct_circle_3(p,sqr,plane);
      const FT ardp = theCompute_area_divided_by_pi_3(circle[j]);
      const FT sldps = theCompute_squared_length_divided_by_pi_square_3(circle[j]);
      assert(ardp == r*r);
      assert(sldps == 4*r*r);
    }

    const double ar1 = theCompute_approximate_area_3(circle[0]);
    const FT ardp1 = theCompute_area_divided_by_pi_3(circle[0]);
    const double asl1 = theCompute_approximate_squared_length_3(circle[0]);
    const FT sldps1 = theCompute_squared_length_divided_by_pi_square_3(circle[0]);
    const double ar2 = theCompute_approximate_area_3(circle[1]);
    const FT ardp2 = theCompute_area_divided_by_pi_3(circle[1]);
    const double asl2 = theCompute_approximate_squared_length_3(circle[1]);
    const FT sldps2 = theCompute_squared_length_divided_by_pi_square_3(circle[1]);
    if(circle[0].squared_radius() > circle[1].squared_radius()) {
      assert(ar1 > ar2); assert(ardp1 > ardp2);
      assert(asl1 > asl2); assert(sldps1 > sldps2);
    } else if(circle[0].squared_radius() == circle[1].squared_radius()) {
      assert(ar1 == ar2); assert(ardp1 == ardp2);
      assert(asl1 == asl2); assert(sldps1 == sldps2);
    } else {
     assert(ar1 < ar2); assert(ardp1 < ardp2);
      assert(asl1 < asl2); assert(sldps1 < sldps2);
    }
  }

  std::cout << "All tests on computations are OK." << std::endl;
}


template <class K>
bool
_test_cls_circle_3(const K& k) {

  std::cout << "Testing class Circle_3" << std::endl;
  _test_circle_construct(k);
  _test_circle_equal(k);
  _test_has_on_predicate(k);
  _test_bounded_side(k);
  _test_construct_radical_plane(k);
  _test_bounding_box_construct(k);
  _test_compute_on_circle_3(k);
 std::cout << "done" << std::endl;
  return true;
}

#endif //CGAL__TEST_CLS_CIRCLE_3_H
