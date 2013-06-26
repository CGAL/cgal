// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#include <CGAL/Random.h>
#include <cassert>

template <class SK>
void _test_circular_arc_point_equal(SK sk) {



  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;


  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;

  typedef typename SK::Equal_3                          Equal_3;




  
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 =
    sk.construct_circular_arc_point_3_object();
  Equal_3 theEqual_3 = sk.equal_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Equal_3 for Circular_arc_point_3..." << std::endl;
  for(int i=0; i<100; i++) {
    int a = theRandom.get_int(random_min,random_max);
    int b = theRandom.get_int(random_min,random_max);
    int c = theRandom.get_int(random_min,random_max);
    Circular_arc_point_3 circ_a_point_test_2 = 
      theConstruct_circular_arc_point_3(a,b,c);
    Circular_arc_point_3 circ_a_point_test_3 = 
      theConstruct_circular_arc_point_3(Point_3(a,b,c));
    assert(theEqual_3(circ_a_point_test_2, circ_a_point_test_3));
  }
}

template <class SK>
void _test_line_arc_equal(SK sk) {



  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;



  typedef typename SK::Line_3                           Line_3;
  typedef typename SK::Line_arc_3                       Line_arc_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Equal_3                          Equal_3;


  typedef typename SK::Construct_line_3                 Construct_line_3;

  typedef typename SK::Construct_line_arc_3             Construct_line_arc_3;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;



  typedef typename AK::Polynomials_for_line_3           Polynomials_for_line_3;


  Equal_3 theEqual_3 = sk.equal_3_object();
  Get_equation theGet_equation = sk.get_equation_object();
  (void)/* Construct_circle_3 theConstruct_circle_3 = */ sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  (void)/* Construct_plane_3 theConstruct_plane_3 = */ sk.construct_plane_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = 
    sk.construct_circular_arc_point_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Equal_3 for Line_arc_3..." << std::endl;
  for(int i=0; i<100; i++) {
    int a,b,c,d,e,f;
    do {
      a = theRandom.get_int(random_min,random_max);
      b = theRandom.get_int(random_min,random_max);
      c = theRandom.get_int(random_min,random_max);
      d = theRandom.get_int(random_min,random_max);
      e = theRandom.get_int(random_min,random_max);
      f = theRandom.get_int(random_min,random_max);
    } while((a == d) && (b == e) && (c == f));
    Point_3 pb = Point_3(d,e,f);
    Point_3 pa = Point_3(a,b,c);
    Polynomials_for_line_3 p = theGet_equation(Line_3(pa,pb));
    Line_3 line = theConstruct_line_3(p);
    Circular_arc_point_3 cpa = theConstruct_circular_arc_point_3(pa);
    Circular_arc_point_3 cpb = theConstruct_circular_arc_point_3(pb);
    Line_arc_3 line_arc[6];
    line_arc[0] = theConstruct_line_arc_3(line, pa, pb);
    line_arc[1] = theConstruct_line_arc_3(line, pb, pa);
    line_arc[2] = theConstruct_line_arc_3(line, cpb, cpa);
    line_arc[3] = theConstruct_line_arc_3(line, cpa, cpb);
    line_arc[4] = theConstruct_line_arc_3(line, cpb, pa);
    line_arc[5] = theConstruct_line_arc_3(line, pa, cpb);
    for(int t1=0;t1<6;t1++) {
      for(int t2=0;t2<6;t2++) {
        assert(theEqual_3(line_arc[t1],line_arc[t2]));
      }
    }
  }
}

template <class SK>
void _test_circular_arc_equal(SK sk) {


  typedef typename SK::FT                               FT;


  typedef typename SK::Circular_arc_3                   Circular_arc_3;
  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Plane_3                          Plane_3;
  typedef typename SK::Circle_3                         Circle_3;

  typedef typename SK::Algebraic_kernel                 AK;

  typedef typename SK::Equal_3                          Equal_3;
  typedef typename SK::Construct_circle_3               Construct_circle_3;

  typedef typename SK::Construct_circular_arc_3         Construct_circular_arc_3;
  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;



  Equal_3 theEqual_3 = sk.equal_3_object();
  (void)/* Get_equation theGet_equation = */ sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Equal_3 for Circular_arc_3..." << std::endl;
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
    const Polynomials_for_circle_3 pfc = 
      std::make_pair(Polynomial_for_spheres_2_3(x,y,z,r),
                     Polynomial_1_3(a,b,c,d));
    Circle_3 circle3 = theConstruct_circle_3(pfc);
    Circle_3 circle = theConstruct_circle_3(p,sqr,plane);
    Circular_arc_3 circle_arc3 = theConstruct_circular_arc_3(circle3);
    Circular_arc_3 circle_arc = theConstruct_circular_arc_3(circle);
    Circle_3 circle2 = theConstruct_circle_3(p,sqr,plane2);
    Circular_arc_3 circle_arc2 = theConstruct_circular_arc_3(circle2);
    assert(theEqual_3(circle_arc,circle_arc2));
    assert(theEqual_3(circle_arc,circle_arc3));
    Circular_arc_3 cother = circle_arc2;
    assert(theEqual_3(circle_arc,cother));
    assert(cother.rep().is_full());
    // When the !if_full(), we use the compare of circular_arc_point, already tested
  }
}

template <class SK>
void _test_has_on_predicate(SK sk) {

  typedef typename SK::FT                               FT;
  typedef typename SK::Root_of_2                        Root_of_2;
  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Plane_3                          Plane_3;
  typedef typename SK::Sphere_3                         Sphere_3;
  typedef typename SK::Circle_3                         Circle_3;
  typedef typename SK::Line_3                           Line_3;
  typedef typename SK::Line_arc_3                       Line_arc_3;
  typedef typename SK::Circular_arc_3                   Circular_arc_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Construct_circle_3               Construct_circle_3;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;
  typedef typename SK::Construct_plane_3                Construct_plane_3;
  typedef typename SK::Construct_line_3                 Construct_line_3;
  typedef typename SK::Construct_line_arc_3             Construct_line_arc_3;
  typedef typename SK::Construct_circular_arc_3         Construct_circular_arc_3;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;
  typedef typename SK::Has_on_3                         Has_on_3;
  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;
  typedef typename AK::Polynomials_for_line_3           Polynomials_for_line_3;
  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;

  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = sk.construct_sphere_3_object();
  Construct_plane_3 theConstruct_plane_3 = sk.construct_plane_3_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = sk.construct_circular_arc_point_3_object();
  Has_on_3 theHas_on_3 = sk.has_on_3_object();

  Sphere_3 s_1 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,1));
  Point_3 p_1_s_1 = Point_3(1,0,0);
  Point_3 p_2_s_1 = Point_3(0,1,0);
  Point_3 p_3_s_1 = Point_3(0,0,1);
  Point_3 p_4_s_1 = Point_3(1,0,1);
  std::cout << "Testing has_on(Sphere,Circular_arc_point)..." << std::endl;
  Root_of_2 sqrt_1_div_3 = CGAL::make_root_of_2(FT(0),FT(1),FT(1) / FT(3));
  Root_of_2 sqrt_1_div_2 = CGAL::make_root_of_2(FT(0),FT(1),FT(FT(1) / FT(2)));
  Root_for_spheres_2_3 r_1_s_1 = Root_for_spheres_2_3(0,sqrt_1_div_2,sqrt_1_div_2);
  Root_for_spheres_2_3 r_2_s_1 = Root_for_spheres_2_3(sqrt_1_div_3,sqrt_1_div_3,sqrt_1_div_3);
  Root_for_spheres_2_3 r_3_s_1 = Root_for_spheres_2_3(sqrt_1_div_3,sqrt_1_div_3,-sqrt_1_div_3);
  Root_for_spheres_2_3 r_4_s_1 = Root_for_spheres_2_3(0,sqrt_1_div_3,-sqrt_1_div_3);
  Circular_arc_point_3 cp_1_s_1 = Circular_arc_point_3(r_1_s_1);
  Circular_arc_point_3 cp_2_s_1 = Circular_arc_point_3(r_2_s_1);
  Circular_arc_point_3 cp_3_s_1 = Circular_arc_point_3(r_3_s_1);
  Circular_arc_point_3 cp_4_s_1 = Circular_arc_point_3(r_4_s_1);
  assert(theHas_on_3(s_1,cp_1_s_1));
  assert(theHas_on_3(s_1,cp_2_s_1));
  assert(theHas_on_3(s_1,cp_3_s_1));
  assert(!theHas_on_3(s_1,cp_4_s_1));

  Plane_3 p_1 = theConstruct_plane_3(Polynomial_1_3(1,2,3,10));
  Point_3 p_1_p_1 = Point_3(-2,-1,-2);
  Point_3 p_2_p_1 = Point_3(-FT(FT(5) / FT(3)),-FT(FT(5) / FT(3)),-FT(FT(5) / FT(3)));
  Point_3 p_3_p_1 = Point_3(-10,0,0);
  Point_3 p_4_p_1 = Point_3(-2,-2,-1);
  std::cout << "Testing has_on(Plane,Circular_arc_point)..." << std::endl;
  Root_of_2 r_1_1_p_1 = CGAL::make_root_of_2(FT(0),FT(4),FT(2));
  Root_of_2 r_1_2_p_1 = CGAL::make_root_of_2(FT(-5),FT(-5),FT(2));
  Root_of_2 r_1_3_p_1 = CGAL::make_root_of_2(FT(0),FT(2),FT(2));
  Root_for_spheres_2_3 r_1_p_1 = Root_for_spheres_2_3(r_1_1_p_1,r_1_2_p_1,r_1_3_p_1);
  Root_for_spheres_2_3 r_2_p_1 = Root_for_spheres_2_3(r_1_2_p_1,r_1_2_p_1,r_1_3_p_1);
  Circular_arc_point_3 cp_1_p_1 = Circular_arc_point_3(r_1_p_1);
  Circular_arc_point_3 cp_2_p_1 = Circular_arc_point_3(r_2_p_1);
  assert(!theHas_on_3(p_1,cp_2_p_1));
  
  Line_3 l_1 = theConstruct_line_3(Polynomials_for_line_3(1,1,-1,3,1,0));
  Point_3 p_1_l_1 = Point_3(1,3,0);
  Point_3 p_2_l_1 = Point_3(0,4,-1);
  Point_3 p_3_l_1 = Point_3(2,2,1);
  Point_3 p_4_l_1 = Point_3(1,1,1);
  std::cout << "Testing has_on(Line,Circular_arc_point)..." << std::endl;
  Root_of_2 r_1_1_l_1 = CGAL::make_root_of_2(FT(1),FT(1),FT(5));
  Root_of_2 r_1_2_l_1 = CGAL::make_root_of_2(FT(3),FT(-1),FT(5));
  Root_of_2 r_1_3_l_1 = CGAL::make_root_of_2(FT(0),FT(1),FT(5));
  Root_for_spheres_2_3 r_1_l_1 = Root_for_spheres_2_3(r_1_1_l_1,r_1_2_l_1,r_1_3_l_1);
  Root_for_spheres_2_3 r_2_l_1 = Root_for_spheres_2_3(r_1_1_l_1,r_1_1_l_1,r_1_1_l_1);
  Circular_arc_point_3 cp_1_l_1 = Circular_arc_point_3(r_1_l_1);
  Circular_arc_point_3 cp_2_l_1 = Circular_arc_point_3(r_2_l_1);
  assert(!theHas_on_3(l_1,cp_2_l_1));

  const Polynomials_for_circle_3 pc1 = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(1,0,0,0));
  Circle_3 c_1 = theConstruct_circle_3(pc1);
  Point_3 p_1_c_1 = Point_3(0,0,1);
  Point_3 p_2_c_1 = Point_3(0,1,0);
  Point_3 p_3_c_1 = Point_3(1,0,0);
  std::cout << "Testing has_on(Circle,Circular_arc_point)..." << std::endl;
  const Polynomials_for_circle_3 pc2 = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(1,1,1,0));
  Circle_3 c_2 = theConstruct_circle_3(pc2);
  Root_of_2 r_1_1_c_2 = FT(FT(1) / FT(2));
  Root_of_2 r_1_2_c_2 = CGAL::make_root_of_2(-FT(FT(1) / FT(4)),-FT(FT(1) / FT(4)),FT(5));
  Root_of_2 r_1_3_c_2 = CGAL::make_root_of_2(-FT(FT(1) / FT(4)),FT(FT(1) / FT(4)),FT(5));
  Root_for_spheres_2_3 r_1_c_2 = Root_for_spheres_2_3(r_1_1_c_2,r_1_2_c_2,r_1_3_c_2);
  Root_for_spheres_2_3 r_2_c_2 = Root_for_spheres_2_3(r_1_2_c_2,r_1_2_c_2,r_1_2_c_2);
  Circular_arc_point_3 cp_1_c_2 = Circular_arc_point_3(r_1_c_2);
  Circular_arc_point_3 cp_2_c_2 = Circular_arc_point_3(r_2_c_2);
  assert(theHas_on_3(c_2,cp_1_c_2));
  assert(!theHas_on_3(c_2,cp_2_c_2));

  // Don't need to test has_on(Plane,Line). It is internal on Cartesian

  Sphere_3 s_2 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,2));
  Sphere_3 s_3 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(5,0,0,26));
  Plane_3 p_2 = theConstruct_plane_3(Polynomial_1_3(1,1,1,0));
  Plane_3 p_3 = theConstruct_plane_3(Polynomial_1_3(3,3,3,0));
  Plane_3 p_4 = theConstruct_plane_3(Polynomial_1_3(1,0,0,0));

  std::cout << "Testing has_on(Line_arc,Circular_arc_point)..." << std::endl;
  for(int vx=0;vx<3;vx++) {
    for(int vy=0;vy<3;vy++) {
      for(int vz=0;vz<3;vz++) { 
        if(vx == 0 && vy == 0 && vz == 0) continue;
        const FT a = FT(vx);
        const FT b = FT(vy);
        const FT c = FT(vz);
        Line_3 l = theConstruct_line_3(Point_3(0,0,0), Point_3(a,b,c));
        for(int t1=-2;t1<3;t1++) {
          Point_3 source = Point_3(a*t1,b*t1,c*t1);
          for(int t2=t1+1;t2<3;t2++) {
            if(t1 == t2) continue;
            Point_3 target = Point_3(a*t2,b*t2,c*t2);
            Line_arc_3 la = theConstruct_line_arc_3(l,source,target);
            for(int t3 = (t1-1); t3 <= (t2+1); t3++) {
              Point_3 pp = Point_3(a*t3,b*t3,c*t3);
              if((t3 < t1) || (t3 > t2)) {
                assert(!theHas_on_3(la,pp));
              } else {
                assert(theHas_on_3(la,pp));
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Testing has_on(Line,Line_arc)..." << std::endl;
  for(int vx=0;vx<3;vx++) {
    for(int vy=0;vy<3;vy++) {
      for(int vz=0;vz<3;vz++) { 
        if(vx == 0 && vy == 0 && vz == 0) continue;
        const FT a = FT(vx);
        const FT b = FT(vy);
        const FT c = FT(vz);
        Line_3 l = theConstruct_line_3(Point_3(0,0,0), Point_3(a,b,c));
        for(int t1=-2;t1<3;t1++) {
          Point_3 source = Point_3(a*t1,b*t1,c*t1);
          for(int t2=t1+1;t2<3;t2++) {
            Point_3 target = Point_3(a*t2,b*t2,c*t2);
            Line_arc_3 la = theConstruct_line_arc_3(l,source,target);
            for(int t3 = (t1-1); t3 <= (t2+1); t3++) {
              Point_3 pp = Point_3(a*t3,b*t3,c*t3);
              theHas_on_3(l,la);
            }
          }
        }
      }
    }
  }

  // both tests depends on other has_on or functions already tested
  std::cout << "Testing has_on(Plane, Line_arc)..." << std::endl;

  std::cout << "Testing has_on(Circular_arc, Circular_arc_point)..." << std::endl;

  // That cover all the cases, since the orientation is setted by default to be the
  // clockwise orientation for a well defined normal vector (read the comments on 
  // include/CGAL/Circular_kernel_3/Circular_arc_3.h)
  Root_for_spheres_2_3 rt[10];

  rt[0] = Root_for_spheres_2_3(0,1,0);
  rt[1] = Root_for_spheres_2_3(-FT(FT(1) / FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(3)), 0);
  rt[2] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)),0);
  rt[3] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(3)), FT(FT(1) / FT(2)), 0);

  rt[4] = Root_for_spheres_2_3(-1,0,0);
  rt[5] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)),0);
  rt[6] = Root_for_spheres_2_3(0,-1,0);
  rt[7] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)),0);
  rt[8] = Root_for_spheres_2_3(1,0,0);
  rt[9] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)),0);

  Circular_arc_point_3 cp[10]; 
  for(int i=0; i<10; i++) {
    cp[i] = theConstruct_circular_arc_point_3(rt[i]);
  }

  const Polynomials_for_circle_3 pcc_test = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(0,0,1,0));
  Circle_3 cc = theConstruct_circle_3(pcc_test);
  for(int i=0; i<10; i++) {
    for(int j=i+1; j<10; j++) {
      Circular_arc_3 ca = theConstruct_circular_arc_3(cc,cp[i],cp[j]);
      Circular_arc_3 cb = theConstruct_circular_arc_3(cc,cp[j],cp[i]);
      for(int t=0;t<10;t++) {
        if((i <= t) && (t <= j)) {
          assert(theHas_on_3(ca,cp[t]));
          if(!(t == i || t == j)) {
            assert(!theHas_on_3(cb,cp[t]));
          }
        }
        else {
          assert(!theHas_on_3(ca,cp[t]));
          assert(theHas_on_3(cb,cp[t]));
        }
      }
    }
  }

  Root_for_spheres_2_3 rt2[8];
  rt2[0] = Root_for_spheres_2_3(0,CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)),CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)));
  rt2[1] = Root_for_spheres_2_3(0,0,1);
  rt2[2] = Root_for_spheres_2_3(0,CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)),CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)));
  rt2[3] = Root_for_spheres_2_3(0,-1,0);
  rt2[4] = Root_for_spheres_2_3(0,CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)),CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)));
  rt2[5] = Root_for_spheres_2_3(0,0,-1);
  rt2[6] = Root_for_spheres_2_3(0,CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)),CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)));
  rt2[7] = Root_for_spheres_2_3(0,1,0);

  for(int i=0; i<8; i++) {
    cp[i] = theConstruct_circular_arc_point_3(rt2[i]);
  }

  const Polynomials_for_circle_3 pcc_test2 = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(1,0,0,0));
  Circle_3 cc2 = theConstruct_circle_3(pcc_test2);

  for(int i=0; i<8; i++) {
    for(int j=i+1; j<8; j++) {
      Circular_arc_3 ca = theConstruct_circular_arc_3(cc2,cp[i],cp[j]);
      Circular_arc_3 cb = theConstruct_circular_arc_3(cc2,cp[j],cp[i]);
      for(int t=0;t<8;t++) {
        if((i <= t) && (t <= j)) {
          assert(theHas_on_3(ca,cp[t]));
          if(!(t == i || t == j)) {
            assert(!theHas_on_3(cb,cp[t]));
          }
        }
        else {
          assert(!theHas_on_3(ca,cp[t]));
          assert(theHas_on_3(cb,cp[t]));
        }
      }
    }
  }

  // Those has_on are already tested before (it is practically the same as the elementar type, say Circle_3)
  std::cout << "Testing has_on(Sphere, Circular_arc)..." << std::endl;
  std::cout << "Testing has_on(Plane, Circular_arc)..." << std::endl;
  std::cout << "Testing has_on(Circle, Circular_arc)..." << std::endl;
  std::cout << "Testing has_on(Circular_arc, Circle)..." << std::endl;
}

// i,j,k must me different
bool simulate_has_on(int i, int j, int k) {
  if(i < j) return i < k && k < j;
  else return !(j < k && k < i);
}

template <class SK>
void _test_do_overlap_predicate(SK sk) {

  typedef typename SK::FT                               FT;

  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;


  typedef typename SK::Circle_3                         Circle_3;
  typedef typename SK::Line_3                           Line_3;
  typedef typename SK::Line_arc_3                       Line_arc_3;
  typedef typename SK::Circular_arc_3                   Circular_arc_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Construct_circle_3               Construct_circle_3;


  typedef typename SK::Construct_line_3                 Construct_line_3;
  typedef typename SK::Construct_line_arc_3             Construct_line_arc_3;
  typedef typename SK::Construct_circular_arc_3         Construct_circular_arc_3;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;

  typedef typename SK::Do_overlap_3                     Do_overlap_3;
  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;

  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;

  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  (void)/* Construct_plane_3 theConstruct_plane_3 = */ sk.construct_plane_3_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = sk.construct_circular_arc_point_3_object();
  (void)/* Has_on_3 theHas_on_3 = */ sk.has_on_3_object();
  Do_overlap_3 theDo_overlap_3 = sk.do_overlap_3_object();

  std::cout << "Testing do_overlap(Line_arc, Line_arc)..." << std::endl;
  for(int vx=0;vx<2;vx++) {
    for(int vy=0;vy<2;vy++) {
      for(int vz=0;vz<2;vz++) { 
        if(vx == 0 && vy == 0 && vz == 0) continue;
        const FT a = FT(vx);
        const FT b = FT(vy);
        const FT c = FT(vz);
        Line_3 l = theConstruct_line_3(Point_3(0,0,0), Point_3(a,b,c));
        for(int t1=-1;t1<2;t1++) {
          Point_3 source = Point_3(a*t1,b*t1,c*t1);
          for(int t2=t1+1;t2<3;t2++) {
            Point_3 target = Point_3(a*t2,b*t2,c*t2);
            Line_arc_3 la = theConstruct_line_arc_3(l,source,target);
            for(int t3=-1;t3<2;t3++) {
              Point_3 sourcel = Point_3(a*t3,b*t3,c*t3);
              for(int t4=t3+1;t4<3;t4++) {
                Point_3 targetl = Point_3(a*t4,b*t4,c*t4);
                Line_arc_3 lb = theConstruct_line_arc_3(l,sourcel,targetl);
                if((t1 == t3) || (t2 == t3)) {
                  assert(theDo_overlap_3(la,lb));
                } else if((t1 == t4) || (t2 == t4)) {
                  assert(theDo_overlap_3(la,lb));
                } else if((t1 < t3) && (t3 < t2 )) {
                  assert(theDo_overlap_3(la,lb));
                } else if((t3 < t1) && (t1 < t4)) {
                  assert(theDo_overlap_3(la,lb));
                } else {
                  assert(!theDo_overlap_3(la,lb));
                } 
              }
            } 
          }
        }
      }
    }
  }

  std::cout << "Testing do_overlap(Circular_arc, Circular_arc)..." << std::endl;
  Root_for_spheres_2_3 rt[8];

  rt[0] = Root_for_spheres_2_3(0,1,0);
  rt[1] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)),0);
  rt[2] = Root_for_spheres_2_3(-1,0,0);
  rt[3] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)),0);
  rt[4] = Root_for_spheres_2_3(0,-1,0);
  rt[5] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)),0);
  rt[6] = Root_for_spheres_2_3(1,0,0);
  rt[7] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)),0);

  Circular_arc_point_3 cp[8]; 
  for(int i=0; i<8; i++) {
    cp[i] = theConstruct_circular_arc_point_3(rt[i]);
  }

  const Polynomials_for_circle_3 pcc_test = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(0,0,1,0));
  Circle_3 cc = theConstruct_circle_3(pcc_test);
  for(int i=0; i<8; i++) {
    for(int j=0; j<8; j++) {
      if(i == j) continue;
      Circular_arc_3 ca = theConstruct_circular_arc_3(cc,cp[i],cp[j]);
      for(int t1=0;t1<8;t1++) {
        for(int t2=0;t2<8;t2++) {
          if(t1 == t2) continue;
          Circular_arc_3 cb = theConstruct_circular_arc_3(cc,cp[t1],cp[t2]);
          if((t1 == i) || (t2 == i) || (t1 == j) || (t2 == j)) {
            assert(theDo_overlap_3(ca,cb));
          } else if(simulate_has_on(i,j,t1) ||
                    simulate_has_on(i,j,t2) ||
                    simulate_has_on(t1,t2,i)) {
            assert(theDo_overlap_3(ca,cb));
          } else {
            assert(!theDo_overlap_3(ca,cb));
          }
        }
      }
    }
  }
}

template <class SK>
void _test_bounded_side(SK sk) {



  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Sphere_3                         Sphere_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;
  typedef typename SK::Bounded_side_3                   Bounded_side_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;




  Get_equation theGet_equation = sk.get_equation_object();
  Construct_sphere_3 theConstruct_sphere_3 = sk.construct_sphere_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 =
    sk.construct_circular_arc_point_3_object();
  Bounded_side_3 theBounded_side_3 = sk.bounded_side_3_object();

  std::cout << "Testing bounded_side(Sphere, Circular_arc_point)..." << std::endl;
  Polynomial_for_spheres_2_3 p = theGet_equation(Sphere_3(Point_3(0,0,0),25));
  Sphere_3 s = theConstruct_sphere_3(p);
  for(int x=-5; x<6; x++) {
    for(int y=-5; y<6; y++) {
      for(int z=-5; z<6; z++) {
        Circular_arc_point_3 cp = theConstruct_circular_arc_point_3(x,y,z);
        CGAL::Bounded_side b = theBounded_side_3(s,cp);
        if((x*x + y*y + z*z) < 25) {
          assert(b == CGAL::ON_BOUNDED_SIDE);
					assert(SK().has_on_bounded_side_3_object()(s,cp));
        } else if((x*x + y*y + z*z) > 25) {
          assert(b == CGAL::ON_UNBOUNDED_SIDE);
					assert(SK().has_on_unbounded_side_3_object()(s,cp));
        } else assert(b == CGAL::ON_BOUNDARY);
      }
    }
  }

  // we dont need to test bounded_side(Circle, Circular_arc_point) because
  // bounded_side(Circle, Circular_arc_point) = bounded_side(Sphere, Circular_arc_point) +
  //         has_on_3(supporting_plane, circular_arc_point) which has already been tested
  std::cout << "Testing bounded_side(Circle, Circular_arc_point)..." << std::endl;

  // Those predicates do not need to be tested because it is an instance of
  // bounded_side(Sphere, Circular_arc_point) and bounded_side(Circle, Circular_arc_point)
  // which have already been tested
  std::cout << "Testing has_on_bounded_side(Sphere, Circular_arc_point)..." << std::endl;
  std::cout << "Testing has_on_bounded_side(Circle, Circular_arc_point)..." << std::endl;
  std::cout << "Testing has_on_unbounded_side(Sphere, Circular_arc_point)..." << std::endl;
  std::cout << "Testing has_on_unbounded_side(Circle, Circular_arc_point)..." << std::endl;
}

template <class SK>
void _test_lexico_operations(SK sk) {



  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;











  (void)// Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 =
    sk.construct_circular_arc_point_3_object();

	Circular_arc_point_3 p[3];
	p[0] = Point_3(1,0,0);
  p[1] = Point_3(1,0,0); 
  p[2] = Point_3(0,1,0);
  std::cout << "Testing lexico_operations(Circular_arc_point, Circular_arc_point)..." << std::endl;
	assert(p[0] > p[2]);
	assert(p[0] >= p[1]);
	assert(p[0] <= p[1]);
	assert(p[2] < p[0]);
}

template <class SK>
void _test_compare(SK /*sk*/) {
	
	typedef CGAL::Circular_arc_point_3<SK>  Circular_arc_point_3;
  typedef CGAL::Point_3<SK>               Point_3;

	Circular_arc_point_3 p[8];
	p[0] = Point_3(1,0,0);
  p[1] = Point_3(1,0,0); 
  p[2] = Point_3(0,1,0);
  p[3] = Point_3(0,1,0);
  p[4] = Point_3(0,0,1);
  p[5] = Point_3(0,0,1);
  p[6] = Point_3(1,1,0);
  p[7] = Point_3(1,1,1);

  std::cout << "Testing compare..." << std::endl;

	assert(compare_x(p[0], p[2]) == CGAL::LARGER);
	assert(compare_x(p[2], p[4]) == CGAL::EQUAL);
	assert(compare_x(p[2], p[0]) == CGAL::SMALLER);
	
	assert(compare_y(p[2], p[0]) == CGAL::LARGER);
	assert(compare_y(p[2], p[3]) == CGAL::EQUAL);
	assert(compare_y(p[0], p[2]) == CGAL::SMALLER);
	
	assert(compare_z(p[4], p[2]) == CGAL::LARGER);
	assert(compare_z(p[4], p[5]) == CGAL::EQUAL);
	assert(compare_z(p[2], p[4]) == CGAL::SMALLER);
	
	assert(compare_xy(p[6], p[0]) == CGAL::LARGER);
	assert(compare_xy(p[6], p[6]) == CGAL::EQUAL);
	assert(compare_xy(p[0], p[6]) == CGAL::SMALLER);
	
	assert(compare_xyz(p[7], p[6]) == CGAL::LARGER);
	assert(compare_xyz(p[7], p[7]) == CGAL::EQUAL);
	assert(compare_xyz(p[6], p[7]) == CGAL::SMALLER);
}

template <class SK>
void _test_spherical_kernel_predicates(SK sk)
{
  std::cout << "TESTING PREDICATES" << std::endl;
  _test_circular_arc_point_equal(sk);
  _test_line_arc_equal(sk);
  _test_circular_arc_equal(sk);
  _test_has_on_predicate(sk);
  _test_do_overlap_predicate(sk);
  _test_bounded_side(sk);
  _test_lexico_operations(sk);
  _test_compare(sk);
  std::cout << "All tests on predicates are OK." << std::endl;
}
