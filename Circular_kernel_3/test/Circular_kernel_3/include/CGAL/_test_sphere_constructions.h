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
void _test_circular_arc_point_construct(SK sk) {

  typedef typename SK::FT                               FT;
  typedef typename SK::Root_of_2                        Root_of_2;
  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;

  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;

  typedef typename SK::Equal_3                          Equal_3;



  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;
  
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 =
    sk.construct_circular_arc_point_3_object();
  Equal_3 theEqual_3 = sk.equal_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  // test the constructor with 3 root_of_2
  Root_of_2 r1 = CGAL::make_root_of_2(FT(3),FT(-1),FT(2));
  Root_of_2 r2 = CGAL::make_root_of_2(FT(4),FT(1),FT(2));
  Root_of_2 r3 = CGAL::make_root_of_2(FT(10),FT(2),FT(2));
  Root_for_spheres_2_3 rs = Root_for_spheres_2_3(r1,r2,r3);
  Circular_arc_point_3 circ_a_point_test_1 = 
    theConstruct_circular_arc_point_3(r1,r2,r3);
  Circular_arc_point_3 circ_a_point_test_4 = 
    theConstruct_circular_arc_point_3(rs);
  std::cout << "Testing Construct_circular_arc_point_3 with 3 Root_of_2..." << std::endl;
  assert(circ_a_point_test_1.x() == r1);
  assert(circ_a_point_test_1.y() == r2);
  assert(circ_a_point_test_1.z() == r3);

  std::cout << "Testing Construct_circular_arc_point_3 with a Root_for_spheres_2_3..." << std::endl;
  assert(circ_a_point_test_4.x() == r1);
  assert(circ_a_point_test_4.y() == r2);
  assert(circ_a_point_test_4.z() == r3);

  std::cout << "Testing many random cases of Construct_circular_arc_point_3..." << std::endl;
  for(int i=0; i<100; i++) {
    int a = theRandom.get_int(random_min,random_max);
    int b = theRandom.get_int(random_min,random_max);
    int c = theRandom.get_int(random_min,random_max);
    Circular_arc_point_3 circ_a_point_test_2 = 
      theConstruct_circular_arc_point_3(a,b,c);
    Circular_arc_point_3 circ_a_point_test_3 = 
      theConstruct_circular_arc_point_3(Point_3(a,b,c));
    assert(circ_a_point_test_2.x() == a);
    assert(circ_a_point_test_2.y() == b);
    assert(circ_a_point_test_2.z() == c);
    assert(circ_a_point_test_3.x() == a);
    assert(circ_a_point_test_3.y() == b);
    assert(circ_a_point_test_3.z() == c);
    assert(theEqual_3(circ_a_point_test_2, circ_a_point_test_3));
  }

  // No need to test the constructors based on intersection
  // _test_intersect_construct will test it

}

template <class SK>
void _test_sphere_construct(SK sk) {




  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Sphere_3                         Sphere_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;




  Get_equation theGet_equation = sk.get_equation_object();
  Construct_sphere_3 theConstruct_sphere_3 = sk.construct_sphere_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Construct_sphere_3 and Get_equation(Sphere_3)......" << std::endl;
  for(int i=0; i<100; i++) {
    int a,b,c,r;
    do {
      a = theRandom.get_int(random_min,random_max);
      b = theRandom.get_int(random_min,random_max);
      c = theRandom.get_int(random_min,random_max);
      r = theRandom.get_int(random_min,random_max);
    } while(r <= 0);
    Polynomial_for_spheres_2_3 p = theGet_equation(Sphere_3(Point_3(a,b,c),r));
    assert(p.a() == a);
    assert(p.b() == b);
    assert(p.c() == c);
    assert(p.r_sq() == r);
    Sphere_3 s = theConstruct_sphere_3(p);
    assert(s.center().x() == a);
    assert(s.center().y() == b);
    assert(s.center().z() == c);
    assert(s.squared_radius() == r);
  }

}

template <class SK>
void _test_plane_construct(SK sk) {





  typedef typename SK::Plane_3                          Plane_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Construct_plane_3                Construct_plane_3;

  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;



  Get_equation theGet_equation = sk.get_equation_object();
  Construct_plane_3 theConstruct_plane_3 = sk.construct_plane_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Construct_plane_3 and Get_equation(Plane_3)..." << std::endl;
  for(int i=0; i<100; i++) {
    int a,b,c,d;
    do {
      a = theRandom.get_int(random_min,random_max);
      b = theRandom.get_int(random_min,random_max);
      c = theRandom.get_int(random_min,random_max);
      d = theRandom.get_int(random_min,random_max);
    } while((a == 0) && (b == 0) && (c == 0));
    Polynomial_1_3 p = theGet_equation(Plane_3(a,b,c,d));
    assert(p.a() == a);
    assert(p.b() == b);
    assert(p.c() == c);
    assert(p.d() == d);
    Plane_3 plane = theConstruct_plane_3(p);
    assert(plane.a() == a);
    assert(plane.b() == b);
    assert(plane.c() == c);
    assert(plane.d() == d);
  }

}

template <class SK>
void _test_line_construct(SK sk) {




  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Line_3                           Line_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Construct_line_3                 Construct_line_3;


  typedef typename AK::Polynomials_for_line_3           Polynomials_for_line_3;
  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;

  Get_equation theGet_equation = sk.get_equation_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Construct_line_3 and Get_equation(Line_3)..." << std::endl;
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
    Root_for_spheres_2_3 r1 = Root_for_spheres_2_3(pb.x(),pb.y(),pb.z());
    Root_for_spheres_2_3 r2 = Root_for_spheres_2_3(pa.x(),pa.y(),pa.z());
    assert(r1.is_on_line(p));
    assert(r2.is_on_line(p));
    Line_3 line = theConstruct_line_3(p);
    assert(line.to_vector().x() == (d-a));
    assert(line.to_vector().y() == (e-b));
    assert(line.to_vector().z() == (f-c)); 
  }

}

template <class SK>
void _test_circle_construct(SK sk) {

  typedef typename SK::FT                               FT;


  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Plane_3                          Plane_3;
  typedef typename SK::Circle_3                         Circle_3;

  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Equal_3                          Equal_3;
  typedef typename SK::Construct_circle_3               Construct_circle_3;

  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;



  Equal_3 theEqual_3 = sk.equal_3_object();
  Get_equation theGet_equation = sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Construct_circle_3 and Get_equation(Circle_3)..." << std::endl;
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
      x = FT(-(b*u + c*v + d)) / FT(a);
      y = FT(u);
      z = FT(v);
    } else if(b != 0) {
      x = FT(u);
      y = FT(-(a*u + c*v + d)) / FT(b);
      z = FT(v);
    } else {
      x = FT(u);
      y = FT(v);
      z = FT(-(a*u + b*v + d)) / FT(c);
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
    assert(circle.supporting_plane().a() == a);
    assert(circle.supporting_plane().b() == b);
    assert(circle.supporting_plane().c() == c);
    assert(circle.supporting_plane().d() == d);
    assert(circle.center().x() == x);
    assert(circle.center().y() == y);
    assert(circle.center().z() == z);
    assert(circle.squared_radius() == sqr);
    Circle_3 circle2 = theConstruct_circle_3(p,sqr,plane2);
    assert(theEqual_3(circle,circle2));
    assert(theEqual_3(circle,circle3));
    Polynomials_for_circle_3 pcircle = theGet_equation(circle);
    assert(pcircle.second.a() == a);
    assert(pcircle.second.b() == b);
    assert(pcircle.second.c() == c);
    assert(pcircle.second.d() == d);
    assert(pcircle.first.a() == x);
    assert(pcircle.first.b() == y);
    assert(pcircle.first.c() == z);
    assert(pcircle.first.r_sq() == sqr);
  }

  // No need to test the constructors based on intersection
  // _test_intersect_construct will test it
}

template <class SK>
void _test_line_arc_construct(SK sk) {



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


  Get_equation theGet_equation = sk.get_equation_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = 
    sk.construct_circular_arc_point_3_object();
  Equal_3 theEqual_3 = sk.equal_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Construct_line_arc_3..." << std::endl;
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
		Line_arc_3 line_arc[7];
    line_arc[0] = theConstruct_line_arc_3(line, pa, pb);
    line_arc[1] = theConstruct_line_arc_3(line, pb, pa);
    line_arc[2] = theConstruct_line_arc_3(line, cpb, cpa);
    line_arc[3] = theConstruct_line_arc_3(line, cpa, cpb);
    line_arc[4] = theConstruct_line_arc_3(line, cpb, pa);
    line_arc[5] = theConstruct_line_arc_3(line, pa, cpb);
    line_arc[6] = theConstruct_line_arc_3(pa, pb);
    for(int t1=0;t1<7;t1++) {
      for(int t2=0;t2<7;t2++) {
        assert(theEqual_3(line_arc[t1],line_arc[t2]));
      }
    }
  }
}

template <class SK>
void _test_circular_arc_construct(SK sk) {


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
  (void)/* Get_equation theGet_equation = */sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127; 

  std::cout << "Testing Construct_circular_arc_3..." << std::endl;
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
      x = FT(-(b*u + c*v + d)) / FT(a);
      y = FT(u);
      z = FT(v);
    } else if(b != 0) {
      x = FT(u);
      y = FT(-(a*u + c*v + d)) / FT(b);
      z = FT(v);
    } else {
      x = FT(u);
      y = FT(v);
      z = FT(-(a*u + b*v + d)) / FT(c);
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
    assert(circle_arc.supporting_plane().a() == a);
    assert(circle_arc.supporting_plane().b() == b);
    assert(circle_arc.supporting_plane().c() == c);
    assert(circle_arc.supporting_plane().d() == d);
    assert(circle_arc.center().x() == x);
    assert(circle_arc.center().y() == y);
    assert(circle_arc.center().z() == z);
    assert(circle_arc.squared_radius() == sqr);
    assert(theEqual_3(circle_arc.source(),circle_arc.target()));
    Circle_3 circle2 = theConstruct_circle_3(p,sqr,plane2);
    Circular_arc_3 circle_arc2 = theConstruct_circular_arc_3(circle2);
    assert(theEqual_3(circle_arc,circle_arc2));
    assert(theEqual_3(circle_arc,circle_arc3));
  }

  // No need to test the constructors based on intersection
  // _test_intersect_construct will test it
}

template <class SK>
void _test_intersection_construct(SK sk) {

  typedef typename SK::FT                               FT;

  typedef CGAL::Circular_arc_point_3<SK>                Circular_arc_point_3;
  typedef CGAL::Point_3<SK>                             Point_3;
  typedef CGAL::Line_3<SK>                              Line_3;
  typedef CGAL::Line_arc_3<SK>                          Line_arc_3;
  typedef CGAL::Circular_arc_3<SK>                      Circular_arc_3;
  typedef CGAL::Plane_3<SK>                             Plane_3;
  typedef CGAL::Sphere_3<SK>                            Sphere_3;
  typedef CGAL::Circle_3<SK>                            Circle_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Get_equation                     Get_equation;
  typedef typename SK::Equal_3                          Equal_3;
  typedef typename SK::Has_on_3                         Has_on_3;
  typedef typename SK::Intersect_3                      Intersect_3;
  typedef typename SK::Do_intersect_3                   Do_intersect_3;
  typedef typename SK::Construct_circle_3               Construct_circle_3;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;
  typedef typename SK::Construct_plane_3                Construct_plane_3;
  typedef typename SK::Construct_line_3                 Construct_line_3;
  typedef typename SK::Construct_line_arc_3             Construct_line_arc_3;
  typedef typename SK::Construct_circular_arc_3         Construct_circular_arc_3;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;
  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;

  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;

  Has_on_3 theHas_on_3 = sk.has_on_3_object();
  Equal_3 theEqual_3 = sk.equal_3_object();
  Intersect_3 theIntersect_3 = sk.intersect_3_object();
  Get_equation theGet_equation = sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = sk.construct_sphere_3_object();
  Construct_plane_3 theConstruct_plane_3 = sk.construct_plane_3_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = sk.construct_circular_arc_point_3_object();
	Do_intersect_3 theDo_intersect_3 = sk.do_intersect_3_object();

  Sphere_3 s = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,1));
  Sphere_3 s_t10 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(10,10,10,1));
  Plane_3 p = theConstruct_plane_3(Polynomial_1_3(1,1,1,-1));
  
  std::cout << "Testing intersection(Sphere,Line)..." << std::endl;
  for(int vx=-2;vx<3;vx++) {
    for(int vy=-2;vy<3;vy++) {
      for(int vz=-2;vz<3;vz++) {
        if(vx == -1 && vy == 0 && vz == 0) continue;
        const FT x = FT(vx);
        const FT y = FT(vy);
        const FT z = FT(vz);
        Line_3 l1 = theConstruct_line_3(Point_3(-1,0,0), Point_3(x,y,z));
        Line_3 l2 = theConstruct_line_3(Point_3(FT(-1)-FT(FT(1) / FT(1000000)),FT(0),FT(0)), Point_3(x,y,z));
        std::vector< CGAL::Object > intersection_1, intersection_2;
        theIntersect_3(s, l1, std::back_inserter(intersection_1));
        theIntersect_3(s, l2, std::back_inserter(intersection_2));

        // tangent case for line 1
        // non-intersection for line 2
        if(vx == -1) {
          assert(theDo_intersect_3(s, l1));
          assert(intersection_1.size() == 1);
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          assert(assign(the_pair1, intersection_1[0]));
          assert(theHas_on_3(s,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));

          assert(intersection_2.size() == 0);
          assert(!theDo_intersect_3(s, l2));
        }

        // 2 intersections for line 1
        // 2 intersections for line 2
        else {
          assert(intersection_1.size() == 2);
          assert(theDo_intersect_3(s, l1));
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          std::pair<Circular_arc_point_3, unsigned > the_pair2;
          assert(assign(the_pair1, intersection_1[0]));
          assert(assign(the_pair2, intersection_1[1]));
          assert(theHas_on_3(s,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));
          assert(theHas_on_3(s,the_pair2.first));
          assert(theHas_on_3(l1,the_pair2.first));

          assert(intersection_2.size() == 2);
          assert(theDo_intersect_3(s, l2));
          std::pair<Circular_arc_point_3, unsigned > the_pair3;
          std::pair<Circular_arc_point_3, unsigned > the_pair4;
          assert(assign(the_pair3, intersection_2[0]));
          assert(assign(the_pair4, intersection_2[1]));
          assert(theHas_on_3(s,the_pair3.first));
          assert(theHas_on_3(l2,the_pair3.first));
          assert(theHas_on_3(s,the_pair4.first));
          assert(theHas_on_3(l2,the_pair4.first));
        }
      }
    }
  }

  std::cout << "Testing global version intersection(Sphere,Line)..." << std::endl;
  for(int vx=-2;vx<3;vx++) {
    for(int vy=-2;vy<3;vy++) {
      for(int vz=-2;vz<3;vz++) {
        if(vx == -1 && vy == 0 && vz == 0) continue;
        const FT x = FT(vx);
        const FT y = FT(vy);
        const FT z = FT(vz);
        Line_3 l1 = theConstruct_line_3(Point_3(-1,0,0), Point_3(x,y,z));
        Line_3 l2 = theConstruct_line_3(Point_3(FT(-1)-FT(FT(1) / FT(1000000)),FT(0),FT(0)), Point_3(x,y,z));
        std::vector< CGAL::Object > intersection_1, intersection_2;
        intersection(s, l1, std::back_inserter(intersection_1));
        intersection(s, l2, std::back_inserter(intersection_2));

        // tangent case for line 1
        // non-intersection for line 2
        if(vx == -1) {
					assert(do_intersect(s, l1));
          assert(intersection_1.size() == 1);
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          assert(assign(the_pair1, intersection_1[0]));
          assert(theHas_on_3(s,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));

          assert(intersection_2.size() == 0);
          assert(!do_intersect(s, l2));
        }

        // 2 intersections for line 1
        // 2 intersections for line 2
        else {
          assert(intersection_1.size() == 2);
          assert(do_intersect(s, l1));
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          std::pair<Circular_arc_point_3, unsigned > the_pair2;
          assert(assign(the_pair1, intersection_1[0]));
          assert(assign(the_pair2, intersection_1[1]));
          assert(theHas_on_3(s,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));
          assert(theHas_on_3(s,the_pair2.first));
          assert(theHas_on_3(l1,the_pair2.first));

          assert(intersection_2.size() == 2);
          assert(do_intersect(s, l2));
          std::pair<Circular_arc_point_3, unsigned > the_pair3;
          std::pair<Circular_arc_point_3, unsigned > the_pair4;
          assert(assign(the_pair3, intersection_2[0]));
          assert(assign(the_pair4, intersection_2[1]));
          assert(theHas_on_3(s,the_pair3.first));
          assert(theHas_on_3(l2,the_pair3.first));
          assert(theHas_on_3(s,the_pair4.first));
          assert(theHas_on_3(l2,the_pair4.first));
        }
      }
    }
  }

  std::cout << "Testing intersection(Sphere,Sphere,Sphere)..." << std::endl;
  Sphere_3 s1 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,16));
  Sphere_3 s2 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(4,0,0,16));
  Sphere_3 s3 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(8,0,0,16));
  for(int vx=-1;vx<5;vx++) {
    for(int vy=-1;vy<5;vy++) {
      for(int vz=-1;vz<5;vz++) {
        for(int vr=1;vr<6;vr++) {
          const FT x = 4*FT(vx);
          const FT y = 4*FT(vy);
          const FT z = 4*FT(vz);
          const FT r = 4*FT(vr) / FT(2);
          Sphere_3 sl = theConstruct_sphere_3(
            Polynomial_for_spheres_2_3(x,y,z,r*r));
          std::vector< CGAL::Object > intersection_1;
          std::vector< CGAL::Object > intersection_2;
          theIntersect_3(s1, s2, sl, std::back_inserter(intersection_1));
          theIntersect_3(s1, s3, sl, std::back_inserter(intersection_2));
          if(intersection_1.size() == 1) {
	          assert(theDo_intersect_3(s1, s2, sl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_1[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s2,circle));
              assert(theHas_on_3(sl,circle));
            }
            if(assign(cap,intersection_1[0])) {
              // This case must never happen
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s2,cap.first));
              assert(theHas_on_3(sl,cap.first));
            }
          }
          if(intersection_1.size() == 2) {
	          assert(theDo_intersect_3(s1, s2, sl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_1[0]));
            assert(assign(cap2,intersection_1[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s2,cap1.first));
            assert(theHas_on_3(sl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s2,cap2.first));
            assert(theHas_on_3(sl,cap2.first));
          }

          if(intersection_2.size() == 1) {
	          assert(theDo_intersect_3(s1, s3, sl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_2[0])) {
              // This case must never happen
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s3,circle));
              assert(theHas_on_3(sl,circle));
            }
            if(assign(cap,intersection_2[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s3,cap.first));
              assert(theHas_on_3(sl,cap.first));
            }
          }
          if(intersection_2.size() == 2) {
	          assert(theDo_intersect_3(s1, s3, sl));
            // This case must never happen
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_2[0]));
            assert(assign(cap2,intersection_2[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s3,cap1.first));
            assert(theHas_on_3(sl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s3,cap2.first));
            assert(theHas_on_3(sl,cap2.first));
          }
        }
      }
    }
  }

  std::cout << "Testing global version intersection(Sphere,Sphere,Sphere)..." << std::endl;
  for(int vx=-1;vx<5;vx++) {
    for(int vy=-1;vy<5;vy++) {
      for(int vz=-1;vz<5;vz++) {
        for(int vr=1;vr<6;vr++) {
          const FT x = 4*FT(vx);
          const FT y = 4*FT(vy);
          const FT z = 4*FT(vz);
          const FT r = 4*FT(vr) / FT(2);
          Sphere_3 sl = theConstruct_sphere_3(
            Polynomial_for_spheres_2_3(x,y,z,r*r));
          std::vector< CGAL::Object > intersection_1;
          std::vector< CGAL::Object > intersection_2;
          intersection(s1, s2, sl, std::back_inserter(intersection_1));
          intersection(s1, s3, sl, std::back_inserter(intersection_2));
          if(intersection_1.size() == 1) {
            assert(CGAL::do_intersect(s1, s2, sl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_1[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s2,circle));
              assert(theHas_on_3(sl,circle));
            }
            if(assign(cap,intersection_1[0])) {
              // This case must never happen
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s2,cap.first));
              assert(theHas_on_3(sl,cap.first));
            }
          }
          if(intersection_1.size() == 2) {
            assert(CGAL::do_intersect(s1, s2, sl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_1[0]));
            assert(assign(cap2,intersection_1[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s2,cap1.first));
            assert(theHas_on_3(sl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s2,cap2.first));
            assert(theHas_on_3(sl,cap2.first));
          }

          if(intersection_2.size() == 1) {
            assert(CGAL::do_intersect(s1, s3, sl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_2[0])) {
              // This case must never happen
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s3,circle));
              assert(theHas_on_3(sl,circle));
            }
            if(assign(cap,intersection_2[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s3,cap.first));
              assert(theHas_on_3(sl,cap.first));
            }
          }
          if(intersection_2.size() == 2) {
            assert(CGAL::do_intersect(s1, s3, sl));
            // This case must never happen
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_2[0]));
            assert(assign(cap2,intersection_2[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s3,cap1.first));
            assert(theHas_on_3(sl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s3,cap2.first));
            assert(theHas_on_3(sl,cap2.first));
          }
        }
      }
    }
  }


  std::cout << "Testing intersection(Sphere,Sphere,Plane)..." << std::endl;
  for(int va=-1;va<3;va++) {
    for(int vb=-1;vb<3;vb++) {
      for(int vc=-1;vc<3;vc++) {
        for(int vd=-8;vd<9;vd++) {
          const FT a = FT(va);
          const FT b = FT(vb);
          const FT c = FT(vc);
          const FT d = -FT(vd) / FT(2);
          if(a == 0 && b == 0 && c == 0) continue;
          Plane_3 pl = theConstruct_plane_3(
            Polynomial_1_3(a,b,c,d));
          std::vector< CGAL::Object > intersection_1;
          std::vector< CGAL::Object > intersection_2;
          theIntersect_3(s1, s2, pl, std::back_inserter(intersection_1));
          theIntersect_3(s1, s3, pl, std::back_inserter(intersection_2));
          if(intersection_1.size() == 1) {
	          assert(theDo_intersect_3(s1, s2, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_1[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s2,circle));
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_1[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s2,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_1.size() == 2) {
	          assert(theDo_intersect_3(s1, s2, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_1[0]));
            assert(assign(cap2,intersection_1[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s2,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s2,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
          if(intersection_2.size() == 1) {
	          assert(theDo_intersect_3(s1, s3, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_2[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s3,circle));
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_2[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s3,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_2.size() == 2) {
	          assert(theDo_intersect_3(s1, s3, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_2[0]));
            assert(assign(cap2,intersection_2[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s3,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s3,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
        }
      }
    }
  }

  std::cout << "Testing global version intersection(Sphere,Sphere,Plane)..." << std::endl;
  for(int va=-1;va<3;va++) {
    for(int vb=-1;vb<3;vb++) {
      for(int vc=-1;vc<3;vc++) {
        for(int vd=-8;vd<9;vd++) {
          const FT a = FT(va);
          const FT b = FT(vb);
          const FT c = FT(vc);
          const FT d = -FT(vd) / FT(2);
          if(a == 0 && b == 0 && c == 0) continue;
          Plane_3 pl = theConstruct_plane_3(
            Polynomial_1_3(a,b,c,d));
          std::vector< CGAL::Object > intersection_1;
          std::vector< CGAL::Object > intersection_2;
          intersection(s1, s2, pl, std::back_inserter(intersection_1));
          intersection(s1, s3, pl, std::back_inserter(intersection_2));
          if(intersection_1.size() == 1) {
	          assert(CGAL::do_intersect(s1, s2, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_1[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s2,circle));
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_1[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s2,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_1.size() == 2) {
	          assert(CGAL::do_intersect(s1, s2, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_1[0]));
            assert(assign(cap2,intersection_1[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s2,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s2,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
          if(intersection_2.size() == 1) {
	          assert(CGAL::do_intersect(s1, s3, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_2[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(s3,circle));
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_2[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(s3,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_2.size() == 2) {
	          assert(CGAL::do_intersect(s1, s3, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_2[0]));
            assert(assign(cap2,intersection_2[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(s3,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(s3,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
        }
      }
    }
  }

  // This test becomes almost Line vs Sphere, is not that interesting to test
  Plane_3 p1 = theConstruct_plane_3(Polynomial_1_3(1,1,1,0));
  Plane_3 p2 = theConstruct_plane_3(Polynomial_1_3(1,0,0,-4));
  std::cout << "Testing intersection(Plane,Plane,Sphere)..." << std::endl;
  for(int va=-1;va<3;va++) {
    for(int vb=-1;vb<3;vb++) {
      for(int vc=-1;vc<3;vc++) {
        for(int vd=-8;vd<9;vd++) {
          const FT a = FT(va);
          const FT b = FT(vb);
          const FT c = FT(vc);
          const FT d = -FT(vd) / FT(2);
          if(a == 0 && b == 0 && c == 0) continue;
          Plane_3 pl = theConstruct_plane_3(
            Polynomial_1_3(a,b,c,d));
          std::vector< CGAL::Object > intersection_1;
          std::vector< CGAL::Object > intersection_2;
          theIntersect_3(s1, p1, pl, std::back_inserter(intersection_1));
          theIntersect_3(s1, p2, pl, std::back_inserter(intersection_2));
          if(intersection_1.size() == 1) {
	          assert(theDo_intersect_3(s1, p1, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_1[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(p1,circle)); 
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_1[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(p1,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_1.size() == 2) {
	          assert(theDo_intersect_3(s1, p1, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_1[0]));
            assert(assign(cap2,intersection_1[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(p1,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(p1,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
          if(intersection_2.size() == 1) {
	          assert(theDo_intersect_3(s1, p2, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_2[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(p2,circle));
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_2[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(p2,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_2.size() == 2) {
	          assert(theDo_intersect_3(s1, p2, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_2[0]));
            assert(assign(cap2,intersection_2[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(p2,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(p2,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
        }
      }
    }
  }

  std::cout << "Testing global version intersection(Plane,Plane,Sphere)..." << std::endl;
  for(int va=-1;va<3;va++) {
    for(int vb=-1;vb<3;vb++) {
      for(int vc=-1;vc<3;vc++) {
        for(int vd=-8;vd<9;vd++) {
          const FT a = FT(va);
          const FT b = FT(vb);
          const FT c = FT(vc);
          const FT d = -FT(vd) / FT(2);
          if(a == 0 && b == 0 && c == 0) continue;
          Plane_3 pl = theConstruct_plane_3(
            Polynomial_1_3(a,b,c,d));
          std::vector< CGAL::Object > intersection_1;
          std::vector< CGAL::Object > intersection_2;
          intersection(s1, p1, pl, std::back_inserter(intersection_1));
          intersection(s1, p2, pl, std::back_inserter(intersection_2));
          if(intersection_1.size() == 1) {
	          assert(CGAL::do_intersect(s1, p1, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_1[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(p1,circle)); 
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_1[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(p1,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_1.size() == 2) {
	          assert(CGAL::do_intersect(s1, p1, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_1[0]));
            assert(assign(cap2,intersection_1[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(p1,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(p1,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
          if(intersection_2.size() == 1) {
	          assert(CGAL::do_intersect(s1, p2, pl));
            Circle_3 circle;
            std::pair<Circular_arc_point_3, unsigned> cap;
            if(assign(circle,intersection_2[0])) {
              assert(theHas_on_3(s1,circle));
              assert(theHas_on_3(p2,circle));
              assert(theHas_on_3(pl,circle));
            }
            if(assign(cap,intersection_2[0])) {
              assert(theHas_on_3(s1,cap.first));
              assert(theHas_on_3(p2,cap.first));
              assert(theHas_on_3(pl,cap.first));
            }
          }
          if(intersection_2.size() == 2) {
	          assert(CGAL::do_intersect(s1, p2, pl));
            std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
            assert(assign(cap1,intersection_2[0]));
            assert(assign(cap2,intersection_2[1]));
            assert(theHas_on_3(s1,cap1.first));
            assert(theHas_on_3(p2,cap1.first));
            assert(theHas_on_3(pl,cap1.first));
            assert(theHas_on_3(s1,cap2.first));
            assert(theHas_on_3(p2,cap2.first));
            assert(theHas_on_3(pl,cap2.first));
          }
        }
      }
    }
  }

  std::cout << "Testing intersection(Circle,Sphere)..." << std::endl;
  std::cout << "Testing intersection(Circle,Plane)..." << std::endl;
  // both tests above are equivalent respectively to 
  // intersection(Sphere,Sphere,Plane) and
  // intersection(Plane,Plane,Sphere)

  std::cout << "Testing intersection(Circle,Circle)..." << std::endl;
  Polynomial_for_spheres_2_3 es1 = Polynomial_for_spheres_2_3(0,0,0,1);
  Polynomial_for_spheres_2_3 es2 = Polynomial_for_spheres_2_3(1,0,0,1);
  Polynomial_for_spheres_2_3 es3 = Polynomial_for_spheres_2_3(2,0,0,1);
  for(int va=-5;va<6;va++) {
    const FT a = -FT(va) / FT(10);
    const FT b = 1;
    const FT c = 0;
    const FT d = 0;
    Polynomial_1_3 pol = Polynomial_1_3(a,b,c,d);
    Circle_3 c1 = theConstruct_circle_3(std::make_pair(es1, pol));
    Circle_3 c2 = theConstruct_circle_3(std::make_pair(es2, pol));
    Circle_3 c3 = theConstruct_circle_3(std::make_pair(es3, pol));

    std::vector< CGAL::Object > intersection_1;
    theIntersect_3(c1, c1, std::back_inserter(intersection_1));
    assert(intersection_1.size() == 1);
    assert(theDo_intersect_3(c1, c1));
    Circle_3 circle;
    assert(assign(circle,intersection_1[0]));
    assert(circle == c1);

    std::vector< CGAL::Object > intersection_2;
    theIntersect_3(c1, c2, std::back_inserter(intersection_2));
    assert(intersection_2.size() == 2);
    assert(theDo_intersect_3(c1, c2));
    std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
    assert(assign(cap1,intersection_2[0]));
    assert(assign(cap2,intersection_2[1]));
    assert(theHas_on_3(c1,cap1.first));
    assert(theHas_on_3(c2,cap1.first));
    assert(theHas_on_3(c1,cap2.first));
    assert(theHas_on_3(c2,cap2.first));

    std::vector< CGAL::Object > intersection_3;
    theIntersect_3(c1, c3, std::back_inserter(intersection_3));
    if(a != 0) {
	    assert(!theDo_intersect_3(c1, c3));
      assert(intersection_3.size() == 0);
    } else {
	    assert(theDo_intersect_3(c1, c3));
      assert(intersection_3.size() == 1);
      std::pair<Circular_arc_point_3, unsigned> cap;
      assert(assign(cap,intersection_3[0]));
      assert(theHas_on_3(c1,cap.first));
      assert(theHas_on_3(c3,cap.first)); 
    }

    for(int vb=-5;vb<6;vb++) {
      const FT al = 1;
      const FT bl = 0;
      const FT cl = -FT(vb) / FT(10);
      const FT dl = 0;
      Polynomial_1_3 pol2 = Polynomial_1_3(al,bl,cl,dl);
      Circle_3 c4 = theConstruct_circle_3(std::make_pair(es1, pol2));
      std::vector< CGAL::Object > intersection_4;
      theIntersect_3(c1, c4, std::back_inserter(intersection_4));
      assert(theDo_intersect_3(c1, c4));
      assert(intersection_4.size() == 2);
      assert(assign(cap1,intersection_4[0]));
      assert(assign(cap2,intersection_4[1]));
      assert(theHas_on_3(c1,cap1.first));
      assert(theHas_on_3(c4,cap1.first));
      assert(theHas_on_3(c1,cap2.first));
      assert(theHas_on_3(c4,cap2.first));
    }

    const FT a_c = FT(va) / FT(10);
    const FT b_c = 1;
    const FT c_c = 0;
    const FT d_c = -FT(va) / FT(10);
    Polynomial_1_3 pol2 = Polynomial_1_3(a_c,b_c,c_c,d_c);
    Circle_3 c5 = theConstruct_circle_3(std::make_pair(es2, pol2));
    std::vector< CGAL::Object > intersection_5;
    theIntersect_3(c1, c5, std::back_inserter(intersection_5));
    assert(theDo_intersect_3(c1, c5));
    assert(intersection_5.size() == 2);
    assert(assign(cap1,intersection_5[0]));
    assert(assign(cap2,intersection_5[1]));
    assert(theHas_on_3(c1,cap1.first));
    assert(theHas_on_3(c5,cap1.first));
    assert(theHas_on_3(c1,cap2.first));
    assert(theHas_on_3(c5,cap2.first));

  }

  std::cout << "Testing global version intersection(Circle,Circle)..." << std::endl;
  for(int va=-5;va<6;va++) {
    const FT a = -FT(va) / FT(10);
    const FT b = 1;
    const FT c = 0;
    const FT d = 0;
    Polynomial_1_3 pol = Polynomial_1_3(a,b,c,d);
    Circle_3 c1 = theConstruct_circle_3(std::make_pair(es1, pol));
    Circle_3 c2 = theConstruct_circle_3(std::make_pair(es2, pol));
    Circle_3 c3 = theConstruct_circle_3(std::make_pair(es3, pol));

    std::vector< CGAL::Object > intersection_1;
    intersection(c1, c1, std::back_inserter(intersection_1));
    assert(intersection_1.size() == 1);
    assert(CGAL::do_intersect(c1, c1));
    Circle_3 circle;
    assert(assign(circle,intersection_1[0]));
    assert(circle == c1);

    std::vector< CGAL::Object > intersection_2;
    intersection(c1, c2, std::back_inserter(intersection_2));
    assert(intersection_2.size() == 2);
    assert(CGAL::do_intersect(c1, c2));
    std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
    assert(assign(cap1,intersection_2[0]));
    assert(assign(cap2,intersection_2[1]));
    assert(theHas_on_3(c1,cap1.first));
    assert(theHas_on_3(c2,cap1.first));
    assert(theHas_on_3(c1,cap2.first));
    assert(theHas_on_3(c2,cap2.first));

    std::vector< CGAL::Object > intersection_3;
    intersection(c1, c3, std::back_inserter(intersection_3));
    if(a != 0) {
	    assert(!do_intersect(c1, c3));
      assert(intersection_3.size() == 0);
    } else {
	    assert(CGAL::do_intersect(c1, c3));
      assert(intersection_3.size() == 1);
      std::pair<Circular_arc_point_3, unsigned> cap;
      assert(assign(cap,intersection_3[0]));
      assert(theHas_on_3(c1,cap.first));
      assert(theHas_on_3(c3,cap.first)); 
    }

    for(int vb=-5;vb<6;vb++) {
      const FT al = 1;
      const FT bl = 0;
      const FT cl = -FT(vb) / FT(10);
      const FT dl = 0;
      Polynomial_1_3 pol2 = Polynomial_1_3(al,bl,cl,dl);
      Circle_3 c4 = theConstruct_circle_3(std::make_pair(es1, pol2));
      std::vector< CGAL::Object > intersection_4;
      intersection(c1, c4, std::back_inserter(intersection_4));
      assert(CGAL::do_intersect(c1, c4));
      assert(intersection_4.size() == 2);
      assert(assign(cap1,intersection_4[0]));
      assert(assign(cap2,intersection_4[1]));
      assert(theHas_on_3(c1,cap1.first));
      assert(theHas_on_3(c4,cap1.first));
      assert(theHas_on_3(c1,cap2.first));
      assert(theHas_on_3(c4,cap2.first));
    }

    const FT a_c = FT(va) / FT(10);
    const FT b_c = 1;
    const FT c_c = 0;
    const FT d_c = -FT(va) / FT(10);
    Polynomial_1_3 pol2 = Polynomial_1_3(a_c,b_c,c_c,d_c);
    Circle_3 c5 = theConstruct_circle_3(std::make_pair(es2, pol2));
    std::vector< CGAL::Object > intersection_5;
    intersection(c1, c5, std::back_inserter(intersection_5));
    assert(CGAL::do_intersect(c1, c5));
    assert(intersection_5.size() == 2);
    assert(assign(cap1,intersection_5[0]));
    assert(assign(cap2,intersection_5[1]));
    assert(theHas_on_3(c1,cap1.first));
    assert(theHas_on_3(c5,cap1.first));
    assert(theHas_on_3(c1,cap2.first));
    assert(theHas_on_3(c5,cap2.first));

  }

  std::cout << "Testing intersection(Circle,Line)..." << std::endl;
  Polynomial_for_spheres_2_3 pol_s = Polynomial_for_spheres_2_3(0,0,0,1);
  for(int vx=-2;vx<1;vx++) {
    for(int vy=-2;vy<1;vy++) {
      for(int vz=-2;vz<1;vz++) {
        if(vx == -1 && vy == 0 && vz == 0) continue;
        const FT x = FT(vx);
        const FT y = FT(vy);
        const FT z = FT(vz);
        Line_3 l1 = theConstruct_line_3(Point_3(-1,0,0), Point_3(x,y,z));
        Plane_3 pl1 = theConstruct_plane_3(Point_3(-1,0,0), Point_3(x,y,z), 
                                           Point_3(3,4,5));
        Line_3 l2 = theConstruct_line_3(Point_3(-FT(1)-FT(FT(1) / FT(1000000)),FT(0),FT(0)), 
                                        Point_3(x,y,z));
        Plane_3 pl2 = theConstruct_plane_3(Point_3(-FT(1)-FT(FT(1) / FT(1000000)),FT(0),FT(0)), 
                                           Point_3(x,y,z), Point_3(3,4,5));
        Polynomial_1_3 pol_pl1 = theGet_equation(pl1);
        Polynomial_1_3 pol_pl2 = theGet_equation(pl2);
        Circle_3 c1 = theConstruct_circle_3(std::make_pair(pol_s, pol_pl1));
        Circle_3 c2 = theConstruct_circle_3(std::make_pair(pol_s, pol_pl2));

        std::vector< CGAL::Object > intersection_1, intersection_2;
        theIntersect_3(c1, l1, std::back_inserter(intersection_1));
        theIntersect_3(c2, l2, std::back_inserter(intersection_2));

        // tangent case for line 1
        // non-intersection for line 2
        if(vx == -1) {
	        assert(theDo_intersect_3(c1, l1));
          assert(intersection_1.size() == 1);
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          assert(assign(the_pair1, intersection_1[0]));
          assert(theHas_on_3(c1,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));

          assert(intersection_2.size() == 0);
          assert(!theDo_intersect_3(c2, l2));
        }

        // 2 intersections for line 1
        // 2 intersections for line 2
        else {
	        assert(theDo_intersect_3(c1, l1));
          assert(intersection_1.size() == 2);
          assert(theDo_intersect_3(c1, l1));
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          std::pair<Circular_arc_point_3, unsigned > the_pair2;
          assert(assign(the_pair1, intersection_1[0]));
          assert(assign(the_pair2, intersection_1[1]));
          assert(theHas_on_3(c1,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));
          assert(theHas_on_3(c1,the_pair2.first));
          assert(theHas_on_3(l1,the_pair2.first));

          assert(intersection_2.size() == 2);
          assert(theDo_intersect_3(c2, l2));
          std::pair<Circular_arc_point_3, unsigned > the_pair3;
          std::pair<Circular_arc_point_3, unsigned > the_pair4;
          assert(assign(the_pair3, intersection_2[0]));
          assert(assign(the_pair4, intersection_2[1]));
          assert(theHas_on_3(c2,the_pair3.first));
          assert(theHas_on_3(l2,the_pair3.first));
          assert(theHas_on_3(c2,the_pair4.first));
          assert(theHas_on_3(l2,the_pair4.first));
        }
      }
    }
  }

  std::cout << "Testing global version intersection(Circle,Line)..." << std::endl;
  for(int vx=-2;vx<1;vx++) {
    for(int vy=-2;vy<1;vy++) {
      for(int vz=-2;vz<1;vz++) {
        if(vx == -1 && vy == 0 && vz == 0) continue;
        const FT x = FT(vx);
        const FT y = FT(vy);
        const FT z = FT(vz);
        Line_3 l1 = theConstruct_line_3(Point_3(-1,0,0), Point_3(x,y,z));
        Plane_3 pl1 = theConstruct_plane_3(Point_3(-1,0,0), Point_3(x,y,z), 
                                           Point_3(3,4,5));
        Line_3 l2 = theConstruct_line_3(Point_3(-FT(1)-FT(FT(1) / FT(1000000)),FT(0),FT(0)), 
                                        Point_3(x,y,z));
        Plane_3 pl2 = theConstruct_plane_3(Point_3(-FT(1)-FT(FT(1) / FT(1000000)),FT(0),FT(0)), 
                                           Point_3(x,y,z), Point_3(3,4,5));
        Polynomial_1_3 pol_pl1 = theGet_equation(pl1);
        Polynomial_1_3 pol_pl2 = theGet_equation(pl2);
        Circle_3 c1 = theConstruct_circle_3(std::make_pair(pol_s, pol_pl1));
        Circle_3 c2 = theConstruct_circle_3(std::make_pair(pol_s, pol_pl2));

        std::vector< CGAL::Object > intersection_1, intersection_2;
        intersection(c1, l1, std::back_inserter(intersection_1));
        intersection(c2, l2, std::back_inserter(intersection_2));

        // tangent case for line 1
        // non-intersection for line 2
        if(vx == -1) {
	        assert(CGAL::do_intersect(c1, l1));
          assert(intersection_1.size() == 1);
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          assert(assign(the_pair1, intersection_1[0]));
          assert(theHas_on_3(c1,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));

          assert(intersection_2.size() == 0);
          assert(!theDo_intersect_3(c2, l2));
        }

        // 2 intersections for line 1
        // 2 intersections for line 2
        else {
	        assert(CGAL::do_intersect(c1, l1));
          assert(intersection_1.size() == 2);
          assert(CGAL::do_intersect(c1, l1));
          std::pair<Circular_arc_point_3, unsigned > the_pair1;
          std::pair<Circular_arc_point_3, unsigned > the_pair2;
          assert(assign(the_pair1, intersection_1[0]));
          assert(assign(the_pair2, intersection_1[1]));
          assert(theHas_on_3(c1,the_pair1.first));
          assert(theHas_on_3(l1,the_pair1.first));
          assert(theHas_on_3(c1,the_pair2.first));
          assert(theHas_on_3(l1,the_pair2.first));

          assert(intersection_2.size() == 2);
          assert(CGAL::do_intersect(c2, l2));
          std::pair<Circular_arc_point_3, unsigned > the_pair3;
          std::pair<Circular_arc_point_3, unsigned > the_pair4;
          assert(assign(the_pair3, intersection_2[0]));
          assert(assign(the_pair4, intersection_2[1]));
          assert(theHas_on_3(c2,the_pair3.first));
          assert(theHas_on_3(l2,the_pair3.first));
          assert(theHas_on_3(c2,the_pair4.first));
          assert(theHas_on_3(l2,the_pair4.first));
        }
      }
    }
  }

  std::cout << "Testing intersection(Line_arc, Line_arc)..." << std::endl;
  // Testing the case where it overlaps, or do not intersect
	int vx=1, vy=1, vz=1;
	
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
          std::vector< CGAL::Object > intersection_1;
          theIntersect_3(la, lb, std::back_inserter(intersection_1));
          if(t1 == t3) {
            Line_arc_3 line_a;
            assert(theDo_intersect_3(la, lb));
            assert(intersection_1.size() == 1);
            assert(assign(line_a, intersection_1[0]));
            if(t2 <= t4) {
              assert(theEqual_3(line_a, la));
            } else {
              assert(theEqual_3(line_a, lb));
            } 
          } else if(t2 == t4) {
            Line_arc_3 line_a;
            assert(theDo_intersect_3(la, lb));
            assert(intersection_1.size() == 1);
            assert(assign(line_a, intersection_1[0]));
            if(t1 > t3) {
              assert(theEqual_3(line_a, la));
            } else {
              assert(theEqual_3(line_a, lb));
            }
          }
          else if((t1 == t4) || (t2 == t3)) {
	          assert(theDo_intersect_3(la, lb));
            std::pair<Circular_arc_point_3, unsigned> pair;
            assert(intersection_1.size() == 1);
            assert(assign(pair, intersection_1[0]));
            if(t2 == t3) assert(theEqual_3(pair.first,target)); 
            if(t1 == t4) assert(theEqual_3(pair.first,source));
          } else if((t1 < t3) && (t3 < t2 )) {
            Line_arc_3 line_a;
            assert(theDo_intersect_3(la, lb));
            assert(intersection_1.size() == 1);
            assert(assign(line_a, intersection_1[0]));
            if(t2 < t4) { 
              Line_arc_3 line_b;
              line_b = theConstruct_line_arc_3(l,sourcel,target);
              assert(theEqual_3(line_a, line_b));
            } else {
              assert(theEqual_3(line_a, lb));
            } 
          } else if((t3 < t1) && (t1 < t4)) {
            Line_arc_3 line_a;
            assert(intersection_1.size() == 1);
            assert(theDo_intersect_3(la, lb));
            assert(assign(line_a, intersection_1[0]));
            if(t4 < t2) {
              Line_arc_3 line_b;
              line_b = theConstruct_line_arc_3(l,source,targetl);
              assert(theEqual_3(line_a, line_b));
            } else {
              assert(theEqual_3(line_a, la));
            }
          } else {
	          assert(!theDo_intersect_3(la, lb));
            assert(intersection_1.size() == 0);
          } 
        }
      } 
    }
  }

  // testing for intersections of Line Arcs without the same line support
  for(int u1=0;u1<=1;u1++) {
    for(int v1=0;v1<=1;v1++) {
      Point_3 p1 = Point_3(3*u1 + 2*v1,
                           u1 + v1 + 3,
                           5 + 2*u1 + v1);
      for(int u2=-1;u2<=0;u2++) {
        for(int v2=-1;v2<=0;v2++) {
          Point_3 p2 = Point_3(3*u2 + 2*v2,
                                 u2 + v2 + 3,
                               5 + 2*u2 + v2);
          if(theEqual_3(p1,p2)) continue;
          for(int u3=0;u3<=1;u3++) {
            for(int v3=0;v3<=1;v3++) {
              Point_3 p3 = Point_3(3*u3 + 2*v3,
                                   u3 + v3 + 3,
                                   5 + 2*u3 + v3);
              if(theEqual_3(p1,p3) || 
                 theEqual_3(p2,p3)) continue;
              for(int u4=-1;u4<=0;u4++) {
                for(int v4=-1;v4<=0;v4++) {
                  Point_3 p4 = Point_3(3*u4 + 2*v4,
                                       u4 + v4 + 3,
                                       5 + 2*u4 + v4);
                  if(theEqual_3(p1,p4) || 
                     theEqual_3(p2,p4) ||
                     theEqual_3(p3,p4)) continue;
                  Line_3 line[6];
                  line[0] = Line_3(p1, p2);
                  line[1] = Line_3(p1, p3);
                  line[2] = Line_3(p1, p4);
                  line[3] = Line_3(p2, p3);
                  line[4] = Line_3(p2, p4);
                  line[5] = Line_3(p3, p4);
                  bool n_equal = true;
                  for(int i=0; i<6; i++) {
                    for(int j=i+1;j<6;j++) { 
                      n_equal = (n_equal && (!theEqual_3(line[i],line[j])));
                      if(!n_equal) break;
                    }
                    if(!n_equal) break;
                  }
                  if(!n_equal) continue;
                  Line_arc_3 l[6];
                  l[0] = theConstruct_line_arc_3(line[0],p1,p2);
                  l[1] = theConstruct_line_arc_3(line[1],p1,p3);
                  l[2] = theConstruct_line_arc_3(line[2],p1,p4);
                  l[3] = theConstruct_line_arc_3(line[3],p2,p3);
                  l[4] = theConstruct_line_arc_3(line[4],p2,p4);
                  l[5] = theConstruct_line_arc_3(line[5],p3,p4);
                  int n_of_diagonals = 0;
                  for(int i=0;i<6;i++) {
                    int n_of_intersection = 0;
                    for(int j=0;j<6;j++) {
                      if(i == j) continue;
                      std::vector< CGAL::Object > intersection_1;
                      theIntersect_3(l[i], l[j], std::back_inserter(intersection_1));
                      assert((intersection_1.size() == 0) || 
                             (intersection_1.size() == 1));
                      if(intersection_1.size() == 1) {
	                      assert(theDo_intersect_3(l[i], l[j]));
                        n_of_intersection++;
                        std::pair<Circular_arc_point_3, unsigned> pair;
                        assert(assign(pair, intersection_1[0]));
                        assert(theHas_on_3(l[i],pair.first));
                        assert(theHas_on_3(l[j],pair.first));
                        int n_of_edges = 2;
                        for(int k=0; k<6; k++) {
                          if(k == i) continue;
                          if(k == j) continue;
                          if(theHas_on_3(l[k],pair.first)) n_of_edges++;
                        }
                        assert((n_of_edges == 2) || (n_of_edges == 3));
                      } 
                    }
                    assert((n_of_intersection == 4) ||
                           (n_of_intersection == 5));
                    if(n_of_intersection == 5) n_of_diagonals++;
                  }
                  // Non-Convex || Convex Polygon
                  assert((n_of_diagonals == 0) || (n_of_diagonals == 2)); 
                }
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Testing global version intersection(Line_arc, Line_arc)..." << std::endl;
  // Testing the case where it overlaps, or do not intersect

  l = theConstruct_line_3(Point_3(0,0,0), Point_3(a,b,c));
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
          std::vector< CGAL::Object > intersection_1;
          intersection(la, lb, std::back_inserter(intersection_1));
          if(t1 == t3) {
            Line_arc_3 line_a;
            assert(CGAL::do_intersect(la, lb));
            assert(intersection_1.size() == 1);
            assert(assign(line_a, intersection_1[0]));
            if(t2 <= t4) {
              assert(theEqual_3(line_a, la));
            } else {
              assert(theEqual_3(line_a, lb));
            } 
          } else if(t2 == t4) {
            Line_arc_3 line_a;
            assert(CGAL::do_intersect(la, lb));
            assert(intersection_1.size() == 1);
            assert(assign(line_a, intersection_1[0]));
            if(t1 > t3) {
              assert(theEqual_3(line_a, la));
            } else {
              assert(theEqual_3(line_a, lb));
            }
          }
          else if((t1 == t4) || (t2 == t3)) {
            assert(CGAL::do_intersect(la, lb));
            std::pair<Circular_arc_point_3, unsigned> pair;
            assert(intersection_1.size() == 1);
            assert(assign(pair, intersection_1[0]));
            if(t2 == t3) assert(theEqual_3(pair.first,target)); 
            if(t1 == t4) assert(theEqual_3(pair.first,source));
          } else if((t1 < t3) && (t3 < t2 )) {
            Line_arc_3 line_a;
            assert(CGAL::do_intersect(la, lb));
            assert(intersection_1.size() == 1);
            assert(assign(line_a, intersection_1[0]));
            if(t2 < t4) { 
              Line_arc_3 line_b;
              line_b = theConstruct_line_arc_3(l,sourcel,target);
              assert(theEqual_3(line_a, line_b));
            } else {
              assert(theEqual_3(line_a, lb));
            } 
          } else if((t3 < t1) && (t1 < t4)) {
            Line_arc_3 line_a;
            assert(intersection_1.size() == 1);
            assert(CGAL::do_intersect(la, lb));
            assert(assign(line_a, intersection_1[0]));
            if(t4 < t2) {
              Line_arc_3 line_b;
              line_b = theConstruct_line_arc_3(l,source,targetl);
              assert(theEqual_3(line_a, line_b));
            } else {
              assert(theEqual_3(line_a, la));
            }
          } else {
            assert(!do_intersect(la, lb));
            assert(intersection_1.size() == 0);
          } 
        }
      } 
    }
  }

  // testing for intersections of Line Arcs without the same line support
  for(int u1=0;u1<=1;u1++) {
    for(int v1=0;v1<=1;v1++) {
      Point_3 p1 = Point_3(3*u1 + 2*v1,
                           u1 + v1 + 3,
                           5 + 2*u1 + v1);
      for(int u2=-1;u2<=0;u2++) {
        for(int v2=-1;v2<=0;v2++) {
          Point_3 p2 = Point_3(3*u2 + 2*v2,
                                 u2 + v2 + 3,
                               5 + 2*u2 + v2);
          if(theEqual_3(p1,p2)) continue;
          for(int u3=0;u3<=1;u3++) {
            for(int v3=0;v3<=1;v3++) {
              Point_3 p3 = Point_3(3*u3 + 2*v3,
                                   u3 + v3 + 3,
                                   5 + 2*u3 + v3);
              if(theEqual_3(p1,p3) || 
                 theEqual_3(p2,p3)) continue;
              for(int u4=-1;u4<=0;u4++) {
                for(int v4=-1;v4<=0;v4++) {
                  Point_3 p4 = Point_3(3*u4 + 2*v4,
                                       u4 + v4 + 3,
                                       5 + 2*u4 + v4);
                  if(theEqual_3(p1,p4) || 
                     theEqual_3(p2,p4) ||
                     theEqual_3(p3,p4)) continue;
                  Line_3 line[6];
                  line[0] = Line_3(p1, p2);
                  line[1] = Line_3(p1, p3);
                  line[2] = Line_3(p1, p4);
                  line[3] = Line_3(p2, p3);
                  line[4] = Line_3(p2, p4);
                  line[5] = Line_3(p3, p4);
                  bool n_equal = true;
                  for(int i=0; i<6; i++) {
                    for(int j=i+1;j<6;j++) { 
                      n_equal = (n_equal && (!theEqual_3(line[i],line[j])));
                      if(!n_equal) break;
                    }
                    if(!n_equal) break;
                  }
                  if(!n_equal) continue;
                  Line_arc_3 l[6];
                  l[0] = theConstruct_line_arc_3(line[0],p1,p2);
                  l[1] = theConstruct_line_arc_3(line[1],p1,p3);
                  l[2] = theConstruct_line_arc_3(line[2],p1,p4);
                  l[3] = theConstruct_line_arc_3(line[3],p2,p3);
                  l[4] = theConstruct_line_arc_3(line[4],p2,p4);
                  l[5] = theConstruct_line_arc_3(line[5],p3,p4);
                  int n_of_diagonals = 0;
                  for(int i=0;i<6;i++) {
                    int n_of_intersection = 0;
                    for(int j=0;j<6;j++) {
                      if(i == j) continue;
                      std::vector< CGAL::Object > intersection_1;
                      intersection(l[i], l[j], std::back_inserter(intersection_1));
                      assert((intersection_1.size() == 0) || 
                             (intersection_1.size() == 1));
                      if(intersection_1.size() == 1) {
	                      assert(CGAL::do_intersect(l[i], l[j]));
                        n_of_intersection++;
                        std::pair<Circular_arc_point_3, unsigned> pair;
                        assert(assign(pair, intersection_1[0]));
                        assert(theHas_on_3(l[i],pair.first));
                        assert(theHas_on_3(l[j],pair.first));
                        int n_of_edges = 2;
                        for(int k=0; k<6; k++) {
                          if(k == i) continue;
                          if(k == j) continue;
                          if(theHas_on_3(l[k],pair.first)) n_of_edges++;
                        }
                        assert((n_of_edges == 2) || (n_of_edges == 3));
                      } 
                    }
                    assert((n_of_intersection == 4) ||
                           (n_of_intersection == 5));
                    if(n_of_intersection == 5) n_of_diagonals++;
                  }
                  // Non-Convex || Convex Polygon
                  assert((n_of_diagonals == 0) || (n_of_diagonals == 2)); 
                }
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Testing intersection(Line_arc,Line)..." << std::endl;
  std::cout << "Testing intersection(Line_arc,Sphere)..." << std::endl;
  std::cout << "Testing intersection(Line_arc,Plane)..." << std::endl;
  std::cout << "Testing intersection(Line_arc,Circle)..." << std::endl;
  // all those tests above are equivalent to
  // an intersection of line with a has_on (already tested)

   std::cout << "Testing intersection(Circular_arc, Circular_arc)..." << std::endl;
  // Testing the case where it overlaps, or do not intersect
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
          std::vector< CGAL::Object > intersection_1;
          theIntersect_3(ca, cb, std::back_inserter(intersection_1));
          Circular_arc_3 cres, cres2;
          std::pair< Circular_arc_point_3, unsigned > cp1, cp2;
          if(t1 == i) {
            if(t2 == j) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(ca,cres));
            } else if(simulate_has_on(i,j,t2)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(ca, cres));
            } 
          } else if(t2 == j) {
            if(simulate_has_on(i,j,t1)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(ca, cres));
            }
          } else if(t1 == j) {
            if(t2 == i) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 2);
              assert(assign(cp1,intersection_1[0]));
              assert(assign(cp2,intersection_1[1]));
              assert(theEqual_3(cp1.first, cp[i]) || 
                     theEqual_3(cp1.first, cp[j]));
              assert(theEqual_3(cp2.first, cp[i]) || 
                     theEqual_3(cp2.first, cp[j]));
              assert(!theEqual_3(cp1.first, cp2.first));
            } else if(simulate_has_on(t1,i,t2)) {
              assert(intersection_1.size() == 1);
              assert(theDo_intersect_3(ca, cb));
              assert(assign(cp1,intersection_1[0]));
              assert(theEqual_3(cp1.first, cp[t1]));
            } else {
              assert(intersection_1.size() == 2);
              assert(theDo_intersect_3(ca, cb));
              if(assign(cp1,intersection_1[0]) && assign(cres,intersection_1[1]));
              else if(assign(cres,intersection_1[0]) && assign(cp1,intersection_1[1]));
              else assert(false);
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert(theEqual_3(conf, cres));
              assert(theEqual_3(cp1.first, cp[t1]));
            } 
          } else if(t2 == i) {
            if(t1 == j) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 2);
              assert(assign(cp1,intersection_1[0]));
              assert(assign(cp2,intersection_1[1]));
              assert(theEqual_3(cp1.first, cp[i]) || 
                     theEqual_3(cp1.first, cp[j]));
              assert(theEqual_3(cp2.first, cp[i]) || 
                     theEqual_3(cp2.first, cp[j]));
              assert(!theEqual_3(cp1.first, cp2.first));
            } else if(simulate_has_on(j,i,t1)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cp1,intersection_1[0]));
              assert(theEqual_3(cp1.first, cp[i]));
            } else {
              assert(intersection_1.size() == 2);
              assert(theDo_intersect_3(ca, cb));
              if(assign(cp1,intersection_1[0]) && assign(cres,intersection_1[1]));
              else if(assign(cres,intersection_1[0]) && assign(cp1,intersection_1[1]));
              else assert(false);
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              assert(theEqual_3(conf, cres));
              assert(theEqual_3(cp1.first, cp[i]));
            }
          } else if(simulate_has_on(i,j,t1)) {
            if(simulate_has_on(t1,j,t2)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else if(simulate_has_on(t2,j,i)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              assert(theEqual_3(cres, conf));
            } else {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 2);
              assert(assign(cres,intersection_1[0]));
              assert(assign(cres2,intersection_1[1]));
              Circular_arc_3 conf1 = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              Circular_arc_3 conf2 = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert((theEqual_3(cres, conf1) && theEqual_3(cres2, conf2)) ||
                     (theEqual_3(cres2, conf1) && theEqual_3(cres, conf2)));
            }
          } else if(simulate_has_on(i,j,t2)) {
            if(simulate_has_on(i,t2,t1)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else if(simulate_has_on(j,i,t1)) {
	            assert(theDo_intersect_3(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert(theEqual_3(cres, conf));
            } else {
              // This case sould never happen, because it already happen before
              assert(intersection_1.size() == 2);
              assert(theDo_intersect_3(ca, cb));
              assert(assign(cres,intersection_1[0]));
              assert(assign(cres2,intersection_1[1]));
              Circular_arc_3 conf1 = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              Circular_arc_3 conf2 = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert((theEqual_3(cres, conf1) && theEqual_3(cres2, conf2)) ||
                     (theEqual_3(cres2, conf1) && theEqual_3(cres, conf2)));
            }
          } else if(simulate_has_on(t1,t2,i)) {
            // the case whether (i,j) contains (t1,t2) is handled before
            assert(intersection_1.size() == 1);
            assert(theDo_intersect_3(ca, cb));
            assert(assign(cres,intersection_1[0]));
            assert(theEqual_3(ca, cres));
          } else {
            assert(intersection_1.size() == 0);
            assert(!theDo_intersect_3(ca, cb));
          }
        }
      }
    }
  }

   std::cout << "Testing global version intersection(Circular_arc, Circular_arc)..." << std::endl;
  // Testing the case where it overlaps, or do not intersect
  for(int i=0; i<8; i++) {
    for(int j=0; j<8; j++) {
      if(i == j) continue;
      Circular_arc_3 ca = theConstruct_circular_arc_3(cc,cp[i],cp[j]);
      for(int t1=0;t1<8;t1++) {
        for(int t2=0;t2<8;t2++) {
          if(t1 == t2) continue;
          Circular_arc_3 cb = theConstruct_circular_arc_3(cc,cp[t1],cp[t2]);
          std::vector< CGAL::Object > intersection_1;
          intersection(ca, cb, std::back_inserter(intersection_1));
          Circular_arc_3 cres, cres2;
          std::pair< Circular_arc_point_3, unsigned > cp1, cp2;
          if(t1 == i) {
            if(t2 == j) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(ca,cres));
            } else if(simulate_has_on(i,j,t2)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(ca, cres));
            } 
          } else if(t2 == j) {
            if(simulate_has_on(i,j,t1)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(ca, cres));
            }
          } else if(t1 == j) {
            if(t2 == i) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 2);
              assert(assign(cp1,intersection_1[0]));
              assert(assign(cp2,intersection_1[1]));
              assert(theEqual_3(cp1.first, cp[i]) || 
                     theEqual_3(cp1.first, cp[j]));
              assert(theEqual_3(cp2.first, cp[i]) || 
                     theEqual_3(cp2.first, cp[j]));
              assert(!theEqual_3(cp1.first, cp2.first));
            } else if(simulate_has_on(t1,i,t2)) {
              assert(intersection_1.size() == 1);
              assert(CGAL::do_intersect(ca, cb));
              assert(assign(cp1,intersection_1[0]));
              assert(theEqual_3(cp1.first, cp[t1]));
            } else {
              assert(intersection_1.size() == 2);
              assert(CGAL::do_intersect(ca, cb));
              if(assign(cp1,intersection_1[0]) && assign(cres,intersection_1[1]));
              else if(assign(cres,intersection_1[0]) && assign(cp1,intersection_1[1]));
              else assert(false);
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert(theEqual_3(conf, cres));
              assert(theEqual_3(cp1.first, cp[t1]));
            } 
          } else if(t2 == i) {
            if(t1 == j) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 2);
              assert(assign(cp1,intersection_1[0]));
              assert(assign(cp2,intersection_1[1]));
              assert(theEqual_3(cp1.first, cp[i]) || 
                     theEqual_3(cp1.first, cp[j]));
              assert(theEqual_3(cp2.first, cp[i]) || 
                     theEqual_3(cp2.first, cp[j]));
              assert(!theEqual_3(cp1.first, cp2.first));
            } else if(simulate_has_on(j,i,t1)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cp1,intersection_1[0]));
              assert(theEqual_3(cp1.first, cp[i]));
            } else {
              assert(intersection_1.size() == 2);
              assert(CGAL::do_intersect(ca, cb));
              if(assign(cp1,intersection_1[0]) && assign(cres,intersection_1[1]));
              else if(assign(cres,intersection_1[0]) && assign(cp1,intersection_1[1]));
              else assert(false);
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              assert(theEqual_3(conf, cres));
              assert(theEqual_3(cp1.first, cp[i]));
            }
          } else if(simulate_has_on(i,j,t1)) {
            if(simulate_has_on(t1,j,t2)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else if(simulate_has_on(t2,j,i)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              assert(theEqual_3(cres, conf));
            } else {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 2);
              assert(assign(cres,intersection_1[0]));
              assert(assign(cres2,intersection_1[1]));
              Circular_arc_3 conf1 = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              Circular_arc_3 conf2 = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert((theEqual_3(cres, conf1) && theEqual_3(cres2, conf2)) ||
                     (theEqual_3(cres2, conf1) && theEqual_3(cres, conf2)));
            }
          } else if(simulate_has_on(i,j,t2)) {
            if(simulate_has_on(i,t2,t1)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              assert(theEqual_3(cb, cres));
            } else if(simulate_has_on(j,i,t1)) {
	            assert(CGAL::do_intersect(ca, cb));
              assert(intersection_1.size() == 1);
              assert(assign(cres,intersection_1[0]));
              Circular_arc_3 conf = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert(theEqual_3(cres, conf));
            } else {
              // This case sould never happen, because it already happen before
              assert(intersection_1.size() == 2);
              assert(CGAL::do_intersect(ca, cb));
              assert(assign(cres,intersection_1[0]));
              assert(assign(cres2,intersection_1[1]));
              Circular_arc_3 conf1 = theConstruct_circular_arc_3(cc,cp[t1],cp[j]);
              Circular_arc_3 conf2 = theConstruct_circular_arc_3(cc,cp[i],cp[t2]);
              assert((theEqual_3(cres, conf1) && theEqual_3(cres2, conf2)) ||
                     (theEqual_3(cres2, conf1) && theEqual_3(cres, conf2)));
            }
          } else if(simulate_has_on(t1,t2,i)) {
            // the case whether (i,j) contains (t1,t2) is handled before
            assert(intersection_1.size() == 1);
            assert(CGAL::do_intersect(ca, cb));
            assert(assign(cres,intersection_1[0]));
            assert(theEqual_3(ca, cres));
          } else {
            assert(intersection_1.size() == 0);
            assert(!do_intersect(ca, cb));
          }
        }
      }
    }
  }

  // The case whether the supporting circle is not the same is
  // the intersection_3(Circle_3, Circle_3) + has_on_3(Circle_3. Circule_arc_3)
  // already tested

  std::cout << "Testing intersection(Circular_arc,Line)..." << std::endl;
  std::cout << "Testing intersection(Circular_arc,Sphere)..." << std::endl;
  std::cout << "Testing intersection(Circular_arc,Plane)..." << std::endl;
  std::cout << "Testing intersection(Circular_arc,Circle)..." << std::endl;
  std::cout << "Testing intersection(Circular_arc,Line_arc)..." << std::endl;
  // all those tests above are equivalent to
  // an intersection of circle with a has_on (already tested)

}

template <class SK>
void _test_bbox(const typename SK::Circular_arc_point_3 &p)
{
  typedef typename SK::FT                               FT;
  typedef typename SK::Construct_bbox_3                 Construct_bbox_3;
  typedef CGAL::Bbox_3                                  Bbox_3;
  Construct_bbox_3 theConstruct_bbox_3 = SK().construct_bbox_3_object();
  Bbox_3 b = theConstruct_bbox_3(p);
  assert(FT(b.xmin()) <= p.x());
  assert(p.x() <= FT(b.xmax()));
  
  if ( FT(b.ymin()) > p.y() ){
    std::cout << "hello" << std::endl;
    assert( typename CGAL::Real_embeddable_traits< typename SK::Root_of_2 >::To_interval()(p.y()).first==b.ymin() );
  }
  
  assert(FT(b.ymin()) <= p.y());
  
  
  assert(p.y() <= FT(b.ymax()));
  assert(FT(b.zmin()) <= p.z());
  assert(p.z() <= FT(b.zmax()));
}

template <class SK>
void _test_bbox(const typename SK::Circle_3 &c)
{

  typedef typename SK::FT                               FT;

  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;









  typedef typename SK::Construct_bbox_3                 Construct_bbox_3;










  typedef CGAL::Bbox_3                                  Bbox_3;

  (void)/* Has_on_3 theHas_on_3 = */ SK().has_on_3_object();
  (void)/* Equal_3 theEqual_3 = */ SK().equal_3_object();
  (void)/* Intersect_3 theIntersect_3 = */ SK().intersect_3_object();
  (void)/* Get_equation theGet_equation = */ SK().get_equation_object();
  (void)/* Construct_circle_3 theConstruct_circle_3 = */ SK().construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ SK().construct_sphere_3_object();
  (void)/* Construct_plane_3 theConstruct_plane_3 = */ SK().construct_plane_3_object();
  (void)/* Construct_line_3 theConstruct_line_3 = */ SK().construct_line_3_object();
  Construct_bbox_3 theConstruct_bbox_3 = SK().construct_bbox_3_object();

  Bbox_3 b = theConstruct_bbox_3(c);
  const FT x1 = FT(b.xmin());
  const FT x2 = FT(b.xmax());
  const FT y1 = FT(b.ymin());
  const FT y2 = FT(b.ymax());
  const FT z1 = FT(b.zmin());
  const FT z2 = FT(b.zmax());

  assert(x1 <= x2);
  assert(y1 <= y2);
  assert(z1 <= z2);

	Circular_arc_point_3 ex1, ex2, ey1, ey2, ez1, ez2;

  if((c.supporting_plane().b() == 0) &&
		 (c.supporting_plane().c() == 0)) ex2 = ex1 = c.center();
	else {
    ex1 = x_extremal_point(c, true);
    ex2 = x_extremal_point(c, false);
  }

  if((c.supporting_plane().a() == 0) &&
		(c.supporting_plane().c() == 0)) ey2 = ey1 = c.center();
	else {
    ey1 = y_extremal_point(c, true);
    ey2 = y_extremal_point(c, false);
  }

  if((c.supporting_plane().a() == 0) &&
		(c.supporting_plane().b() == 0)) ez2 = ez1 = c.center();
  else {
	  ez1 = z_extremal_point(c, true);
    ez2 = z_extremal_point(c, false);
  }

	if(ex1.x() > ex2.x()) std::swap(ex1, ex2);
	if(ey1.y() > ey2.y()) std::swap(ey1, ey2);
	if(ez1.z() > ez2.z()) std::swap(ez1, ez2);
	
  assert(x1 <= ex1.x());
  assert(x2 >= ex2.x());
  assert(y1 <= ey1.y());
	assert(y2 >= ey2.y());
  assert(z1 <= ez1.z());
	assert(z2 >= ez2.z());	
}

template <class SK>
void _test_bounding_box_construct(SK sk)
{

  typedef typename SK::FT                               FT;

  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Line_3                           Line_3;

  typedef typename SK::Sphere_3                         Sphere_3;
  typedef typename SK::Circle_3                         Circle_3;
  typedef typename SK::Line_arc_3                       Line_arc_3;
  typedef typename SK::Algebraic_kernel                 AK;



  typedef typename SK::Construct_bbox_3                 Construct_bbox_3;
  typedef typename SK::Intersect_3                      Intersect_3;
  typedef typename SK::Construct_circle_3               Construct_circle_3;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;

  typedef typename SK::Construct_line_3                 Construct_line_3;
  typedef typename SK::Construct_line_arc_3             Construct_line_arc_3;

  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;


  typedef CGAL::Bbox_3                                  Bbox_3;

  (void)/* Has_on_3 theHas_on_3 = */ sk.has_on_3_object();
  (void)/* Equal_3 theEqual_3 = */ sk.equal_3_object();
  Intersect_3 theIntersect_3 = sk.intersect_3_object();
  (void)/* Get_equation theGet_equation = */ sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = sk.construct_sphere_3_object();
  (void)/* Construct_plane_3 theConstruct_plane_3 = */ sk.construct_plane_3_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_bbox_3 theConstruct_bbox_3 = sk.construct_bbox_3_object();

  std::cout << "Testing the bbox of Circular_arc_point_3..." << std::endl;
  Polynomial_for_spheres_2_3 es1 = Polynomial_for_spheres_2_3(0,0,0,1);
  Polynomial_for_spheres_2_3 es2 = Polynomial_for_spheres_2_3(1,0,0,1);
  Polynomial_for_spheres_2_3 es3 = Polynomial_for_spheres_2_3(2,0,0,1);
  for(int va=-5;va<6;va++) {
    const FT a = -FT(va) / FT(10);
    const FT b = 1;
    const FT c = 0;
    const FT d = 0;
    Polynomial_1_3 pol = Polynomial_1_3(a,b,c,d);
    Circle_3 c1 = theConstruct_circle_3(std::make_pair(es1, pol));
    Circle_3 c2 = theConstruct_circle_3(std::make_pair(es2, pol));

    std::vector< CGAL::Object > intersection_2;
    theIntersect_3(c1, c2, std::back_inserter(intersection_2));
    std::pair<Circular_arc_point_3, unsigned> cap1, cap2;
    assign(cap1,intersection_2[0]);
    assign(cap2,intersection_2[1]);
    _test_bbox<SK>(cap1.first);
    _test_bbox<SK>(cap2.first);

    for(int vb=-5;vb<6;vb++) {
      const FT al = 1;
      const FT bl = 0;
      const FT cl = -FT(vb) / FT(10);
      const FT dl = 0;
      Polynomial_1_3 pol2 = Polynomial_1_3(al,bl,cl,dl);
      Circle_3 c4 = theConstruct_circle_3(std::make_pair(es1, pol2));
      std::vector< CGAL::Object > intersection_4;
      theIntersect_3(c1, c4, std::back_inserter(intersection_4));
      assign(cap1,intersection_4[0]);
      assign(cap2,intersection_4[1]);
      _test_bbox<SK>(cap1.first);
      _test_bbox<SK>(cap2.first);
    }

    const FT a_c = FT(va) / FT(10);
    const FT b_c = 1;
    const FT c_c = 0;
    const FT d_c = -FT(va) / FT(10);
    Polynomial_1_3 pol2 = Polynomial_1_3(a_c,b_c,c_c,d_c);
    Circle_3 c5 = theConstruct_circle_3(std::make_pair(es2, pol2));
    std::vector< CGAL::Object > intersection_5;
    theIntersect_3(c1, c5, std::back_inserter(intersection_5));
    assign(cap1,intersection_5[0]);
    assign(cap2,intersection_5[1]);
    _test_bbox<SK>(cap1.first);
    _test_bbox<SK>(cap2.first);
  }

  std::cout << "Testing the bbox of Circle_3..." << std::endl;

  Sphere_3 s = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,1));
  Sphere_3 s_t10 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(10,10,10,1));
  for(int vx=-3;vx<4;vx++) {
    for(int vy=-3;vy<4;vy++) {
      for(int vz=-3;vz<4;vz++) {
        for(int vr=1;vr<6;vr++) {
          const FT x = FT(vx);
          const FT y = FT(vy);
          const FT z = FT(vz);
          const FT r = FT(vr) / FT(2);
          Sphere_3 sl_1 = theConstruct_sphere_3(
            Polynomial_for_spheres_2_3(x,y,z,r*r));
          Sphere_3 sl_2 = theConstruct_sphere_3(
            Polynomial_for_spheres_2_3(x+10,y+10,z+10,r*r));
          int d2 = (vx*vx + vy*vy + vz*vz);
          CGAL::Object intersection_1 = theIntersect_3(s, sl_1);
          CGAL::Object intersection_2 = theIntersect_3(s_t10, sl_2);
          // No intersection
          if((d2 > (r+1)*(r+1)) || (d2 < (r-1)*(r-1))) continue; 
          if(x == 0 && y == 0 && z == 0) continue;
          if((d2 == (r+1)*(r+1)) || (d2 == (r-1)*(r-1))) continue; 
          Circle_3 circle1, circle2;
          assign(circle1, intersection_1);
          assign(circle2, intersection_2);
          _test_bbox<SK>(circle1);
          _test_bbox<SK>(circle2);
        }
      }
    }
  }

  std::cout << "Testing the bbox of Line_arc_3..." << std::endl;
  for(int vx=0;vx<6;vx++) {
    for(int vy=0;vy<6;vy++) {
      for(int vz=0;vz<6;vz++) { 
        if(vx == 0 && vy == 0 && vz == 0) continue;
        const FT a = FT(vx);
        const FT b = FT(vy);
        const FT c = FT(vz);
        Line_3 l = theConstruct_line_3(Point_3(0,0,0), Point_3(a,b,c));
        for(int t1=-2;t1<3;t1++) {
          Point_3 source = Point_3(a*t1,b*t1,c*t1);
          for(int t2=t1+1;t2<5;t2++) {
            Point_3 target = Point_3(a*t2,b*t2,c*t2);
            Line_arc_3 la = theConstruct_line_arc_3(l,source,target);
            Bbox_3 b = theConstruct_bbox_3(la);
            if(vx != 0) {
              assert(FT(b.xmin()) > (t1-1)*vx);
              assert(FT(b.xmin()) <= t1*vx);
              assert(FT(b.xmax()) >= t2*vx);
              assert(FT(b.xmax()) < (t2+1)*vx);
            }
            if(vy != 0) {
              assert(FT(b.ymin()) > (t1-1)*vy);
              assert(FT(b.ymin()) <= t1*vy);
              assert(FT(b.ymax()) >= t2*vy);
              assert(FT(b.ymax()) < (t2+1)*vy);
            }
            if(vz != 0) {
              assert(FT(b.zmin()) > (t1-1)*vz);
              assert(FT(b.zmin()) <= t1*vz);
              assert(FT(b.zmax()) >= t2*vz);
              assert(FT(b.zmax()) < (t2+1)*vz);
            }
          }
        }
      }
    }
  }

  // At the moment this bbox is the same as the Circle_3 one
  // the bbox of the circular arc is the bbox of its supporting circle
  std::cout << "Testing the bbox of Circular_arc_3..." << std::endl;
}

template <class SK>
void _test_split_construct(SK sk) {

  typedef typename SK::FT                               FT;

  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Line_3                           Line_3;


  typedef typename SK::Circle_3                         Circle_3;
  typedef typename SK::Line_arc_3                       Line_arc_3;
  typedef typename SK::Circular_arc_3                   Circular_arc_3;
  typedef typename SK::Algebraic_kernel                 AK;

  typedef typename SK::Equal_3                          Equal_3;

  typedef typename SK::Split_3                          Split_3;

  typedef typename SK::Construct_circular_arc_3         Construct_circular_arc_3;
  typedef typename SK::Construct_circular_arc_point_3   Construct_circular_arc_point_3;
  typedef typename SK::Construct_circle_3               Construct_circle_3;


  typedef typename SK::Construct_line_3                 Construct_line_3;
  typedef typename SK::Construct_line_arc_3             Construct_line_arc_3;
  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;

  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;

  (void)/* Has_on_3 theHas_on_3 = */ sk.has_on_3_object();
  Equal_3 theEqual_3 = sk.equal_3_object();
  Split_3 theSplit_3 = sk.split_3_object();
  (void)/* Intersect_3 theIntersect_3 = */ sk.intersect_3_object();
  (void)/* Get_equation theGet_equation = */ sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  (void)/* Construct_plane_3 theConstruct_plane_3 = */ sk.construct_plane_3_object();
  Construct_line_3 theConstruct_line_3 = sk.construct_line_3_object();
  Construct_line_arc_3 theConstruct_line_arc_3 = sk.construct_line_arc_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = sk.construct_circular_arc_point_3_object();

  std::cout << "Testing Split a Line_arc_3..." << std::endl;
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
          for(int t2=t1+1;t2<5;t2++) {
            Point_3 target = Point_3(a*t2,b*t2,c*t2);
            Line_arc_3 la = theConstruct_line_arc_3(l,source,target);
            const FT tm = FT(t1+t2)/2;
            Point_3 mdl = Point_3(a*tm,b*tm,c*tm);
            Line_arc_3 l1, l2;
            theSplit_3(la, mdl, l1, l2);
            assert(theEqual_3(l1.source(), la.source()));
            assert(theEqual_3(l1.target(), mdl));
            assert(theEqual_3(l2.source(), mdl));
            assert(theEqual_3(l2.target(), la.target()));
          }
        }
      }
    }
  }

  std::cout << "Testing Split a Circular_arc_3..." << std::endl;
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
      for(int t=0;t<8;t++) {
        if(i == t || j == t) continue;
        if(simulate_has_on(i,j,t)) {
          Circular_arc_3 c1, c2;
          theSplit_3(ca, cp[t], c1, c2);
          if(SK().compare_xyz_3_object()(ca.source(), cp[t]) == CGAL::NEGATIVE) {
            assert(theEqual_3(c1.source(), ca.source()));
            assert(theEqual_3(c1.target(), cp[t]));
            assert(theEqual_3(c2.source(), cp[t]));
            assert(theEqual_3(c2.target(), ca.target()));
          } else {
            assert(theEqual_3(c1.source(), cp[t]));
            assert(theEqual_3(c1.target(), ca.target()));
            assert(theEqual_3(c2.source(), ca.source()));
            assert(theEqual_3(c2.target(), cp[t]));
          }
        }
      }
    }
  }
}

template <class SK>
void _test_extremal_points_construct(SK sk) {

  typedef typename SK::FT                               FT;

  typedef CGAL::Circular_arc_point_3<SK>                Circular_arc_point_3;



  typedef CGAL::Sphere_3<SK>                            Sphere_3;
  typedef CGAL::Circle_3<SK>                            Circle_3;


  typedef typename SK::Algebraic_kernel                 AK;







  typedef typename SK::Construct_circle_3               Construct_circle_3;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;



  typedef typename SK::Polynomials_for_circle_3         Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3       Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                   Polynomial_1_3;

  typedef typename AK::Root_for_spheres_2_3             Root_for_spheres_2_3;

  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = sk.construct_sphere_3_object();

  std::cout << "Testing {x,y,z}_extremal_points..." << std::endl;

  const Polynomials_for_circle_3 pcc_test = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(1,0,1,0));
  Circle_3 c = theConstruct_circle_3(pcc_test);
	Sphere_3 s = theConstruct_sphere_3(Polynomial_for_spheres_2_3(1,1,1,2));
	
	Circular_arc_point_3 pc[4], ps[6], res[6];

	pc[0] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)), 
	                             0,
	                             CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)));
	pc[1] = Root_for_spheres_2_3(0,1,0);
	pc[2] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1) / FT(2)),FT(2)), 
	                             0,
	                             CGAL::make_root_of_2(FT(0),-FT(FT(1) / FT(2)),FT(2)));
	pc[3] = Root_for_spheres_2_3(0,-1,0);
	
	res[0] = x_extremal_point(c, true);
	res[1] = x_extremal_point(c, false);
	res[2] = y_extremal_point(c, true);
	res[3] = y_extremal_point(c, false);
	res[4] = z_extremal_point(c, true);
	res[5] = z_extremal_point(c, false);
	
	assert(res[0] == pc[0]);
	assert(res[1] == pc[2]);
	assert(res[2] == pc[3]);
	assert(res[3] == pc[1]);
	assert(res[4] == pc[0]);
	assert(res[5] == pc[2]);
	
	ps[0] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(1),-FT(1),FT(2)), 1, 1);
	ps[1] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(1),FT(1),FT(2)), 1, 1);
	ps[2] = Root_for_spheres_2_3(1, CGAL::make_root_of_2(FT(1),-FT(1),FT(2)), 1);
	ps[3] = Root_for_spheres_2_3(1, CGAL::make_root_of_2(FT(1),FT(1),FT(2)), 1);
	ps[4] = Root_for_spheres_2_3(1, 1, CGAL::make_root_of_2(FT(1),-FT(1),FT(2)));
	ps[5] = Root_for_spheres_2_3(1, 1, CGAL::make_root_of_2(FT(1),FT(1),FT(2)));

	res[0] = x_extremal_point(s, true);
	res[1] = x_extremal_point(s, false);
	res[2] = y_extremal_point(s, true);
	res[3] = y_extremal_point(s, false);
	res[4] = z_extremal_point(s, true);
	res[5] = z_extremal_point(s, false);

	assert(res[0] == ps[0]);
	assert(res[1] == ps[1]);
	assert(res[2] == ps[2]);
	assert(res[3] == ps[3]);
	assert(res[4] == ps[4]);
	assert(res[5] == ps[5]);

}

template <class SK>
void _test_calls(SK /*sk*/) {
	CGAL::Line_arc_3<SK> la =
		SK().construct_line_arc_3_object()(CGAL::Point_3<SK>(1,1,1), CGAL::Point_3<SK>(1,2,1));
	CGAL::Circular_arc_point_3<SK> cap = SK().construct_circular_max_vertex_3_object()(la);
	cap = SK().construct_circular_min_vertex_3_object()(la);
	cap = SK().construct_circular_source_vertex_3_object()(la);
	cap = SK().construct_circular_target_vertex_3_object()(la);
}

template <class SK>
void _test_spherical_kernel_construct(SK sk)
{
  std::cout << "TESTING CONSTRUCTIONS" << std::endl;
  _test_circular_arc_point_construct(sk);
  _test_sphere_construct(sk);
  _test_plane_construct(sk);
  _test_line_construct(sk);
  _test_circle_construct(sk);
  _test_line_arc_construct(sk);
  _test_circular_arc_construct(sk);
  _test_intersection_construct(sk);
  _test_split_construct(sk);
  _test_bounding_box_construct(sk);
  _test_extremal_points_construct(sk);
	_test_calls(sk);
  std::cout << "All tests on construction are OK." << std::endl;
}
