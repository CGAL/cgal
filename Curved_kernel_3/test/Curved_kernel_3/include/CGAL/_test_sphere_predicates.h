#include <CGAL/Random.h>

template <class SK>
void _test_has_on_predicate(SK sk) {
  typedef typename SK::RT                               RT;
  typedef typename SK::FT                               FT;
  typedef typename SK::Root_of_2                        Root_of_2;
  typedef typename SK::Circular_arc_point_3             Circular_arc_point_3;
  typedef typename SK::Point_3                          Point_3;
  typedef typename SK::Plane_3                          Plane_3;
  typedef typename SK::Sphere_3                         Sphere_3;
  typedef typename SK::Circle_3                         Circle_3;
  typedef typename SK::Line_3                           Line_3;
  typedef typename SK::Algebraic_kernel                 AK;
  typedef typename SK::Construct_circle_3               Construct_circle_3;
  typedef typename SK::Construct_sphere_3               Construct_sphere_3;
  typedef typename SK::Construct_plane_3                Construct_plane_3;
  typedef typename SK::Construct_line_3                 Construct_line_3;
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
  Has_on_3 theHas_on_3 = sk.has_on_3_object();

  Sphere_3 s_1 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,1));
  std::cout << "Testing has_on(Sphere,Point)..." << std::endl;
  Point_3 p_1_s_1 = Point_3(1,0,0);
  Point_3 p_2_s_1 = Point_3(0,1,0);
  Point_3 p_3_s_1 = Point_3(0,0,1);
  Point_3 p_4_s_1 = Point_3(1,0,1);
  assert(theHas_on_3(s_1,p_1_s_1));
  assert(theHas_on_3(s_1,p_2_s_1));
  assert(theHas_on_3(s_1,p_3_s_1));
  assert(!theHas_on_3(s_1,p_4_s_1));
  std::cout << "Testing has_on(Sphere,Circular_arc_point)..." << std::endl;
  Root_of_2 sqrt_1_div_3 = make_root_of_2(FT(0),FT(1),FT(1,3));
  Root_of_2 sqrt_1_div_2 = make_root_of_2(FT(0),FT(1),FT(1,2));
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
  std::cout << "Testing has_on(Plane,Point)..." << std::endl;
  Point_3 p_1_p_1 = Point_3(-2,-1,-2);
  Point_3 p_2_p_1 = Point_3(-FT(5,3),-FT(5,3),-FT(5,3));
  Point_3 p_3_p_1 = Point_3(-10,0,0);
  Point_3 p_4_p_1 = Point_3(-2,-2,-1);
  assert(theHas_on_3(p_1,p_1_p_1));
  assert(theHas_on_3(p_1,p_2_p_1));
  assert(theHas_on_3(p_1,p_3_p_1));
  assert(!theHas_on_3(p_1,p_4_p_1));
  std::cout << "Testing has_on(Plane,Circular_arc_point)..." << std::endl;
  Root_of_2 r_1_1_p_1 = make_root_of_2(FT(0),FT(4),FT(2));
  Root_of_2 r_1_2_p_1 = make_root_of_2(FT(-5),FT(-5),FT(2));
  Root_of_2 r_1_3_p_1 = make_root_of_2(FT(0),FT(2),FT(2));
  Root_for_spheres_2_3 r_1_p_1 = Root_for_spheres_2_3(r_1_1_p_1,r_1_2_p_1,r_1_3_p_1);
  Root_for_spheres_2_3 r_2_p_1 = Root_for_spheres_2_3(r_1_2_p_1,r_1_2_p_1,r_1_3_p_1);
  Circular_arc_point_3 cp_1_p_1 = Circular_arc_point_3(r_1_p_1);
  Circular_arc_point_3 cp_2_p_1 = Circular_arc_point_3(r_2_p_1);
  assert(!theHas_on_3(p_1,cp_2_p_1));
  
  Line_3 l_1 = theConstruct_line_3(Polynomials_for_line_3(1,1,-1,3,1,0));
  std::cout << "Testing has_on(Line,Point)..." << std::endl;
  Point_3 p_1_l_1 = Point_3(1,3,0);
  Point_3 p_2_l_1 = Point_3(0,4,-1);
  Point_3 p_3_l_1 = Point_3(2,2,1);
  Point_3 p_4_l_1 = Point_3(1,1,1);
  assert(theHas_on_3(l_1,p_1_l_1));
  assert(theHas_on_3(l_1,p_2_l_1));
  assert(theHas_on_3(l_1,p_3_l_1));
  assert(!theHas_on_3(l_1,p_4_l_1));
  std::cout << "Testing has_on(Line,Circular_arc_point)..." << std::endl;
  Root_of_2 r_1_1_l_1 = make_root_of_2(FT(1),FT(1),FT(5));
  Root_of_2 r_1_2_l_1 = make_root_of_2(FT(3),FT(-1),FT(5));
  Root_of_2 r_1_3_l_1 = make_root_of_2(FT(0),FT(1),FT(5));
  Root_for_spheres_2_3 r_1_l_1 = Root_for_spheres_2_3(r_1_1_l_1,r_1_2_l_1,r_1_3_l_1);
  Root_for_spheres_2_3 r_2_l_1 = Root_for_spheres_2_3(r_1_1_l_1,r_1_1_l_1,r_1_1_l_1);
  Circular_arc_point_3 cp_1_l_1 = Circular_arc_point_3(r_1_l_1);
  Circular_arc_point_3 cp_2_l_1 = Circular_arc_point_3(r_2_l_1);
  assert(!theHas_on_3(l_1,cp_2_l_1));

  const Polynomials_for_circle_3 pc1 = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(1,0,0,0));
  Circle_3 c_1 = theConstruct_circle_3(pc1);
  std::cout << "Testing has_on(Circle,Point)..." << std::endl;
  Point_3 p_1_c_1 = Point_3(0,0,1);
  Point_3 p_2_c_1 = Point_3(0,1,0);
  Point_3 p_3_c_1 = Point_3(1,0,0);
  assert(theHas_on_3(c_1,p_1_c_1));
  assert(theHas_on_3(c_1,p_2_c_1));
  assert(!theHas_on_3(c_1,p_3_c_1));
  std::cout << "Testing has_on(Circle,Circular_arc_point)..." << std::endl;
  const Polynomials_for_circle_3 pc2 = 
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(1,1,1,0));
  Circle_3 c_2 = theConstruct_circle_3(pc2);
  Root_of_2 r_1_1_c_2 = Root_of_2(FT(1,2));
  Root_of_2 r_1_2_c_2 = make_root_of_2(-FT(1,4),-FT(1,4),FT(5));
  Root_of_2 r_1_3_c_2 = make_root_of_2(-FT(1,4),FT(1,4),FT(5));
  Root_for_spheres_2_3 r_1_c_2 = Root_for_spheres_2_3(r_1_1_c_2,r_1_2_c_2,r_1_3_c_2);
  Root_for_spheres_2_3 r_2_c_2 = Root_for_spheres_2_3(r_1_2_c_2,r_1_2_c_2,r_1_2_c_2);
  Circular_arc_point_3 cp_1_c_2 = Circular_arc_point_3(r_1_c_2);
  Circular_arc_point_3 cp_2_c_2 = Circular_arc_point_3(r_2_c_2);
  assert(theHas_on_3(c_2,cp_1_c_2));
  assert(!theHas_on_3(c_2,cp_2_c_2));

  // Don't need to test has_on(Plane,Line). It is internal on Cartesian

  std::cout << "Testing has_on(Sphere,Circle)..." << std::endl;
  assert(theHas_on_3(s_1,c_1));
  assert(theHas_on_3(s_1,c_2));
  Sphere_3 s_2 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(0,0,0,2));
  assert(!theHas_on_3(s_2,c_1));
  assert(!theHas_on_3(s_2,c_2));
  Sphere_3 s_3 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(5,0,0,26));
  assert(theHas_on_3(s_3,c_1));
  assert(!theHas_on_3(s_3,c_2));

  std::cout << "Testing has_on(Plane,Circle)..." << std::endl;
  Plane_3 p_2 = theConstruct_plane_3(Polynomial_1_3(1,1,1,0));
  Plane_3 p_3 = theConstruct_plane_3(Polynomial_1_3(3,3,3,0));
  Plane_3 p_4 = theConstruct_plane_3(Polynomial_1_3(1,0,0,0));
  assert(theHas_on_3(p_2,c_2));
  assert(theHas_on_3(p_3,c_2));
  assert(!theHas_on_3(p_4,c_2));
  assert(!theHas_on_3(p_2,c_1));
  assert(!theHas_on_3(p_3,c_1));
  assert(theHas_on_3(p_4,c_1));

}

template <class SK>
void _test_spherical_kernel_predicates(SK sk)
{
  std::cout << "TESTING PREDICATES" << std::endl;
  _test_has_on_predicate(sk);
  std::cout << "All tests on predicates are OK." << std::endl;
}
