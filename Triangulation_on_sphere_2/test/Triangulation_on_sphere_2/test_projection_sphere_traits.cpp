#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>

#include <CGAL/enum.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point;

typedef CGAL::Projection_sphere_traits_3<K>                      Gt;
typedef Gt::Point_2                                              Point_2;
typedef Gt::Compare_xyz_3                                        Compare_xyz_3;
typedef Gt::Compute_squared_distance_2                           Compute_squared_distance_2;
typedef Gt::Construct_projected_point_3                          Construct_projected_point_3;
typedef Gt::Coradial_sphere_2                                    Coradial_sphere_2;
typedef Gt::Inside_cone_2                                        Inside_cone_2;
typedef Gt::Orientation_2                                        Orientation_2;
typedef Gt::Power_test_2                                         Power_test_2;

int main(int, char**)
{
  Gt traits(Point(0, 0, 0), 10);
  traits.set_radius(1);

  Construct_projected_point_3 cst = traits.construct_projected_point_3_object();
  Coradial_sphere_2 coradial = traits.coradial_sphere_2_object();
  Compute_squared_distance_2 squared_distance= traits.compute_squared_distance_2_object();
  Compare_xyz_3 compare_xyz = traits.compare_xyz_3_object();
  Inside_cone_2 inside_cone = traits.inside_cone_2_object();
  Orientation_2 orientation = traits.orientation_2_object();
  Power_test_2 power_test = traits.power_test_2_object();

  Point c(0, 0, 0);

  // Testing with Points projected on unit sphere
  Point p11(1.5,  1.5, 1.5);
  Point p12( -1,    1,   1);
  Point p13(  0,  0.7,   0);
  Point p14(  0,   -2,   0);

  Point p21(   2, 0,   -1);
  Point p22(-1.8, 0, -0.7);
  Point p23(   0, 0,  2.6);

  // Points with same coordinates
  Point p31(0.6, 0.3, -1);
  Point p32(0.6, 0.3, -1);
  Point p33(0.6, 0.3,  1);

  // Coradial points
  Point p41(  1,   1,    1);
  Point p42(0.7, 0.7,  0.7);
  Point p43(0.3, 0.6, -2.4);
  Point p44(0.6, 1.2, -4.8);

  // Inside Cone
  Point p51(  1,   1,   1);
  Point p52( -1,   1,   1);
  Point p53(  1,   3,   3); // inside
  Point p54(0.8, 0.8, 0.8); // boundary
  Point p55(  0,  -2,  -2); // outside

  //distance
  Point p61(  1,   1,   0);
  Point p62( -1,   1,   0);

  // transform to projected points
  Point_2 pp11(p11, c);
  Point_2 pp12(p12, c);
  Point_2 pp13(p13, c);
  Point_2 pp14(p14, c);
  Point_2 pp21(p21, c);
  Point_2 pp22(p22, c);
  Point_2 pp23(p23, c);
  Point_2 pp31(p31, c);
  Point_2 pp32(p32, c);
  Point_2 pp33(p33, c);
  Point_2 pp41(p41, c);
  Point_2 pp42(p42, c);
  Point_2 pp43(p43, c);
  Point_2 pp44(p44, c);
  Point_2 pp51(p51, c);
  Point_2 pp52(p52, c);
  Point_2 pp53(p53, c);
  Point_2 pp54(p54, c);
  Point_2 pp55(p55, c);
  Point_2 pp61(p61, c);
  Point_2 pp62(p62, c);

  std::cout << "Test Orientation" << std::endl;
  assert(orientation(pp11, pp12, pp13) == CGAL::NEGATIVE);
  assert(orientation(pp13, pp12, pp11) == CGAL::POSITIVE);
  assert(orientation(pp21, pp22, pp23) == CGAL::ON_ORIENTED_BOUNDARY);
  assert(orientation(pp11, pp12, pp13, pp14) == CGAL::POSITIVE);

  std::cout << "Test Power_Test_2" << std::endl;
  assert(power_test(pp11, pp12, pp13, pp14) == CGAL::POSITIVE);
  assert(power_test(pp14, pp11, pp12, pp13) == CGAL::NEGATIVE);
  assert(power_test(pp21, pp22, pp23, pp11) == CGAL::POSITIVE);

  std::cout << "Test Coradial_sphere_2" << std::endl;
  assert(coradial(pp41, pp41));
  assert(coradial(pp41, pp42));
  assert(coradial(pp43, pp44));

  std::cout << "Test Inside_Cone_sphere_2" << std::endl;
  assert(inside_cone(pp51, pp52, pp53));
  assert(!inside_cone(pp51, pp52, pp54));
  assert(!inside_cone(pp51, pp52, pp55));

  std::cout << "Test Squared_Distance_sphere_2" << std::endl;
  double dist1 = squared_distance(pp41, pp42);
  assert(dist1 == 0);

  // double dist2 = squared_distance(p61, p62);
  // assert( dist2 == 2);

  std::cout << "Test compare_xyz_3" << std::endl;

  assert(compare_xyz(pp31, pp31) == CGAL::EQUAL);
  assert(compare_xyz(pp31, pp32) == CGAL::EQUAL);
  assert(compare_xyz(pp31, pp33) == CGAL::SMALLER);
  assert(compare_xyz(pp33, pp31) == CGAL::LARGER);

  return EXIT_SUCCESS;
}
