#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Projection_on_sphere_traits_3.h>

#include <CGAL/enum.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;

typedef CGAL::Projection_on_sphere_traits_3<K>                   Gt;

typedef Gt::Point_3                                              Point_3;
typedef Gt::Point_on_sphere_2                                    Point_on_sphere_2;
typedef Gt::Construct_point_on_sphere_2                          Construct_point_on_sphere_2;

typedef Gt::Compare_on_sphere_2                                  Compare_on_sphere_2;
typedef Gt::Equal_on_sphere_2                                    Equal_on_sphere_2;
typedef Gt::Collinear_are_strictly_ordered_on_great_circle_2     Collinear_are_strictly_ordered_on_great_circle_2;
typedef Gt::Orientation_on_sphere_2                              Orientation_on_sphere_2;
typedef Gt::Side_of_oriented_circle_on_sphere_2                  Side_of_oriented_circle_on_sphere_2;

typedef Gt::Compare_on_sphere_2                                  Compare_on_sphere_2;

int main(int, char**)
{
  Point_3 c(0, 0, 0);

  Gt traits(c, 10);

  Equal_on_sphere_2 equal_on_sphere = traits.equal_on_sphere_2_object();
  Compare_on_sphere_2 compare_on_sphere = traits.compare_on_sphere_2_object();
  Collinear_are_strictly_ordered_on_great_circle_2 inside_cone = traits.collinear_are_strictly_ordered_on_great_circle_2_object();
  Orientation_on_sphere_2 orientation = traits.orientation_on_sphere_2_object();
  Side_of_oriented_circle_on_sphere_2 side_of_oriented_circle = traits.side_of_oriented_circle_on_sphere_2_object();

  Construct_point_on_sphere_2 cst = traits.construct_point_on_sphere_2_object();

  // Testing with Points projected on unit sphere
  Point_3 p11(1.5,  1.5, 1.5);
  Point_3 p12( -1,    1,   1);
  Point_3 p13(  0,  0.7,   0);
  Point_3 p14(  0,   -2,   0);

  Point_3 p21(   2, 0,   -1);
  Point_3 p22(-1.8, 0, -0.7);
  Point_3 p23(   0, 0,  2.6);

  // Points with same coordinates
  Point_3 p31(0.6, 0.3, -1);
  Point_3 p32(0.6, 0.3, -1);
  Point_3 p33(0.6, 0.3,  1);

  // Equal points
  Point_3 p41(  1,   1,    1);
  Point_3 p42(0.7, 0.7,  0.7);
  Point_3 p43(0.3, 0.6, -2.4);
  Point_3 p44(0.6, 1.2, -4.8);

  // Inside Cone
  Point_3 p51(  1,   1,   1);
  Point_3 p52( -1,   1,   1);
  Point_3 p53(  1,   3,   3); // inside
  Point_3 p54(0.8, 0.8, 0.8); // boundary
  Point_3 p55(  0,  -2,  -2); // outside

  //distance
  Point_3 p61(  1,   1,   0);
  Point_3 p62( -1,   1,   0);

  // transform to projected points
  Point_on_sphere_2 pp11 = cst(p11);
  Point_on_sphere_2 pp12 = cst(p12);
  Point_on_sphere_2 pp13 = cst(p13);
  Point_on_sphere_2 pp14 = cst(p14);
  Point_on_sphere_2 pp21 = cst(p21);
  Point_on_sphere_2 pp22 = cst(p22);
  Point_on_sphere_2 pp23 = cst(p23);
  Point_on_sphere_2 pp31 = cst(p31);
  Point_on_sphere_2 pp32 = cst(p32);
  Point_on_sphere_2 pp33 = cst(p33);
  Point_on_sphere_2 pp41 = cst(p41);
  Point_on_sphere_2 pp42 = cst(p42);
  Point_on_sphere_2 pp43 = cst(p43);
  Point_on_sphere_2 pp44 = cst(p44);
  Point_on_sphere_2 pp51 = cst(p51);
  Point_on_sphere_2 pp52 = cst(p52);
  Point_on_sphere_2 pp53 = cst(p53);
  Point_on_sphere_2 pp54 = cst(p54);
  Point_on_sphere_2 pp55 = cst(p55);

  std::cout << "Test Orientation" << std::endl;
  assert(orientation(pp11, pp12, pp13) == CGAL::NEGATIVE);
  assert(orientation(pp13, pp12, pp11) == CGAL::POSITIVE);
  assert(orientation(pp21, pp22, pp23) == CGAL::ON_ORIENTED_BOUNDARY);

  std::cout << "Test Side_of_oriented_circle_on_sphere_2" << std::endl;
  assert(side_of_oriented_circle(pp11, pp12, pp13, pp14) == CGAL::POSITIVE);
  assert(side_of_oriented_circle(pp14, pp11, pp12, pp13) == CGAL::NEGATIVE);
  assert(side_of_oriented_circle(pp21, pp22, pp23, pp11) == CGAL::POSITIVE);

  std::cout << "Test Equal_on_sphere_2" << std::endl;
  assert(equal_on_sphere(pp41, pp41));
  assert(equal_on_sphere(pp41, pp42));
  assert(equal_on_sphere(pp43, pp44));

  std::cout << "Test Inside_cone_sphere_2" << std::endl;
  assert(inside_cone(pp51, pp52, pp53));
  assert(!inside_cone(pp51, pp52, pp54));
  assert(!inside_cone(pp51, pp52, pp55));

  std::cout << "Test Compare_on_sphere_2" << std::endl;
  assert(compare_on_sphere(pp31, pp31) == CGAL::EQUAL);
  assert(compare_on_sphere(pp31, pp32) == CGAL::EQUAL);
  assert(compare_on_sphere(pp31, pp33) == CGAL::SMALLER);
  assert(compare_on_sphere(pp33, pp31) == CGAL::LARGER);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
