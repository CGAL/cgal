#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>

#include <CGAL/enum.h>
#include <CGAL/Origin.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_3                                                  Point;

typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<K>          Gt;

typedef Gt::Orientation_on_sphere_2                                 Orientation_on_sphere_2;
typedef Gt::Side_of_oriented_circle_on_sphere_2                     Side_of_oriented_circle_on_sphere_2;

// This only tests non-trivial functions (a lot of functors are just calls to some functors
// of the 3D linear kernel, such as:
// typedef typename LK::Compare_xyz_3 Compare_on_sphere_2;

int main(int, char**)
{
  Point P0 (   0,    0,    0);
  Point P1 (   1,    0,    0);
  Point P2 (   0,    1,    0);
  Point P3 (   0,    0,    1);
  Point P4 (   1,    1,    1);
  Point P5 (  -1,   -1,   -1);
  Point P6 ( 0.5,  0.5,    0);
  Point P7 (   2,    0,    0);
  Point P8 ( 0.5,    0,    0);
  Point P9 (   1,    1,    0);
  Point P10(0.25, 0.25,    0);
  Point P11( 0.5,    0,  0.5);
  Point P12( 0.5,    0, -0.5);
  Point P13(   0,    0,  0.5);
  Point P14(  -1,    0,    0);

  Gt traits;
  Gt traits_2(CGAL::ORIGIN);
  Gt traits_3(Point(1, 2, 3), 15);

  Gt traits_4(traits_3);
  assert(traits_4.radius() == traits_3.radius());
  assert(traits_4.center() == traits_3.center());

  Gt traits_5 = traits_4;
  assert(traits_4.radius() == traits_5.radius());
  assert(traits_4.center() == traits_5.center());

  std::cout << "Test Side_of_oriented_circle_on_sphere_2" << std::endl;
  CGAL::Oriented_side result;
  result = traits.side_of_oriented_circle_on_sphere_2_object()(P1, P2, P3, P4);
  assert(result == CGAL::ON_POSITIVE_SIDE);
  result = traits.side_of_oriented_circle_on_sphere_2_object()(P1, P2, P3, P5);
  assert(result == CGAL::ON_NEGATIVE_SIDE);
  result = traits.side_of_oriented_circle_on_sphere_2_object()(P1, P2, P3, P6);
  assert(result == CGAL::ON_ORIENTED_BOUNDARY);

  // Former tests for the degenerate cases,
  // https://github.com/MaelRL/cgal/blob/Triangulation_on_sphere_2-GF-old/
  //   Triangulation_on_sphere_2/include/CGAL/Delaunay_triangulation_sphere_traits_2.h#L29
  //
  // power_test_2(p,q,r) where p,q and r are coplanar points with the center of the sphere
//  result = traits.power_test_2_object()(P1, P2, P9);
//  assert(result == CGAL::ON_POSITIVE_SIDE);
//  result = traits.power_test_2_object()(P2, P1, P9);
//  assert(result == CGAL::ON_POSITIVE_SIDE);
//  result = traits.power_test_2_object()(P1, P2, P6);
//  assert(result == CGAL::ON_ORIENTED_BOUNDARY);
//  result = traits.power_test_2_object()(P1, P2, P10);
//  assert(result == CGAL::ON_NEGATIVE_SIDE);
//  result = traits.power_test_2_object()(P2, P1, P10);
//  assert(result == CGAL::ON_NEGATIVE_SIDE);

  // power_test_2(p,q) where p, q, and the center of the sphere are colinear
//  result = traits.power_test_2_object()(P1, P7);
//  assert(result == CGAL::ON_POSITIVE_SIDE);
//  result = traits.power_test_2_object()(P1, P1);
//  assert(result == CGAL::ON_ORIENTED_BOUNDARY);
//  result = traits.power_test_2_object()(P1, P8);
//  assert(result == CGAL::ON_NEGATIVE_SIDE);

  // With a center different from the origin
  Gt traits_6(Point(0,0,1));

  // power_test_2(p,q,r,s)
  result = traits_6.side_of_oriented_circle_on_sphere_2_object()(P0, P2, P1, P5);
  assert(result == CGAL::ON_POSITIVE_SIDE);
  result = traits_6.side_of_oriented_circle_on_sphere_2_object()(P0, P2, P1, P4);
  assert(result == CGAL::ON_NEGATIVE_SIDE);
  result = traits_6.side_of_oriented_circle_on_sphere_2_object()(P0, P1, P1, P10);
  assert(result == CGAL::ON_ORIENTED_BOUNDARY);

  // power_test_2(p,q,r) where p,q and r are coplanar points
//  result = traits_6.power_test_2_object()(P14, P1, P12);
//  assert(result == CGAL::ON_POSITIVE_SIDE);
//  result = traits_6.power_test_2_object()(P1, P2, P6);
//  assert(result == CGAL::ON_ORIENTED_BOUNDARY);
//  result = traits_6.power_test_2_object()(P14, P1, P11);
//  assert(result == CGAL::ON_NEGATIVE_SIDE);

  // power_test_2(p,q) where p, q and sphere are colinear
//  result = traits_6.power_test_2_object()(P13, P0);
//  assert(result == CGAL::ON_POSITIVE_SIDE);
//  result = traits_6.power_test_2_object()(P13, P13);
//  assert(result == CGAL::ON_ORIENTED_BOUNDARY);
//  result = traits_6.power_test_2_object()(P0, P13);
//  assert(result == CGAL::ON_NEGATIVE_SIDE);

  std::cout << "Test Orientation_on_sphere_2" << std::endl;
  Point p21( 0.5, 0.5,   sqrt(0.75));
  Point p22(-0.5, 0.5,   sqrt(0.75));
  Point p23(   0,  -1,   sqrt(0.75));
  Point p24(   0,   0,           -1);
  Point p25(   0,   0,         -1.5);
  Point p26(   1,   1, 2*sqrt(0.75));

  result = traits.orientation_on_sphere_2_object()(p21, p22, p23);
  assert(result == CGAL::ON_POSITIVE_SIDE);
  result = traits.orientation_on_sphere_2_object()(p23, p22, p21);
  assert(result == CGAL::ON_NEGATIVE_SIDE);

  result = traits.side_of_oriented_circle_on_sphere_2_object()(p24, p21, p22, p23);
  assert(result == CGAL::ON_POSITIVE_SIDE);

  std::cout << "Test Equal_on_sphere_2" << std::endl;
  bool coradial;
  coradial = traits.equal_on_sphere_2_object()(p24, p25);
  assert(coradial);
  coradial = traits.equal_on_sphere_2_object()(p25, p24);
  assert(coradial);
  coradial = traits.equal_on_sphere_2_object()(p21, p26);
  assert(coradial);
  coradial = traits.equal_on_sphere_2_object()(p22, p26);
  assert(!coradial);

  std::cout << "Test Collinear_are_strictly_ordered_on_great_circle_2" << std::endl;
  Point p31(1,1,1);
  Point p32(-1,1,1);
  Point p33(0, sqrt(1.5),sqrt(1.5));
  Point p34(0, -sqrt(2),-1);
  Point p35(0.9,0.9,0.9);
  Point p36(0, sqrt(2.5),0);

  bool inside;
  // cone defined by 0, p31 and p32
  inside = traits.collinear_are_strictly_ordered_on_great_circle_2_object()(p31, p32, p33);
  assert(inside);
  inside = traits.collinear_are_strictly_ordered_on_great_circle_2_object()(p31, p32, p34);
  assert(!inside);
  inside = traits.collinear_are_strictly_ordered_on_great_circle_2_object()(p31, p32, p35);
  assert(!inside);
  //not coplanar
  inside = traits.collinear_are_strictly_ordered_on_great_circle_2_object()(p31, p32, p36);
  assert(!inside);

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
