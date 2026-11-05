// Homogeneous kernel with an FT that is not CGAL::Quotient, see issue #5967:
// the functors used below called FT::numerator() and FT::denominator().

#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Exact_rational.h>
#include <cassert>

typedef CGAL::Homogeneous<CGAL::Exact_integer, CGAL::Exact_rational> K;
typedef K::FT       FT;
typedef K::Point_2  Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Circle_2 Circle_2;
typedef K::Point_3  Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Sphere_3 Sphere_3;

int main()
{
  const FT f(2, 3);

  // Construct_divided_vector_2 and Construct_scaled_vector_2
  const Vector_2 v2(1, 2);
  assert(v2 / f == Vector_2(3, 6, 2));
  assert(v2 * f == Vector_2(2, 4, 3));

  // Construct_divided_vector_3 and Construct_scaled_vector_3
  const Vector_3 v3(1, 2, 3);
  assert(v3 / f == Vector_3(3, 6, 9, 2));
  assert(v3 * f == Vector_3(2, 4, 6, 3));

  // Construct_radical_line_2: x^2 + y^2 - 1/4 = (x-1)^2 + y^2 - 1/9 gives x = 41/72
  const Circle_2 c1(Point_2(0, 0), FT(1, 4));
  const Circle_2 c2(Point_2(1, 0), FT(1, 9));
  const K::Line_2 l = CGAL::radical_line(c1, c2);
  assert(l.has_on(Point_2(41, 0, 72)));
  assert(l.has_on(Point_2(41, 72, 72)));

  // Construct_radical_plane_3: same equation, the plane x = 41/72
  const Sphere_3 s1(Point_3(0, 0, 0), FT(1, 4));
  const Sphere_3 s2(Point_3(1, 0, 0), FT(1, 9));
  const K::Plane_3 p = CGAL::radical_plane(s1, s2);
  assert(p.has_on(Point_3(41, 0, 0, 72)));
  assert(p.has_on(Point_3(41, 72, 144, 72)));

  return 0;
}
