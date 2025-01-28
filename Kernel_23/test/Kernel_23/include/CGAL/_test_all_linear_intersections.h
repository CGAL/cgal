#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>
#include <optional>

template <class T>
bool is_intersection_empty(const std::optional<T>& t)
{
  return bool(t);
}

bool is_intersection_empty(CGAL::Object o)
{
  return o.empty();
}

template <class K1, class K2, class T1, class T2>
void test_2d(const T1& t1, const T2& t2)
{
  CGAL::Cartesian_converter<K1,K2> convert;
  const auto res1 = CGAL::intersection(t1, t2);
  assert( is_intersection_empty( convert(res1) ) == is_intersection_empty(res1) );
}

template <class K1, class K2, class T1, class T2>
void test_3d(const T1& t1, const T2& t2)
{
  CGAL::Cartesian_converter<K1,K2> convert;
  const auto res1 = CGAL::intersection(t1, t2);
  assert( is_intersection_empty( convert(res1) ) == is_intersection_empty(res1) );
}

template <class K1, class K2>
void test_linear_intersections()
{
  typedef typename K1::Iso_rectangle_2 Iso_rectangle_2;
  typedef typename K1::Line_2 Line_2;
  typedef typename K1::Ray_2 Ray_2;
  typedef typename K1::Segment_2 Segment_2;
  typedef typename K1::Triangle_2 Triangle_2;
  typedef typename K1::Line_3 Line_3;
  typedef typename K1::Plane_3 Plane_3;
  typedef typename K1::Ray_3 Ray_3;
  typedef typename K1::Segment_3 Segment_3;
  typedef typename K1::Sphere_3 Sphere_3;
  typedef typename K1::Triangle_3 Triangle_3;
  typedef typename K1::Point_2 Point_2;
  typedef typename K1::Point_3 Point_3;


  Point_3 points_3D[4] = {
    Point_3(0,0,0),
    Point_3(0,0,1),
    Point_3(1,0,0),
    Point_3(0,1,0)
  };

  Point_2 points_2D[4]={
    Point_2(0,0),
    Point_2(0,1),
    Point_2(1,1),
    Point_2(-1,1)
  };

  Iso_rectangle_2 iso_rectangle_2(points_2D[0], points_2D[1]);
  Line_2 line_2(points_2D[0], points_2D[1]);
  Ray_2 ray_2(points_2D[0], points_2D[1]);
  Segment_2 segment_2(points_2D[0], points_2D[1]);
  Triangle_2 triangle_2(points_2D[0], points_2D[1], points_2D[2]);;
  Line_3 line_3(points_3D[0], points_3D[1]);
  Plane_3 plane_3(points_3D[0], points_3D[1], points_3D[2]);
  Ray_3 ray_3(points_3D[0], points_3D[1]);
  Segment_3 segment_3(points_3D[0], points_3D[1]);
  Sphere_3 sphere_3(points_3D[0], 1);
  Triangle_3 triangle_3(points_3D[0], points_3D[1], points_3D[2]);

  // Test 2D intersections from the linear kernel
  test_2d<K1,K2>(iso_rectangle_2, iso_rectangle_2);
  test_2d<K1,K2>(iso_rectangle_2, line_2);
  test_2d<K1,K2>(iso_rectangle_2, ray_2);
  test_2d<K1,K2>(iso_rectangle_2, segment_2);
  test_2d<K1,K2>(iso_rectangle_2, triangle_2);
  test_2d<K1,K2>(line_2, line_2);
  test_2d<K1,K2>(line_2, ray_2);
  test_2d<K1,K2>(line_2, segment_2);
  test_2d<K1,K2>(line_2, triangle_2);
  test_2d<K1,K2>(ray_2, ray_2);
  test_2d<K1,K2>(ray_2, segment_2);
  test_2d<K1,K2>(ray_2, triangle_2);
  test_2d<K1,K2>(segment_2, segment_2);
  test_2d<K1,K2>(segment_2, triangle_2);
  test_2d<K1,K2>(triangle_2, triangle_2);

  // Test 3D intersections from the linear kernel
  test_3d<K1,K2>(line_3, line_3);
  test_3d<K1,K2>(line_3, plane_3);
  test_3d<K1,K2>(line_3, ray_3);
  test_3d<K1,K2>(line_3, segment_3);
  test_3d<K1,K2>(line_3, triangle_3);
  test_3d<K1,K2>(plane_3, plane_3);
  test_3d<K1,K2>(plane_3, ray_3);
  test_3d<K1,K2>(plane_3, segment_3);
  test_3d<K1,K2>(plane_3, sphere_3);
  test_3d<K1,K2>(plane_3, triangle_3);
  test_3d<K1,K2>(ray_3, ray_3);
  test_3d<K1,K2>(ray_3, segment_3);
  test_3d<K1,K2>(ray_3, triangle_3);
  test_3d<K1,K2>(segment_3, segment_3);
  test_3d<K1,K2>(segment_3, triangle_3);
  test_3d<K1,K2>(sphere_3, sphere_3);
  test_3d<K1,K2>(triangle_3, triangle_3);

  CGAL::Cartesian_converter<K1,K2> convert;
  const auto res1 = CGAL::intersection(plane_3, plane_3, plane_3);
  assert( is_intersection_empty( convert(res1) ) == is_intersection_empty(res1) );
}
