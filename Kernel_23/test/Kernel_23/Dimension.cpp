// Test program for Dimension<>.
// Sylvain Pion, 2005.

#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>

#include <CGAL/Kernel/Dimension.h>

template < typename K >
void test_2()
{
  assert( 2 == CGAL::Dimension<typename K::Point_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Vector_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Direction_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Line_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Ray_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Segment_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Triangle_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Iso_rectangle_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Circle_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Conic_2>::value );
  assert( 2 == CGAL::Dimension<typename K::Aff_transformation_2>::value );
}

template < typename K >
void test_3()
{
  assert( 3 == CGAL::Dimension<typename K::Point_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Plane_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Vector_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Direction_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Line_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Ray_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Segment_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Triangle_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Tetrahedron_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Iso_cuboid_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Sphere_3>::value );
  assert( 3 == CGAL::Dimension<typename K::Aff_transformation_3>::value );
}

template < typename K >
void test_d()
{
  assert( 0 == CGAL::Dimension<typename K::Point_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Hyperplane_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Vector_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Direction_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Line_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Ray_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Segment_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Iso_box_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Sphere_d>::value );
  assert( 0 == CGAL::Dimension<typename K::Aff_transformation_d>::value );
}

int main()
{
  test_2<CGAL::Cartesian<int> >();
  test_3<CGAL::Cartesian<int> >();
  test_d<CGAL::Cartesian_d<int> >();

  return 0;
}
