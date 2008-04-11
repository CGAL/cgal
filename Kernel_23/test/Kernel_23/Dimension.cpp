// Test program for Dimension<>.
// Sylvain Pion, 2005.

#include <CGAL/Cartesian.h>
#include <cassert>
#include <CGAL/Cartesian_d.h>

#include <CGAL/Kernel/Dimension.h>
#include <CGAL/Dimension_tag.h>

template < typename K >
void test(CGAL::Dimension_tag<2>)
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
void test(CGAL::Dimension_tag<3>)
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
void test(CGAL::Dimension_tag<0>)
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
  test<CGAL::Cartesian<int> >(CGAL::Dimension_tag<2>());
  test<CGAL::Cartesian<int> >(CGAL::Dimension_tag<3>());
  test<CGAL::Cartesian_d<int> >(CGAL::Dimension_tag<0>());
  
  return 0;
}
