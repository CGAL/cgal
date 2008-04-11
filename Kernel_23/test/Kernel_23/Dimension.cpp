// Test program for Dimension<>.
// Sylvain Pion, 2005.

#include <cassert>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Dimension.h>
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
void test(CGAL::Dimension_tag<CGAL::Dynamic_dimension>)
{
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Point_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Hyperplane_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Vector_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Direction_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Line_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Ray_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Segment_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Iso_box_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Sphere_d>::value );
  assert( CGAL::Dynamic_dimension == CGAL::Dimension<typename K::Aff_transformation_d>::value );
}

int main()
{
  test<CGAL::Cartesian<int> >(CGAL::Dimension_tag<2>());
  test<CGAL::Cartesian<int> >(CGAL::Dimension_tag<3>());
  test<CGAL::Cartesian_d<int> >(CGAL::Dimension_tag<CGAL::Dynamic_dimension>());

  test<CGAL::Homogeneous<int> >(CGAL::Dimension_tag<2>());
  test<CGAL::Homogeneous<int> >(CGAL::Dimension_tag<3>());
  test<CGAL::Homogeneous_d<int> >(CGAL::Dimension_tag<CGAL::Dynamic_dimension>());
  
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(CGAL::Dimension_tag<2>());
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(CGAL::Dimension_tag<3>());

  test<CGAL::Exact_predicates_exact_constructions_kernel>(CGAL::Dimension_tag<2>());
  test<CGAL::Exact_predicates_exact_constructions_kernel>(CGAL::Dimension_tag<3>());

  return 0;
}
