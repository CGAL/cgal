// Test program for Ambient_dimension<> and Feature_dimension<>.
// Sylvain Pion, 2005, 2008.

#include <cassert>
#include <CGAL/Cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Dimension.h>
#include <CGAL/Dimension.h>

template < typename K >
void test(CGAL::Dimension_tag<2>)
{
  assert( 2 == CGAL::Ambient_dimension<typename K::Point_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Vector_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Direction_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Line_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Ray_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Segment_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Triangle_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Iso_rectangle_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Circle_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Conic_2>::value );
  assert( 2 == CGAL::Ambient_dimension<typename K::Aff_transformation_2>::value );
  assert( 2 == CGAL::Ambient_dimension<CGAL::Bbox_2>::value );

  assert( 0 == CGAL::Feature_dimension<typename K::Point_2>::value );
  assert( 0 == CGAL::Feature_dimension<typename K::Vector_2>::value );
  assert( 0 == CGAL::Feature_dimension<typename K::Direction_2>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Line_2>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Ray_2>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Segment_2>::value );
  assert( 2 == CGAL::Feature_dimension<typename K::Triangle_2>::value );
  assert( 2 == CGAL::Feature_dimension<typename K::Iso_rectangle_2>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Circle_2>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Conic_2>::value );
  // assert( ? == CGAL::Feature_dimension<typename K::Aff_transformation_2>::value );
  assert( 2 == CGAL::Feature_dimension<CGAL::Bbox_2>::value );
}

template < typename K >
void test(CGAL::Dimension_tag<3>)
{
  assert( 3 == CGAL::Ambient_dimension<typename K::Point_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Plane_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Vector_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Direction_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Line_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Ray_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Segment_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Triangle_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Tetrahedron_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Iso_cuboid_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Sphere_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Circle_3>::value );
  assert( 3 == CGAL::Ambient_dimension<typename K::Aff_transformation_3>::value );
  assert( 3 == CGAL::Ambient_dimension<CGAL::Bbox_3>::value );

  assert( 0 == CGAL::Feature_dimension<typename K::Point_3>::value );
  assert( 2 == CGAL::Feature_dimension<typename K::Plane_3>::value );
  assert( 0 == CGAL::Feature_dimension<typename K::Vector_3>::value );
  assert( 0 == CGAL::Feature_dimension<typename K::Direction_3>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Line_3>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Ray_3>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Segment_3>::value );
  assert( 2 == CGAL::Feature_dimension<typename K::Triangle_3>::value );
  assert( 3 == CGAL::Feature_dimension<typename K::Tetrahedron_3>::value );
  assert( 3 == CGAL::Feature_dimension<typename K::Iso_cuboid_3>::value );
  assert( 2 == CGAL::Feature_dimension<typename K::Sphere_3>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Circle_3>::value );
  // assert( ? == CGAL::Feature_dimension<typename K::Aff_transformation_3>::value );
  assert( 3 == CGAL::Feature_dimension<CGAL::Bbox_3>::value );
}

void check_dyn_dim(CGAL::Dynamic_dimension_tag) {}

template < typename K >
void test(CGAL::Dynamic_dimension_tag)
{
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Point_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Hyperplane_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Vector_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Direction_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Line_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Ray_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Segment_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Iso_box_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Sphere_d>::type() );
  check_dyn_dim(typename CGAL::Ambient_dimension<typename K::Aff_transformation_d>::type() );

  assert( 0 == CGAL::Feature_dimension<typename K::Point_d>::value );
  assert( 0 == CGAL::Feature_dimension<typename K::Vector_d>::value );
  assert( 0 == CGAL::Feature_dimension<typename K::Direction_d>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Line_d>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Ray_d>::value );
  assert( 1 == CGAL::Feature_dimension<typename K::Segment_d>::value );
  check_dyn_dim(typename CGAL::Feature_dimension<typename K::Hyperplane_d>::type() );
  check_dyn_dim(typename CGAL::Feature_dimension<typename K::Iso_box_d>::type() );
  check_dyn_dim(typename CGAL::Feature_dimension<typename K::Sphere_d>::type() );
  // check_dyn_dim(typename CGAL::Feature_dimension<typename K::Aff_transformation_d>::type() );
}

int main()
{
  test<CGAL::Cartesian<int> >(CGAL::Dimension_tag<2>());
  test<CGAL::Cartesian<int> >(CGAL::Dimension_tag<3>());
  test<CGAL::Cartesian_d<int> >(CGAL::Dynamic_dimension_tag());

  test<CGAL::Homogeneous<int> >(CGAL::Dimension_tag<2>());
  test<CGAL::Homogeneous<int> >(CGAL::Dimension_tag<3>());
  test<CGAL::Homogeneous_d<int> >(CGAL::Dynamic_dimension_tag());

  test<CGAL::Exact_predicates_inexact_constructions_kernel>(CGAL::Dimension_tag<2>());
  test<CGAL::Exact_predicates_inexact_constructions_kernel>(CGAL::Dimension_tag<3>());

  test<CGAL::Exact_predicates_exact_constructions_kernel>(CGAL::Dimension_tag<2>());
  test<CGAL::Exact_predicates_exact_constructions_kernel>(CGAL::Dimension_tag<3>());

  return 0;
}
