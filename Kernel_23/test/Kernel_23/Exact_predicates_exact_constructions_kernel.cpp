#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_kth_root.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>
#endif
#include <CGAL/use.h>

#include "CGAL/_test_cls_kernel.h"
#define TEST_FILENAME "Test-Cartesian-IO.out"
#include "CGAL/_test_io.h"
#include "CGAL/_test_2.h"
#include "CGAL/_test_3.h"

#include "CGAL/_test_new_2.h"
#include "CGAL/_test_new_3.h"

#include "CGAL/_test_fct_points_implicit_sphere.h"
#include "CGAL/_test_orientation_and_bounded_side.h"
#include "CGAL/_test_fct_constructions_2.h"
#include "CGAL/_test_fct_constructions_3.h"
#include "CGAL/_test_fct_point_3.h"
#include "CGAL/_test_fct_coplanar_3.h"
#include "CGAL/_test_cls_iso_cuboid_3.h"
#include "CGAL/_test_angle.h"
#include "CGAL/_test_cls_circle_3.h"

#include "CGAL/_test_mf_plane_3_to_2d.h"

template <typename K>
void test()
{
  CGAL_USE_TYPE(K);

  std::cout << "Testing nested types with:" << typeid(K).name() << std::endl;
  _test_kernel( K() );

  std::cout << "Testing IO with:" << typeid(K).name() << std::endl;
  _test_io( K() );

  std::cout << "Testing 2d with:" << typeid(K).name() << std::endl;
  _test_2( K() );

  std::cout << "Testing 3d with:" << typeid(K).name() << std::endl;
  _test_3( K() );
  _test_cls_circle_3( K() );

  std::cout << "Testing new 2d with:" << typeid(K).name() << std::endl;
  test_new_2( K() );
  _test_cls_new_2( K() );

  std::cout << "Testing new 3d with:" << typeid(K).name() << std::endl;
  test_new_3( K() );

  std::cout << "Testing new parts with:" << typeid(K).name() << std::endl;
  _test_orientation_and_bounded_side( K() );
  _test_fct_points_implicit_sphere( K() );
  _test_fct_constructions_2( K() );
  _test_fct_constructions_3( K() );
  _test_fct_point_3( K() );
  _test_fct_coplanar_3( K() );
  _test_cls_iso_cuboid_3( K() );
  _test_angle( K() );

  std::cout << "Testing 3d-2d with:" << typeid(K).name() << std::endl;
  _test_mf_plane_3_to_2d( K() );

  std::cout << "All tests done" << std::endl;
}

int main()
{
  test<CGAL::Exact_predicates_exact_constructions_kernel>();

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
  test<CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt>();
  test<CGAL::Exact_predicates_exact_constructions_kernel_with_kth_root>();
  test<CGAL::Exact_predicates_exact_constructions_kernel_with_root_of>();
#endif

  return 0;
}
