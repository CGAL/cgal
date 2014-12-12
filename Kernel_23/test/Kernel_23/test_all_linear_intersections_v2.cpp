#ifdef CGAL_INTERSECTION_VERSION
#undef CGAL_INTERSECTION_VERSION
#endif
#define CGAL_INTERSECTION_VERSION 2

#include <CGAL/_test_all_linear_intersections.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// test that there is no conflict with overloads defined in BSO_2 package
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K1;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K2;
typedef CGAL::Simple_cartesian<double> K3;

int main()
{
  test_linear_intersections<K1,K2>();
  test_linear_intersections<K2,K1>();
  test_linear_intersections<K1,K3>();
  test_linear_intersections<K3,K1>();
  test_linear_intersections<K3,K2>();
  test_linear_intersections<K2,K3>();
}
