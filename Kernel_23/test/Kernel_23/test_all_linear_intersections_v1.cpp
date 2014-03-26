#define CGAL_INTERSECTION_VERSION 1

#include <CGAL/_test_all_linear_intersections.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


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
