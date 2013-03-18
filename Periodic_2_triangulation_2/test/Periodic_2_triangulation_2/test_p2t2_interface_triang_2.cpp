#include "types.h"
#include "interface_test.h"
#include <CGAL/Periodic_2_triangulation_2.h>

int main() {
  typedef Periodic_2_triangulation_2<Gt>              P2T2;
  typedef Periodic_2_Delaunay_triangulation_2<Gt>     DP2T2;

  test<P2T2>();
  test<DP2T2>();
  test<PTH_Dt>();
  test<Delaunay_triangulation_hierarchy>();

  test_nearest<DP2T2>();
  test_nearest<PTH_Dt>();
  test_nearest<Delaunay_triangulation_hierarchy>();

  return 0;
}
