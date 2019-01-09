#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/_test_cls_periodic_3_triangulation_3.h>

#include <CGAL/Periodic_3_regular_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel         Epeck;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<Epeck>    PRTT_Exact;

typedef CGAL::Exact_predicates_inexact_constructions_kernel       Epick;
typedef CGAL::Periodic_3_regular_triangulation_traits_3<Epick>    PRTT_Inexact;
int main(int, char**)
{
  CGAL::Timer t;
  t.start();
  typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Exact>    P3RT3_Exact;
  _test_cls_periodic_3_triangulation_3(P3RT3_Exact(),
                                       PRTT_Exact::Weighted_point_3(0.816982, 0.161518, -0.0942375),
                                       "data/P3RT3_covering_test_HOM.tri",
                                       "data/P3RT3_covering_test.tri",
                                       true);

  typedef CGAL::Periodic_3_regular_triangulation_3<PRTT_Inexact>  P3RT3_Inexact;
  _test_cls_periodic_3_triangulation_3(P3RT3_Inexact(),
                                       PRTT_Inexact::Weighted_point_3(0.816982, 0.161518, -0.0942375),
                                       "data/P3RT3_covering_test_HOM.tri",
                                       "data/P3RT3_covering_test.tri");

  std::cout << t.time() << " sec." << std::endl;
  return 0;
}
