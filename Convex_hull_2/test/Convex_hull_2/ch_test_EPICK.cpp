#include <CGAL/_test_fct_ch_I_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/convex_hull_constructive_traits_2.h>

int main()
{
  CGAL::Exact_predicates_inexact_constructions_kernel ch_EPICK;
  std::cout << "EPICK:" << std::endl;
  CGAL::ch__batch_test(ch_EPICK);

  CGAL::Convex_hull_constructive_traits_2<CGAL::Exact_predicates_inexact_constructions_kernel> cch_EPICK;
  std::cout << "Constructive EPICK:" << std::endl;
  CGAL::ch__batch_test(cch_EPICK);

  return EXIT_SUCCESS;
}
