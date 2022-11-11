#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "include/utils.h"
#include "include/wrappers.h"

// Typedefs.
using SCKER = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename Kernel>
bool test_kernel() {
  const wrappers::Tangent_wrapper<Kernel> tan;
  const wrappers::Mean_value_wrapper<Kernel> mvw;
  return tests::test_analytic_weight<Kernel>(tan, mvw);
}

int main() {
  assert(test_kernel<SCKER>());
  assert(test_kernel<EPICK>());
  assert(test_kernel<EPECK>());
  std::cout << "* test_tangent_weights: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
