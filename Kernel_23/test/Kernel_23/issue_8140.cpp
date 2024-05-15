#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

int main() { // For testing, we just print the combinations of types
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Vec_3 = Kernel::Vector_3;
  using Cp_fnc3 = Vec_3(*)(const Vec_3&, const Vec_3&);

  Cp_fnc3 f = static_cast<Cp_fnc3>(&CGAL::cross_product<Kernel>);
  Vec_3 v1, v2;
  Vec_3 v = f(v1, v2);
  return 0;
}
