#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Ft_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Geom_traits = CGAL::Arrangement_on_curve_1::Ft_traits_1<Kernel>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Geom_traits>;

int main() {
  auto traits_ptr = std::make_shared<const Geom_traits>();
  Arrangement arr(traits_ptr);
  auto v1 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(10.5));
  auto v2 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(2.0));
  auto v3 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(5.25));
  return 0;
}
