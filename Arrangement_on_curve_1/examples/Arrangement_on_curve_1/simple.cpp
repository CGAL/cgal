#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1_functions.h>

// Minimal geometry traits model utilizing direct coordinate FT parameters
template <typename Kernel>
class Aoc_ft_traits_1 {
public:
  using Point_1 = typename Kernel::FT;

  class Compare_x_1 {
  public:
    CGAL::Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      if (p1 < p2) return CGAL::SMALLER;
      if (p1 > p2) return CGAL::LARGER;
      return CGAL::EQUAL;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(); }
};

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = Aoc_ft_traits_1<Kernel>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits, std::string, void>;

int main() {
  Traits traits;
  Arrangement arr(traits);

  std::cout << "Inserting raw coordinate fields into 1D track...\n";
  auto v1 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(10.5));
  auto v2 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(2.0));
  auto v3 = CGAL::Arrangement_on_curve_1::insert(arr, Kernel::FT(5.25));

  // Modify user-extended data fields on vertices
  v1->data() = "Right Node";
  v2->data() = "Left Node";
  v3->data() = "Center Node";

  std::cout << "\nTraversing generated topological components:\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "  Vertex coordinate: " << vit->point() << " | Label: " << vit->data() << "\n";
  }

  return 0;
}
