/* Test the aggregate insertion (the surface sweep) of overlapping vertical
 * linear curves together with a line that passes through the endpoint of one
 * of them:
 *
 *   the line y = 2x + 6, which passes through (-4,-2);
 *   the vertical ray that emanates upwards from (-4,-2);
 *   the vertical line x = -4; and
 *   the vertical segment from (-4,-5) to (-4,5).
 *
 * The sweep computes the overlap of curves with common ancestors, and asserts
 * that the overlap starts at the current event. The event comparer obtains the
 * location of the end of the overlap from boundary conditions that must be set
 * beforehand; they used to be left over from a previous comparison (here, of
 * the end of the vertical line on the bottom boundary), and the comparison
 * failed with an error in Compare_x_on_boundary_2.
 */

#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits_2 = CGAL::Arr_linear_traits_2<Kernel>;
using Point_2 = Kernel::Point_2;
using X_monotone_curve_2 = Traits_2::X_monotone_curve_2;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2>;

/*! checks the sizes of an arrangement.
 */
bool check(const char* name, const Arrangement_2& arr) {
  const char* names[] = {"number of vertices", "number of vertices at infinity", "number of edges",
                         "number of faces", "number of unbounded faces"};
  const Arrangement_2::Size expected[] = {3, 4, 6, 4, 4};
  const Arrangement_2::Size sizes[] = {arr.number_of_vertices(), arr.number_of_vertices_at_infinity(),
                                       arr.number_of_edges(), arr.number_of_faces(),
                                       arr.number_of_unbounded_faces()};
  bool ok = arr.is_valid();
  if (! ok) std::cerr << name << ": the arrangement is not valid!" << std::endl;
  for (int i = 0; i < 5; ++i) {
    if (sizes[i] == expected[i]) continue;
    std::cerr << name << ": " << names[i] << ": " << sizes[i] << ", expected: " << expected[i] << std::endl;
    ok = false;
  }
  return ok;
}

int main() {
  std::vector<X_monotone_curve_2> curves;
  curves.push_back(Kernel::Line_2(Point_2(-5, -4), Point_2(-1, 4)));
  curves.push_back(Kernel::Ray_2(Point_2(-4, -2), Point_2(-4, 0)));
  curves.push_back(Kernel::Line_2(Point_2(-4, -1), Point_2(-4, 1)));
  curves.push_back(Kernel::Segment_2(Point_2(-4, -5), Point_2(-4, 5)));

  // aggregate insertion (the surface sweep)
  Arrangement_2 arr1;
  CGAL::insert(arr1, curves.begin(), curves.end());
  bool ok = check("aggregate insertion", arr1);

  // incremental insertion, for reference
  Arrangement_2 arr2;
  for (const auto& cv : curves) CGAL::insert(arr2, cv);
  ok = check("incremental insertion", arr2) && ok;

  std::cout << (ok ? "passed" : "failed") << std::endl;
  return ok ? 0 : 1;
}
