/* Test the aggregate insertion (the surface sweep) of a curve that starts on
 * the left boundary and ends on the top boundary, and of a bounded curve below
 * it that ends before it.
 *
 *          -1
 * o---------o------o
 * |        /       |
 * |-------/        |   y = 1/(-1-x), x < -1
 * |                |
 * |  ----          |   y = -2, -3 <= x <= -2
 * |                |
 * o----------------o
 *
 * The first curve has no curve above it at its left end, so the index of its
 * subcurve is recorded in the list of the top face; this index is mapped to a
 * halfedge only after the curve is inserted. The second curve gets a larger
 * index, which is mapped first. When the first curve is inserted, at the top
 * boundary, it closes a new face, and the relocation of the holes and the
 * isolated vertices into the new face meets the unmapped index of the first
 * curve. This index must not be mistaken for an index of an isolated vertex.
 */

#include <CGAL/config.h>

#include <iostream>

#ifndef CGAL_USE_CORE

int main() {
  std::cout << "NOTE: Core is not installed; the test is skipped." << std::endl;
  return 0;
}

#else

#include <vector>

#include <CGAL/CORE_BigInt.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Arrangement_2.h>

using Number_type = CORE::BigInt;
using AK1 = CGAL::Algebraic_kernel_d_1<Number_type>;
using Traits_2 = CGAL::Arr_rational_function_traits_2<AK1>;
using Curve_2 = Traits_2::Curve_2;
using Polynomial_1 = Traits_2::Polynomial_1;
using Alg_real_1 = Traits_2::Algebraic_real_1;
using Arrangement_2 = CGAL::Arrangement_2<Traits_2>;

/*! checks the sizes of an arrangement.
 */
bool check(const char* name, const Arrangement_2& arr) {
  const char* names[] = {"number of vertices", "number of vertices at infinity", "number of edges",
                         "number of faces", "number of unbounded faces"};
  // y = 1/(-1-x) consists of 2 branches, each with 2 ends at infinity; the segment has 2 vertices.
  const Arrangement_2::Size expected[] = {2, 4, 3, 3, 3};
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
  AK1 ak1;
  Traits_2 traits(&ak1);
  Traits_2::Construct_curve_2 construct = traits.construct_curve_2_object();
  Polynomial_1 x = CGAL::shift(Polynomial_1(1), 1);

  std::vector<Curve_2> curves;
  curves.push_back(construct(Polynomial_1(1), Polynomial_1(-1) - x));       // y = 1/(-1-x)
  curves.push_back(construct(Polynomial_1(-2), Alg_real_1(-3), Alg_real_1(-2)));   // y = -2, -3 <= x <= -2

  // aggregate insertion (the surface sweep)
  Arrangement_2 arr1(&traits);
  CGAL::insert(arr1, curves.begin(), curves.end());
  bool ok = check("aggregate insertion", arr1);

  // incremental insertion, for reference
  Arrangement_2 arr2(&traits);
  for (const auto& cv : curves) CGAL::insert(arr2, cv);
  ok = check("incremental insertion", arr2) && ok;

  std::cout << (ok ? "passed" : "failed") << std::endl;
  return ok ? 0 : 1;
}

#endif
