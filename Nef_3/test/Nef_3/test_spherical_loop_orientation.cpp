// Regression for facet orientation when SNC_simplify converts a spherical
// edge into a spherical loop. All coordinates and operations are exact.
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <array>
#include <iostream>
#include <utility>
#include <vector>

using Kernel = CGAL::Cartesian<CGAL::Gmpq>;
using Point = Kernel::Point_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using Nef = CGAL::Nef_polyhedron_3<Kernel>;

Nef tetrahedron(Point a, Point b, Point c, Point d) {
  if (CGAL::orientation(a, b, c, d) == CGAL::NEGATIVE)
    std::swap(a, b);
  Polyhedron p;
  p.make_tetrahedron(a, b, c, d);
  return Nef(p);
}

Nef convex_hull(const std::vector<Point>& points) {
  Polyhedron p;
  CGAL::convex_hull_3(points.begin(), points.end(), p);
  return Nef(p);
}

Point transform(const Point& p, int variant) {
  if (variant == 1) return Point(p.z() + 2, p.x() - 3, p.y() + 5);
  if (variant == 2) return Point(-p.x() + 2, p.y() - 3, p.z() + 5);
  return p;
}

int run_case(int variant) {
  const Point apex = transform(Point(0, 0, 0), variant);
  std::array<Point, 4> bottom = {
    Point(-4, -4, 0), Point(4, -4, 0),
    Point(4, 4, 0), Point(-4, 4, 0)};
  std::array<Point, 4> top = {
    Point(-2, -2, 2), Point(2, -2, 2),
    Point(2, 2, 2), Point(-2, 2, 2)};
  for (auto& p : bottom) p = transform(p, variant);
  for (auto& p : top) p = transform(p, variant);

  // Independent construction is validated first, including on the unpatched
  // library. The oracle has no internal fan partitions to simplify away.
  std::vector<Point> outer(bottom.begin(), bottom.end());
  outer.insert(outer.end(), top.begin(), top.end());
  std::vector<Point> cavity(top.begin(), top.end());
  cavity.push_back(apex);
  Nef expected = (convex_hull(outer) - convex_hull(cavity)).regularization();
  if (!expected.is_valid()) return 6;
  std::cout << "case " << variant << ": exact oracle valid\n" << std::flush;

  // Four side pyramids of a square frustum, each decomposed into two
  // tetrahedra with disjoint interiors and the same apex. The exact solid
  // is the frustum minus the pyramid from apex to the top square.
  std::array<Nef, 4> side;
  for (int i = 0; i < 4; ++i) {
    const int j = (i + 1) % 4;
    side[i] = (tetrahedron(apex, bottom[i], bottom[j], top[j]) ^
               tetrahedron(apex, bottom[i], top[j], top[i])).regularization();
  }
  const Nef a = (side[0] ^ side[1]).regularization();
  const Nef b = (side[2] ^ side[3]).regularization();
  Nef result = a ^ b;
  if (!result.is_valid()) return 1;
  result.extract_interior();
  if (!result.is_valid()) return 2;
  std::cout << "XOR and interior valid\n" << std::flush;

  // Closure removes the remaining internal faces. At apex the spherical
  // map retains the cavity edges and gains a loop for the bottom plane.
  result.extract_closure();
  if (!result.is_valid()) {
    std::cerr << "closure invalid\n";
    return 3;
  }

  // Independent construction of the same closed point set. In particular,
  // the bottom-plane loop must not connect the material and exterior shells.
  if (result != expected) return 4;
  if (result.number_of_vertices() != 9 || result.number_of_facets() != 9 ||
      result.number_of_volumes() != 2) return 5;
  std::cout << "closure valid; equals exact frustum minus pyramid; "
               "9 vertices, 9 facets, 2 volumes\n";
  return 0;
}

int run_five_tetra_case() {
  const Point apex(0, 0, 0), tip(0, 0, 4);
  const std::array<Point, 4> bottom = {
    Point(-4, -4, 0), Point(4, -4, 0),
    Point(4, 4, 0), Point(-4, 4, 0)};
  const Nef cavity = tetrahedron(apex, Point(-1, -1, 1),
                                 Point(1, -1, 1), Point(0, 1, 1));
  std::array<Nef, 4> fan;
  for (int i = 0; i < 4; ++i)
    fan[i] = tetrahedron(apex, bottom[i], bottom[(i + 1) % 4], tip);
  const Nef a = (fan[0] ^ fan[1]).regularization();
  const Nef b = ((fan[2] ^ fan[3]).regularization() ^ cavity).regularization();
  Nef result = a ^ b;
  if (!result.is_valid()) return 10;
  result.extract_interior();
  if (!result.is_valid()) return 11;
  result.extract_closure();
  if (!result.is_valid()) return 12;

  // The cavity touches the outside at apex. Loop conversion during closure
  // must replace an incident SFace that the facet-removal union-find merged.
  std::vector<Point> outer(bottom.begin(), bottom.end());
  outer.push_back(tip);
  Nef expected = (convex_hull(outer) - cavity).regularization();
  if (!expected.is_valid() || result != expected) return 13;
  if (result.number_of_vertices() != 9 || result.number_of_facets() != 9 ||
      result.number_of_volumes() != 3) return 14;
  std::cout << "five tetrahedra: valid closure; equals exact pyramid minus "
               "tetrahedron; 9 vertices, 9 facets, 3 volumes\n";
  return 0;
}

int main() {
  for (int variant = 0; variant < 3; ++variant)
    if (int error = run_case(variant)) return error;
  return run_five_tetra_case();
}
