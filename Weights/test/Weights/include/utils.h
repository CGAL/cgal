#ifndef CGAL_WEIGHTS_TESTS_UTILS_H
#define CGAL_WEIGHTS_TESTS_UTILS_H

// STL includes.
#include <array>
#include <vector>
#include <cassert>
#include <string>

// CGAL includes.
#include <CGAL/assertions.h>

namespace tests {

template<typename FT>
FT get_tolerance() {
  return FT(1) / FT(10000000000);
}

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_default_triangles() {

  using Point_2 = typename Kernel::Point_2;
  const std::array<Point_2, 3> triangle0 = {
    Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0)
  };
  const std::array<Point_2, 3> triangle1 = {
    Point_2(-2, 0), Point_2(0, -1), Point_2(2, 0)
  };
  const std::array<Point_2, 3> triangle2 = {
    Point_2(-1, 0), Point_2(-1, -2), Point_2(1, 0)
  };
  const std::array<Point_2, 3> triangle3 = {
    Point_2(-1, 0), Point_2(1, -2), Point_2(1, 0)
  };
  return { triangle0, triangle1, triangle2, triangle3 };
}

template<typename Kernel>
std::vector< std::vector<typename Kernel::Point_2> >
get_default_polygons() {

  using Point_2 = typename Kernel::Point_2;
  const std::vector<Point_2> polygon0 = {
    Point_2(-1, 0), Point_2(1, 0), Point_2(0, 1)
  };
  const std::vector<Point_2> polygon1 = {
    Point_2(0, 0), Point_2(1, 0), Point_2(1, 1), Point_2(0, 1)
  };
  const std::vector<Point_2> polygon2 = {
    Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0), Point_2(0, 1)
  };
  const std::vector<Point_2> polygon3 = {
    Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0), Point_2(1, 1), Point_2(-1, 1)
  };
  return { polygon0, polygon1, polygon2, polygon3 };
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_query(
  const Weight_wrapper& wrapper,
  const typename Kernel::Point_2& query,
  const std::array<typename Kernel::Point_2, 3>& neighbors) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const FT tol = get_tolerance<FT>();

  // 2D configuration.
  const Point_2& t2 = neighbors[0];
  const Point_2& r2 = neighbors[1];
  const Point_2& p2 = neighbors[2];
  const Point_2& q2 = query;

  // 3D configuration.
  const Point_3 t3(t2.x(), t2.y(), 1);
  const Point_3 r3(r2.x(), r2.y(), 1);
  const Point_3 p3(p2.x(), p2.y(), 1);
  const Point_3 q3(q2.x(), q2.y(), 1);

  const auto a2 = wrapper.weight_a(t2, r2, p2, q2);
  const auto b2 = wrapper.weight_b(t2, r2, p2, q2);
  CGAL_assertion(a2 >= FT(0) && b2 >= FT(0));
  if (a2 < FT(0) || b2 < FT(0)) return false;
  CGAL_assertion(CGAL::abs(a2 - b2) < tol);
  if (CGAL::abs(a2 - b2) >= tol) return false;

  if (wrapper.supports_3d()) {
    const auto a3 = wrapper.weight_a(t3, r3, p3, q3);
    const auto b3 = wrapper.weight_b(t3, r3, p3, q3);
    CGAL_assertion(a3 >= FT(0) && b3 >= FT(0));
    if (a3 < FT(0) || b3 < FT(0)) return false;
    CGAL_assertion(CGAL::abs(a3 - b3) < tol);
    if (CGAL::abs(a3 - b3) >= tol) return false;
    CGAL_assertion(CGAL::abs(a2 - a3) < tol);
    CGAL_assertion(CGAL::abs(b2 - b3) < tol);
    if (CGAL::abs(a2 - a3) >= tol) return false;
    if (CGAL::abs(b2 - b3) >= tol) return false;
  }
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_symmetry_x(
  const Weight_wrapper& wrapper,
  const std::array<typename Kernel::Point_2, 3>& neighbors,
  const typename Kernel::FT& x) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const FT tol = get_tolerance<FT>();

  // 2D configuration.
  const Point_2& t2 = neighbors[0];
  const Point_2& r2 = neighbors[1];
  const Point_2& p2 = neighbors[2];

  // 3D configuration.
  const Point_3 t3(t2.x(), t2.y(), 1);
  const Point_3 r3(r2.x(), r2.y(), 1);
  const Point_3 p3(p2.x(), p2.y(), 1);

  const auto a2 = wrapper.weight_a(t2, r2, p2, Point_2(-x, 0));
  const auto b2 = wrapper.weight_a(t2, r2, p2, Point_2(+x, 0));
  CGAL_assertion(a2 >= FT(0) && b2 >= FT(0));
  if (a2 < FT(0) || b2 < FT(0)) return false;
  CGAL_assertion(CGAL::abs(a2 - b2) < tol);
  if (CGAL::abs(a2 - b2) >= tol) return false;

  if (wrapper.supports_3d()) {
    const auto a3 = wrapper.weight_a(t3, r3, p3, Point_3(-x, 0, 1));
    const auto b3 = wrapper.weight_a(t3, r3, p3, Point_3(+x, 0, 1));
    CGAL_assertion(a3 >= FT(0) && b3 >= FT(0));
    if (a3 < FT(0) || b3 < FT(0)) return false;
    CGAL_assertion(CGAL::abs(a3 - b3) < tol);
    if (CGAL::abs(a3 - b3) >= tol) return false;
    CGAL_assertion(CGAL::abs(a2 - a3) < tol);
    CGAL_assertion(CGAL::abs(b2 - b3) < tol);
    if (CGAL::abs(a2 - a3) >= tol) return false;
    if (CGAL::abs(b2 - b3) >= tol) return false;
  }
  return true;
}

template<
typename Kernel,
typename Weight_wrapper_1,
typename Weight_wrapper_2>
bool test_compare(
  const Weight_wrapper_1& wrapper1,
  const Weight_wrapper_2& wrapper2,
  const typename Kernel::Point_2& query,
  const std::array<typename Kernel::Point_2, 3>& neighbors) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const FT tol = get_tolerance<FT>();

  // 2D configuration.
  const Point_2& t2 = neighbors[0];
  const Point_2& r2 = neighbors[1];
  const Point_2& p2 = neighbors[2];
  const Point_2& q2 = query;

  // 3D configuration.
  const Point_3 t3(t2.x(), t2.y(), 1);
  const Point_3 r3(r2.x(), r2.y(), 1);
  const Point_3 p3(p2.x(), p2.y(), 1);
  const Point_3 q3(q2.x(), q2.y(), 1);

  const auto a2 = wrapper1.weight_a(t2, r2, p2, q2);
  const auto b2 = wrapper2.weight_a(t2, r2, p2, q2);
  CGAL_assertion(a2 >= FT(0) && b2 >= FT(0));
  if (a2 < FT(0) || b2 < FT(0)) return false;
  CGAL_assertion(CGAL::abs(a2 - b2) < tol);
  if (CGAL::abs(a2 - b2) >= tol) return false;

  if (wrapper1.supports_3d() && wrapper2.supports_3d()) {
    const auto a3 = wrapper1.weight_a(t3, r3, p3, q3);
    const auto b3 = wrapper2.weight_a(t3, r3, p3, q3);
    CGAL_assertion(a3 >= FT(0) && b3 >= FT(0));
    if (a3 < FT(0) || b3 < FT(0)) return false;
    CGAL_assertion(CGAL::abs(a3 - b3) < tol);
    if (CGAL::abs(a3 - b3) >= tol) return false;
  }
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_neighbors(
  const Weight_wrapper& wrapper,
  const std::array<typename Kernel::Point_2, 3>& neighbors) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const FT tol = get_tolerance<FT>();

  // 2D configuration.
  const Point_2& p2 = neighbors[0];
  const Point_2& q2 = neighbors[1];
  const Point_2& r2 = neighbors[2];

  // 3D configuration.
  const Point_3 p3(p2.x(), p2.y(), 1);
  const Point_3 q3(q2.x(), q2.y(), 1);
  const Point_3 r3(r2.x(), r2.y(), 1);

  const auto a2 = wrapper.weight_a(p2, q2, r2);
  const auto a3 = wrapper.weight_a(p3, q3, r3);
  CGAL_assertion(a2 <= CGAL::area(p2, q2, r2));
  if (a2 > CGAL::area(p2, q2, r2)) return false;
  CGAL_assertion(a2 >= FT(0) && a3 >= FT(0));
  if (a2 < FT(0) || a3 < FT(0)) return false;
  CGAL_assertion(CGAL::abs(a2 - a3) < tol);
  if (CGAL::abs(a2 - a3) >= tol) return false;
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_on_polygon(
  const Weight_wrapper& wrapper,
  const typename Kernel::Point_2& query,
  const std::vector<typename Kernel::Point_2>& polygon) {

  // Get weights.
  using FT = typename Kernel::FT;
  CGAL_assertion(polygon.size() >= 3);
  if (polygon.size() < 3) return false;
  const FT tol = get_tolerance<FT>();

  std::vector<FT> weights;
  weights.reserve(polygon.size());
  wrapper.compute_on_polygon(
    polygon, query, std::back_inserter(weights));
  CGAL_assertion(weights.size() == polygon.size());
  if (weights.size() != polygon.size()) return false;

  // Compute the sum of weights.
  FT sum = FT(0);
  for (const auto& weight : weights) {
    sum += weight;
  }
  CGAL_assertion(sum >= tol);
  if (sum < tol) return false;

  // Compute coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(weights.size());
  for (const auto& weight : weights) {
    coordinates.push_back(weight / sum);
  }
  CGAL_assertion(coordinates.size() == weights.size());
  if (coordinates.size() != weights.size()) return false;

  // Test partition of unity.
  sum = FT(0);
  for (const auto& coordinate : coordinates) {
    sum += coordinate;
  }
  CGAL_assertion(CGAL::abs(FT(1) - sum) < tol);
  if (CGAL::abs(FT(1) - sum) >= tol) return false;

  // Test linear precision.
  FT x = FT(0), y = FT(0);
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    x += coordinates[i] * polygon[i].x();
    y += coordinates[i] * polygon[i].y();
  }
  CGAL_assertion(CGAL::abs(query.x() - x) < tol);
  CGAL_assertion(CGAL::abs(query.y() - y) < tol);
  if (CGAL::abs(query.x() - x) >= tol) return false;
  if (CGAL::abs(query.y() - y) >= tol) return false;
  return true;
}

template<
typename Kernel,
typename Weight_wrapper_1,
typename Weight_wrapper_2>
bool test_analytic_weight(
  const Weight_wrapper_1& weight,
  const Weight_wrapper_2& alternative) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  // Data.
  const FT q = FT(1) / FT(4);
  const FT h = FT(1) / FT(2);
  const FT t = FT(3) / FT(4);

  const Point_2 zero(0, 0);
  const std::vector<Point_2> queries = {
    Point_2(-q, 0), Point_2(+q, 0), Point_2(0, -q), Point_2(0, +q),
    Point_2(-h, 0), Point_2(+h, 0), Point_2(0, -h), Point_2(0, +h),
    Point_2(-t, 0), Point_2(+t, 0), Point_2(0, -t), Point_2(0, +t)
  };
  const auto configs = get_default_triangles<Kernel>();

  // Test query points.
  for (const auto& config : configs) {
    if (!test_query<Kernel>(weight, zero, config)) return false;
    for (const auto& query : queries) {
      if (!test_query<Kernel>(weight, query, config)) return false;
    }
  }

  // Test symmetry along x axis.
  for (const auto& config : configs) {
    if (!test_symmetry_x<Kernel>(weight, config, q)) return false;
    if (!test_symmetry_x<Kernel>(weight, config, h)) return false;
    if (!test_symmetry_x<Kernel>(weight, config, t)) return false;
  }

  // Test alternative formulations.
  for (const auto& config : configs) {
    if (!test_compare<Kernel>(weight, alternative, zero, config)) {
      return false;
    }
    for (const auto& query : queries) {
      if (!test_compare<Kernel>(weight, alternative, query, config)) {
        return false;
      }
    }
  }
  return true;
}

template<
typename Kernel,
typename Weight_wrapper_1,
typename Weight_wrapper_2>
bool test_barycentric_weight(
  const Weight_wrapper_1& weight,
  const Weight_wrapper_2& alternative) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  // Data.
  const FT q = FT(1) / FT(4);
  const FT h = FT(1) / FT(2);
  const Point_2 zero(0, 0);
  const std::vector<Point_2> queries = {
    Point_2(-h,  0), Point_2(+h,  0), Point_2(-q,  0), Point_2(+q,  0),
    Point_2( 0, -h), Point_2( 0, +h), Point_2( 0, -q), Point_2( 0, +q),
    Point_2(-h, -h), Point_2(+h, +h), Point_2(-q, -q), Point_2(+q, +q),
    Point_2(-h, +q), Point_2(+h, -q), Point_2(-q, +h), Point_2(+q, -h),
    Point_2(-1,  0), Point_2( 0, -1), Point_2( 1,  0), Point_2( 0,  1),
    Point_2(-1, -1), Point_2( 1, -1), Point_2( 1,  1), Point_2(-1,  1)
  };

  // Test analytic formulations.
  if (!test_analytic_weight<Kernel>(weight, alternative)) {
    return false;
  }

  // Test on polygons.
  const auto polygons = get_default_polygons<Kernel>();
  for (const auto& polygon : polygons) {
    if (!test_on_polygon<Kernel>(weight, zero, polygon)) return false;
    for (const auto& query : queries) {
      if (!test_on_polygon<Kernel>(weight, query, polygon)) {
        return false;
      }
    }
  }
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_region_weight(const Weight_wrapper& weight) {

  // Test neighborhoods.
  const auto configs = get_default_triangles<Kernel>();
  for (const auto& config : configs) {
    if (!test_neighbors<Kernel>(weight, config)) return false;
  }
  return true;
}

} // namespace tests

#endif // CGAL_WEIGHTS_TESTS_UTILS_H
