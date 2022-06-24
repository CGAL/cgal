#ifndef CGAL_WEIGHTS_TESTS_UTILS_H
#define CGAL_WEIGHTS_TESTS_UTILS_H

// STL includes.
#include <array>
#include <vector>
#include <cassert>
#include <algorithm>
#include <string>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Weights/utils.h>

namespace tests {

template<typename FT>
FT get_tolerance() {
  return FT(1) / FT(10000000000);
}

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_all_triangles() {

  using Point_2 = typename Kernel::Point_2;
  const std::array<Point_2, 3> triangle0 = {
    Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0)
  };
  const std::array<Point_2, 3> triangle1 = {
    Point_2(-2, 0), Point_2(0, -1), Point_2(2, 0)
  };
  const std::array<Point_2, 3> triangle2 = {
    Point_2(-2, 0), Point_2(-2, -2), Point_2(2, 0)
  };
  const std::array<Point_2, 3> triangle3 = {
    Point_2(-2, 0), Point_2(2, -2), Point_2(2, 0)
  };
  return { triangle0, triangle1, triangle2, triangle3 };
}

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_symmetric_triangles() {

  using Point_2 = typename Kernel::Point_2;
  const std::array<Point_2, 3> triangle0 = {
    Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0)
  };
  const std::array<Point_2, 3> triangle1 = {
    Point_2(-2, 0), Point_2(0, -1), Point_2(2, 0)
  };
  const std::array<Point_2, 3> triangle2 = {
    Point_2(-3, 0), Point_2(0, -1), Point_2(3, 0)
  };
  return { triangle0, triangle1, triangle2 };
}

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_uniform_triangles() {

  using Point_2 = typename Kernel::Point_2;
  const std::array<Point_2, 3> triangle0 = {
    Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0)
  };
  const std::array<Point_2, 3> triangle1 = {
    Point_2(-2, 0), Point_2(0, -2), Point_2(2, 0)
  };
  const std::array<Point_2, 3> triangle2 = {
    Point_2(1, 0), Point_2(-1, 0), Point_2(-1, -2)
  };
  const std::array<Point_2, 3> triangle3 = {
    Point_2(1, -2), Point_2(1, 0), Point_2(-1, 0)
  };
  return { triangle0, triangle1, triangle2, triangle3 };
}

template<typename Kernel>
std::vector< std::vector<typename Kernel::Point_2> >
get_all_polygons() {

  using Point_2 = typename Kernel::Point_2;
  const std::vector<Point_2> polygon0 = {
    Point_2(-2, -2), Point_2(2, -2), Point_2(0, 2)
  };
  const std::vector<Point_2> polygon1 = {
    Point_2(-1, -1), Point_2(1, -1), Point_2(1, 1), Point_2(-1, 1)
  };
  const std::vector<Point_2> polygon2 = {
    Point_2(-2, 0), Point_2(0, -2), Point_2(2, 0), Point_2(0, 2)
  };
  const std::vector<Point_2> polygon3 = {
    Point_2(-2, -2), Point_2(2, -2), Point_2(2, 0), Point_2(0, 2), Point_2(-2, 0)
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
  assert(a2 >= FT(0) && b2 >= FT(0));
  if (a2 < FT(0) || b2 < FT(0)) return false;
  assert(CGAL::abs(a2 - b2) < tol);
  if (CGAL::abs(a2 - b2) >= tol) return false;

  if (wrapper.supports_3d()) {
    const auto a3 = wrapper.weight_a(t3, r3, p3, q3);
    const auto b3 = wrapper.weight_b(t3, r3, p3, q3);
    assert(a3 >= FT(0) && b3 >= FT(0));
    if (a3 < FT(0) || b3 < FT(0)) return false;
    assert(CGAL::abs(a3 - b3) < tol);
    if (CGAL::abs(a3 - b3) >= tol) return false;
    assert(CGAL::abs(a2 - a3) < tol);
    assert(CGAL::abs(b2 - b3) < tol);
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
  assert(a2 >= FT(0) && b2 >= FT(0));
  if (a2 < FT(0) || b2 < FT(0)) return false;
  assert(CGAL::abs(a2 - b2) < tol);
  if (CGAL::abs(a2 - b2) >= tol) return false;

  if (wrapper.supports_3d()) {
    const auto a3 = wrapper.weight_a(t3, r3, p3, Point_3(-x, 0, 1));
    const auto b3 = wrapper.weight_a(t3, r3, p3, Point_3(+x, 0, 1));
    assert(a3 >= FT(0) && b3 >= FT(0));
    if (a3 < FT(0) || b3 < FT(0)) return false;
    assert(CGAL::abs(a3 - b3) < tol);
    if (CGAL::abs(a3 - b3) >= tol) return false;
    assert(CGAL::abs(a2 - a3) < tol);
    assert(CGAL::abs(b2 - b3) < tol);
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
  assert(a2 >= FT(0) && b2 >= FT(0));
  if (a2 < FT(0) || b2 < FT(0)) return false;
  assert(CGAL::abs(a2 - b2) < tol);
  if (CGAL::abs(a2 - b2) >= tol) return false;

  if (wrapper1.supports_3d() && wrapper2.supports_3d()) {
    const auto a3 = wrapper1.weight_a(t3, r3, p3, q3);
    const auto b3 = wrapper2.weight_a(t3, r3, p3, q3);
    assert(a3 >= FT(0) && b3 >= FT(0));
    if (a3 < FT(0) || b3 < FT(0)) return false;
    assert(CGAL::abs(a3 - b3) < tol);
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

  const auto a2 = wrapper.weight(p2, q2, r2);
  const auto a3 = wrapper.weight(p3, q3, r3);
  assert(a2 >= FT(0) && a3 >= FT(0));
  if (a2 < FT(0) || a3 < FT(0)) return false;
  assert(CGAL::abs(a2 - a3) < tol);
  if (CGAL::abs(a2 - a3) >= tol) return false;
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_area(
  const Weight_wrapper& wrapper,
  const std::array<typename Kernel::Point_2, 3>& neighbors) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  // 2D configuration.
  const Point_2& p2 = neighbors[0];
  const Point_2& q2 = neighbors[1];
  const Point_2& r2 = neighbors[2];

  // 3D configuration.
  const Point_3 p3(p2.x(), p2.y(), 1);
  const Point_3 q3(q2.x(), q2.y(), 1);
  const Point_3 r3(r2.x(), r2.y(), 1);

  const auto a2 = wrapper.weight(p2, q2, r2);
  const auto a3 = wrapper.weight(p3, q3, r3);
  assert(a2 <= CGAL::Weights::area(p2, q2, r2));
  assert(a3 <= CGAL::Weights::area(p3, q3, r3));
  if (a2 > CGAL::Weights::area(p2, q2, r2)) return false;
  if (a3 > CGAL::Weights::area(p3, q3, r3)) return false;
  assert(a2 >= FT(0));
  assert(a3 >= FT(0));
  if (a2 < FT(0)) return false;
  if (a3 < FT(0)) return false;
  return true;
}

template<typename FT, typename Point>
bool test_coordinates(
  const Point& query,
  const std::vector<Point>& polygon,
  const std::vector<FT>& weights) {

  assert(weights.size() > 0);
  if (weights.size() == 0) return false;

  // Compute the sum of weights.
  const FT tol = get_tolerance<FT>();
  FT sum = FT(0);
  for (const FT& weight : weights) {
    sum += weight;
  }
  assert(sum >= tol);
  if (sum < tol) return false;

  // Compute coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(weights.size());
  for (const FT& weight : weights) {
    coordinates.push_back(weight / sum);
  }
  assert(coordinates.size() == weights.size());
  if (coordinates.size() != weights.size()) return false;

  // Test partition of unity.
  sum = FT(0);
  for (const FT& coordinate : coordinates) {
    sum += coordinate;
  }
  assert(CGAL::abs(FT(1) - sum) < tol);
  if (CGAL::abs(FT(1) - sum) >= tol) return false;

  // Test linear precision.
  FT x = FT(0), y = FT(0);
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    x += coordinates[i] * polygon[i].x();
    y += coordinates[i] * polygon[i].y();
  }
  assert(CGAL::abs(query.x() - x) < tol);
  assert(CGAL::abs(query.y() - y) < tol);
  if (CGAL::abs(query.x() - x) >= tol) return false;
  if (CGAL::abs(query.y() - y) >= tol) return false;
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_on_polygon(
  const Weight_wrapper& wrapper,
  const typename Kernel::Point_2& query_2,
  const std::vector<typename Kernel::Point_2>& polygon_2) {

  // Get weights.
  using FT = typename Kernel::FT;
  assert(polygon_2.size() >= 3);
  if (polygon_2.size() < 3) return false;

  // 2D version.
  std::vector<FT> weights_2;
  weights_2.reserve(polygon_2.size());
  wrapper.compute_on_polygon(
    polygon_2, query_2, Kernel(), std::back_inserter(weights_2));
  assert(weights_2.size() == polygon_2.size());
  if (weights_2.size() != polygon_2.size()) return false;
  if (!test_coordinates(query_2, polygon_2, weights_2)) return false;

  // 3D version.
  using Point_3 = typename Kernel::Point_3;
  const Point_3 query_3(query_2.x(), query_2.y(), 1);
  std::vector<Point_3> polygon_3;
  polygon_3.reserve(polygon_2.size());
  for (const auto& vertex_2 : polygon_2) {
    polygon_3.push_back(Point_3(vertex_2.x(), vertex_2.y(), 1));
  }
  assert(polygon_3.size() == polygon_2.size());
  if (polygon_3.size() != polygon_2.size()) return false;
  const CGAL::Projection_traits_xy_3<Kernel> ptraits;

  std::vector<FT> weights_3;
  weights_3.reserve(polygon_3.size());
  wrapper.compute_on_polygon(
    polygon_3, query_3, ptraits, std::back_inserter(weights_3));
  assert(weights_3.size() == polygon_3.size());
  if (weights_3.size() != polygon_3.size()) return false;
  if (!test_coordinates(query_3, polygon_3, weights_3)) return false;
  return true;
}

template<
typename Kernel,
typename Weight_wrapper>
bool test_barycentric_properties(
  const Weight_wrapper& wrapper,
  const typename Kernel::Point_2& query,
  const std::vector<typename Kernel::Point_2>& polygon) {

  // Get weights.
  using FT = typename Kernel::FT;
  const std::size_t n = polygon.size();
  assert(n >= 3);
  if (n < 3) return false;

  // Check properties.
  std::vector<FT> weights;
  weights.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t im = (i + n - 1) % n;
    const std::size_t ip = (i + 1) % n;
    const auto& t = polygon[im];
    const auto& r = polygon[i];
    const auto& p = polygon[ip];
    const auto& q = query;
    const FT weight = wrapper.weight_a(t, r, p, q);
    weights.push_back(weight);
  }
  assert(weights.size() == n);
  if (weights.size() != n) return false;
  if (!test_coordinates(query, polygon, weights)) return false;
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

  // Test query points.
  auto configs = get_all_triangles<Kernel>();
  for (const auto& config : configs) {
    if (!test_query<Kernel>(weight, zero, config)) return false;
    for (const auto& query : queries) {
      if (!test_query<Kernel>(weight, query, config)) return false;
    }
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

  // Test symmetry along x axis.
  configs = get_symmetric_triangles<Kernel>();
  for (const auto& config : configs) {
    if (!test_symmetry_x<Kernel>(weight, config, q)) return false;
    if (!test_symmetry_x<Kernel>(weight, config, h)) return false;
    if (!test_symmetry_x<Kernel>(weight, config, t)) return false;
  }

  // Test barycentric properties.
  if (weight.is_barycentric()) {
    const auto polygons = get_all_polygons<Kernel>();
    for (const auto& polygon : polygons) {
      if (!test_barycentric_properties<Kernel>(weight, zero, polygon)) {
        return false;
      }
      for (const auto& query : queries) {
        if (!test_barycentric_properties<Kernel>(weight, query, polygon)) {
          return false;
        }
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
    Point_2(-h, +q), Point_2(+h, -q), Point_2(-q, +h), Point_2(+q, -h)
  };

  // Test analytic formulations.
  if (!test_analytic_weight<Kernel>(weight, alternative)) {
    return false;
  }

  // Test on polygons.
  const auto polygons = get_all_polygons<Kernel>();
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
  auto configs = get_all_triangles<Kernel>();
  for (const auto& config : configs) {
    if (!test_neighbors<Kernel>(weight, config)) return false;
  }

  // Test areas.
  configs = get_uniform_triangles<Kernel>();
  for (const auto& config : configs) {
    if (!test_area<Kernel>(weight, config)) return false;
  }
  return true;
}

} // namespace tests

#endif // CGAL_WEIGHTS_TESTS_UTILS_H
