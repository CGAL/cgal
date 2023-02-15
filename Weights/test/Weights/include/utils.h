#ifndef CGAL_WEIGHTS_TESTS_UTILS_H
#define CGAL_WEIGHTS_TESTS_UTILS_H

#include <CGAL/Weights/utils.h>

#include <CGAL/assertions.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <array>
#include <vector>
#include <cassert>
#include <algorithm>
#include <string>

namespace CGAL {
namespace Weights {
namespace internal {

template<typename Kernel>
typename Kernel::FT squared_distance(const CGAL::Point_2<Kernel>& p,
                                     const CGAL::Point_2<Kernel>& q)
{
  const Kernel traits;
  auto squared_distance_2 = traits.compute_squared_distance_2_object();
  return squared_distance_2(p, q);
}

template<typename Kernel>
typename Kernel::FT squared_distance(const CGAL::Point_3<Kernel>& p,
                                     const CGAL::Point_3<Kernel>& q)
{
  const Kernel traits;
  auto squared_distance_3 = traits.compute_squared_distance_3_object();
  return squared_distance_3(p, q);
}

template<typename Kernel>
typename Kernel::FT distance(const CGAL::Point_2<Kernel>& p,
                             const CGAL::Point_2<Kernel>& q)
{
  const Kernel traits;
  return CGAL::Weights::internal::distance_2(p, q, traits);
}

template<typename Kernel>
typename Kernel::FT distance(const CGAL::Point_3<Kernel>& p,
                             const CGAL::Point_3<Kernel>& q)
{
  const Kernel traits;
  return CGAL::Weights::internal::distance_3(p, q, traits);
}

template<typename Kernel>
typename Kernel::FT area(const CGAL::Point_2<Kernel>& p,
                         const CGAL::Point_2<Kernel>& q,
                         const CGAL::Point_2<Kernel>& r)
{
  const Kernel traits;
  return CGAL::Weights::internal::positive_area_2(p, q, r, traits);
}

template<typename Kernel>
typename Kernel::FT area(const CGAL::Point_3<Kernel>& p,
                         const CGAL::Point_3<Kernel>& q,
                         const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  return CGAL::Weights::internal::positive_area_3(p, q, r, traits);
}

template<typename Kernel>
typename Kernel::FT scalar_product(const CGAL::Point_2<Kernel>& p,
                                   const CGAL::Point_2<Kernel>& q,
                                   const CGAL::Point_2<Kernel>& r)
{
  const Kernel traits;
  auto scalar_product_2 = traits.compute_scalar_product_2_object();
  auto vector_2 = traits.construct_vector_2_object();

  const auto v1 = vector_2(q, r);
  const auto v2 = vector_2(q, p);
  return scalar_product_2(v1, v2);
}

template<typename Kernel>
typename Kernel::FT scalar_product(const CGAL::Point_3<Kernel>& p,
                                   const CGAL::Point_3<Kernel>& q,
                                   const CGAL::Point_3<Kernel>& r)
{
  const Kernel traits;
  auto scalar_product_3 = traits.compute_scalar_product_3_object();
  auto vector_3 = traits.construct_vector_3_object();

  const auto v1 = vector_3(q, r);
  const auto v2 = vector_3(q, p);

  return scalar_product_3(v1, v2);
}

} // namespace internal
} // namespace Weights
} // namespace CGAL

namespace tests {

template<typename FT>
FT get_tolerance() { return FT(1e-10); }

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_all_triangles()
{
  using Point_2 = typename Kernel::Point_2;

  const std::array<Point_2, 3> triangle0 = { Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0) };
  const std::array<Point_2, 3> triangle1 = { Point_2(-2, 0), Point_2(0, -1), Point_2(2, 0) };
  const std::array<Point_2, 3> triangle2 = { Point_2(-2, 0), Point_2(-2, -2), Point_2(2, 0) };
  const std::array<Point_2, 3> triangle3 = { Point_2(-2, 0), Point_2(2, -2), Point_2(2, 0) };

  return { triangle0, triangle1, triangle2, triangle3 };
}

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_symmetric_triangles()
{
  using Point_2 = typename Kernel::Point_2;

  const std::array<Point_2, 3> triangle0 = { Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0) };
  const std::array<Point_2, 3> triangle1 = { Point_2(-2, 0), Point_2(0, -1), Point_2(2, 0) };
  const std::array<Point_2, 3> triangle2 = { Point_2(-3, 0), Point_2(0, -1), Point_2(3, 0) };

  return { triangle0, triangle1, triangle2 };
}

template<typename Kernel>
std::vector< std::array<typename Kernel::Point_2, 3> >
get_uniform_triangles()
{
  using Point_2 = typename Kernel::Point_2;

  const std::array<Point_2, 3> triangle0 = { Point_2(-1, 0), Point_2(0, -1), Point_2(1, 0) };
  const std::array<Point_2, 3> triangle1 = { Point_2(-2, 0), Point_2(0, -2), Point_2(2, 0) };
  const std::array<Point_2, 3> triangle2 = { Point_2(1, 0), Point_2(-1, 0), Point_2(-1, -2) };
  const std::array<Point_2, 3> triangle3 = { Point_2(1, -2), Point_2(1, 0), Point_2(-1, 0) };

  return { triangle0, triangle1, triangle2, triangle3 };
}

template<typename Kernel>
std::vector< std::vector<typename Kernel::Point_2> >
get_all_polygons()
{
  using Point_2 = typename Kernel::Point_2;

  const std::vector<Point_2> polygon0 = { Point_2(-2, -2), Point_2(2, -2), Point_2(0, 2) };
  const std::vector<Point_2> polygon1 = { Point_2(-1, -1), Point_2(1, -1), Point_2(1, 1), Point_2(-1, 1) };
  const std::vector<Point_2> polygon2 = { Point_2(-2, 0), Point_2(0, -2), Point_2(2, 0), Point_2(0, 2) };
  const std::vector<Point_2> polygon3 = { Point_2(-2, -2), Point_2(2, -2), Point_2(2, 0),
                                          Point_2(0, 2), Point_2(-2, 0) };

  return { polygon0, polygon1, polygon2, polygon3 };
}

template<typename Kernel,
         typename Weight_wrapper>
void test_query(const Weight_wrapper& wrapper,
                const typename Kernel::Point_2& query,
                const std::array<typename Kernel::Point_2, 3>& neighbors)
{
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

  const FT a2 = wrapper.weight_a(t2, r2, p2, q2);
  const FT b2 = wrapper.weight_b(t2, r2, p2, q2);
  assert(a2 >= FT(0) && b2 >= FT(0));
  assert(CGAL::abs(a2 - b2) < tol);

  if (wrapper.supports_3d())
  {
    const FT a3 = wrapper.weight_a(t3, r3, p3, q3);
    const FT b3 = wrapper.weight_b(t3, r3, p3, q3);
    assert(a3 >= FT(0) && b3 >= FT(0));
    assert(CGAL::abs(a3 - b3) < tol);
    assert(CGAL::abs(a2 - a3) < tol);
    assert(CGAL::abs(b2 - b3) < tol);
  }
}

template<typename Kernel,
         typename Weight_wrapper>
void test_symmetry_x(const Weight_wrapper& wrapper,
                     const std::array<typename Kernel::Point_2, 3>& neighbors,
                     const typename Kernel::FT& x)
{
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

  const FT a2 = wrapper.weight_a(t2, r2, p2, Point_2(-x, 0));
  const FT b2 = wrapper.weight_a(t2, r2, p2, Point_2(+x, 0));
  assert(a2 >= FT(0) && b2 >= FT(0));
  assert(CGAL::abs(a2 - b2) < tol);

  if (wrapper.supports_3d())
  {
    const FT a3 = wrapper.weight_a(t3, r3, p3, Point_3(-x, 0, 1));
    const FT b3 = wrapper.weight_a(t3, r3, p3, Point_3(+x, 0, 1));
    assert(a3 >= FT(0) && b3 >= FT(0));
    assert(CGAL::abs(a3 - b3) < tol);
    assert(CGAL::abs(a2 - a3) < tol);
    assert(CGAL::abs(b2 - b3) < tol);
  }
}

template<typename Kernel,
         typename Weight_wrapper_1,
         typename Weight_wrapper_2>
void test_compare(const Weight_wrapper_1& wrapper1,
                  const Weight_wrapper_2& wrapper2,
                  const typename Kernel::Point_2& query,
                  const std::array<typename Kernel::Point_2, 3>& neighbors)
{
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

  const FT a2 = wrapper1.weight_a(t2, r2, p2, q2);
  const FT b2 = wrapper2.weight_a(t2, r2, p2, q2);
  assert(a2 >= FT(0) && b2 >= FT(0));
  assert(CGAL::abs(a2 - b2) < tol);

  if (wrapper1.supports_3d() && wrapper2.supports_3d())
  {
    const FT a3 = wrapper1.weight_a(t3, r3, p3, q3);
    const FT b3 = wrapper2.weight_a(t3, r3, p3, q3);
    assert(a3 >= FT(0) && b3 >= FT(0));
    assert(CGAL::abs(a3 - b3) < tol);
  }
}

template<typename Kernel,
         typename Weight_wrapper>
void test_neighbors(const Weight_wrapper& wrapper,
                    const std::array<typename Kernel::Point_2, 3>& neighbors)
{
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

  const FT a2 = wrapper.weight(p2, q2, r2);
  const FT a3 = wrapper.weight(p3, q3, r3);
  assert(a2 >= FT(0) && a3 >= FT(0));
  assert(CGAL::abs(a2 - a3) < tol);
}

template<typename Kernel,
         typename Weight_wrapper>
void test_area(const Weight_wrapper& wrapper,
               const std::array<typename Kernel::Point_2, 3>& neighbors)
{
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

  const FT a2 = wrapper.weight(p2, q2, r2);
  const FT a3 = wrapper.weight(p3, q3, r3);
  assert(a2 <= CGAL::Weights::internal::area(p2, q2, r2));
  assert(a3 <= CGAL::Weights::internal::area(p3, q3, r3));

  assert(a2 >= FT(0));
  assert(a3 >= FT(0));
}

template<typename FT, typename Point>
void test_coordinates(const Point& query,
                      const std::vector<Point>& polygon,
                      const std::vector<FT>& weights)
{
  assert(weights.size() > 0);

  // Compute the sum of weights.
  const FT tol = get_tolerance<FT>();
  FT sum = FT(0);
  for (const FT& weight : weights)
    sum += weight;
  assert(sum >= tol);

  // Compute coordinates.
  std::vector<FT> coordinates;
  coordinates.reserve(weights.size());
  for (const FT& weight : weights)
    coordinates.push_back(weight / sum);

  assert(coordinates.size() == weights.size());

  // Test partition of unity.
  sum = FT(0);
  for (const FT& coordinate : coordinates)
    sum += coordinate;
  assert(CGAL::abs(FT(1) - sum) < tol);

  // Test linear precision.
  FT x = FT(0), y = FT(0);
  for (std::size_t i = 0; i < polygon.size(); ++i)
  {
    x += coordinates[i] * polygon[i].x();
    y += coordinates[i] * polygon[i].y();
  }
  assert(CGAL::abs(query.x() - x) < tol);
  assert(CGAL::abs(query.y() - y) < tol);
}

template<typename Kernel,
         typename Weight_wrapper>
void test_on_polygon(const Weight_wrapper& wrapper,
                     const typename Kernel::Point_2& query_2,
                     const std::vector<typename Kernel::Point_2>& polygon_2)
{
  // Get weights.
  using FT = typename Kernel::FT;
  assert(polygon_2.size() >= 3);

  // 2D version.
  std::vector<FT> weights_2;
  weights_2.reserve(polygon_2.size());
  wrapper.compute_on_polygon(polygon_2, query_2, Kernel(), std::back_inserter(weights_2));
  assert(weights_2.size() == polygon_2.size());
  test_coordinates(query_2, polygon_2, weights_2);

  // 3D version.
  using Point_3 = typename Kernel::Point_3;
  const Point_3 query_3(query_2.x(), query_2.y(), 1);
  std::vector<Point_3> polygon_3;
  polygon_3.reserve(polygon_2.size());
  for (const auto& vertex_2 : polygon_2)
    polygon_3.emplace_back(vertex_2.x(), vertex_2.y(), 1);
  assert(polygon_3.size() == polygon_2.size());

  const CGAL::Projection_traits_xy_3<Kernel> ptraits;
  std::vector<FT> weights_3;
  weights_3.reserve(polygon_3.size());
  wrapper.compute_on_polygon(polygon_3, query_3, ptraits, std::back_inserter(weights_3));
  assert(weights_3.size() == polygon_3.size());

  test_coordinates(query_3, polygon_3, weights_3);
}

template<typename Kernel,
         typename Weight_wrapper>
void test_barycentric_properties(const Weight_wrapper& wrapper,
                                 const typename Kernel::Point_2& query,
                                 const std::vector<typename Kernel::Point_2>& polygon)
{
  // Get weights.
  using FT = typename Kernel::FT;
  const std::size_t n = polygon.size();
  assert(n >= 3);

  // Check properties.
  std::vector<FT> weights;
  weights.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
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

  test_coordinates(query, polygon, weights);
}

template<typename Kernel,
         typename Weight_wrapper_1,
         typename Weight_wrapper_2>
void test_analytic_weight(const Weight_wrapper_1& weight,
                          const Weight_wrapper_2& alternative)
{
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

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
  for (const auto& config : configs)
  {
    test_query<Kernel>(weight, zero, config);
    for (const auto& query : queries)
      test_query<Kernel>(weight, query, config);
  }

  // Test alternative formulations.
  for (const auto& config : configs)
  {
    test_compare<Kernel>(weight, alternative, zero, config);
    for (const auto& query : queries)
      test_compare<Kernel>(weight, alternative, query, config);
  }

  // Test symmetry along x axis.
  configs = get_symmetric_triangles<Kernel>();
  for (const auto& config : configs)
  {
    test_symmetry_x<Kernel>(weight, config, q);
    test_symmetry_x<Kernel>(weight, config, h);
    test_symmetry_x<Kernel>(weight, config, t);
  }

  // Test barycentric properties.
  if (weight.is_barycentric())
  {
    const auto polygons = get_all_polygons<Kernel>();
    for (const auto& polygon : polygons)
    {
      test_barycentric_properties<Kernel>(weight, zero, polygon);
      for (const auto& query : queries)
        test_barycentric_properties<Kernel>(weight, query, polygon);
    }
  }

  return;
}

template<typename Kernel,
         typename Weight_wrapper_1,
         typename Weight_wrapper_2>
void test_barycentric_weight(const Weight_wrapper_1& weight,
                             const Weight_wrapper_2& alternative)
{
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  const FT q = FT(1) / FT(4);
  const FT h = FT(1) / FT(2);
  const Point_2 zero(0, 0);
  const std::vector<Point_2> queries =
  {
    Point_2(-h,  0), Point_2(+h,  0), Point_2(-q,  0), Point_2(+q,  0),
    Point_2( 0, -h), Point_2( 0, +h), Point_2( 0, -q), Point_2( 0, +q),
    Point_2(-h, -h), Point_2(+h, +h), Point_2(-q, -q), Point_2(+q, +q),
    Point_2(-h, +q), Point_2(+h, -q), Point_2(-q, +h), Point_2(+q, -h)
  };

  // Test analytic formulations.
  test_analytic_weight<Kernel>(weight, alternative);

  // Test on polygons.
  const auto polygons = get_all_polygons<Kernel>();
  for (const auto& polygon : polygons)
  {
    test_on_polygon<Kernel>(weight, zero, polygon);
    for (const auto& query : queries)
      test_on_polygon<Kernel>(weight, query, polygon);
  }
}

template<typename Kernel,
         typename Weight_wrapper>
void test_region_weight(const Weight_wrapper& weight)
{
  auto configs = get_all_triangles<Kernel>();
  for (const auto& config : configs)
    test_neighbors<Kernel>(weight, config);

  configs = get_uniform_triangles<Kernel>();
  for (const auto& config : configs)
    test_area<Kernel>(weight, config);
}

} // namespace tests

#endif // CGAL_WEIGHTS_TESTS_UTILS_H
