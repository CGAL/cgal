#ifndef CGAL_WEIGHTS_TESTS_WRAPPERS_H
#define CGAL_WEIGHTS_TESTS_WRAPPERS_H

// STL includes.
#include <vector>
#include <cassert>
#include <string>

// CGAL includes.
#include <CGAL/Weights.h>

namespace wrappers {

template<typename Kernel>
struct Authalic_wrapper {
  using FT = typename Kernel::FT;
  template<typename Point>
  FT weight_a(
    const Point& t, const Point& r, const Point& p, const Point& q) const {
    return CGAL::Weights::authalic_weight(t, r, p, q);
  }
  template<typename Point>
  FT weight_b(
    const Point& t, const Point& r, const Point& p, const Point& q) const {
    return
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(t, r, q),
      CGAL::Weights::squared_distance(q, r)) +
    CGAL::Weights::half_authalic_weight(
      CGAL::Weights::cotangent(q, r, p),
      CGAL::Weights::squared_distance(q, r));
  }
  bool supports_3d() const { return true; }
};

template<typename Kernel>
struct Cotangent_wrapper {
  using FT = typename Kernel::FT;
  template<typename Point>
  FT weight_a(
    const Point& t, const Point& r, const Point& p, const Point& q) const {
    return CGAL::Weights::cotangent_weight(t, r, p, q);
  }
  template<typename Point>
  FT weight_b(
    const Point& t, const Point& r, const Point& p, const Point& q) const {
    return
    CGAL::Weights::half_cotangent_weight(
      CGAL::Weights::cotangent(q, t, r)) +
    CGAL::Weights::half_cotangent_weight(
      CGAL::Weights::cotangent(r, p, q));
  }
  bool supports_3d() const { return true; }
};

template<typename Kernel>
struct Tangent_wrapper {
  using FT = typename Kernel::FT;
  template<typename Point>
  FT weight_a(
    const Point& t, const Point& r, const Point& p, const Point& q) const {
    return CGAL::Weights::tangent_weight(t, r, p, q);
  }
  template<typename Point>
  FT weight_b(
    const Point& t, const Point& r, const Point& p, const Point& q) const {
    return
    CGAL::Weights::half_tangent_weight(
      CGAL::Weights::distance(r, q),
      CGAL::Weights::distance(t, q),
      CGAL::Weights::area(r, q, t),
      CGAL::Weights::scalar_product(r, q, t)) +
    CGAL::Weights::half_tangent_weight(
      CGAL::Weights::distance(r, q),
      CGAL::Weights::distance(p, q),
      CGAL::Weights::area(p, q, r),
      CGAL::Weights::scalar_product(p, q, r));
  }
  bool supports_3d() const { return true; }
};

template<typename Kernel>
struct Wachspress_wrapper {
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  FT weight_a(
    const Point_2& t, const Point_2& r, const Point_2& p, const Point_2& q) const {
    return CGAL::Weights::wachspress_weight(t, r, p, q);
  }
  FT weight_a(
    const Point_3&, const Point_3&, const Point_3&, const Point_3&) const {
    return FT(-1);
  }
  template<typename Point>
  FT weight_b(
    const Point&, const Point&, const Point&, const Point&) const {
    return FT(-1);
  }
  bool supports_3d() const { return false; }
};

template<typename Kernel>
struct Discrete_harmonic_wrapper {
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  FT weight_a(
    const Point_2& t, const Point_2& r, const Point_2& p, const Point_2& q) const {
    return CGAL::Weights::discrete_harmonic_weight(t, r, p, q);
  }
  FT weight_a(
    const Point_3&, const Point_3&, const Point_3&, const Point_3&) const {
    return FT(-1);
  }
  template<typename Point>
  FT weight_b(
    const Point&, const Point&, const Point&, const Point&) const {
    return FT(-1);
  }
  bool supports_3d() const { return false; }
};

template<typename Kernel>
struct Mean_value_wrapper {
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  FT weight_a(
    const Point_2& t, const Point_2& r, const Point_2& p, const Point_2& q) const {
    return CGAL::Weights::mean_value_weight(t, r, p, q);
  }
  FT weight_a(
    const Point_3&, const Point_3&, const Point_3&, const Point_3&) const {
    return FT(-1);
  }
  template<typename Point>
  FT weight_b(
    const Point&, const Point&, const Point&, const Point&) const {
    return FT(-1);
  }
  bool supports_3d() const { return false; }
};

template<typename Kernel>
struct Three_point_family_wrapper {
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  const FT a;
  Three_point_family_wrapper(const FT a) : a(a) { }
  FT weight_a(
    const Point_2& t, const Point_2& r, const Point_2& p, const Point_2& q) const {
    return CGAL::Weights::three_point_family_weight(t, r, p, q, a);
  }
  FT weight_a(
    const Point_3&, const Point_3&, const Point_3&, const Point_3&) const {
    return FT(-1);
  }
  template<typename Point>
  FT weight_b(
    const Point&, const Point&, const Point&, const Point&) const {
    return FT(-1);
  }
  bool supports_3d() const { return false; }
};

template<typename Kernel>
struct Barycentric_region_wrapper {
  using FT = typename Kernel::FT;
  template<typename Point>
  FT weight_a(
    const Point& p, const Point& q, const Point& r) const {
    return CGAL::Weights::barycentric_area(p, q, r);
  }
  template<typename Point>
  FT weight_b(
    const Point&, const Point&, const Point&) const {
    return FT(-1);
  }
  bool supports_3d() const { return true; }
};

} // namespace wrappers

#endif // CGAL_WEIGHTS_TESTS_WRAPPERS_H
