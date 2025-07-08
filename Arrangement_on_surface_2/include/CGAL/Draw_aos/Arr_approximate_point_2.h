#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_H

#include "CGAL/Arr_rational_function_traits_2.h"
#include "CGAL/number_utils.h"
#include "CGAL/Draw_aos/type_utils.h"

namespace CGAL {
namespace draw_aos {
namespace internal {

template <typename GeomTraits, bool Has_approximate_2>
class Arr_approximate_point_2_impl;

template <typename GeomTraits>
class Arr_approximate_point_2_impl<GeomTraits, true>
{
  using Approx_traits = Arr_approximation_geometry_traits<GeomTraits>;
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using Approx_point = typename Approx_traits::Approx_point;

public:
  Arr_approximate_point_2_impl(const GeomTraits& traits)
      : m_approx(traits.approximate_2_object()) {}

  /**
   * @brief Approximate a point.
   * TODO: make it work for spherical traits.
   * @param pt
   * @return Point_geom
   */
  Approx_point operator()(const Point_2& pt) const { return {m_approx(pt, 0), m_approx(pt, 1)}; }

  /**
   * @brief Approximate a specific dimension of a point.
   *
   * @param pt
   * @param dim 0 for x, 1 for y, 2 for z if applicable. An exception is thrown if the dimension is invalid.
   * @return double
   */
  double operator()(const Point_2& pt, int dim) const { return m_approx(pt, dim); }

private:
  const typename GeomTraits::Approximate_2 m_approx;
};

// Specialized for Arr_rational_function_traits_2.
template <typename Kernel>
class Arr_approximate_point_2_impl<Arr_rational_function_traits_2<Kernel>, true>
{
  using Geom_traits = Arr_rational_function_traits_2<Kernel>;
  using Approx_traits = Arr_approximation_geometry_traits<Geom_traits>;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approximate_2 = typename Geom_traits::Approximate_2;

public:
  Arr_approximate_point_2_impl(const Geom_traits& traits) {}

  Approx_point operator()(const Point_2& pt) const { return {pt.x().to_double(), pt.y().to_double()}; }

  double operator()(const Point_2& pt, int dim) const { return (dim == 0) ? pt.x().to_double() : pt.y().to_double(); }
};

// Fallback to use CGAL::to_double for traits that do not have an approximate_2_object.
template <typename GeomTraits>
class Arr_approximate_point_2_impl<GeomTraits, false>
{
  using Approx_traits = Arr_approximation_geometry_traits<GeomTraits>;
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using Approx_point = typename Approx_traits::Approx_point;

public:
  // traits object is not used in the fallback implementation, but we keep it for consistency.
  Arr_approximate_point_2_impl(const GeomTraits& traits) {}

  /**
   * @brief Approximate a point.
   *
   * @note this functions does not check if the point is within the bounding box.
   * @param pt
   * @return Point_geom
   */
  Approx_point operator()(const Point_2& pt) const { return {CGAL::to_double(pt.x()), CGAL::to_double(pt.y())}; }

  /**
   * @brief Approximate a specific dimension of a point.
   *
   * @param pt
   * @param dim 0 for x, 1 for y, 2 for z if applicable. An exception is thrown if the dimension is invalid.
   * @return double
   */
  double operator()(const Point_2& pt, int dim) const {
    if(dim == 0) {
      return CGAL::to_double(pt.x());
    } else if(dim == 1) {
      return CGAL::to_double(pt.y());
    } else if(dim == 2 && pt.dimension() == 3) {
      return CGAL::to_double(pt.z());
    } else {
      throw std::invalid_argument("Invalid dimension for approximate_point");
    }
  }
};

} // namespace internal

template <typename GeomTraits>
using Arr_approximate_point_2 =
    internal::Arr_approximate_point_2_impl<GeomTraits,
                                           has_approximate_2_object_v<GeomTraits> &&
                                               has_operator_point_v<GeomTraits, typename GeomTraits::Approximate_2>>;

} // namespace draw_aos
} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_H