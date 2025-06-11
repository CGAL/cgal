#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_H

#include "CGAL/Arr_has.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include "CGAL/number_utils.h"

namespace CGAL {

namespace internal {

template <typename Geom_traits, bool Has_approximate_2>
class Arr_approximate_point_2_impl;

template <typename Geom_traits>
class Arr_approximate_point_2_impl<Geom_traits, true>
{
  using Approx_kernel = Arr_approximation_geometry_traits::Approximation_kernel;
  using Point_2 = typename Geom_traits::Point_2;
  using Approx_point = Arr_approximation_geometry_traits::Approx_point;

public:
  Arr_approximate_point_2_impl(const Geom_traits& traits)
      : m_approx(traits.approximate_2_object()) {}

  /**
   * @brief Approximate a point.
   * TODO: make it work for spherical traits.
   *
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
  typename Geom_traits::Approximate_2 m_approx;
};

// Fallback to use CGAL::to_double for traits that do not have an approximate_2_object.
template <typename Geom_traits>
class Arr_approximate_point_2_impl<Geom_traits, false>
{
  using Approx_kernel = Arr_approximation_geometry_traits::Approximation_kernel;
  using Point_2 = typename Geom_traits::Point_2;
  using Approx_point = Arr_approximation_geometry_traits::Approx_point;

public:
  // traits object  is not used in the fallback implementation, but we keep it for consistency.
  Arr_approximate_point_2_impl(const Geom_traits& traits) {}

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

template <typename Geom_traits>
using Arr_approximate_point_2 =
    internal::Arr_approximate_point_2_impl<Geom_traits, has_approximate_2<Geom_traits>::value>;

} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_H