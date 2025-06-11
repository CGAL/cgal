#ifndef CGAL_ARR_CONSTRUCT_CURVE_END_H
#define CGAL_ARR_CONSTRUCT_CURVE_END_H

#include "CGAL/Arr_enums.h"
#include "CGAL/Arr_has.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include <optional>

namespace CGAL {

namespace internal {

template <typename Geom_traits>
class Arr_construct_curve_end_base
{
  using Construct_min_vertex_2 = typename Geom_traits::Construct_min_vertex_2;
  using Construct_max_vertex_2 = typename Geom_traits::Construct_max_vertex_2;

protected:
  Arr_construct_curve_end_base(const Geom_traits& traits)
      : m_cst_min_vertex(traits.construct_min_vertex_2_object())
      , m_cst_max_vertex(traits.construct_max_vertex_2_object()) {}

  Construct_max_vertex_2 m_cst_max_vertex;
  Construct_min_vertex_2 m_cst_min_vertex;
};

template <typename Geom_traits, bool Has_parameter_space_in_x>
class Arr_construct_curve_end_impl;

template <typename Geom_traits>
class Arr_construct_curve_end_impl<Geom_traits, true> : public Arr_construct_curve_end_base<Geom_traits>
{
  using Approx_geom_traits = Arr_approximation_geometry_traits;
  using Point_2 = typename Geom_traits::Point_2;
  using X_monotone_curve = typename Geom_traits::X_monotone_curve_2;
  using Parameter_space_in_x_2 = typename Geom_traits::Parameter_space_in_x_2;
  using Parameter_space_in_y_2 = typename Geom_traits::Parameter_space_in_y_2;

public:
  Arr_construct_curve_end_impl(const Geom_traits& traits)
      : Arr_construct_curve_end_base<Geom_traits>(traits)
      , m_param_space_in_x(traits.parameter_space_in_x_2_object())
      , m_param_space_in_y(traits.parameter_space_in_y_2_object()) {}

  /**
   * @brief Returns an optional type contains the approximated point the given end of the curve is bounded.
   * Otherwise it's std::nullopt.
   *
   * @param curve
   * @param ce
   * @return std::optional<Point_2>
   */
  std::optional<Point_2> operator()(const X_monotone_curve& curve, Arr_curve_end ce) const {
    auto ps_x = m_param_space_in_x(curve, ce);
    auto ps_y = m_param_space_in_y(curve, ce);

    if(ps_x == ARR_INTERIOR && ps_y == ARR_INTERIOR) {
      return std::make_optional(ce == ARR_MIN_END ? this->m_cst_min_vertex(curve) : this->m_cst_max_vertex(curve));
    }

    return std::nullopt;
  }

private:
  Parameter_space_in_x_2 m_param_space_in_x;
  Parameter_space_in_y_2 m_param_space_in_y;
};

template <typename Geom_traits>
class Arr_construct_curve_end_impl<Geom_traits, false> : public Arr_construct_curve_end_base<Geom_traits>
{
  using Approx_geom_traits = Arr_approximation_geometry_traits;
  using Point_2 = typename Geom_traits::Point_2;
  using X_monotone_curve = typename Geom_traits::X_monotone_curve_2;

public:
  Arr_construct_curve_end_impl(const Geom_traits& traits)
      : Arr_construct_curve_end_base<Geom_traits>(traits) {}

  /**
   * @brief Construct the min or max vertex of the x-monotone curve based on the curve end enum.
   *
   * @param curve
   * @param ce
   * @return std::optional<Point_2>
   */
  std::optional<Point_2> operator()(const X_monotone_curve& curve, Arr_curve_end ce) const {
    return ce == ARR_MIN_END ? this->m_cst_min_vertex(curve) : this->m_cst_max_vertex(curve);
  }
};

} // namespace internal

template <typename Geom_traits>
using Arr_construct_curve_end =
    internal::Arr_construct_curve_end_impl<Geom_traits, has_parameter_space_in_x_2<Geom_traits>::value>;

} // namespace CGAL

#endif // CGAL_ARR_CONSTRUCT_CURVE_END_H