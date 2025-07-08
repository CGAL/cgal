#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_AT_X_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATE_POINT_2_AT_X_H
#include <boost/iterator/function_output_iterator.hpp>

#include "CGAL/number_utils.h"
#include "CGAL/Arr_rational_function_traits_2.h"
#include "CGAL/Draw_aos/Arr_construct_curve_end.h"
#include "CGAL/Draw_aos/Arr_construct_segments.h"
#include "CGAL/Draw_aos/Arr_approximate_point_2.h"
#include "CGAL/Draw_aos/type_utils.h"

namespace CGAL {
namespace draw_aos {

/**
 * @brief Functor to compute the point at a given x-coordinate on an x-monotone curve within a bounding box.
 */
template <typename GeomTraits>
class Arr_approximate_point_2_at_x
{
  using Point_2 = typename Traits_adaptor<GeomTraits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<GeomTraits>::X_monotone_curve_2;
  using Intersect_2 = typename Traits_adaptor<GeomTraits>::Intersect_2;
  using FT = typename Traits_adaptor<GeomTraits>::FT;
  using Approx_traits = Arr_approximation_geometry_traits<GeomTraits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Is_vertical_2 = typename Traits_adaptor<GeomTraits>::Is_vertical_2;

public:
  Arr_approximate_point_2_at_x(const GeomTraits& traits)
      : m_approx_pt(traits)
      , m_cst_vertical_segment(traits)
      , m_is_vertical_2(traits.is_vertical_2_object())
      , m_intersect_2(traits.intersect_2_object()) {}

  /**
   * @brief Computes the point at a given x-coordinate on an x-monotone curve
   *
   * @precondition: The curve is not verical
   * @param curve
   * @param x
   * @return true if there is an intersection at given x,
   * @return false otherwise.
   */
  std::optional<Approx_point> operator()(const X_monotone_curve_2& curve, FT x) const {
    CGAL_assertion(!m_is_vertical_2(curve));

    using Multiplicity = typename GeomTraits::Multiplicity;
    using Intersect_point = std::pair<Point_2, Multiplicity>;
    using Intersect_curve = X_monotone_curve_2;
    using Intersect_type = std::variant<Intersect_point, Intersect_curve>;

    auto vertical_line = m_cst_vertical_segment(x, m_ymin, m_ymax);
    std::optional<Approx_point> pt;
    auto func_out_iter = boost::make_function_output_iterator([&pt, this](const Intersect_type& res) {
      CGAL_assertion_msg(std::holds_alternative<Intersect_point>(res),
                         "Unexpected intersection type, expected Intersect_point");
      pt = this->m_approx_pt(std::get<Intersect_point>(res).first);
    });
    m_intersect_2(vertical_line, curve, func_out_iter);

    if(pt.has_value()) {
      pt = Approx_point(CGAL::to_double(x), pt->y());
    }
    return pt;
  }

private:
  const Arr_approximate_point_2<GeomTraits> m_approx_pt;
  const Arr_construct_vertical_segment<GeomTraits> m_cst_vertical_segment;
  const Is_vertical_2 m_is_vertical_2;
  const Intersect_2 m_intersect_2;

  // Should be enough for visualization purposes.
  // constexpr static double m_ymin = std::numeric_limits<double>::lowest();
  // constexpr static double m_ymax = std::numeric_limits<double>::max();
  // maximum of double is too large for CORE number types, no idea why
  constexpr static double m_ymin = -1e8;
  constexpr static double m_ymax = 1e8;
};

// Specialization for rational function traits, which has no vertical curves but provides
// evaluation at x directly.
template <typename Kernel>
class Arr_approximate_point_2_at_x<Arr_rational_function_traits_2<Kernel>>
{
  using Geom_traits = CGAL::Arr_rational_function_traits_2<Kernel>;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<Geom_traits>::X_monotone_curve_2;
  using FT = typename Traits_adaptor<Geom_traits>::FT;
  using Approx_point = typename Arr_approximation_geometry_traits<Geom_traits>::Approx_point;
  using Construct_curve_end = Arr_construct_curve_end<Geom_traits>;

public:
  Arr_approximate_point_2_at_x(const Geom_traits& traits)
      : m_cst_curve_end(traits) {}

  std::optional<Approx_point> operator()(const X_monotone_curve_2& curve, FT x) const {
    std::optional<Point_2> min_end = m_cst_curve_end(curve, ARR_MIN_END);
    if(min_end.has_value() && x < min_end->x()) {
      return std::nullopt;
    }
    std::optional<Point_2> max_end = m_cst_curve_end(curve, ARR_MAX_END);
    if(max_end.has_value() && x > max_end->x()) {
      return std::nullopt;
    }

    using Bound = typename Geom_traits::Bound;
    const auto& numerator = curve.numerator();
    const auto& denominator = curve.denominator();
    Bound approx_x = (x.lower() + x.upper()) / 2.0;
    double enum_at_x = CGAL::to_double(numerator.evaluate(approx_x));
    double denom_at_x = CGAL::to_double(denominator.evaluate(approx_x));
    return std::make_optional(Approx_point(CGAL::to_double(x), enum_at_x / denom_at_x));
  }

private:
  const Construct_curve_end m_cst_curve_end;
};

} // namespace draw_aos
} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_COMPUTE_Y_AT_X_H