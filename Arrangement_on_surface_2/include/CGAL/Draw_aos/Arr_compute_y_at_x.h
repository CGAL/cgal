#ifndef CGAL_DRAW_AOS_ARR_COMPUTE_Y_AT_X_H
#define CGAL_DRAW_AOS_ARR_COMPUTE_Y_AT_X_H
#include "CGAL/Bbox_2.h"
#include "CGAL/Draw_aos/Arr_render_context.h"
#include "CGAL/Draw_aos/helpers.h"
#include <boost/iterator/function_output_iterator.hpp>

namespace CGAL {

/**
 * @brief Functor to compute the y-coordinate at a given x-coordinate for an x-monotone curve within a bounding box.
 */
class Arr_compute_y_at_x
{
public:
  using Point_2 = Geom_traits::Point_2;
  using X_monotone_curve_2 = Geom_traits::X_monotone_curve_2;
  using Intersect_2 = Geom_traits::Intersect_2;
  using Construct_min_vertex_2 = Geom_traits::Construct_min_vertex_2;
  using FT = Geom_traits::FT;
  using Approximate_2 = Geom_traits::Approximate_2;
  using Is_vertical_2 = Geom_traits::Is_vertical_2;

  Arr_compute_y_at_x(const Arr_render_context& ctx)
      : m_ctx(ctx)
      // TODO: some traits does not have approximate_2_object. we'll need a specialization for them.
      , m_approx(ctx.traits.approximate_2_object()) {}

  /**
   * @brief Computes the y-coordinate at a given x-coordinate for an x-monotone curve trimmed
   * to the bounding box.
   *
   * The bounding box here is considered as closed.
   *
   * @precondition The curve is not verical
   * @param curve
   * @param x
   * @return true if there is an intersection at given x within the bounding box,
   * @return false otherwise.
   */
  std::optional<FT> operator()(const X_monotone_curve_2& curve, const FT& x) const {
    CGAL_assertion(!m_ctx.is_vertical_2(curve));

    auto min_pt = m_ctx.cst_curve_end(curve, ARR_MIN_END);
    auto max_pt = m_ctx.cst_curve_end(curve, ARR_MAX_END);
    if(min_pt.has_value() && min_pt->x() == x) {
      return min_pt->y();
    }
    if(max_pt.has_value() && max_pt->x() == x) {
      return max_pt->y();
    }

    using Multiplicity = Geom_traits::Multiplicity;
    using Intersect_point = std::pair<Point_2, Multiplicity>;
    using Intersect_curve = X_monotone_curve_2;
    using Intersect_type = std::variant<Intersect_point, Intersect_curve>;

    auto vertical_line = m_ctx.cst_vertical_segment(x, m_ymin, m_ymax);
    std::optional<FT> y;
    auto func_out_iter = boost::make_function_output_iterator(
        [&y, this](const Intersect_type& res) { y = std::get<Intersect_point>(res).first.y(); });
    m_ctx.intersect_2(curve, vertical_line, func_out_iter);
    return y;
  }

private:
  Approximate_2 m_approx;
  const Arr_render_context& m_ctx;

  // Should be enough for visualization purposes.
  constexpr static double m_ymin = std::numeric_limits<double>::lowest();
  constexpr static double m_ymax = std::numeric_limits<double>::max();
};

} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_COMPUTE_Y_AT_X_H