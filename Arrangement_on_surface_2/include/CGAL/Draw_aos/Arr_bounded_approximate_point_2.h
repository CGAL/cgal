
#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_POINT_2_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_POINT_2_H

#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_render_context.h>

namespace CGAL {
namespace draw_aos {

template <typename Arrangement>
class Arr_bounded_approximate_point_2
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Point_2 = typename Geom_traits::Point_2;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Point_geom = typename Arr_approximate_traits<Geom_traits>::Point_geom;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

public:
  Arr_bounded_approximate_point_2(const Bounded_render_context& ctx)
      : m_ctx(ctx)
      , m_approx(ctx.m_traits.approximate_2_object()) {}

  /**
   * @brief Approximate a vertex within the x-bounded range.
   *
   * The function uses cached values if available.
   * @precondition: The vertex must have an associated point.
   *
   * @param vh the vertex handle
   * @return const Point_geom&
   */
  const Point_geom& operator()(const Vertex_const_handle& vh) const {
    auto [point, inserted] = m_ctx.m_cache.try_emplace(vh);
    if(!inserted) return point;
    return point = m_approx(vh->point());
  }

private:
  const Bounded_render_context& m_ctx;
  const typename Geom_traits::Approximate_2 m_approx;
};

} // namespace draw_aos
} // namespace CGAL

#endif
