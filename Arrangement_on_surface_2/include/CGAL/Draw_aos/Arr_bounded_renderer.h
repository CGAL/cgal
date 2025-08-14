#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Draw_aos/Arr_approximation_cache.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_face.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/Arr_portals.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/**
 * @brief Render arrangement on surface within a bounding box.
 */
template <typename Arrangement>
class Arr_bounded_renderer
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Render_context = Arr_render_context<Arrangement>;
  using Approx_cache = Arr_approximation_cache<Arrangement>;

public:
  Arr_bounded_renderer(const Render_context& ctx, Bbox_2 bbox)
      : m_ctx(ctx)
      , m_bbox(bbox) {}

  Approx_cache render() const {
    Approx_cache cache;
    if(m_ctx.is_cancelled()) return cache;
    cache.reserve_faces(m_ctx.m_arr.number_of_faces());
    cache.reserve_halfedges(m_ctx.m_arr.number_of_halfedges());
    cache.reserve_vertices(m_ctx.m_arr.number_of_vertices());

    Arr_bounded_render_context<Arrangement> ctx(m_ctx, m_bbox, cache);
    Arr_bounded_approximate_face<Arrangement> bounded_approx_face(ctx);
    for(Face_const_handle fh = m_ctx.m_arr.faces_begin(); fh != m_ctx.m_arr.faces_end(); ++fh) bounded_approx_face(fh);

    return cache;
  }

private:
  const Render_context& m_ctx;
  const Bbox_2 m_bbox;
};

} // namespace draw_aos
} // namespace CGAL

#endif