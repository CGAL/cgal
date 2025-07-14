#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H

#include <vector>

#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Draw_aos/Arr_approximation_cache.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_curve_2.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_face_2.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_portals.h>

namespace CGAL {
namespace draw_aos {

/**
 * @brief Render arrangement on surface within a bounding box.
 *
 * @note The class is not thread-safe.
 */
template <typename Arrangement>
class Arr_bounded_renderer
{
  using Color = IO::Color;
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using FT = typename Traits_adaptor<Geom_traits>::FT;
  using Point_2 = typename Traits_adaptor<Geom_traits>::Point_2;
  using X_monotone_curve_2 = typename Traits_adaptor<Geom_traits>::X_monotone_curve_2;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Edge_const_handle = typename Arrangement::Edge_const_iterator;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_handle = typename Arrangement::Vertex_handle;
  using Halfedge_handle = typename Arrangement::Halfedge_handle;
  using Face_handle = typename Arrangement::Face_handle;
  using Point_Location = Arr_trapezoid_ric_point_location<Arrangement>;
  using Feature_const = std::variant<Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;
  using Feature_const_vector = std::vector<Feature_const>;
  using Feature_const_vector_inserter = std::back_insert_iterator<Feature_const_vector>;
  using Approx_traits = Arr_approximation_geometry_traits<Geom_traits>;
  using Point_geom = typename Approx_traits::Point_geom;
  using Polyline_geom = typename Approx_traits::Polyline_geom;
  using Feature_portal_map = typename Arr_portals<Arrangement>::Feature_portals_map;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;
  using Render_context = Arr_render_context<Arrangement>;
  using Bounded_approx_point_2 = Arr_bounded_approximate_point_2<Arrangement>;
  using Bounded_approx_curve_2 = Arr_bounded_approximate_curve_2<Arrangement>;
  using Bounded_approx_face_2 = Arr_bounded_approximate_face_2<Arrangement>;
  using Approx_cache = Arr_approximation_cache<Arrangement>;

  // QFlags implement this pattern better, but we try not to reply on Qt classes.
  template <typename E>
  class Arr_flags
  {
  public:
    Arr_flags(const E& flags)
        : m_flags(static_cast<int>(flags)) {}
    Arr_flags(std::initializer_list<E> flags)
        : m_flags(0) {
      for(const auto& flag : flags) {
        m_flags |= static_cast<int>(flag);
      }
    }

    bool is_set(E flag) const { return (m_flags & static_cast<int>(flag)) != 0; }

  private:
    int m_flags;
  };

  enum class Feature_type {
    Vertex = 1 << 0,
    Halfedge = 1 << 1,
    Face = 1 << 2,
  };

  struct Execution_context : Arr_context_delegator<Bounded_render_context>
  {
    Execution_context(const Bounded_render_context& ctx)
        : Arr_context_delegator<Bounded_render_context>(ctx)
        , bounded_approx_pt(ctx)
        , bounded_approx_curve(ctx)
        , bounded_approx_face(ctx, bounded_approx_pt, bounded_approx_curve) {}

    const Bounded_approx_point_2 bounded_approx_pt;
    const Bounded_approx_curve_2 bounded_approx_curve;
    const Bounded_approx_face_2 bounded_approx_face;
  };

private:
  static Face_const_handle inbound_face_of_edge(const Halfedge_const_handle& he, Side_of_boundary side) {
    bool is_left_to_right = he->direction() == ARR_LEFT_TO_RIGHT;
    switch(side) {
    case Side_of_boundary::Top:
    case Side_of_boundary::Left:
      return is_left_to_right ? he->twin()->face() : he->face();
    case Side_of_boundary::Bottom:
    case Side_of_boundary::Right:
      return is_left_to_right ? he->face() : he->twin()->face();
    default:
      CGAL_assertion(false && "Invalid side of boundary");
    }
    return Face_const_handle();
  }
  static void approx_intersecting_features(Execution_context& ctx,
                                           const X_monotone_curve_2& curve,
                                           Side_of_boundary side,
                                           Arr_flags<Feature_type> feats) {
    CGAL_assertion(side != Side_of_boundary::None);
    using Feature = std::variant<Vertex_handle, Halfedge_handle, Face_handle>;
    using Feature_vector = std::vector<Feature>;

    auto func_out_iter = boost::make_function_output_iterator([&ctx, feats, side](const Feature& feature) {
      if(auto* vh = std::get_if<Vertex_handle>(&feature)) {
        if(!feats.is_set(Feature_type::Vertex) || (*vh)->is_at_open_boundary()) {
          return;
        }
        ctx.bounded_approx_pt(*vh);
      } else if(auto* he = std::get_if<Halfedge_handle>(&feature)) {
        if(!feats.is_set(Feature_type::Halfedge) || (*he)->is_fictitious()) {
          return;
        }
        ctx.bounded_approx_curve(*he);

        // There's a rare case that some faces' outer ccb overlaps with the bounding box edges, causing
        // no face to be discovered. We need to approximate the inbound face incident to the halfedge
        discover_faces(ctx, inbound_face_of_edge(*he, side));
      } else if(auto* fh = std::get_if<Face_handle>(&feature)) {
        if(!feats.is_set(Feature_type::Face)) {
          return;
        }
        discover_faces(ctx, *fh);
      }
    });

    using Zone_visitor = Arr_compute_zone_visitor<Arrangement, decltype(func_out_iter)>;
    using Zone_2 = Arrangement_zone_2<Arrangement, Zone_visitor>;

    Zone_visitor zone_visitor(func_out_iter);
    // use const_cast to fit in the interface but we will not modify the arrangement
    Zone_2 zone(const_cast<Arrangement&>(ctx->arr), &zone_visitor);
    zone.init(curve, ctx->point_location);
    zone.compute_zone();
  }

  static void discover_faces(Execution_context& ctx, const Face_const_handle& fh) {
    if(ctx->cache.has(fh)) {
      // The face has already been discovered and approximated.
      return;
    }
    ctx.bounded_approx_face(fh);

    for(auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
      auto circ = *inner_ccb;

      // We have to traverse the entire inner ccb instead of using arbitary vertex,
      // because the inner ccb may contain degenerate edges.
      do {
        auto inner_face = circ->twin()->face();
        if(inner_face != fh && ctx->strictly_contains(circ->source()->point())) {
          discover_faces(ctx, inner_face);
        }
      } while(++circ != *inner_ccb);
    }

    // We don't need to handle isolated vertices.
    // for(...)

    if(!fh->has_outer_ccb()) {
      // The unbounded face of bounded arrangements has no outer CCB.
      return;
    }
    // find fully contained adjacent faces
    auto outer_ccb = fh->outer_ccb();
    auto circ = outer_ccb;
    do {
      auto adj_face = circ->twin()->face();
      if(adj_face->is_fictitious() || adj_face->is_unbounded()) {
        // Unbounded faces cannot be contained in out bounding box.
        continue;
      }

      if(!ctx->strictly_contains(adj_face->outer_ccb()->source()->point())) {
        continue;
      }
      // For a face that is not one of the seeding faces,
      // the face must be contained iff one of its vertices is contained.
      // For a face that is, we'll touch it sooner or later.
      discover_faces(ctx, adj_face);
    } while(++circ != outer_ccb);
  }

public:
  Arr_bounded_renderer(const Render_context& ctx, Bbox_2 bbox)
      : m_ctx(ctx)
      , to_ft()
      , m_bbox(bbox) {}

  Approx_cache render() const {
    Approx_cache cache;
    cache.reserve_face_cache(m_ctx.arr.number_of_faces());
    cache.reserve_halfedge_cache(m_ctx.arr.number_of_halfedges());
    cache.reserve_vertex_cache(m_ctx.arr.number_of_vertices());

    if(m_ctx.is_cancelled()) {
      return cache;
    }

    Bounded_render_context bounded_ctx(m_ctx, m_bbox, cache);
    Execution_context ctx(bounded_ctx);

    auto top = ctx->cst_horizontal_segment(to_ft(ctx->ymax()), to_ft(ctx->xmin()), to_ft(ctx->xmax()));
    auto bottom = ctx->cst_horizontal_segment(to_ft(ctx->ymin()), to_ft(ctx->xmin()), to_ft(ctx->xmax()));
    auto left = ctx->cst_vertical_segment(to_ft(ctx->xmin()), to_ft(ctx->ymin()), to_ft(ctx->ymax()));
    auto right = ctx->cst_vertical_segment(to_ft(ctx->xmax()), to_ft(ctx->ymin()), to_ft(ctx->ymax()));

    // top, left are open edges while bottom, right are closed.
    approx_intersecting_features(ctx, top, Side_of_boundary::Top, Feature_type::Face);
    approx_intersecting_features(ctx, left, Side_of_boundary::Left, Feature_type::Face);
    approx_intersecting_features(ctx, bottom, Side_of_boundary::Bottom, {Feature_type::Face, Feature_type::Halfedge});
    approx_intersecting_features(ctx, right, Side_of_boundary::Right,
                                 {Feature_type::Face, Feature_type::Halfedge, Feature_type::Vertex});

    return cache;
  }

private:
  const Render_context& m_ctx;
  const Construct_coordinate<Geom_traits> to_ft;
  const Bbox_2 m_bbox;
};

} // namespace draw_aos
} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H