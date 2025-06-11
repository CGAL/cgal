#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H

#include "CGAL/Arr_trapezoid_ric_point_location.h"
#include "CGAL/Bbox_2.h"
#include "CGAL/Draw_aos/Arr_approximation_cache.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include "CGAL/Draw_aos/Arr_bounded_approximate_curve_2.h"
#include "CGAL/Draw_aos/Arr_bounded_approximate_face_2.h"
#include "CGAL/Draw_aos/Arr_render_context.h"
#include <CGAL/Draw_aos/helpers.h>
#include <CGAL/Draw_aos/Arr_portals.h>
#include <CGAL/IO/Color.h>
#include <vector>

namespace CGAL {
/**
 * @brief Render arrangement on surface within a bounding box.
 *
 * @note The class is not thread-safe.
 */
class Arr_bounded_renderer
{
  using Color = IO::Color;
  using Point_2 = Geom_traits::Point_2;
  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Edge_const_handle = Arrangement::Edge_const_iterator;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Vertex_handle = Arrangement::Vertex_handle;
  using Halfedge_handle = Arrangement::Halfedge_handle;
  using Face_handle = Arrangement::Face_handle;
  using X_monotone_curve_2 = Geom_traits::X_monotone_curve_2;
  using Point_Location = Arr_trapezoid_ric_point_location<Arrangement>;
  using Feature_const = std::variant<Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;
  using Feature_const_vector = std::vector<Feature_const>;
  using Feature_const_vector_inserter = std::back_insert_iterator<Feature_const_vector>;
  using Arr_approx_geom_traits = Arr_approximation_geometry_traits;
  using Point_geom = Arr_approx_geom_traits::Point_geom;
  using Polyline_geom = Arr_approx_geom_traits::Polyline_geom;
  using Feature_portal_map = Arr_portals::Feature_portals_map;
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

  struct Execution_context : Arr_context_delegator<Arr_bounded_render_context>
  {
    Execution_context(const Arr_bounded_render_context& ctx)
        : Arr_context_delegator(ctx)
        , bounded_approx_pt(ctx)
        , bounded_approx_curve(ctx, bounded_approx_pt)
        , bounded_approx_face(ctx, bounded_approx_pt, bounded_approx_curve) {}

    const Arr_bounded_approximate_point_2 bounded_approx_pt;
    const Arr_bounded_approximate_curve_2 bounded_approx_curve;
    const Arr_bounded_approximate_face_2 bounded_approx_face;
  };

private:
  template <typename OutputIterator>
  static void locate_intersecting_features(const Execution_context& ctx,
                                           const X_monotone_curve_2& curve,
                                           OutputIterator out_iter,
                                           Arr_flags<Feature_type> feats) {
    using Feature = std::variant<Vertex_handle, Halfedge_handle, Face_handle>;
    using Feature_vector = std::vector<Feature>;

    auto func_out_iter = boost::make_function_output_iterator([&out_iter, &feats](const Feature& feature) {
      if(auto* vh = std::get_if<Vertex_handle>(&feature)) {
        if(!feats.is_set(Feature_type::Vertex) || (*vh)->is_at_open_boundary()) {
          return;
        }
        *out_iter++ = Vertex_const_handle(*vh);
      } else if(auto* he = std::get_if<Halfedge_handle>(&feature)) {
        if(!feats.is_set(Feature_type::Halfedge) || (*he)->is_fictitious()) {
          return;
        }
        *out_iter++ = Halfedge_const_handle(*he);
      } else if(auto* fh = std::get_if<Face_handle>(&feature)) {
        if(!feats.is_set(Feature_type::Face)) {
          return;
        }
        *out_iter++ = Face_const_handle(*fh);
      } else {
        CGAL_assertion(false && "Unexpected feature type");
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
      auto he = *inner_ccb;
      auto inner_face = he->twin()->face();

      bool is_degenerate = inner_face == fh;
      bool within_bounds = ctx->strictly_contains(inner_face->outer_ccb()->source()->point());
      if(is_degenerate || !within_bounds) {
        continue;
      }
      discover_faces(ctx, inner_face);
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

      bool within_bounds = ctx->strictly_contains(adj_face->outer_ccb()->source()->point());
      if(!within_bounds) {
        continue;
      }
      // For a face that is not one of the seeding faces,
      // the face must be contained iff one of its vertices is contained.
      // For a face that is, we'll touch it sooner or later.
      discover_faces(ctx, adj_face);
    } while(++circ != outer_ccb);
  }

public:
  Arr_bounded_renderer(Arr_render_context& ctx, Bbox_2 bbox)
      : m_ctx(ctx)
      , m_bbox(bbox) {}

  Arr_approximation_cache render() const {
    Arr_approximation_cache cache;

    if(m_ctx.is_cancelled()) {
      return cache;
    }

    Execution_context ctx(Arr_bounded_render_context(m_ctx, m_bbox, cache));

    auto insert_features = boost::make_function_output_iterator([&ctx](const Feature_const& feature) {
      if(auto* vh = std::get_if<Vertex_const_handle>(&feature)) {
        ctx.bounded_approx_pt(*vh);
      } else if(auto* he = std::get_if<Halfedge_const_handle>(&feature)) {
        ctx.bounded_approx_curve(*he);
      } else if(auto* fh = std::get_if<Face_const_handle>(&feature)) {
        discover_faces(ctx, *fh);
      } else {
        CGAL_assertion(false && "Unexpected feature type");
      }
    });

    auto top = ctx->cst_horizontal_segment(ctx->ymax(), ctx->xmin(), ctx->xmax());
    auto bottom = ctx->cst_horizontal_segment(ctx->ymin(), ctx->xmin(), ctx->xmax());
    auto left = ctx->cst_vertical_segment(ctx->xmin(), ctx->ymin(), ctx->ymax());
    auto right = ctx->cst_vertical_segment(ctx->xmax(), ctx->ymin(), ctx->ymax());

    // top, left are open edges while bottom, right are closed.
    locate_intersecting_features(ctx, top, insert_features, Feature_type::Face);
    locate_intersecting_features(ctx, left, insert_features, Feature_type::Face);
    locate_intersecting_features(ctx, bottom, insert_features, {Feature_type::Face, Feature_type::Halfedge});
    locate_intersecting_features(ctx, right, insert_features,
                                 {Feature_type::Face, Feature_type::Halfedge, Feature_type::Vertex});

    return cache;
  }

private:
  Arr_render_context m_ctx;
  const Bbox_2 m_bbox;
};

} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_BOUNDED_RENDERER_H