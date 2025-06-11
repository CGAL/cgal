#ifndef CGAL_DRAW_AOS_ARR_CREATE_PORTALS_H
#define CGAL_DRAW_AOS_ARR_CREATE_PORTALS_H

#include "CGAL/Arr_vertical_decomposition_2.h"
#include "CGAL/Draw_aos/Arr_approximate_point_2.h"
#include "CGAL/Draw_aos/Arr_construct_segments.h"
#include "CGAL/Draw_aos/helpers.h"
#include "CGAL/Draw_aos/Arr_graph_conn.h"
#include "CGAL/Object.h"
#include "CGAL/basic.h"
#include "CGAL/unordered_flat_map.h"
#include <boost/iterator/function_output_iterator.hpp>
#include <CGAL/Draw_aos/helpers.h>
#include <CGAL/Draw_aos/Arr_approximation_geometry_traits.h>
#include <limits>
#include <utility>

namespace CGAL {
/**
 * @brief Portals are virtual vertical segments that connect the outer
 * connected component boundary (ccb) of a face with its inner ccbs.
 *
 * We use vertical decomposition to cast upward rays in a batched manner.
 * For each inner ccb, only one "portal" is created at a vertex whose upward ray
 * hits another ccb. Eventually, faces of the arrangement become hole-free and can
 * be drawn with Graphics_scene.
 */
class Arr_portals
{

  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Point_2 = Geom_traits::Point_2;
  using Approx_point = Arr_approximation_geometry_traits::Approx_point;
  using X_monotone_curve_2 = Geom_traits::X_monotone_curve_2;
  using Feature_const = std::variant<Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;

public:
  // Pair composed of a portal entry point and the portaled vertex on the inner ccb.
  using Portal = std::pair<std::optional<Approx_point>, Vertex_const_handle>;
  using Portal_vector = std::vector<Portal>;
  // Map from a feature to its portals sorted by the x coordinate of the virtual vertical segments.
  using Feature_portals_map = unordered_flat_map<Feature_const, Portal_vector>;

private:
  // Use this function to locate the intersection point of a vertical ray shooted from a point
  // to it's upper x-monotone curve(it can't be vertical, or we should've got its min vertex from vertical
  // decomposition).
  static Point_2 upper_intersection(const Geom_traits& traits, const Point_2& pt, const X_monotone_curve_2& curve) {
    Arr_construct_vertical_segment cst_vertical_segment(traits);
    Point_2 intersection_point;
    auto intersect = traits.intersect_2_object();
    auto vertical_line = cst_vertical_segment(pt.x(), pt.y(), std::numeric_limits<double>::max());
    bool found_intersection = false;

    using Multiplicity = Geom_traits::Multiplicity;
    using Intersect_point = std::pair<Point_2, Multiplicity>;
    using Intersect_curve = X_monotone_curve_2;
    using Intersect_type = std::variant<Intersect_point, Intersect_curve>;

    intersect(curve, vertical_line, boost::make_function_output_iterator([&](const Intersect_type& res) {
                found_intersection = true;
                if(std::holds_alternative<Intersect_point>(res)) {
                  intersection_point = std::get<Intersect_point>(res).first;
                  return;
                }
                CGAL_assertion(false && "Unexpected intersection type");
              }));

    return found_intersection ? intersection_point : Point_2(pt.x(), std::numeric_limits<double>::max());
  }

public:
  Feature_portals_map create(const Arrangement& arr) const {
    using Object_pair = std::pair<Object, Object>;
    using Vert_decomp_entry = std::pair<Vertex_const_handle, Object_pair>;

    Arr_graph_conn conn(arr);
    auto visited_ccbs = std::unordered_set<Vertex_const_handle>();
    Feature_portals_map feature_portals;
    auto intersect = arr.traits()->intersect_2_object();
    auto approx_pt = Arr_approximate_point_2<Geom_traits>(*arr.traits());

    auto func_out_iter = boost::make_function_output_iterator([&](const Vert_decomp_entry& entry) {
      const auto& [vh, obj_pair] = entry;

      const auto& ccb_main_vertex = conn.ccb_representative_vertex(vh);
      if(visited_ccbs.find(ccb_main_vertex) != visited_ccbs.end()) {
        // This ccb has already been processed
        return;
      }

      const auto& above_feat = obj_pair.second;
      if(Vertex_const_handle above_vh; CGAL::assign(above_vh, above_feat)) {
        if(conn.is_connected(above_vh, vh)) {
          // This upper vertex is connected to vh, skip it
          return;
        }

        const auto& [it, _] = feature_portals.try_emplace(above_vh, std::vector<Portal>{});
        if((above_vh)->is_at_open_boundary()) {
          it->second.emplace_back(std::nullopt, vh);
        } else {
          it->second.emplace_back(approx_pt((above_vh)->point()), vh);
        }
      } else if(Halfedge_const_handle above_he; CGAL::assign(above_he, above_feat)) {
        if(conn.is_connected((above_he)->source(), vh)) {
          return;
        }

        const auto& [it, _] = feature_portals.try_emplace(above_he, std::vector<Portal>{});
        if((above_he)->is_fictitious()) {
          it->second.emplace_back(std::nullopt, vh);
          return;
        }
        it->second.emplace_back(approx_pt(upper_intersection(*arr.traits(), vh->point(), (above_he)->curve())), vh);
      } else if(Face_const_handle above_fh; CGAL::assign(above_fh, above_feat)) {
        // We don't create portals for the unbounded face in bounded arrangements.
        CGAL_assertion(above_fh->is_unbounded() && !above_fh->has_outer_ccb());
      } else {
        return;
      }

      visited_ccbs.insert(ccb_main_vertex);
    });

    decompose(arr, func_out_iter);
    return feature_portals;
  }
};
} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_CREATE_PORTALS_H