#ifndef CGAL_DRAW_AOS_ARR_CREATE_PORTALS_H
#define CGAL_DRAW_AOS_ARR_CREATE_PORTALS_H
#include <utility>
#include <variant>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Object.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/Draw_aos/Arr_graph_conn.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/**
 * @brief Portals are virtual vertical segments that connect the outer
 * connected component boundary (ccb) of a face with its inner ccbs.
 *
 * We use vertical decomposition to cast upward rays in a batched manner.
 * For each inner ccb, only one "portal" is created at a vertex whose upward ray
 * hits another ccb. Eventually, faces of the arrangement become hole-free and can
 * be drawn with Graphics_scene.
 */
template <typename Arrangement>
class Arr_portals
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Point_2 = typename Geom_traits::Point_2;
  using X_monotone_curve_2 = typename Geom_traits::X_monotone_curve_2;
  using Approx_point = typename Arr_approximate_traits<Geom_traits>::Approx_point;
  using Feature_const = std::variant<Vertex_const_handle, Halfedge_const_handle, Face_const_handle>;

public:
  using Portal_exit = Vertex_const_handle;
  using Portal_exit_vector = std::vector<Portal_exit>;
  // Map from a feature to its portals sorted by the x coordinate of the virtual vertical segments.
  using Feature_portals_map = unordered_flat_map<Feature_const, Portal_exit_vector>;

public:
  Arr_portals(const Arrangement& arr) { init(arr); }

private:
  void init(const Arrangement& arr) {
    using Object_pair = std::pair<Object, Object>;
    using Vert_decomp_entry = std::pair<Vertex_const_handle, Object_pair>;

    Arr_graph_conn conn(arr);
    auto visited_ccbs = std::unordered_set<Vertex_const_handle>();

    auto func_out_iter = boost::make_function_output_iterator([&](const Vert_decomp_entry& entry) {
      const auto& [vh, obj_pair] = entry;
      const auto& above_feat = obj_pair.second;
      Vertex_const_handle ccb_main_vertex = conn.ccb_representative_vertex(vh);
      // Skip processed ccb.
      if(visited_ccbs.find(ccb_main_vertex) != visited_ccbs.end()) return;

      if(Vertex_const_handle above_vh; CGAL::assign(above_vh, above_feat)) {
        // Skip vertex connected to vh
        if(conn.is_connected(above_vh, vh)) return;
        const auto& [it, _] = m_map.try_emplace(above_vh, Portal_exit_vector{});
        it->second.emplace_back(vh);
      } else if(Halfedge_const_handle above_he; CGAL::assign(above_he, above_feat)) {
        // The given halfedge is always directed from right to left (exactly what we need).
        if(conn.is_connected((above_he)->source(), vh)) return;
        const auto& [it, _] = m_map.try_emplace(above_he, Portal_exit_vector{});
        it->second.emplace_back(vh);
      } else // Don't handle faces.
        return;

      visited_ccbs.insert(ccb_main_vertex);
    });

    decompose(arr, func_out_iter);

    // reverse portals on each halfedge, as they are stored in left-to-right order
    for(auto& [feat, portals] : m_map) {
      if(!std::holds_alternative<Halfedge_const_handle>(feat)) continue;
      std::reverse(portals.begin(), portals.end());
    }
  }

public:
  const Feature_portals_map& get() const { return m_map; }

private:
  Feature_portals_map m_map;
};

} // namespace draw_aos
} // namespace CGAL
#endif