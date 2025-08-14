#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H
#include <cstddef>

#include <boost/range/iterator_range.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/**
 * @brief Cache class for approximating arrangement on surface.
 *
 * When iterating over the arrangement dcel, a feature(vertex, halfedge, face) might be visited multiple times.
 * This cache stores the approximated geometry for each feature to avoid redundant calculations.
 * @tparam Arrangement
 */
template <typename Arrangement>
class Arr_approximation_cache
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;

public:
  using Vertex_cache_obj = typename Approx_traits::Point_geom;
  using Halfedge_cache_obj = typename Approx_traits::Polyline_geom;
  using Face_cache_obj = typename Approx_traits::Triangulated_face;

private:
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Edge_const_handle = typename Arrangement::Edge_const_iterator;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_iterator;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_cache = unordered_flat_map<Vertex_const_handle, Vertex_cache_obj>;
  using Halfedge_cache = unordered_flat_map<Halfedge_const_handle, Halfedge_cache_obj>;
  using Face_cache = unordered_flat_map<Face_const_handle, Face_cache_obj>;

  using Vertex_cache_const_iterator = typename Vertex_cache::const_iterator;
  using Halfedge_cache_const_iterator = typename Halfedge_cache::const_iterator;
  using Face_cache_const_iterator = typename Face_cache::const_iterator;
  using Vertex_cache_range = boost::iterator_range<Vertex_cache_const_iterator>;
  using Halfedge_cache_range = boost::iterator_range<Halfedge_cache_const_iterator>;
  using Face_cache_range = boost::iterator_range<Face_cache_const_iterator>;

public:
  Arr_approximation_cache() = default;

  void reserve_vertices(std::size_t size) { m_vertices.reserve(size); }
  void reserve_halfedges(std::size_t size) { m_halfedges.reserve(size); }
  void reserve_faces(std::size_t size) { m_faces.reserve(size); }

  std::pair<Vertex_cache_obj&, bool> try_emplace(const Vertex_const_handle& vh) {
    const auto& [it, inserted] = m_vertices.try_emplace(vh, Vertex_cache_obj());
    return {it->second, inserted};
  }
  std::pair<Halfedge_cache_obj&, bool> try_emplace(const Halfedge_const_handle& he) {
    const auto& [it, inserted] = m_halfedges.try_emplace(he, Halfedge_cache_obj());
    return {it->second, inserted};
  }
  std::pair<Face_cache_obj&, bool> try_emplace(const Face_const_handle& fh) {
    const auto& [it, inserted] = m_faces.try_emplace(fh, Face_cache_obj());
    return {it->second, inserted};
  }

  std::pair<const Vertex_cache_obj&, bool> get(const Vertex_const_handle& vh) const {
    auto it = m_vertices.find(vh);
    if(it != m_vertices.end()) {
      return {it->second, true};
    }
    return {Vertex_cache_obj(), false};
  }
  std::pair<const Halfedge_cache_obj&, bool> get(const Halfedge_const_handle& he) const {
    auto it = m_halfedges.find(he);
    if(it != m_halfedges.end()) {
      return {it->second, true};
    }
    return {Halfedge_cache_obj(), false};
  }
  std::pair<const Face_cache_obj&, bool> get(const Face_const_handle& fh) const {
    auto it = m_faces.find(fh);
    if(it != m_faces.end()) {
      return {it->second, true};
    }
    return {Face_cache_obj(), false};
  }

  bool has(const Vertex_const_handle& vh) const { return m_vertices.find(vh) != m_vertices.end(); }
  bool has(const Halfedge_const_handle& he) const { return m_halfedges.find(he) != m_halfedges.end(); }
  bool has(const Face_const_handle& fh) const { return m_faces.find(fh) != m_faces.end(); }

  Vertex_cache_const_iterator vertex_begin() const { return m_vertices.begin(); }
  Vertex_cache_const_iterator vertex_end() const { return m_vertices.end(); }
  Vertex_cache_range vertices() const { return boost::make_iterator_range(vertex_begin(), vertex_end()); }
  std::size_t number_of_vertices() const { return m_vertices.size(); }

  Halfedge_cache_const_iterator halfedge_begin() const { return m_halfedges.begin(); }
  Halfedge_cache_const_iterator halfedge_end() const { return m_halfedges.end(); }
  Halfedge_cache_range halfedges() const { return boost::make_iterator_range(halfedge_begin(), halfedge_end()); }
  std::size_t number_of_halfedges() const { return m_halfedges.size(); }

  Face_cache_const_iterator face_begin() const { return m_faces.begin(); }
  Face_cache_const_iterator face_end() const { return m_faces.end(); }
  Face_cache_range faces() const { return boost::make_iterator_range(face_begin(), face_end()); }
  std::size_t number_of_faces() const { return m_faces.size(); }

private:
  Vertex_cache m_vertices;
  Halfedge_cache m_halfedges;
  Face_cache m_faces;
};

} // namespace draw_aos
} // namespace CGAL
#endif