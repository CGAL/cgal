#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H
#include "CGAL/Arr_enums.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include "CGAL/Draw_aos/helpers.h"
#include "CGAL/unordered_flat_map.h"
#include <boost/range/iterator_range.hpp>

namespace CGAL {
class Arr_approximation_cache
{
  using Approx_geom_traits = Arr_approximation_geometry_traits;

public:
  using Vertex_cache_obj = Approx_geom_traits::Point_geom;
  using Halfedge_cache_obj = Approx_geom_traits::Polyline_geom;
  using Face_cache_obj = Approx_geom_traits::Triangulated_face;

private:
  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Edge_const_handle = Arrangement::Edge_const_iterator;
  using Halfedge_const_handle = Arrangement::Halfedge_const_iterator;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Vertex_cache = unordered_flat_map<Vertex_const_handle, Vertex_cache_obj>;
  using Halfedge_cache = unordered_flat_map<Halfedge_const_handle, Halfedge_cache_obj>;
  using Face_cache = unordered_flat_map<Face_const_handle, Face_cache_obj>;

  using Vertex_cache_const_iterator = typename Vertex_cache::const_iterator;
  using Halfedge_cache_const_iterator = typename Halfedge_cache::const_iterator;
  using Face_cache_const_iterator = typename Face_cache::const_iterator;
  using Vertex_cache_range = boost::iterator_range<Vertex_cache_const_iterator>;
  using Halfedge_cache_range = boost::iterator_range<Halfedge_cache_const_iterator>;
  using Face_cache_range = boost::iterator_range<Face_cache_const_iterator>;

private:
  Halfedge_const_handle identify_halfedge(const Halfedge_const_handle& he) const {
    return he->direction() == ARR_RIGHT_TO_LEFT ? he : he->twin();
  }

public:
  std::pair<Vertex_cache_obj&, bool> try_emplace(const Vertex_const_handle& vh) {
    const auto& [it, inserted] = m_vertex_cache.try_emplace(vh, Vertex_cache_obj());
    return {it->second, inserted};
  }
  std::pair<Halfedge_cache_obj&, bool> try_emplace(const Halfedge_const_handle& he) {
    const auto& [it, inserted] = m_halfedge_cache.try_emplace(identify_halfedge(he), Halfedge_cache_obj());
    return {it->second, inserted};
  }
  std::pair<Face_cache_obj&, bool> try_emplace(const Face_const_handle& fh) {
    const auto& [it, inserted] = m_face_cache.try_emplace(fh, Face_cache_obj());
    return {it->second, inserted};
  }

  std::pair<const Vertex_cache_obj&, bool> get(const Vertex_const_handle& vh) const {
    auto it = m_vertex_cache.find(vh);
    if(it != m_vertex_cache.end()) {
      return {it->second, true};
    }
    return {Vertex_cache_obj(), false};
  }
  std::pair<const Halfedge_cache_obj&, bool> get(const Halfedge_const_handle& he) const {
    auto it = m_halfedge_cache.find(identify_halfedge(he));
    if(it != m_halfedge_cache.end()) {
      return {it->second, true};
    }
    return {Halfedge_cache_obj(), false};
  }
  std::pair<const Face_cache_obj&, bool> get(const Face_const_handle& fh) const {
    auto it = m_face_cache.find(fh);
    if(it != m_face_cache.end()) {
      return {it->second, true};
    }
    return {Face_cache_obj(), false};
  }

  bool has(const Vertex_const_handle& vh) const { return m_vertex_cache.find(vh) != m_vertex_cache.end(); }
  bool has(const Halfedge_const_handle& he) const {
    return m_halfedge_cache.find(identify_halfedge(he)) != m_halfedge_cache.end();
  }
  bool has(const Face_const_handle& fh) const { return m_face_cache.find(fh) != m_face_cache.end(); }

  Vertex_cache_const_iterator vertex_cache_begin() const { return m_vertex_cache.begin(); }
  Vertex_cache_const_iterator vertex_cache_end() const { return m_vertex_cache.end(); }
  Vertex_cache_range vertex_cache() const {
    return boost::make_iterator_range(vertex_cache_begin(), vertex_cache_end());
  }
  std::size_t vertex_cache_size() const { return m_vertex_cache.size(); }

  Halfedge_cache_const_iterator halfedge_cache_begin() const { return m_halfedge_cache.begin(); }
  Halfedge_cache_const_iterator halfedge_cache_end() const { return m_halfedge_cache.end(); }
  Halfedge_cache_range halfedge_cache() const {
    return boost::make_iterator_range(halfedge_cache_begin(), halfedge_cache_end());
  }
  std::size_t halfedge_cache_size() const { return m_halfedge_cache.size(); }

  Face_cache_const_iterator face_cache_begin() const { return m_face_cache.begin(); }
  Face_cache_const_iterator face_cache_end() const { return m_face_cache.end(); }
  Face_cache_range face_cache() const { return boost::make_iterator_range(face_cache_begin(), face_cache_end()); }
  std::size_t face_cache_size() const { return m_face_cache.size(); }

private:
  Vertex_cache m_vertex_cache;
  Halfedge_cache m_halfedge_cache;
  Face_cache m_face_cache;
};

} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H