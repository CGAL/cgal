
#ifndef CGAL_DRAW_AOS_ARR_GRAPH_CONN_H
#define CGAL_DRAW_AOS_ARR_GRAPH_CONN_H
#include <boost/range/iterator_range.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/Union_find.h>
#include <CGAL/unordered_flat_map.h>

namespace CGAL {
/**
 * @brief Arr_graph_conn provides fast connectivity queries for arrangement vertices
 * based on union-find data structure.
 */
template <typename Arr>
class Arr_graph_conn
{
  using Vertex_const_handle = typename Arr::Vertex_const_handle;
  using Edge_const_handle = typename Arr::Edge_const_iterator;
  using Halfedge_const_handle = typename Arr::Halfedge_const_iterator;
  using Union_find_handle = typename Union_find<Vertex_const_handle>::handle;

private:
  void insert_halfedge(Halfedge_const_handle he) {
    const auto& source = he->source();
    const auto& target = he->target();

    auto [source_it, source_inserted] = m_lookup.try_emplace(source, Union_find_handle());
    if(source_inserted) {
      source_it->second = m_uf.make_set(source);
    }
    auto [target_it, target_inserted] = m_lookup.try_emplace(target, Union_find_handle());
    if(target_inserted) {
      target_it->second = m_uf.make_set(target);
    }

    m_uf.unify_sets(source_it->second, target_it->second);
  }

  void insert_isolated_vertex(const Vertex_const_handle& vh) { m_lookup[vh] = m_uf.make_set(vh); }

public:
  Arr_graph_conn(const Arr& arr) {
    m_lookup.reserve(arr.number_of_vertices());

    for(const auto& he : arr.halfedge_handles()) {
      if(he->direction() != ARR_LEFT_TO_RIGHT) {
        continue;
      }
      insert_halfedge(he);
    }

    for(const auto& vh : arr.vertex_handles()) {
      if(!vh->is_isolated()) {
        continue;
      }
      insert_isolated_vertex(vh);
    }

    // add fictitious edges and open vertices
    for(auto fh = arr.unbounded_faces_begin(); fh != arr.unbounded_faces_end(); ++fh) {
      if(!fh->has_outer_ccb()) {
        continue;
      }
      auto outer_ccb = fh->outer_ccb();
      auto curr = outer_ccb;
      do {
        if(!curr->is_fictitious()) {
          continue;
        }
        insert_halfedge(curr);
      } while(++curr != outer_ccb);
    }
  }

  bool is_connected(const Vertex_const_handle& v1, const Vertex_const_handle& v2) const {
    auto it1 = m_lookup.find(v1);
    auto it2 = m_lookup.find(v2);
    // This can't happen
    CGAL_assertion(it1 != m_lookup.end() && it2 != m_lookup.end());
    return m_uf.same_set(it1->second, it2->second);
  }

  /**
   * @brief Returns the representative vertex of the connected component containing the given vertex.
   * For each connected component boundary (CCB), the same representative vertex is consistently returned.
   * @param vh a vertex handle in the ccb
   * @return Vertex_const_handle
   */
  Vertex_const_handle ccb_representative_vertex(const Vertex_const_handle& vh) const {
    // path compression typically does not mutate the representative member of a set.
    auto it = m_lookup.find(vh);
    CGAL_assertion(it != m_lookup.end());
    return *m_uf.find(it->second);
  }

private:
  Union_find<Vertex_const_handle> m_uf;
  unordered_flat_map<Vertex_const_handle, Union_find_handle> m_lookup;
};
} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_GRAPH_CONN_H