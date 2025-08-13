#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_COLLAPSE_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_COLLAPSE_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/collapse_short_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <boost/container/small_vector.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/functional/hash.hpp>

#include <utility>
#include <optional>
#include <vector>
#include <iostream>
#include <tbb/concurrent_unordered_map.h>
#include <unordered_set>
#include <iterator>
#include <algorithm>


namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template<typename C3t3, typename SizingFunction, typename CellSelector>
class EdgeCollapseOperation :
    public ElementaryOperation<C3t3,
                          std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                          std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,

template<typename C3t3, typename SizingFunction, typename CellSelector, typename Visitor>
class EdgeCollapseOperation 
  : public ElementaryOperation<C3t3,
                              typename C3t3::Triangulation::Edge,
                              std::vector<typename C3t3::Triangulation::Edge>,
                              typename C3t3::Triangulation::Cell_handle>
{
public:
  using Base = ElementaryOperation<C3t3,
                              typename C3t3::Triangulation::Edge,
                              std::vector<typename C3t3::Triangulation::Edge>,
                              typename C3t3::Triangulation::Cell_handle>;
  using ElementType = typename Base::ElementType;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;
  using Triangulation = typename C3t3::Triangulation;
  using Edge = typename Triangulation::Edge;
  using Cell_handle = typename Triangulation::Cell_handle;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Point = typename Triangulation::Point;
  using FT = typename Triangulation::Geom_traits::FT;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  const Visitor& m_visitor;
  tbb::concurrent_unordered_map<Edge,bool> should_skip_edge;
  
  // Debug data member to store locked cells
  //thread_local static std::set<std::size_t> m_locked_cells_timestamps;

  // Helper to get incident cells for a vertex
  boost::container::small_vector<Cell_handle, 64> get_incident_cells(Vertex_handle vh, const C3t3& c3t3) const {
    boost::container::small_vector<Cell_handle, 64> inc_cells;
#ifdef USE_THREADSAFE_INCIDENT_CELLS
    c3t3.triangulation().incident_cells_threadsafe(vh, std::back_inserter(inc_cells));
#else
    c3t3.triangulation().incident_cells(vh, std::back_inserter(inc_cells));
#endif
    return inc_cells;
  }

public:
  EdgeCollapseOperation(const SizingFunction& sizing,
                       const CellSelector& cell_selector,
                       const bool protect_boundaries,
                       const Visitor& visitor)
    : m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
    , m_visitor(visitor)
  {}

  std::vector<ElementType> get_short_edges(const C3t3& c3t3) const {
    std::vector<std::pair<Edge,FT>> short_edges_with_length;
    const auto& tr = c3t3.triangulation();
    
    for (const Edge& e : tr.finite_edges()) {
      auto [collapsible, boundary] = can_be_collapsed(e, c3t3, m_protect_boundaries, m_cell_selector);
      if (!collapsible)
        continue;

      const auto sqlen = is_too_short(e, boundary, m_sizing, c3t3, m_cell_selector);
      if (sqlen != std::nullopt) {
        short_edges_with_length.push_back(std::make_pair(e, sqlen.value()));
      }
    }
    // Sort ascending: shortest first
    auto length_comp = [&](const std::pair<Edge, FT>& a, const std::pair<Edge, FT>& b) { return a.second < b.second; };
    std::stable_sort(short_edges_with_length.begin(), short_edges_with_length.end(), length_comp);

    std::vector<ElementType> short_edges;
    short_edges.reserve(short_edges_with_length.size());
    for (const auto& ef : short_edges_with_length) short_edges.push_back(ef.first);
    return short_edges;
  }

  ElementSource get_element_source(const C3t3& c3t3) const override {
    // Collect short edges
    return get_short_edges(c3t3);
  }

  bool lock_zone(const ElementType& edge, const C3t3& c3t3) const override {
    auto& tr = c3t3.triangulation();
    if(should_skip_edge.contains(edge)) {
      return true;
    }
    // Clear previous locked cells
    //m_locked_cells_timestamps.clear();
    const Vertex_handle v0 = edge.first->vertex(edge.second);
    const Vertex_handle v1 = edge.first->vertex(edge.third);
    if(!(tr.try_lock_vertex(v0) && tr.try_lock_vertex(v1))) {
      return false;
   }
    bool locked =true ;
    //Cell_handle c;
    //int i,j;
    //if(!tr.is_edge(v0, v1, c, i, j)) {

    //  return false;
    //}
    //// Lock all incident cells to both vertices (similar to edge split)
    #if 0
    std::vector<Cell_handle> inc_cells_0,inc_cells_1;
    std::vector<Vertex_handle> adj_vertices_0,adj_vertices_1;
    if(locked) {
      bool ok_adj0 = tr.try_lock_and_get_adjacent_vertices_and_cells_3(v0, std::back_inserter(adj_vertices_0), inc_cells_0);
      bool ok_adj1 = tr.try_lock_and_get_adjacent_vertices_and_cells_3(v1, std::back_inserter(adj_vertices_1), inc_cells_1);
      locked = ok_adj0 && ok_adj1;

  Lock_zone get_lock_zone(const ElementType& vp, const C3t3& c3t3) const override {
    Lock_zone zone;
    std::cout<<"EdgeCollapseOperation::get_lock_zone"<<std::endl;
    return zone;
      if(locked){
        for(const auto& adj_vertex : adj_vertices_0) {
            std::vector<Cell_handle> adj_inc_cells;
            locked = tr.try_lock_and_get_incident_cells(adj_vertex, adj_inc_cells);
            if(!locked){
              break;
            }
        }
        if(locked){
          for(const auto& adj_vertex : adj_vertices_1) {
              std::vector<Cell_handle> adj_inc_cells;
              locked = tr.try_lock_and_get_incident_cells(adj_vertex, adj_inc_cells);
              if(!locked){
                break;
              }
          }
        }
      }
    }
    
    return locked;
    
    #else
    if(!tr.is_vertex(v0) || !tr.is_vertex(v1)) {
      return false;
    }
    std::vector<Cell_handle> inc_cells_0,inc_cells_1;
    return tr.try_lock_and_get_incident_cells(v0, inc_cells_0) && tr.try_lock_and_get_incident_cells(v1, inc_cells_1);
    // We need to lock v individually first, to be sure v->cell() is valid
    //if(!tr.try_lock_vertex(v0) || !tr.try_lock_vertex(v1))
    //  return false;

    //// Cell_handle ch;
    //// int i0, i1;
    //// boost::container::small_vector<Cell_handle, 64> inc_vh = get_incident_cells(v0, c3t3);
    //// is_edge_uv(v0, v1, inc_vh, ch, i0, i1);
    ////   Edge edge_to_flip(ch, i0, i1);

    //C3t3::Triangulation::Cell_circulator ccirc(edge);
    //C3t3::Triangulation::Cell_circulator cdone = ccirc;

    //do {
    //  if(!tr.is_cell(ccirc) || !tr.try_lock_cell(ccirc)) // LOCK
    //    return false;
    //  ++ccirc;
    //} while(ccirc != cdone);
    //
    //return true;
          #endif
  }

  bool execute_pre_operation(const ElementType& vp, C3t3& c3t3) override {
    std::cout<<"EdgeCollapseOperation::execute_pre_operation"<<std::endl;
    return true;
    }

  }


  bool execute_operation(const ElementType& edge, C3t3& c3t3) override {
    auto& tr = c3t3.triangulation();
    if (should_skip_edge.contains(edge)) {
      return true;
    }
    // Get vertices for debug logging
    const Vertex_handle v0 = edge.first->vertex(edge.second);
    const Vertex_handle v1 = edge.first->vertex(edge.third);


    // Use existing collapse_edge function from collapse_short_edges.h
    Vertex_handle result = collapse_edge(edge, c3t3, m_sizing, m_protect_boundaries, 
                                        m_cell_selector, should_skip_edge, m_visitor);
    
    bool success = (result != Vertex_handle());


    return success;
  }

  std::string operation_name() const override {
    return "Edge Collapse";
  }
  bool requires_ordered_processing() const override {
    return true;
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_COLLAPSE_OPERATION_H 