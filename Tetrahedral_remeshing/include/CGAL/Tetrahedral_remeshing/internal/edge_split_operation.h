#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/split_long_edges.h>

#include <optional>
#include <iostream>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {



template<typename C3t3, typename SizingFunction, typename CellSelector>
class EdgeSplitOperation 
  : public ElementaryOperation<C3t3,
                          std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                          std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using Base = ElementaryOperation<C3t3,
                              std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                              std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,
                              typename C3t3::Triangulation::Cell_handle>;
  using Complex = C3t3;
  using Triangulation = typename C3t3::Triangulation;
  using Edge = typename Triangulation::Edge;
  using ElementType = typename Base::ElementType;
  using ElementSource = typename Base::ElementSource;
  using VertexPair = std::pair<typename C3t3::Triangulation::Vertex_handle, typename C3t3::Triangulation::Vertex_handle>;
  using Cell_handle = typename Triangulation::Cell_handle;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Point = typename Triangulation::Point;
  using Lock_zone = typename Base::Lock_zone;
  using Facet = typename Triangulation::Facet;
  using Subdomain_index = typename C3t3::Subdomain_index;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Curve_index = typename C3t3::Curve_index;
  using FT = typename Triangulation::Geom_traits::FT;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;

public:
  EdgeSplitOperation(const SizingFunction& sizing,
                    const CellSelector& cell_selector,
                    const bool protect_boundaries)
    : m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
  {}

  bool should_process_element(const ElementType& vp, const Complex& c3t3) const override {
    // Convert vertex pair to edge for existing logic
    Cell_handle cell;
    int i, j;
    if (!is_edge_uv(vp.first, vp.second, get_incident_cells(vp.first, c3t3), cell, i, j)) {
      return false; // Edge doesn't exist
    }
    
    Edge edge(cell, i, j);
    
    // Direct wrapper around existing logic from split_long_edges.h
    auto [splittable, boundary] = can_be_split(edge, c3t3, m_protect_boundaries, m_cell_selector);
    if (!splittable)
      return false;

    const auto sqlen = is_too_long(edge, boundary, m_sizing, c3t3, m_cell_selector);
    if (sqlen.has_value()) {
      return true;
    }
    
    return false;
  }

  ElementSource get_element_source(const C3t3& c3t3) const override {
    // Collect long edges as vertex pairs
    std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>> long_edges;
    get_long_edges(c3t3, m_sizing, m_cell_selector, std::back_inserter(long_edges));
    return long_edges;
  }

  bool can_apply_operation(const ElementType& vp, const Complex& c3t3) const override {
    // This corresponds to lines 409-416 in split_long_edges.h
    // The algorithm needs to re-validate that the edge can still be split
    // because the triangulation may have been modified by previous operations
    
    const Triangulation& tr = c3t3.triangulation();
    Cell_handle cell;
    int i1, i2;
    
    // Check if the edge still exists in the triangulation
    if (!tr.tds().is_edge(vp.first, vp.second, cell, i1, i2)) {
      // Edge no longer exists (was modified by previous operations)
      return false;
    }
    
    // Create the actual edge from the found cell and indices
    Edge actual_edge(cell, i1, i2);
    
    // Re-check splittability (equivalent to the can_be_split call in line 414)
    auto [splittable, boundary] = can_be_split(actual_edge, c3t3, m_protect_boundaries, m_cell_selector);
    
    return splittable;
  }

  Lock_zone get_lock_zone(const ElementType& vp, const Complex& c3t3) const override {
    Lock_zone zone;
    return zone;
  }

  bool execute_operation(const ElementType& vp, Complex& c3t3) override {
    // Combines all three phases: pre-operation + operation + post-operation
    // Corresponds to the complete sequence in split_long_edges.h lines 421-425:
    //   visitor.before_split(tr, edge);
    //   Vertex_handle vh = split_edge(edge, cell_selector, c3t3);
    //   if(vh != Vertex_handle()) visitor.after_split(tr, vh);
    
    // Find the actual edge in the triangulation
    const Triangulation& tr = c3t3.triangulation();
    Cell_handle cell;
    int i1, i2;
    
    // Since can_apply_operation already validated the edge exists, we can safely call is_edge
    if (!tr.tds().is_edge(vp.first, vp.second, cell, i1, i2)) {
      // Edge was modified between can_apply_operation and execute_operation
      return false;
    }
    
    Edge actual_edge(cell, i1, i2);
    
    // PRE-OPERATION: Call visitor.before_split
    // m_visitor.before_split(tr, actual_edge);
    
    // MAIN OPERATION: Split the edge
    Vertex_handle created_vertex = split_edge(actual_edge, m_cell_selector, c3t3);
    
    // POST-OPERATION: Call visitor.after_split if successful
    // if (created_vertex != Vertex_handle()) {
    //   m_visitor.after_split(tr, created_vertex);
    // }
    
    return (created_vertex != Vertex_handle());
  }

  std::string operation_name() const override {
    return "Edge Split";
  }

};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_OPERATION_H 