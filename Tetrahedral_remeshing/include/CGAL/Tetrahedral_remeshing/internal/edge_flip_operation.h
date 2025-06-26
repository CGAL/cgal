#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H

#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Tetrahedral_remeshing/internal/flip_edges.h>

#include <boost/unordered_map.hpp>

#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <mutex>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

// Base class for shared functionality between internal and boundary edge flipping
template<typename C3t3, typename CellSelector, typename Visitor>
class EdgeFlipOperationBase {
protected:
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Cell_handle           Cell_handle;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename Tr::Edge                  Edge;
  typedef typename Tr::Facet                 Facet;
  typedef typename Tr::Geom_traits           Gt;
  typedef typename Gt::Point_3               Point_3;
  typedef typename Gt::FT                    FT;

  C3t3& m_c3t3;
  CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Visitor m_visitor;
  
  // Shared inc_cells data structure (used by both internal and boundary operations)
  static thread_local std::unordered_map<Vertex_handle, boost::container::small_vector<Cell_handle, 64>> s_inc_cells;

public:
  EdgeFlipOperationBase(C3t3& c3t3,
                       CellSelector& cell_selector,
                       const bool protect_boundaries,
                       const Visitor& visitor)
    : m_c3t3(c3t3)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
    , m_visitor(visitor)
  {}

  // Provide access to member variables for execution classes
  CellSelector& get_cell_selector() { return m_cell_selector; }
  bool get_protect_boundaries() const { return m_protect_boundaries; }

protected:

  // Helper to get incident cells for a vertex
  boost::container::small_vector<Cell_handle, 64> get_incident_cells(Vertex_handle vh, const C3t3& c3t3) const {
    boost::container::small_vector<Cell_handle, 64> inc_cells;
    c3t3.triangulation().incident_cells(vh, std::back_inserter(inc_cells));
    return inc_cells;
  }
};

// Internal Edge Flip Operation - processes vertex pairs like original get_internal_edges
template<typename C3t3, typename CellSelector, typename Visitor>
class InternalEdgeFlipOperation 
  : public EdgeFlipOperationBase<C3t3, CellSelector, Visitor>,
    public ElementaryOperation<C3t3, 
                          std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using BaseClass = EdgeFlipOperationBase<C3t3, CellSelector, Visitor>;
  using VertexPair = std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>;
  
  using Base = ElementaryOperation<C3t3, 
                              VertexPair,
                              typename C3t3::Triangulation::Cell_handle>;
  using ElementType = typename Base::ElementType;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_c3t3;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_visitor;
  using BaseClass::s_inc_cells;
  using BaseClass::get_incident_cells;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;


public:
  InternalEdgeFlipOperation(C3t3& c3t3,
                           CellSelector& cell_selector,
                           const bool protect_boundaries,
                           const Visitor& visitor)
    : BaseClass(c3t3, cell_selector, protect_boundaries, visitor)
  {}

  void perform_global_preprocessing(const C3t3& c3t3) const {
    // Reset cache validity for all cells (matches original flip_edges.h line 1925)
    for (auto c : c3t3.cells_in_complex())
      c->reset_cache_validity();
    
    // Clear the shared inc_cells to match original behavior where inc_cells starts empty
    s_inc_cells.clear();
  }

  virtual bool should_process_element(const ElementType& e, const C3t3& c3t3) const override {
    // For internal edges, we already filtered during collection, so always process
    return true;
  }

  std::vector<ElementType> get_element_source(const C3t3& c3t3) const override {
    // Perform global preprocessing (cache validity reset, inc_cells clear)
    perform_global_preprocessing(c3t3);
    
    // Collect internal vertex pairs fresh each time
    std::vector<ElementType> internal_vertex_pairs;
    get_internal_edges(c3t3, m_cell_selector, std::back_inserter(internal_vertex_pairs));
    
    // Return container by value
    return internal_vertex_pairs;
  }

  bool can_apply_operation(const ElementType& e, const C3t3& c3t3) const override {
    // Detailed validation will be done in the actual implementation
    return true;
  }

  Lock_zone get_lock_zone(const ElementType& e, const C3t3& c3t3) const override {
    Lock_zone zone;
    
    // For vertex pairs, we need to find the corresponding edge to get incident cells
    Cell_handle cell;
    int i0, i1;
    if (is_edge_uv(e.first, e.second, get_incident_cells(e.first, c3t3), cell, i0, i1)) {
      Edge edge(cell, i0, i1);
      auto circ = c3t3.triangulation().incident_cells(edge);
      auto done = circ;
      do {
        if (!c3t3.triangulation().is_infinite(circ)) {
          zone.push_back(circ);
        }
      } while (++circ != done);
    }
    
    return zone;
  }

  bool execute_operation(const ElementType& e, C3t3& c3t3) override {
    return execute_internal_edge_flip(e, c3t3);
  }

  std::string operation_name() const override {
    return "Edge Flip (Internal Edges)";
  }

  bool execute_pre_operation(const ElementType& e, C3t3& c3t3) override {
    
    return true;
  } 

  bool execute_post_operation(const ElementType& e, C3t3& c3t3) override {
    // Element-specific cleanup if needed
    return true;
  }

private:
  bool execute_internal_edge_flip(const ElementType& e, C3t3& c3t3) {
    // For internal edges, e is a vertex pair
    const auto& vp = e;
    
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    // DEBUG: Log edge processing order with edge info (first 5000 only)
    static thread_local std::size_t internal_edge_counter = 0;
    static thread_local bool show_debug = true;
    
    if (show_debug && internal_edge_counter < 5000) {
      std::cout << "REFACTORED: Processing internal edge #" << internal_edge_counter << " from vertex " << vp.first->point() << " to vertex " << vp.second->point() << std::endl;
    }
    
    if (internal_edge_counter == 4999) show_debug = false; // Stop after first 5000 edges (0-4999)
    internal_edge_counter++;
#endif
    
    auto& tr = c3t3.triangulation();
    
    // Use shared inc_cells cache (matches original flip_all_edges behavior)
    boost::container::small_vector<Cell_handle, 64>& o_inc_vh = s_inc_cells[vp.first];
    if (o_inc_vh.empty())
      tr.incident_cells(vp.first, std::back_inserter(o_inc_vh));
    
    Cell_handle ch;
    int i0, i1;
    if (is_edge_uv(vp.first, vp.second, o_inc_vh, ch, i0, i1))
    {
      Edge edge_to_flip(ch, i0, i1);
      
      Sliver_removal_result res = find_best_flip(edge_to_flip, c3t3, MIN_ANGLE_BASED, 
                                                 s_inc_cells, m_cell_selector, m_visitor);
      
      if (res == INVALID_CELL || res == INVALID_VERTEX || res == INVALID_ORIENTATION)
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
              if (show_debug && internal_edge_counter-1 < 5000) {
        std::cout << "REFACTORED: Edge flip problem for edge #" << internal_edge_counter-1 << " from vertex " << vp.first->point() << " to vertex " << vp.second->point() << std::endl;
      }
#endif
        return false; // Flip problem
      }
      
      if (res == VALID_FLIP) {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        if (show_debug && internal_edge_counter-1 < 5000) {
          std::cout << "REFACTORED: Successfully flipped edge #" << internal_edge_counter-1 << " from vertex " << vp.first->point() << " to vertex " << vp.second->point() << std::endl;
        }
#endif
      }
      
      return (res == VALID_FLIP);
    }
    
    return false; // Edge not found
  }
};

// Boundary Edge Flip Operation - processes vertex pairs like original flipBoundaryEdges
template<typename C3t3, typename CellSelector, typename Visitor>
class BoundaryEdgeFlipOperation 
  : public EdgeFlipOperationBase<C3t3, CellSelector, Visitor>,
    public ElementaryOperation<C3t3, 
                          std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using BaseClass = EdgeFlipOperationBase<C3t3, CellSelector, Visitor>;
  using VertexPair = std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>;
  
  using Base = ElementaryOperation<C3t3, 
                              VertexPair,
                              typename C3t3::Triangulation::Cell_handle>;
  using ElementType = typename Base::ElementType;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_c3t3;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_visitor;
  using BaseClass::s_inc_cells;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Facet;
  using typename BaseClass::Surface_patch_index;

private:
  // Boundary-specific preprocessing data - recreated every call like original
  static thread_local boost::unordered_map<Vertex_handle, 
                               boost::unordered_map<Surface_patch_index, unsigned int>> s_boundary_vertices_valences;
  static thread_local boost::unordered_map<Vertex_handle, std::unordered_set<typename C3t3::Subdomain_index>> s_vertices_subdomain_indices;
  
public:
  BoundaryEdgeFlipOperation(C3t3& c3t3,
                           CellSelector& cell_selector,
                           const bool protect_boundaries,
                           const Visitor& visitor)
    : BaseClass(c3t3, cell_selector, protect_boundaries, visitor)
  {}

  void perform_global_preprocessing(const C3t3& c3t3) const {
    // Clear boundary data structures every time (like original flip_edges)
    s_boundary_vertices_valences.clear();
    s_vertices_subdomain_indices.clear();
    
    // Collect boundary edges and compute vertices valences (needed for boundary flipping)
    std::vector<typename C3t3::Edge> boundary_edges; // We don't need to store this
    collectBoundaryEdgesAndComputeVerticesValences(c3t3, m_cell_selector, boundary_edges, 
                                                   s_boundary_vertices_valences, s_vertices_subdomain_indices);
    
  }

  virtual bool should_process_element(const ElementType& e, const C3t3& c3t3) const override {
    // For boundary edges, we already filtered during collection, so always process
    return true;
  }

  std::vector<ElementType> get_element_source(const C3t3& c3t3) const override {
    // Perform global preprocessing (valences computation, subdomain collection)
    perform_global_preprocessing(c3t3);
    
    // Collect boundary vertex pairs fresh each time
    std::vector<ElementType> boundary_vertex_pairs;
    
    // Collect boundary edges as vertex pairs to match original behavior
    for (auto eit = c3t3.triangulation().finite_edges_begin();
         eit != c3t3.triangulation().finite_edges_end(); ++eit) {
      const auto& e = *eit;
      if (is_boundary(c3t3, e, m_cell_selector) && !c3t3.is_in_complex(e)) {
        boundary_vertex_pairs.push_back(make_vertex_pair(e));
      }
    }
    
    // Return container by value
    return boundary_vertex_pairs;
  }

  bool can_apply_operation(const ElementType& e, const C3t3& c3t3) const override {
    return true;
  }

  Lock_zone get_lock_zone(const ElementType& e, const C3t3& c3t3) const override {
    Lock_zone zone;
    
    // For vertex pairs, we need to find the corresponding edge to get incident cells
    Cell_handle cell;
    int i0, i1;
    if (is_edge_uv(e.first, e.second, get_incident_cells(e.first, c3t3), cell, i0, i1)) {
      Edge edge(cell, i0, i1);
      auto circ = c3t3.triangulation().incident_cells(edge);
      auto done = circ;
      do {
        if (!c3t3.triangulation().is_infinite(circ)) {
          zone.push_back(circ);
        }
      } while (++circ != done);
    }
    
    return zone;
  }

  bool execute_operation(const ElementType& e, C3t3& c3t3) override {
    return execute_boundary_edge_flip(e, c3t3);
  }

  std::string operation_name() const override {
    return "Edge Flip (Boundary Edges)";
  }

  bool execute_pre_operation(const ElementType& e, C3t3& c3t3) override {
    
    return true;
  }

  bool execute_post_operation(const ElementType& e, C3t3& c3t3) override {
    // Element-specific cleanup if needed
    return true;
  }

private:
  bool execute_boundary_edge_flip(const ElementType& e, C3t3& c3t3) {
    // For boundary edges, e is a vertex pair
    const auto& vp = e;
    const auto& vh0 = vp.first;
    const auto& vh1 = vp.second;
    
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    // DEBUG: Log boundary edge processing order with edge info (first 5000 only)
    static thread_local std::size_t boundary_edge_counter = 0;
    static thread_local bool show_boundary_debug = true;
    
    if (show_boundary_debug && boundary_edge_counter < 5000) {
      std::cout << "REFACTORED: Processing boundary edge #" << boundary_edge_counter << " from vertex " << vh0->point() << " to vertex " << vh1->point() << std::endl;
    }
    
    if (boundary_edge_counter == 4999) show_boundary_debug = false; // Stop after first 5000 edges (0-4999)
    boundary_edge_counter++;
#endif
    
    auto& tr = c3t3.triangulation();
    
    // Use shared inc_cells cache (matches original behavior)
    boost::container::small_vector<Cell_handle, 64>& inc_vh0 = s_inc_cells[vh0];
    if (inc_vh0.empty())
      tr.incident_cells(vh0, std::back_inserter(inc_vh0));
    
    Cell_handle c;
    int i, j;
    if (!is_edge_uv(vh0, vh1, inc_vh0, c, i, j))
      return false;
    
    Edge edge_to_flip(c, i, j);
    std::vector<typename C3t3::Triangulation::Facet> boundary_facets;
    const bool on_boundary = is_boundary_edge(edge_to_flip, c3t3, m_cell_selector, boundary_facets);
    
    if (!on_boundary || boundary_facets.size() != 2)
      return false;
    
    const auto& f0 = boundary_facets[0];
    const auto& f1 = boundary_facets[1];
    
    // Find 3rd and 4th vertices to flip on surface
    const Vertex_handle vh2 = third_vertex(f0, vh0, vh1, tr);
    const Vertex_handle vh3 = third_vertex(f1, vh0, vh1, tr);
    
    if (!tr.tds().is_edge(vh2, vh3)) // Most-likely early exit
    {      
      // Implement valence-based cost calculation from lines 1810-1850
      const Surface_patch_index surfi = c3t3.surface_patch_index(boundary_facets[0]);
      
      // Get current valences
      int v0 = s_boundary_vertices_valences.at(vh0).at(surfi);
      int v1 = s_boundary_vertices_valences.at(vh1).at(surfi);
      int v2 = s_boundary_vertices_valences.at(vh2).at(surfi);
      int v3 = s_boundary_vertices_valences.at(vh3).at(surfi);
      
      if(v0 < 2 || v1 < 2 || v2 < 2 || v3 < 2)
        return false;
      
      // Calculate target valences
      int m0 = (s_boundary_vertices_valences.at(vh0).size() > 1 ? 4 : 6);
      int m1 = (s_boundary_vertices_valences.at(vh1).size() > 1 ? 4 : 6);
      int m2 = (s_boundary_vertices_valences.at(vh2).size() > 1 ? 4 : 6);
      int m3 = (s_boundary_vertices_valences.at(vh3).size() > 1 ? 4 : 6);
      
      // Calculate cost before flip
      int initial_cost = (v0 - m0)*(v0 - m0)
                       + (v1 - m1)*(v1 - m1)
                       + (v2 - m2)*(v2 - m2)
                       + (v3 - m3)*(v3 - m3);
      
      // Calculate cost after flip (vh0,vh1 decrease by 1, vh2,vh3 increase by 1)
      v0--; v1--; v2++; v3++;
      int final_cost = (v0 - m0)*(v0 - m0)
                     + (v1 - m1)*(v1 - m1)
                     + (v2 - m2)*(v2 - m2)
                     + (v3 - m3)*(v3 - m3);
      
      if (initial_cost > final_cost) {
        // Perform the flip (using shared inc_cells for consistency with original)
        Sliver_removal_result db = flip_on_surface(c3t3, edge_to_flip, vh2, vh3,
                                                   s_inc_cells, MIN_ANGLE_BASED, m_visitor);
        
        if (db == VALID_FLIP) {
          // Add new facets to complex
          Cell_handle c_new;
          int li, lj, lk;
          if (tr.tds().is_facet(vh2, vh3, vh0, c_new, li, lj, lk)) {
            c3t3.add_to_complex(c_new, (6 - li - lj - lk), surfi);
          }
          if (tr.tds().is_facet(vh2, vh3, vh1, c_new, li, lj, lk)) {
            c3t3.add_to_complex(c_new, (6 - li - lj - lk), surfi);
          }
          
          // Update valences
          s_boundary_vertices_valences[vh0][surfi]--;
          s_boundary_vertices_valences[vh1][surfi]--;
          s_boundary_vertices_valences[vh2][surfi]++;
          s_boundary_vertices_valences[vh3][surfi]++;
          
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
          if (show_boundary_debug && boundary_edge_counter-1 < 5000) {
            std::cout << "REFACTORED: Successfully flipped boundary edge #" << boundary_edge_counter-1 << " from vertex " << vh0->point() << " to vertex " << vh1->point() << std::endl;
          }
#endif
          return true;
        }
      }
    }
    
    return false;
  }
};

// Static member definitions for EdgeFlipOperationBase (shared)
template<typename C3t3, typename CellSelector, typename Visitor>
thread_local std::unordered_map<typename C3t3::Vertex_handle, boost::container::small_vector<typename C3t3::Cell_handle, 64>> 
EdgeFlipOperationBase<C3t3, CellSelector, Visitor>::s_inc_cells;

// Static member definitions for InternalEdgeFlipOperation
// (No longer needed since we removed std::once_flag)

// Static member definitions for BoundaryEdgeFlipOperation  
// (No longer needed since we removed std::once_flag)

template<typename C3t3, typename CellSelector, typename Visitor>
thread_local boost::unordered_map<typename C3t3::Vertex_handle, 
                       boost::unordered_map<typename C3t3::Surface_patch_index, unsigned int>> 
BoundaryEdgeFlipOperation<C3t3, CellSelector, Visitor>::s_boundary_vertices_valences;

template<typename C3t3, typename CellSelector, typename Visitor>
thread_local boost::unordered_map<typename C3t3::Vertex_handle, std::unordered_set<typename C3t3::Subdomain_index>> 
BoundaryEdgeFlipOperation<C3t3, CellSelector, Visitor>::s_vertices_subdomain_indices;

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H 