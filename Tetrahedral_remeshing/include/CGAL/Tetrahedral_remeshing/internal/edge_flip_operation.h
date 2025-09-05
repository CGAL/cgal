#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H

#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Tetrahedral_remeshing/internal/flip_edges.h>
#include <CGAL/tss.h>

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
  typedef typename boost::container::small_vector<Cell_handle, 64> Cells_vector;

  C3t3& m_c3t3;
  CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  Visitor m_visitor;
  mutable tbb::concurrent_unordered_map<Vertex_handle, Cells_vector> inc_cells;

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

};

// Internal Edge Flip Operation - processes vertex pairs like original get_internal_edges
template<typename C3t3, typename CellSelector, typename Visitor>
class InternalEdgeFlipOperation
  : public EdgeFlipOperationBase<C3t3, CellSelector, Visitor>,
    public ElementaryOperation<C3t3,
                          std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                          std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using BaseClass = EdgeFlipOperationBase<C3t3, CellSelector, Visitor>;
  using VertexPair = std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>;

  using Base = ElementaryOperation<C3t3,
                              std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                              std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,
                              typename C3t3::Triangulation::Cell_handle>;
  using ElementType = typename Base::ElementType;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_c3t3;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_visitor;
  using BaseClass::inc_cells;

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
  }


  ElementSource get_element_source(const C3t3& c3t3) const override {
    // Perform global preprocessing (cache validity reset)
    perform_global_preprocessing(c3t3);

    // Collect internal vertex pairs fresh each time
    std::vector<ElementType> internal_vertex_pairs;
    get_internal_edges(c3t3, m_cell_selector, std::back_inserter(internal_vertex_pairs));
    // Return container by value
    return internal_vertex_pairs;
  }

  bool lock_zone(const ElementType& e, const C3t3& c3t3) const override {

      auto& tr = c3t3.triangulation();
      #if 1
    boost::container::small_vector<Cell_handle, 64> inc_cells_first,inc_cells_second;
    bool successfully_locked=tr.try_lock_and_get_incident_cells(e.first, inc_cells_first) &&tr.try_lock_and_get_incident_cells(e.second, inc_cells_second);
    //Cache the incident cells
    inc_cells[e.first] = inc_cells_first;
    inc_cells[e.second] = inc_cells_second;
return successfully_locked;
#else

    const auto& vp = e;
    // We need to lock v individually first, to be sure v->cell() is valid
    if(!tr.try_lock_vertex(vp.first) || !tr.try_lock_vertex(vp.second))
      return false;

    Cell_handle ch;
    int i0, i1;
    boost::container::small_vector<Cell_handle, 64> inc_vh = get_incident_cells(vp.first, c3t3);
    is_edge_uv(vp.first, vp.second, inc_vh, ch, i0, i1);
      Edge edge_to_flip(ch, i0, i1);

    C3t3::Triangulation::Cell_circulator ccirc(edge_to_flip);
    C3t3::Triangulation::Cell_circulator cdone = ccirc;

    do {
      if(!tr.try_lock_cell(ccirc)) // LOCK
        return false;
      ++ccirc;
    } while(ccirc != cdone);

    return true;
    #endif
  }

  bool execute_operation(const ElementType& e, C3t3& c3t3) override {
    return execute_internal_edge_flip(e, c3t3);
  }

  bool requires_ordered_processing() const override {
    return false; // InternalEdgeFlip can use unordered parallel processing
  }

  std::string operation_name() const override {
    return "Edge Flip (Internal Edges)";
  }

private:
  bool execute_internal_edge_flip(const ElementType& e, C3t3& c3t3) {
    const auto& vp = e;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    static thread_local std::size_t internal_edge_counter = 0;
    static thread_local bool show_debug = true;

    if (show_debug && internal_edge_counter < 5000) {
      std::cout << "REFACTORED: Processing internal edge #" << internal_edge_counter << " from vertex " << vp.first->point() << " to vertex " << vp.second->point() << std::endl;
    }

    if (internal_edge_counter == 4999) show_debug = false;
    internal_edge_counter++;
#endif

    auto& tr = c3t3.triangulation();

    // Get incident cells for first vertex only
    //auto inc_vh = get_incident_cells(vp.first, c3t3);
    boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vp.first];
    if (o_inc_vh.empty())
      o_inc_vh=get_incident_cells(vp.first, c3t3);

    Cell_handle ch;
    int i0, i1;
    if (is_edge_uv(vp.first, vp.second, o_inc_vh, ch, i0, i1))
    {
      Edge edge_to_flip(ch, i0, i1);

      // Create a map starting with the cells we already computed

      Sliver_removal_result res = find_best_flip(edge_to_flip, c3t3, MIN_ANGLE_BASED,
                                               inc_cells, m_cell_selector, m_visitor);

      if (res == INVALID_CELL || res == INVALID_VERTEX || res == INVALID_ORIENTATION)
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        if (show_debug && internal_edge_counter-1 < 5000) {
          std::cout << "REFACTORED: Edge flip problem for edge #" << internal_edge_counter-1 << " from vertex " << vp.first->point() << " to vertex " << vp.second->point() << std::endl;
        }
#endif
        return false;
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

    return false;
  }
};

// Boundary Edge Flip Operation - processes vertex pairs like original flipBoundaryEdges
template<typename C3t3, typename CellSelector, typename Visitor>
class BoundaryEdgeFlipOperation
  : public EdgeFlipOperationBase<C3t3, CellSelector, Visitor>,
    public ElementaryOperation<C3t3,
                          std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                          std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,
                          typename C3t3::Triangulation::Cell_handle>
{
public:
  using BaseClass = EdgeFlipOperationBase<C3t3, CellSelector, Visitor>;
  using VertexPair = std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>;

  using Base = ElementaryOperation<C3t3,
                              std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                              std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>,
                              typename C3t3::Triangulation::Cell_handle>;
  using ElementType = typename Base::ElementType;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_c3t3;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_visitor;
  using BaseClass::inc_cells;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Facet;
  using typename BaseClass::Surface_patch_index;

private:

    //NOTE: Is the usage of s_boundary_vertices_valences static variable in a multithreaded enviroment safe?
  //// Boundary-specific preprocessing data - recreated every call like original
  using BVV =boost::unordered_map<Vertex_handle,
                                   boost::unordered_map<Surface_patch_index, unsigned int>>;
  static BVV s_boundary_vertices_valences;
  static BVV& get_static_boundary_vertices_valences() {
    //CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(BVV, s_boundary_vertices_valences);
    return s_boundary_vertices_valences;
  }


public:
  BoundaryEdgeFlipOperation(C3t3& c3t3,
                           CellSelector& cell_selector,
                           const bool protect_boundaries,
                           const Visitor& visitor)
    : BaseClass(c3t3, cell_selector, protect_boundaries, visitor)
  {}

  //TODO: collectBoundaryEdgesAndComputeVerticesValences is quite wasteful and should be refactored
  void perform_global_preprocessing(const C3t3& c3t3) const {
    // Collect boundary edges and compute vertices valences (needed for boundary flipping)
    std::vector<typename C3t3::Edge> boundary_edges; // We don't need to store this
  boost::unordered_map<Vertex_handle,
                                   std::unordered_set<typename C3t3::Subdomain_index>> vertices_subdomain_indices;
    collectBoundaryEdgesAndComputeVerticesValences(c3t3, m_cell_selector, boundary_edges,
                                                   get_static_boundary_vertices_valences(),
                                                   vertices_subdomain_indices);

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

  bool lock_zone(const ElementType& e, const C3t3& c3t3) const override {
    auto& tr = c3t3.triangulation();
  #if 1
    boost::container::small_vector<Cell_handle, 64> inc_cells_first,inc_cells_second;
    bool successfully_locked = tr.try_lock_and_get_incident_cells(e.first, inc_cells_first) &&tr.try_lock_and_get_incident_cells(e.second, inc_cells_second);
    inc_cells[e.first] = inc_cells_first;
    inc_cells[e.second] = inc_cells_second;
    return successfully_locked;
    //      const auto& vp = e;
    //// We need to lock v individually first, to be sure v->cell() is valid
    //if(!tr.try_lock_vertex(vp.first) || !tr.try_lock_vertex(vp.second))
    //  return false;

    //Cell_handle ch;
    //int i0, i1;
    //boost::container::small_vector<Cell_handle, 64> inc_vh = get_incident_cells(vp.first, c3t3);
    //is_edge_uv(vp.first, vp.second, inc_vh, ch, i0, i1);
    //Edge edge_to_flip(ch, i0, i1);

    //C3t3::Triangulation::Cell_circulator ccirc(edge_to_flip);
    //C3t3::Triangulation::Cell_circulator cdone = ccirc;

    //do {
    //  if(!tr.try_lock_cell(ccirc)) // LOCK
    //    return false;
    //  ++ccirc;
    //} while(ccirc != cdone);

    //return true;
    #else
      // QUESTION: why does locking the incident facets cause a data race?
      const auto& vp = e;
    auto& tr = c3t3.triangulation();
    // We need to lock v individually first, to be sure v->cell() is valid
    if(!tr.try_lock_vertex(vp.first) || !tr.try_lock_vertex(vp.second))
      return false;

    Cell_handle ch;
    int i0, i1;
    boost::container::small_vector<Cell_handle, 64> inc_vh = get_incident_cells(vp.first, c3t3);
    is_edge_uv(vp.first, vp.second, inc_vh, ch, i0, i1);
      Edge edge_to_flip(ch, i0, i1);

    C3t3::Triangulation::Facet_circulator ccirc(edge_to_flip);
  C3t3::Triangulation::Facet_circulator cdone = ccirc;

    do
    {
      if(!tr.try_lock_facet(*ccirc))
        ++ccirc;
    } while(ccirc!=cdone);

    return true;

       #endif

  }

  bool execute_operation(const ElementType& e, C3t3& c3t3) override {
    return execute_boundary_edge_flip(e, c3t3);
  }

  bool requires_ordered_processing() const override {
    return false; // BoundaryEdgeFlip can use unordered parallel processing
  }

  std::string operation_name() const override {
    return "Edge Flip (Boundary Edges)";
  }

private:
  bool execute_boundary_edge_flip(const ElementType& e, C3t3& c3t3) {
    assert(get_static_boundary_vertices_valences().size() > 0 &&
           "Boundary vertices valences must be initialized before flipping boundary edges.");
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
    //boost::container::small_vector<Cell_handle, 64> inc_vh0 = get_incident_cells(vh0, c3t3);
    boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vp.first];
    if (o_inc_vh.empty())
      o_inc_vh=get_incident_cells(vp.first, c3t3);

    Cell_handle c;
    int i, j;
    if (!is_edge_uv(vh0, vh1, o_inc_vh, c, i, j))
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
      int v0 = get_static_boundary_vertices_valences().at(vh0).at(surfi);
      int v1 = get_static_boundary_vertices_valences().at(vh1).at(surfi);
      int v2 = get_static_boundary_vertices_valences().at(vh2).at(surfi);
      int v3 = get_static_boundary_vertices_valences().at(vh3).at(surfi);

      if(v0 < 2 || v1 < 2 || v2 < 2 || v3 < 2)
        return false;

      // Calculate target valences
      int m0 = (get_static_boundary_vertices_valences().at(vh0).size() > 1 ? 4 : 6);
      int m1 = (get_static_boundary_vertices_valences().at(vh1).size() > 1 ? 4 : 6);
      int m2 = (get_static_boundary_vertices_valences().at(vh2).size() > 1 ? 4 : 6);
      int m3 = (get_static_boundary_vertices_valences().at(vh3).size() > 1 ? 4 : 6);

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
        // Perform the flip using temporary inc_cells map
        Sliver_removal_result db = flip_on_surface(c3t3, edge_to_flip, vh2, vh3,
                                                 inc_cells, MIN_ANGLE_BASED, m_visitor);

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
          //TODO: In a parallel context we should be sure that this is (thread) safe.
          // Whatever tests we have done there were no crashes.
          // tbb offers a concurrent hash map
          get_static_boundary_vertices_valences()[vh0][surfi]--;
          get_static_boundary_vertices_valences()[vh1][surfi]--;
          get_static_boundary_vertices_valences()[vh2][surfi]++;
          get_static_boundary_vertices_valences()[vh3][surfi]++;

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
  template<typename C3t3, typename CellSelector, typename Visitor>
  typename BoundaryEdgeFlipOperation<C3t3, CellSelector, Visitor>::BVV
  BoundaryEdgeFlipOperation<C3t3, CellSelector, Visitor>::s_boundary_vertices_valences;

  } // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_FLIP_OPERATION_H