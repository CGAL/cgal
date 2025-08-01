#ifndef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/split_long_edges.h>
#include <boost/bimap.hpp>

#include <optional>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <thread>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

//#define CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
// Thread-safe debug logger for edge split operations
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
static ThreadSafeLogger edge_split_logger("debug_logs/edge_split_debug.log");
static ThreadSafeLogger long_edges_logger("debug_logs/long_edges_collection.log");
#endif

template<typename C3t3, typename SizingFunction, typename CellSelector>
class EdgeSplitOperation 
  : public ElementaryOperation<C3t3,
                              std::pair<typename C3t3::Triangulation::Geom_traits::FT,std::pair<typename C3t3::Triangulation::Vertex_handle, typename C3t3::Triangulation::Vertex_handle>>,
                              std::vector<std::pair<typename C3t3::Triangulation::Geom_traits::FT,std::pair<typename C3t3::Triangulation::Vertex_handle, typename C3t3::Triangulation::Vertex_handle>>>,
                              typename C3t3::Triangulation::Cell_handle>
{

public:
  using Complex = C3t3;
  using Tr = typename C3t3::Triangulation;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using FT = typename Tr::Geom_traits::FT;
  using Base = ElementaryOperation<C3t3,
                              std::pair<FT ,std::pair<typename C3t3::Triangulation::Vertex_handle, typename C3t3::Triangulation::Vertex_handle>>,
                              std::vector<std::pair<double,std::pair<typename C3t3::Triangulation::Vertex_handle, typename C3t3::Triangulation::Vertex_handle>>>,
                              typename C3t3::Triangulation::Cell_handle>;
  using Edge = typename Tr::Edge;
  using ElementType = typename Base::ElementType;
  using ElementSource = typename Base::ElementSource;
  using VertexPair = std::pair<typename C3t3::Triangulation::Vertex_handle, typename C3t3::Triangulation::Vertex_handle>;
  using Cell_handle = typename Tr::Cell_handle;
  using Point = typename Tr::Point;
  using Lock_zone = typename Base::Lock_zone;
  using Facet = typename Tr::Facet;
  using Subdomain_index = typename C3t3::Subdomain_index;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Curve_index = typename C3t3::Curve_index;

private:
  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  static int s_operation_counter; // Static counter for operation tracking
  
  boost::container::small_vector<Cell_handle, 64> get_incident_cells(Vertex_handle vh, const C3t3& c3t3) const {
    boost::container::small_vector<Cell_handle, 64> inc_cells;
#ifdef USE_THREADSAFE_INCIDENT_CELLS
    c3t3.triangulation().incident_cells_threadsafe(vh, std::back_inserter(inc_cells));
#else
    c3t3.triangulation().incident_cells(vh, std::back_inserter(inc_cells));
#endif
    return inc_cells;
  }

  // Helper function to create a unique identifier for an edge
  std::string edge_id(const ElementType& vp) const {
    std::ostringstream oss;
    oss << "(" << vp.first->point() << " -> " << vp.second->point() << ")";
    return oss.str();
  }

public:
  EdgeSplitOperation(const SizingFunction& sizing,
                    const CellSelector& cell_selector,
                    const bool protect_boundaries)
    : m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
  {}


  std::vector<std::pair<FT,std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>>
  get_long_edges(const C3t3& c3t3) const {
    std::vector<std::pair<FT, std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>>
        long_edges_with_lengths;
    const Tr& tr = c3t3.triangulation();

    for(Edge e : tr.finite_edges()) {
      auto [splittable, boundary] = can_be_split(e, c3t3, m_protect_boundaries, m_cell_selector);
      if(!splittable)
        continue;

      const std::optional<FT> sqlen = is_too_long(e, boundary, m_sizing, c3t3, m_cell_selector);
      if(sqlen != std::nullopt) {
        auto edge_pair = make_vertex_pair(e);
        long_edges_with_lengths.push_back(make_pair(*sqlen, edge_pair));
      }

      for(const auto& length_edge_pair : long_edges_with_lengths) {
        auto edge_pair = length_edge_pair.second;
      }
    }

    // Custom comparator to match bimap's tie-breaking behavior
    // Bimap: multiset_of<FT, std::greater<FT>> preserves insertion order for equal keys
    // We use stable_sort with simple length comparison to achieve the same behavior
    auto bimap_comparator = [](const std::pair<double, std::pair<Vertex_handle, Vertex_handle>>& a,
                               const std::pair<double, std::pair<Vertex_handle, Vertex_handle>>& b) {
      // Only compare by length - stable_sort will preserve insertion order for equal lengths
      return a.first > b.first; // std::greater<FT> behavior (descending order)
    };

    std::stable_sort(long_edges_with_lengths.begin(), long_edges_with_lengths.end(), bimap_comparator);
    return long_edges_with_lengths;

    //// Extract just the edges (without lengths) into long_edges
    //std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>> long_edges;
    //long_edges.reserve(long_edges_with_lengths.size());
    //for(const auto& length_edge_pair : long_edges_with_lengths) {
    //  long_edges.push_back(length_edge_pair.second);
    //}

    //return long_edges;
  }
  ElementSource get_element_source(const C3t3& c3t3) const override {
    // Collect long edges as vertex pairs
    auto long_edges = get_long_edges(c3t3);
    
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
    // Log the collection of long edges to separate file
    long_edges_logger.log("=== LONG EDGES COLLECTION - Thread: ", std::this_thread::get_id(), " ===");
    long_edges_logger.log("Total long edges found: ", long_edges.size());
    
    for (size_t i = 0; i < long_edges.size(); ++i) {
      const auto& edge_data = long_edges[i];
      const auto& length = edge_data.first;
      
      long_edges_logger.log("Edge[", i, "] length:", length);
    }
    long_edges_logger.log("=== END LONG EDGES COLLECTION ===");
#endif
    
    return long_edges;
  }


  bool lock_zone(const ElementType& el, const Complex& c3t3) const override {
    auto& tr = c3t3.triangulation();
    const auto& vertex_pair =el.second;
    const auto& edge_length =el.first;
      #if 1
    std::vector<Cell_handle> inc_cells_first,inc_cells_second;
    
    bool lock_success = tr.try_lock_and_get_incident_cells(vertex_pair.first, inc_cells_first) && 
                       tr.try_lock_and_get_incident_cells(vertex_pair.second, inc_cells_second);
    
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
    if(!lock_success) {
        edge_split_logger.log("[LOCK_FAILED] length:", edge_length, " - Thread: ", std::this_thread::get_id());
    }
#endif
    
    return lock_success;
    
    #else
    
    //auto& tr = c3t3.triangulation();
    // We need to lock v individually first, to be sure v->cell() is valid
    if(!tr.try_lock_vertex(vertex_pair.first) || !tr.try_lock_vertex(vertex_pair.second))
      return false;

    Cell_handle ch;
    int i0, i1;
    boost::container::small_vector<Cell_handle, 64> inc_vh = get_incident_cells(vertex_pair.first, c3t3);
    is_edge_uv(vertex_pair.first, vertex_pair.second, inc_vh, ch, i0, i1);
      Edge edge_to_flip(ch, i0, i1);

    C3t3::Triangulation::Cell_circulator ccirc(edge_to_flip);
    C3t3::Triangulation::Cell_circulator cdone = ccirc;

    do {
      if(!tr.try_lock_cell(ccirc)) // LOCK
        return false;
      for(int i=0; i < 4; ++i) {
        if(!tr.try_lock_cell(ccirc->neighbor(i))) // LOCK
          return false;
      }
      ++ccirc;
    } while(ccirc != cdone);
    
    return true;
    #endif
  }

  bool execute_operation(const ElementType& vp, Complex& c3t3) override {
    auto& tr = c3t3.triangulation();

    Edge_vv e = vp.second;

    Cell_handle cell;
    int i1, i2;
    if(tr.tds().is_edge(e.first, e.second, cell, i1, i2)) {
      Edge edge(cell, i1, i2);

      // check that splittability has not changed
      auto [splittable, boundary] = can_be_split(edge, c3t3, m_protect_boundaries, m_cell_selector);
      if(!splittable) {
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
        // Debug log for split failure
        edge_split_logger.log("[SPLIT_FAILED] edge not splittable. length:", vp.first, " - Thread: ", std::this_thread::get_id());
#endif
        return false;
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      else
        can_be_split_ofs << "2 " << edge.first->vertex(edge.second)->point() << " "
                         << edge.first->vertex(edge.third)->point() << std::endl;
#endif
      // m_visitor.before_split(tr, edge);
      Vertex_handle vh = split_edge(edge, m_cell_selector, c3t3);

      // Debug output for split result
      if (vh == Vertex_handle()) {
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
        // Debug log for split failure
        edge_split_logger.log("[SPLIT_FAILED] operation failed. length:", vp.first, " - Thread: ", std::this_thread::get_id());
#endif
        return false;
      }
      
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
      // Debug log for successful split
      edge_split_logger.log("[SPLIT_SUCCESS] on edge with length:", vp.first, " - Thread: ", std::this_thread::get_id());
#endif
      
        //if(vh != Vertex_handle())
        //  m_visitor.after_split(tr, vh);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      else
        split_failed_ofs << "2 " << edge.first->vertex(edge.second)->point() << " "
                         << edge.first->vertex(edge.third)->point() << std::endl;
      if (vh != Vertex_handle())
        ofs << vh->point() << std::endl;
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
      std::cout << "\rSplit... (" << long_edges.left.size() << " long edges, "
                << "length  = " << std::sqrt(sqlen) << ", " << nb_splits << " splits)";
      std::cout.flush();
#endif
      return true;
   }
   
#ifdef CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_DEBUG
   // Debug log for edge not found
   edge_split_logger.log("[SPLIT_FAILED] edge not found. length:", vp.first, " - Thread: ", std::this_thread::get_id());
#endif
    return false;
  }

  bool requires_ordered_processing() const override {
    return true; // EdgeSplit requires ordered processing for optimal performance
  }

  std::string operation_name() const override {
    return "Edge Split";
  }

};

// Initialize static counter
template<typename C3t3, typename SizingFunction, typename CellSelector>
int EdgeSplitOperation<C3t3, SizingFunction, CellSelector>::s_operation_counter = 0;

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_EDGE_SPLIT_OPERATION_H 