#ifndef CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTHING_CONTEXT_H
#define CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTHING_CONTEXT_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <boost/container/small_vector.hpp>

#include <unordered_map>
#include <vector>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template<typename C3t3, typename SizingFunction, typename CellSelector>
class VertexSmoothingContext {
public:
  using Tr = typename C3t3::Triangulation;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Cell_handle = typename Tr::Cell_handle;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Edge = typename Tr::Edge;
  using Facet = typename Tr::Facet;

  using Gt = typename Tr::Geom_traits;
  using Vector_3 = typename Gt::Vector_3;
  using Point_3 = typename Gt::Point_3;
  using FT = typename Gt::FT;

  struct Move {
    Vector_3 move;
    int neighbors;
    FT mass;
  };

  // Shared pre-computed data (used by at least 2 operations)
  std::unordered_map<Vertex_handle, std::size_t> m_vertex_id;
  std::vector<bool> m_free_vertices;
  std::vector<Move> m_moves;

  // Incident cells data (used by all 3 operations)
  using Incident_cells_vector = boost::container::small_vector<Cell_handle, 64>;
  std::vector<Incident_cells_vector> m_inc_cells;

public:
  VertexSmoothingContext(C3t3& c3t3,
                         const SizingFunction& sizing,
                         const CellSelector& cell_selector,
                         const bool protect_boundaries)
  {
    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3, protect_boundaries, cell_selector);
    collect_incident_cells(c3t3.triangulation(), cell_selector);
  }

private:
  void reset_vertex_id_map(const Tr& tr) {
    m_vertex_id.clear();
    std::size_t id = 0;
    for (const Vertex_handle v : tr.finite_vertex_handles()) {
      m_vertex_id[v] = id++;
    }
  }

  void reset_free_vertices(const C3t3& c3t3, bool protect_boundaries, const CellSelector& cell_selector) {
      const Tr& tr = c3t3.triangulation();
      m_free_vertices.clear();
      m_free_vertices.resize(tr.number_of_vertices(), false);

      for (const Cell_handle c : tr.finite_cell_handles()) {
          if (!get(cell_selector, c))
              continue;

          for (auto vi : tr.vertices(c)) {
              const std::size_t idi = m_vertex_id.at(vi);
              const int dim = vi->in_dimension();

              switch (dim) {
              case 3:
                  m_free_vertices[idi] = true;
                  break;
              case 2:
                  m_free_vertices[idi] = !protect_boundaries;
                  break;
              case 1:
                  m_free_vertices[idi] = !protect_boundaries; // and smooth_constrained_edges
                  break;
              case 0:
                  m_free_vertices[idi] = false;
                  break;
              default:
                  CGAL_unreachable();
              }
          }
      }
  }

  void collect_incident_cells(const Tr& tr, const CellSelector& cell_selector) {
    const std::size_t nbv = tr.number_of_vertices();
    m_inc_cells.resize(nbv, Incident_cells_vector{});
    
    for (const Cell_handle c : tr.finite_cell_handles()) {
      if (!get(cell_selector, c))
        continue;

      for (auto vi : tr.vertices(c)) {
        const std::size_t idi = m_vertex_id.at(vi);
        if (m_free_vertices[idi])
          m_inc_cells[idi].push_back(c);
      }
    }
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTHING_CONTEXT_H 