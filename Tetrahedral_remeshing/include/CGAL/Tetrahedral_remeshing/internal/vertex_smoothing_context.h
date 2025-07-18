#ifndef CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTHING_CONTEXT_H
#define CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTHING_CONTEXT_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <vector>
#include <list>
#include <optional>

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

  // Surface-related data (used by ComplexEdge and SurfaceVertex operations)
  std::unordered_map<Vertex_handle, std::vector<Surface_patch_index>> m_vertices_surface_indices;
  std::unordered_map<Vertex_handle, std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> m_vertices_normals;

private:
  const CellSelector& m_cell_selector;

public:
  VertexSmoothingContext(C3t3& c3t3,
                         const SizingFunction& sizing,
                         const CellSelector& cell_selector,
                         const bool protect_boundaries)
    : m_cell_selector(cell_selector)
  {
    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3, protect_boundaries, cell_selector);
    collect_incident_cells(c3t3.triangulation(), cell_selector);
    
    #ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    // Initialize surface data if boundaries are not protected
    if (!protect_boundaries) {
      collect_vertices_surface_indices(c3t3);
      compute_vertices_normals(c3t3);
    }
    #endif
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

  void collect_vertices_surface_indices(const C3t3& c3t3) {
    for (const Facet& fit : c3t3.facets_in_complex()) {
      const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);

      for (const Vertex_handle vi : c3t3.triangulation().vertices(fit)) {
        std::vector<Surface_patch_index>& v_surface_indices = m_vertices_surface_indices[vi];
        if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index) == v_surface_indices.end())
          v_surface_indices.push_back(surface_index);
      }
    }
  }

  void compute_vertices_normals(const C3t3& c3t3) {
    typename Tr::Geom_traits gt = c3t3.triangulation().geom_traits();
    typename Tr::Geom_traits::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();

    const Tr& tr = c3t3.triangulation();

    //collect all facet normals
    std::unordered_map<Facet, Vector_3, boost::hash<Facet>> fnormals;
    for (const Facet& f : tr.finite_facets())
    {
      if (is_boundary(c3t3, f, m_cell_selector))
      {
        const Facet cf = canonical_facet(f);
        fnormals[cf] = CGAL::NULL_VECTOR;
      }
    }

    for (const auto& fn : fnormals)
    {
      if(fn.second != CGAL::NULL_VECTOR)
        continue;

      const Facet& f = fn.first;
      const Facet& mf = tr.mirror_facet(f);
      CGAL_expensive_assertion(is_boundary(c3t3, f, m_cell_selector));

      Vector_3 start_ref = CGAL::Tetrahedral_remeshing::normal(f, tr.geom_traits());
      if (c3t3.triangulation().is_infinite(mf.first)
          || c3t3.subdomain_index(mf.first) < c3t3.subdomain_index(f.first))
        start_ref = opp(start_ref);
      fnormals[f] = start_ref;

      std::list<Facet> facets;
      facets.push_back(f);
      while (!facets.empty())
      {
        const Facet ff = facets.front();
        facets.pop_front();

        const Vector_3& ref = fnormals[ff];
        for (const Edge& ei : facet_edges(ff.first, ff.second, tr))
        {
          if (std::optional<Facet> neighbor
              = find_adjacent_facet_on_surface(ff, ei, c3t3))
          {
            const Facet neigh = *neighbor; //already a canonical_facet
            if (fnormals[neigh] == CGAL::NULL_VECTOR) //check it's not already computed
            {
              fnormals[neigh] = compute_normal(neigh, ref, gt);
              facets.push_back(neigh);
            }
          }
        }
      }
    }

    for (const auto& [f, n] : fnormals)
    {
      const Surface_patch_index& surf_i = c3t3.surface_patch_index(f);

      for (const Vertex_handle vi : tr.vertices(f))
      {
        auto patch_vector_it = m_vertices_normals.find(vi);

        if (patch_vector_it == m_vertices_normals.end()
            || patch_vector_it->second.find(surf_i) == patch_vector_it->second.end())
        {
          m_vertices_normals[vi][surf_i] = n;
        }
        else
        {
          m_vertices_normals[vi][surf_i] += n;
        }
      }
    }

    //normalize the computed normals
    for (auto vnm_it = m_vertices_normals.begin(); vnm_it != m_vertices_normals.end(); ++vnm_it)
    {
      //value type is map<Surface_patch_index, Vector_3>
      for (auto it = vnm_it->second.begin(); it != vnm_it->second.end(); ++it)
      {
        Vector_3& n = it->second;
        CGAL::Tetrahedral_remeshing::normalize(n, gt);
      }
    }
  }

  std::optional<Facet>
  find_adjacent_facet_on_surface(const Facet& f,
                                 const Edge& edge,
                                 const C3t3& c3t3)
  {
    CGAL_expensive_assertion(is_boundary(c3t3, f, m_cell_selector));

    typedef typename Tr::Facet_circulator Facet_circulator;

    if (c3t3.is_in_complex(edge))
      return {}; //do not "cross" complex edges
    //they are likely to be sharp and not to follow the > 0 dot product criterion

    const Surface_patch_index& patch = c3t3.surface_patch_index(f);
    const Facet& mf = c3t3.triangulation().mirror_facet(f);

    Facet_circulator fcirc = c3t3.triangulation().incident_facets(edge);
    Facet_circulator fend = fcirc;
    do
    {
      const Facet fi = *fcirc;
      if (f != fi
          && mf != fi
          && is_boundary(c3t3, fi, m_cell_selector)
          && patch == c3t3.surface_patch_index(fi))
      {
        return canonical_facet(fi); //"canonical" is important
      }
    } while (++fcirc != fend);

    return {};
  }

  template<typename Gt>
  Vector_3 compute_normal(const Facet& f,
                          const Vector_3& reference_normal,
                          const Gt& gt)
  {
    typename Gt::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();
    typename Gt::Compute_scalar_product_3
      scalar_product = gt.compute_scalar_product_3_object();

    Vector_3 n = CGAL::Tetrahedral_remeshing::normal(f, gt);
    if (scalar_product(n, reference_normal) < 0.)
      n = opp(n);

    return n;
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTHING_CONTEXT_H 