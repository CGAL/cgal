// Copyright (c) 2025 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Iasonas Manolas, Jane Tournois

#ifndef CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H
#define CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Tetrahedral_remeshing/internal/smooth_vertices.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>


#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_segment_primitive_3.h>
#include <CGAL/use.h>

#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>
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
  bool m_flip_smooth_steps = false; //TODO: when the smooth-steps start this needs to be true
  std::vector<Move> m_moves;

  // Incident cells data (used by all 3 operations)
  using Incident_cells_vector = boost::container::small_vector<Cell_handle, 64>;
  std::vector<Incident_cells_vector> m_inc_cells;

    // Segment AABB tree components (used by both SurfaceVertexSmoothOperation and ComplexEdgeVertexSmoothOperation)
  using Segment_vec = std::vector<typename Gt::Segment_3>;
  using Segment_iter = typename Segment_vec::iterator;
  using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
  using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
  using AABB_segment_tree = CGAL::AABB_tree<AABB_segment_traits>;
  // Triangle AABB tree (only used by SurfaceVertexSmoothOperation)
  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  Triangle_vec m_aabb_triangles;
  AABB_triangle_tree m_triangles_aabb_tree;
  FT m_aabb_epsilon;

  Segment_vec m_aabb_segments;
  AABB_segment_tree m_segments_aabb_tree;

  std::unordered_map<Vertex_handle, std::vector<Surface_patch_index>> m_vertices_surface_indices;
  std::unordered_map<Vertex_handle, std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> m_vertices_normals;

  FT m_total_move = 0;
private:
  const CellSelector& m_cell_selector;
  const bool m_protect_boundaries;
  const bool m_smooth_constrained_edges;

public:
  using FMLS = CGAL::Tetrahedral_remeshing::internal::FMLS<Gt>;
  std::vector<FMLS> subdomain_FMLS;
  std::unordered_map<Surface_patch_index, std::size_t, boost::hash<Surface_patch_index>> subdomain_FMLS_indices;
  VertexSmoothingContext(C3t3& c3t3,
                         const CellSelector& cell_selector,
                         const bool protect_boundaries,
      const bool smooth_constrained_edges)
      : m_cell_selector(cell_selector)
      , m_protect_boundaries(protect_boundaries)
      , m_smooth_constrained_edges(smooth_constrained_edges)
  {
    refresh(c3t3);

    #ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    createMLSSurfaces(subdomain_FMLS,
                      subdomain_FMLS_indices,
                      m_vertices_normals,
                      m_vertices_surface_indices,
                      c3t3);
    #else
    build_triangle_aabb_tree(c3t3);
    build_segment_aabb_tree(c3t3);
    #endif
  }
  // Refresh context data right before use (called from smooth() method)
  void refresh(C3t3& c3t3) {
    if(!m_protect_boundaries) {
      collect_vertices_surface_indices(c3t3);
      compute_vertices_normals(c3t3);
    }

    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3.triangulation());
    collect_incident_cells(c3t3.triangulation());
  }

  void start_flip_smooth_steps(const C3t3& c3t3) {
    CGAL_assertion(!m_flip_smooth_steps);
    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3.triangulation());

    // once this variable is set to true,
    // m_vertex_id becomes constant and
    // m_free_vertices can refer to it safely
    m_flip_smooth_steps = true;
  }

private:
  void build_segment_aabb_tree(const C3t3& c3t3) {
    for(const Edge& e : c3t3.edges_in_complex()) {
      auto segment = c3t3.triangulation().segment(e);
      m_aabb_segments.push_back(segment);
    }

    m_segments_aabb_tree.rebuild(m_aabb_segments.begin(), m_aabb_segments.end());
    m_segments_aabb_tree.accelerate_distance_queries();
  }

  void build_triangle_aabb_tree(const C3t3& c3t3) {
    for (const Facet& f : c3t3.facets_in_complex()) {
      m_aabb_triangles.push_back(c3t3.triangulation().triangle(f));
    }
    m_triangles_aabb_tree.rebuild(m_aabb_triangles.begin(), m_aabb_triangles.end());
    m_triangles_aabb_tree.accelerate_distance_queries();

        const CGAL::Bbox_3& bb = m_triangles_aabb_tree.bbox();
        m_aabb_epsilon = 1e-3 * (std::min)(bb.xmax() - bb.xmin(),
                                (std::min)(bb.ymax() - bb.ymin(),
                                           bb.zmax() - bb.zmin()));
  }


  void reset_vertex_id_map(const Tr& tr) {
    m_vertex_id.clear();
    std::size_t id = 0;
    for (const Vertex_handle v : tr.finite_vertex_handles()) {
      m_vertex_id[v] = id++;
    }
  }

  bool is_selected(const Cell_handle c) const { return get(m_cell_selector, c); }

  std::size_t vertex_id(const Vertex_handle v) const {
    CGAL_expensive_assertion(m_vertex_id.find(v) != m_vertex_id.end());
    return m_vertex_id.at(v);
  }

  void reset_free_vertices(const Tr& tr) {
    // when flip-smooth steps start,
    // m_free_vertices should already be initialized
    // by the last smoothing step.
    // Then, it does not need to be recomputed
    // because no vertices are inserted nor removed anymore
    if(m_flip_smooth_steps) {
      CGAL_assertion(m_free_vertices.size() == tr.number_of_vertices());
      return;
    }
    m_free_vertices.clear();
    m_free_vertices.resize(tr.number_of_vertices(), false);

    for(const Cell_handle c : tr.finite_cell_handles()) {
      if(!is_selected(c))
        continue;

      for(auto vi : tr.vertices(c)) {
        const std::size_t idi =vertex_id(vi);
        const int dim = vi->in_dimension();

        switch(dim) {
        case 3:
          m_free_vertices[idi] = true;
          break;
        case 2:
          m_free_vertices[idi] = !m_protect_boundaries;
          break;
        case 1:
          m_free_vertices[idi] = !m_protect_boundaries && m_smooth_constrained_edges;
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

  void collect_incident_cells(const Tr& tr) {
    m_inc_cells.clear();
    const std::size_t nbv = tr.number_of_vertices();
    m_inc_cells.resize(nbv, Incident_cells_vector{});

    for(const Cell_handle c : tr.finite_cell_handles()) {
      for(auto vi : tr.vertices(c)) {
        const std::size_t idi = m_vertex_id.at(vi);
        if(m_free_vertices[idi])
          m_inc_cells[idi].push_back(c);
      }
    }
  }

  void collect_vertices_surface_indices(const C3t3& c3t3) {
    m_vertices_surface_indices.clear();
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
    m_vertices_normals.clear();
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


// Base class for shared functionality between different vertex smoothing operations
template<typename C3t3, typename SizingFunction, typename CellSelector>
class VertexSmoothOperationBase
{
protected:
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Cell_handle           Cell_handle;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename Tr::Edge                  Edge;
  typedef typename Tr::Facet                 Facet;

  typedef typename Tr::Geom_traits           Gt;
  typedef typename Gt::Vector_3              Vector_3;
  typedef typename Gt::Point_3               Point_3;
  typedef typename Gt::FT                    FT;

  // AABB tree types (commented out until needed)
   using Triangle_vec = std::vector<typename Tr::Triangle>;
   using Triangle_iter = typename Triangle_vec::iterator;
   using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
   using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
   using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  using Context = VertexSmoothingContext<C3t3, SizingFunction, CellSelector>;

  const SizingFunction& m_sizing;
  const CellSelector& m_cell_selector;
  bool m_protect_boundaries;
  std::shared_ptr<Context> m_context{nullptr}; // Pointer to shared context

public:
  bool m_smooth_constrained_edges;
  VertexSmoothOperationBase(
                       const SizingFunction& sizing,
                       const CellSelector& cell_selector,
                       const bool protect_boundaries,
                       const bool smooth_constrained_edges,
                            std::shared_ptr<Context> context)
    :
     m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
    , m_smooth_constrained_edges(smooth_constrained_edges)
    , m_context(context)
  {}

  void set_context(std::shared_ptr<Context> p_context) {
    m_context = p_context;
  }

protected:

  Point_3 project_on_tangent_plane(const Point_3& gi,
                                    const Point_3& pi,
                                    const Vector_3& normal)
  {
    Vector_3 diff(gi, pi);
    Point_3 result = gi + (normal * diff) * normal;
    return result;
  }

  FT density_along_segment(const Edge& e, const C3t3& c3t3, bool boundary_edge = false) const {
    // Simplified density calculation - in practice this would use sizing function
    const auto [pt, dim, index] = midpoint_with_info(e, boundary_edge, c3t3);
    const FT s = sizing_at_midpoint(e, pt, dim, index, m_sizing, c3t3, m_cell_selector);
    const FT density = 1. / s; //density = 1 / size^(dimension)
    return density;
  }

  bool is_selected(const Cell_handle c) const { return get(m_cell_selector, c); }

  template<typename CellRange>
  Dihedral_angle_cosine max_cosine(const Tr& tr,
    const CellRange& cells) const
  {
    Dihedral_angle_cosine max_cos_dh = cosine_of_90_degrees();// = 0.
    for (Cell_handle c : cells)
    {
      if(!is_selected(c))
        continue;
      Dihedral_angle_cosine cos_dh = max_cos_dihedral_angle(tr, c, false);
      if (max_cos_dh < cos_dh)
        max_cos_dh = cos_dh;
    }
    return max_cos_dh;
  }



  template<typename CellRange, typename Tr>
  bool check_inversion_and_move(const typename Tr::Vertex_handle v,
                                const typename Tr::Geom_traits::Point_3& final_pos,
                                const CellRange& inc_cells,
                                const Tr& tr,
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
                                FT& total_move) const
#else
                                FT&) const
#endif
  {
    // Debug logging for check_inversion_and_move (ComplexEdgeVertexSmoothOperation)

    const typename Tr::Point backup = v->point(); // backup v's position
    const typename Tr::Geom_traits::Point_3 pv = point(backup);

    bool valid_orientation = false;
    bool angles_improved = true;
    double frac = 1.0;
    typename Tr::Geom_traits::Vector_3 move(pv, final_pos);

    // TODO: This variable is only used in the smooth-steps, so it should be removed
    const Dihedral_angle_cosine curr_max_cos = m_context->m_flip_smooth_steps
                                                   ? max_cosine(tr, inc_cells)
                                                   : Dihedral_angle_cosine(CGAL::ZERO, 0., 1.); // Dummy unused value

    bool valid_try = true;
    do {
      v->set_point(typename Tr::Point(pv + frac * move));

      valid_try = true;
      valid_orientation = true;
      angles_improved = true;

      for(const typename Tr::Cell_handle& ci : inc_cells) {
        if(CGAL::POSITIVE != CGAL::orientation(point(ci->vertex(0)->point()), point(ci->vertex(1)->point()),
                                               point(ci->vertex(2)->point()), point(ci->vertex(3)->point())))
        {
          frac = 0.5 * frac;
          valid_try = false;
          valid_orientation = false;
          break;
        } else if(m_context->m_flip_smooth_steps) // check that dihedral angles get improved
        {
          if(is_selected(ci)) {
            Dihedral_angle_cosine max_cos_ci = max_cos_dihedral_angle(tr, ci, false);
            if(curr_max_cos < max_cos_ci) {
              // keep move only if new cosine is smaller than previous one
              // i.e. if angle is larger
              frac = 0.5 * frac;
              valid_try = false;
              angles_improved = false;
              break;
            }
          }
        }
      }
    } while(!valid_try && frac > 0.1);

    // if move failed, cancel move
    bool valid_move = valid_orientation && angles_improved;

      if(!valid_move)
        v->set_point(backup);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      else
        total_move += CGAL::approximate_sqrt(CGAL::squared_distance(pv, point(v->point())));
#endif

    return valid_move;
  }

};

// Internal Vertex Smooth Operation - processes internal vertices
template<typename C3t3, typename SizingFunction, typename CellSelector>
class InternalVertexSmoothOperation
  : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
    public ElementaryOperation<C3t3,
                          typename C3t3::Triangulation::Vertex_handle,
                          typename C3t3::Triangulation::Finite_vertex_handles,
                          typename C3t3::Triangulation::Vertex_handle>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Base = ElementaryOperation<C3t3,
                                   typename C3t3::Triangulation::Vertex_handle,
                                   typename C3t3::Triangulation::Finite_vertex_handles,
                                   typename C3t3::Triangulation::Vertex_handle>;
  using ElementType = typename C3t3::Triangulation::Vertex_handle;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_sizing;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_smooth_constrained_edges;
  using BaseClass::m_context;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Point_3;
  using typename BaseClass::Vector_3;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;

public:
  InternalVertexSmoothOperation(const C3t3& c3t3,
                               const SizingFunction& sizing,
                               const CellSelector& cell_selector,
                               const bool protect_boundaries,
                               const bool smooth_constrained_edges,
                               std::shared_ptr<typename BaseClass::Context> context)
    : BaseClass(sizing, cell_selector, protect_boundaries,smooth_constrained_edges, context)
  {}

  void perform_global_preprocessing(const C3t3& c3t3) const
  {
  auto& tr = c3t3.triangulation();

  const std::size_t nbv = tr.number_of_vertices();
  const typename BaseClass::Context::Move default_move{CGAL::NULL_VECTOR, 0 /*neighbors*/, 0. /*mass*/};
  m_context->m_moves.assign(nbv, default_move);

  for (const Edge& e : tr.finite_edges())
  {
    if (is_outside(e, c3t3, m_cell_selector))
      continue;

      const auto [vh0, vh1] =  make_vertex_pair(e);

      const std::size_t& i0 = m_context->m_vertex_id.at(vh0);
      const std::size_t& i1 = m_context->m_vertex_id.at(vh1);

      const bool vh0_moving = (c3t3.in_dimension(vh0) == 3 && m_context->m_free_vertices[i0]);
      const bool vh1_moving = (c3t3.in_dimension(vh1) == 3 && m_context->m_free_vertices[i1]);

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = BaseClass::density_along_segment(e, c3t3);

      if (vh0_moving)
      {
        m_context->m_moves[i0].move += density * Vector_3(p0, p1);
        m_context->m_moves[i0].mass += density;
        ++m_context->m_moves[i0].neighbors;
      }
      if (vh1_moving)
      {
        m_context->m_moves[i1].move += density * Vector_3(p1, p0);
        m_context->m_moves[i1].mass += density;
        ++m_context->m_moves[i1].neighbors;
      }
  }


  }

  ElementSource get_element_source(const C3t3& c3t3) const override {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }


  bool lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    auto& tr = c3t3.triangulation();
    std::vector<Cell_handle> inc_cells;
    return tr.try_lock_and_get_incident_cells(v, inc_cells);
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
    const std::size_t vid = m_context->m_vertex_id.at(v);
    if(!(m_context->m_free_vertices[vid] && c3t3.in_dimension(v) == 3 && m_context->m_moves[vid].neighbors > 1)){
      return false;
    }

    const Vector_3 move = m_context->m_moves[vid].move / m_context->m_moves[vid].mass;
    Point_3 new_pos = point(v->point()) + move;
    const auto& inc_cells = m_context->m_inc_cells[vid];
    bool result = BaseClass::check_inversion_and_move(v, new_pos, inc_cells, c3t3.triangulation(), m_context->m_total_move);
    return result;
  }

  bool requires_ordered_processing() const override {
    return false; // InternalVertexSmooth can use unordered parallel processing
  }

  std::string operation_name() const override {
    return "Vertex Smooth (Internal Vertices)";
  }
};

// Surface Vertex Smooth Operation - processes surface vertices (2-junctions)
template<typename C3t3, typename SizingFunction, typename CellSelector>
class SurfaceVertexSmoothOperation
  : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
    public ElementaryOperation<C3t3,
                          typename C3t3::Triangulation::Vertex_handle,
                          typename C3t3::Triangulation::Finite_vertex_handles,
                          typename C3t3::Triangulation::Vertex_handle>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Base = ElementaryOperation<C3t3,
                              typename C3t3::Triangulation::Vertex_handle,
                              typename C3t3::Triangulation::Finite_vertex_handles,
                              typename C3t3::Triangulation::Vertex_handle>;
  using ElementType = typename C3t3::Triangulation::Vertex_handle;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using BaseClass::m_sizing;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_smooth_constrained_edges;
  using BaseClass::m_context;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Point_3;
  using typename BaseClass::Vector_3;
  using typename BaseClass::Gt;
  using typename BaseClass::Surface_patch_index;
  using typename BaseClass::FT;
  using typename BaseClass::AABB_triangle_tree;

private:

  void perform_global_preprocessing(const C3t3& c3t3) const {
    auto& tr = c3t3.triangulation();
    const std::size_t nbv = tr.number_of_vertices();

    // Initialize moves vector
    const typename BaseClass::Context::Move default_move{CGAL::NULL_VECTOR, 0 /*neighbors*/, 0. /*mass*/};
    m_context->m_moves.assign(nbv, default_move);

    for (const Edge& e : tr.finite_edges()) {
      if (!c3t3.is_in_complex(e) && is_boundary(c3t3, e, m_cell_selector)) {
        const Vertex_handle vh0 = e.first->vertex(e.second);
        const Vertex_handle vh1 = e.first->vertex(e.third);

        const std::size_t& i0 = m_context->m_vertex_id.at(vh0);
        const std::size_t& i1 = m_context->m_vertex_id.at(vh1);

        const bool vh0_moving = !is_on_feature(vh0) && m_context->m_free_vertices[i0];
        const bool vh1_moving = !is_on_feature(vh1) && m_context->m_free_vertices[i1];

        if (!vh0_moving && !vh1_moving)
          continue;

        const Point_3& p0 = point(vh0->point());
        const Point_3& p1 = point(vh1->point());
        const FT density = BaseClass::density_along_segment(e, c3t3, true);

        if (vh0_moving) {
          m_context->m_moves[i0].move += density * Vector_3(p0, p1);
          m_context->m_moves[i0].mass += density;
          ++m_context->m_moves[i0].neighbors;
        }
        if (vh1_moving) {
          m_context->m_moves[i1].move += density * Vector_3(p1, p0);
          m_context->m_moves[i1].mass += density;
          ++m_context->m_moves[i1].neighbors;
        }
      }
    }
  }

  std::optional<Point_3> project(const Surface_patch_index& si,
                                 const Point_3& gi)
  {
    CGAL_expensive_assertion(m_context->subdomain_FMLS_indices.find(si) != m_context->subdomain_FMLS_indices.end());
    CGAL_assertion(!std::isnan(gi.x()) && !std::isnan(gi.y()) && !std::isnan(gi.z()));

    Vector_3 point(gi.x(), gi.y(), gi.z());
    Vector_3 res_normal = CGAL::NULL_VECTOR;
    Vector_3 result(CGAL::ORIGIN, gi);

    const typename BaseClass::Context::FMLS& fmls
      = m_context->subdomain_FMLS[m_context->subdomain_FMLS_indices.at(si)];

    int it_nb = 0;
    const int max_it_nb = 5;
    const double epsilon = fmls.getPNScale() / 1000.;
    const double sq_eps = CGAL::square(epsilon);

    do
    {
      point = result;

      fmls.fastProjectionCPU(point, result, res_normal);

      if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2])) {
        std::cout << "MLS error detected si " //<< si
                  << "\t(size : "       << fmls.getPNSize() << ")"
                  << "\t(point = "      << point      << " )" << std::endl;
        return {};
      }
    } while ((result - point).squared_length() > sq_eps && ++it_nb < max_it_nb);

    return Point_3(result[0], result[1], result[2]);
  }

public:
  SurfaceVertexSmoothOperation(const C3t3& c3t3,
                              const SizingFunction& sizing,
                              const CellSelector& cell_selector,
                              const bool protect_boundaries,
                              const bool smooth_constrained_edges,
                               std::shared_ptr<typename BaseClass::Context> context)
    : BaseClass( sizing, cell_selector, protect_boundaries, smooth_constrained_edges, context)
  {
  }


  ElementSource get_element_source(const C3t3& c3t3) const override {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }


  bool lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    auto& tr = c3t3.triangulation();
    std::vector<Cell_handle> inc_cells;
    return tr.try_lock_and_get_incident_cells(v, inc_cells);
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
    auto& tr = c3t3.triangulation();
    const std::size_t vid = m_context->m_vertex_id.at(v);
    if(!(m_context->m_free_vertices[vid] && v->in_dimension() == 2)){
      return false;
    }

    const std::size_t nb_neighbors = m_context->m_moves[vid].neighbors;

    const Point_3 current_pos = point(v->point());

    CGAL_assertion(m_context->m_vertices_surface_indices.find(v)!=m_context->m_vertices_surface_indices.end());
    const auto& incident_surface_patches = m_context->m_vertices_surface_indices.at(v);

    if (incident_surface_patches.size() > 1) {
      return false;
    }


    const Surface_patch_index si = incident_surface_patches[0];
    CGAL_assertion(si != Surface_patch_index());
    CGAL_expensive_assertion_code(auto siv = surface_patch_index(v, c3t3));
    CGAL_expensive_assertion(si == siv);

    Point_3 new_pos;
    bool result = false;

    if (nb_neighbors > 1) {
      const Vector_3 move = m_context->m_moves[vid].move / m_context->m_moves[vid].mass;
      const Point_3 smoothed_position = point(v->point()) + move;

#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      Point_3 normal_projection = BaseClass::project_on_tangent_plane(smoothed_position,
                                                           current_pos,
                                                           m_context->m_vertices_normals.at(v).at(si));
      std::optional<Point_3> mls_projection = project(si, normal_projection);

      new_pos = (mls_projection != std::nullopt) ? *mls_projection : smoothed_position;

#else // AABB_tree projection
      if (m_context->m_triangles_aabb_tree.squared_distance(smoothed_position) < m_context->m_aabb_epsilon) {
        new_pos = m_context->m_triangles_aabb_tree.closest_point(smoothed_position);
      } else {
        using Ray = typename Tr::Geom_traits::Ray_3;
        using Projection = std::optional<
          typename AABB_triangle_tree::template Intersection_and_primitive_id<Ray>::Type>;

        auto get_intersection_point = [](const Projection& proj) -> std::optional<Point_3> {
          const auto intersection = proj.value().first;
          if (const Point_3* p = std::get_if<Point_3>(&intersection))
            return *p;
          else
            return std::nullopt;
        };

        auto get_intersection_midpoint = [](const Projection& proj) -> std::optional<Point_3> {
          const auto intersection = proj.value().first;
          using Segment = typename Tr::Geom_traits::Segment_3;
          if (const Segment* s = std::get_if<Segment>(&intersection))
            return CGAL::midpoint(s->source(), s->target());
          else{
            CGAL_assertion(false);
            return std::nullopt;
          }
        };

        const auto n = m_context->m_vertices_normals.at(v).at(si);
        const Ray ray = tr.geom_traits().construct_ray_3_object()(current_pos, n);

        const Projection proj = m_context->m_triangles_aabb_tree.first_intersection(ray);
        const Projection proj_opp = m_context->m_triangles_aabb_tree.first_intersection(
          tr.geom_traits().construct_opposite_ray_3_object()(ray));

        if (proj != std::nullopt && proj_opp == std::nullopt) {
          const auto p = get_intersection_point(proj);
          if (p != std::nullopt)
            new_pos = p.value();
          else
            new_pos = get_intersection_midpoint(proj).value();
        } else if (proj == std::nullopt && proj_opp != std::nullopt) {
          const auto p = get_intersection_point(proj_opp);
          if (p != std::nullopt)
            new_pos = p.value();
          else
            new_pos = get_intersection_midpoint(proj_opp).value();
        } else if (proj != std::nullopt && proj_opp != std::nullopt) {
          const auto op1 = get_intersection_point(proj);
          const auto op2 = get_intersection_point(proj_opp);

          const FT sqd1 = (op1 == std::nullopt) ? 0. : CGAL::squared_distance(smoothed_position, op1.value());
          const FT sqd2 = (op2 == std::nullopt) ? 0. : CGAL::squared_distance(smoothed_position, op2.value());

          if (sqd1 != 0. && sqd1 < sqd2)
            new_pos = op1.value();
          else if (sqd2 != 0)
            new_pos = op2.value();
          else
            new_pos = smoothed_position;
        } else {
          new_pos = smoothed_position;
        }
      }
#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      const auto& inc_cells = m_context->m_inc_cells[vid];

      result = BaseClass::check_inversion_and_move(v, new_pos, inc_cells, tr, m_context->m_total_move);

    } else if (nb_neighbors > 0) {
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      std::optional<Point_3> mls_proj = project(si, current_pos);
      if (mls_proj == std::nullopt) {
        return false;
      }

      new_pos = *mls_proj;
#else
      new_pos = m_context->m_segments_aabb_tree.closest_point(current_pos);
#endif
      const auto& inc_cells = m_context->m_inc_cells[vid];

      result = BaseClass::check_inversion_and_move(v, new_pos, inc_cells, tr, m_context->m_total_move);
    }

    return result;
  }

  bool requires_ordered_processing() const override {
    return false; // SurfaceVertexSmooth can use unordered parallel processing
  }

  std::string operation_name() const override {
    return "Vertex Smooth (Surface Vertices)";
  }
};

// Complex Edge Vertex Smooth Operation - processes vertices on complex edges (1-junctions)
template<typename C3t3, typename SizingFunction, typename CellSelector>
class ComplexEdgeVertexSmoothOperation
  : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
    public ElementaryOperation<C3t3,
                          typename C3t3::Triangulation::Vertex_handle,
                          typename C3t3::Triangulation::Finite_vertex_handles,
                          typename C3t3::Triangulation::Vertex_handle>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Vertex_handle = typename C3t3::Triangulation::Vertex_handle;

  using Base = ElementaryOperation<C3t3,
                                   Vertex_handle,
                                   typename C3t3::Triangulation::Finite_vertex_handles,
                                   Vertex_handle>;
  using ElementType = Vertex_handle;
  using ElementSource = typename Base::ElementSource;
  using Lock_zone = typename Base::Lock_zone;

  using Surface_patch_index = typename C3t3::Surface_patch_index;

  using BaseClass::m_sizing;
  using BaseClass::m_cell_selector;
  using BaseClass::m_protect_boundaries;
  using BaseClass::m_smooth_constrained_edges;
  using BaseClass::m_context;

  // Import types from base class
  using typename BaseClass::Tr;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Point_3;
  using typename BaseClass::Vector_3;
  using typename BaseClass::FT;

public:

  ComplexEdgeVertexSmoothOperation(const C3t3& c3t3,
                                  const SizingFunction& sizing,
                                  const CellSelector& cell_selector,
                                  const bool protect_boundaries,
                                  const bool smooth_constrained_edges,
                                  std::shared_ptr<typename BaseClass::Context> context)
    : BaseClass( sizing, cell_selector, protect_boundaries,smooth_constrained_edges, context)
  {}

  ElementSource get_element_source(const C3t3& c3t3) const override {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }

  // NOTE: There is an incident cells orientation check during smoothing by check_inversion_and_move
  bool lock_zone(const ElementType& v, const C3t3& c3t3) const override {
    auto& tr = c3t3.triangulation();
    std::vector<Cell_handle> inc_cells;
    return tr.try_lock_and_get_incident_cells(v, inc_cells);
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override {
    auto& tr = c3t3.triangulation();
    const std::size_t vid = m_context->m_vertex_id.at(v);
    if (!(m_context->m_free_vertices[vid] && is_on_feature(v)))
      return false;

    const Point_3 current_pos = point(v->point());
    const auto& moves = m_context->m_moves;

    if (moves[vid].neighbors == 0)
      return false;

    CGAL_assertion(moves[vid].mass > 0);
    const Vector_3 move = moves[vid].move / moves[vid].mass;

    const Point_3 smoothed_position = current_pos + move;

    //TODO: Test #ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      Vector_3 sum_projections = CGAL::NULL_VECTOR;
      Point_3 tmp_pos = current_pos;

#ifndef CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION
      const std::vector<Surface_patch_index>& v_surface_indices = m_context->m_vertices_surface_indices.at(v);
      for (const Surface_patch_index& si : v_surface_indices)
      {
        Point_3 normal_projection = BaseClass::project_on_tangent_plane(smoothed_position,
                                                             current_pos,
                                                             m_context->m_vertices_normals.at(v).at(si));

        sum_projections += Vector_3(tmp_pos, normal_projection);
        tmp_pos = normal_projection;
      }
#endif //CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION

      const Point_3 new_pos = current_pos + sum_projections;

#else // AABB_tree projection
    // Use AABB tree projection for complex edges
    const Point_3 new_pos = m_context->m_segments_aabb_tree.closest_point(smoothed_position);
#endif

    const auto& inc_cells = m_context->m_inc_cells[vid];

    if (BaseClass::check_inversion_and_move(v, new_pos, inc_cells, tr,m_context->m_total_move)) {
      return true;
    }
    return false;
  }

  bool requires_ordered_processing() const override {
    return false; // ComplexEdgeVertexSmooth can use unordered parallel processing
  }

  std::string operation_name() const override {
    return "Vertex Smooth (Complex Edge Vertices)";
  }

  // Debug method to print final total_move
  void print_final_total_move(const C3t3& c3t3) const {
    // Method content removed - no longer needed without debug logging
  }

  // This method will be called by the execution framework to collect moves
  void perform_global_preprocessing(const C3t3& c3t3) const {
    auto& tr = c3t3.triangulation();
    const std::size_t nbv = tr.number_of_vertices();

    // Initialize moves vector
    m_context->m_moves.assign(nbv, typename BaseClass::Context::Move{CGAL::NULL_VECTOR, 0, 0.});

    // Collect moves from complex edges
    for (const Edge& e : c3t3.edges_in_complex()) {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      CGAL_expensive_assertion(is_on_feature(vh0));
      CGAL_expensive_assertion(is_on_feature(vh1));

      const std::size_t& i0 = m_context->m_vertex_id.at(vh0);
      const std::size_t& i1 = m_context->m_vertex_id.at(vh1);

      const bool vh0_moving = m_context->m_free_vertices[i0];
      const bool vh1_moving = m_context->m_free_vertices[i1];

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = BaseClass::density_along_segment(e, c3t3, true);

      if (vh0_moving) {
        Vector_3 move_vector = density * Vector_3(p0, p1);
        m_context->m_moves[i0].move += move_vector;
        m_context->m_moves[i0].mass += density;
        ++m_context->m_moves[i0].neighbors;
      }
      if (vh1_moving) {
        Vector_3 move_vector = density * Vector_3(p1, p0);
        m_context->m_moves[i1].move += move_vector;
        m_context->m_moves[i1].mass += density;
        ++m_context->m_moves[i1].neighbors;
      }
    }
  }
};

} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_VERTEX_SMOOTH_OPERATION_H