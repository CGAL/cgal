// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_INTERNAL_SMOOTH_VERTICES_H
#define CGAL_INTERNAL_SMOOTH_VERTICES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Vector_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>

#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_segment_primitive_3.h>
#include <CGAL/use.h>

#include <optional>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <unordered_map>
#include <vector>
#include <cmath>
#include <list>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
template <typename C3t3, typename SizingFunction, typename CellSelector>
class Vertex_smoothing_context
{
public:
  using Tr                  = typename C3t3::Triangulation;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Cell_handle         = typename Tr::Cell_handle;
  using Vertex_handle       = typename Tr::Vertex_handle;
  using Edge                = typename Tr::Edge;
  using Facet               = typename Tr::Facet;

  using Gt       = typename Tr::Geom_traits;
  using Vector_3 = typename Gt::Vector_3;
  using Point_3  = typename Gt::Point_3;
  using FT       = typename Gt::FT;

  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  using Segment_vec = std::vector<typename Gt::Segment_3>;
  using Segment_iter = typename Segment_vec::iterator;
  using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
  using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
  using AABB_segment_tree = CGAL::AABB_tree<AABB_segment_traits>;

private:
  using FMLS = CGAL::Tetrahedral_remeshing::internal::FMLS<Gt>;
  std::vector<FMLS> subdomain_FMLS;
  std::unordered_map<Surface_patch_index, std::size_t, boost::hash<Surface_patch_index>> subdomain_FMLS_indices;

  Triangle_vec m_aabb_triangles;
  Segment_vec m_aabb_segments;

public:
  AABB_triangle_tree m_triangles_aabb_tree;
  AABB_segment_tree m_segments_aabb_tree;
  FT m_aabb_epsilon;

  const SizingFunction& m_sizing;

  using Incident_cells_vector = boost::container::small_vector<Cell_handle, 64>;
  std::vector<Incident_cells_vector> m_inc_cells;

  using Vertices_surface_indices_map = std::unordered_map<Vertex_handle, std::vector<Surface_patch_index>>;
  using Vertices_normals_map =
      std::unordered_map<Vertex_handle,
                         std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>>;

  Vertices_surface_indices_map m_vertices_surface_indices;
  Vertices_normals_map m_vertices_normals;

  const CellSelector& m_cell_selector;
  const bool m_protect_boundaries;

  const bool m_smooth_constrained_edges;

  // the 2 following variables become useful and valid
  // just before flip/smooth steps, when no vertices get inserted
  // nor removed anymore
  std::unordered_map<Vertex_handle, std::size_t> m_vertex_id;
  std::vector<bool> m_free_vertices{};
  bool m_flip_smooth_steps{false};

public:
  struct Move
  {
    Vector_3 move;
    int neighbors;
    FT mass;
  };
  std::vector<Move> m_moves{};
  FT m_total_move{0};

  Vertex_smoothing_context(C3t3& c3t3,
                           const SizingFunction& sizing,
                           const CellSelector& cell_selector,
                           const bool protect_boundaries,
                           const bool smooth_constrained_edges)
      : m_cell_selector(cell_selector)
      , m_sizing(sizing)
      , m_protect_boundaries(protect_boundaries)
      , m_smooth_constrained_edges(smooth_constrained_edges)
  {
    refresh(c3t3);
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    if (m_protect_boundaries)
    {
      collect_vertices_surface_indices(c3t3);
      compute_vertices_normals(c3t3);
    }
    createMLSSurfaces(subdomain_FMLS,
                      subdomain_FMLS_indices,
                      m_vertices_normals,
                      m_vertices_surface_indices,
                      c3t3);
#else
    build_aabb_trees(c3t3);
#endif
  }

  void refresh(C3t3& c3t3)
  {
    if (!m_protect_boundaries)
    {
      collect_vertices_surface_indices(c3t3);
      compute_vertices_normals(c3t3);
    }
    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3.triangulation());
    collect_incident_cells(c3t3.triangulation());
  }

  void start_flip_smooth_steps(const C3t3& c3t3)
  {
    CGAL_assertion(!m_flip_smooth_steps);
    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3.triangulation());

    // once this variable is set to true,
    // m_vertex_id becomes constant and
    // m_free_vertices can refer to it safely
    m_flip_smooth_steps = true;
  }

  std::size_t vertex_id(const Vertex_handle v) const
  {
    CGAL_expensive_assertion(m_vertex_id.find(v) != m_vertex_id.end());
    return m_vertex_id.at(v);
  }

  bool is_free(const Vertex_handle v) const  { return m_free_vertices[vertex_id(v)]; }
  bool is_free(const std::size_t& vid) const { return m_free_vertices[vid]; }

  const Incident_cells_vector& incident_cells(const Vertex_handle v) const
  {
    return m_inc_cells[vertex_id(v)];
  }
  const Incident_cells_vector& incident_cells(const std::size_t& vid) const
  {
    return m_inc_cells[vid];
  }

private:
  bool is_selected(const Cell_handle c) const
  {
    return get(m_cell_selector, c);
  }

  // this function can be used iff m_vertex_id
  // has already been initialized
  void reset_free_vertices(const Tr& tr)
  {
    // when flip-smooth steps start,
    // m_free_vertices should already be initialized
    // by the last smoothing step.
    // Then, it does not need to be recomputed
    // because no vertices are inserted nor removed anymore
    if (m_flip_smooth_steps)
    {
      CGAL_assertion(m_free_vertices.size() == tr.number_of_vertices());
      return;
    }
    m_free_vertices.clear();
    m_free_vertices.resize(tr.number_of_vertices(), false);

    for (const Cell_handle c : tr.finite_cell_handles())
    {
      if (!is_selected(c))
        continue;

      for (auto vi : tr.vertices(c))
      {
        const std::size_t idi = vertex_id(vi);
        const int dim = vi->in_dimension();

        switch (dim)
        {
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

  void build_aabb_trees(const C3t3& c3t3)
  {
    // build AABB tree of facets in complex
    for (const Facet& f : c3t3.facets_in_complex())
    {
      m_aabb_triangles.push_back(c3t3.triangulation().triangle(f));
    }
    m_triangles_aabb_tree.rebuild(m_aabb_triangles.begin(), m_aabb_triangles.end());
    m_triangles_aabb_tree.accelerate_distance_queries();

    // build AABB tree of edges in complex
    for (const Edge& e : c3t3.edges_in_complex())
    {
      m_aabb_segments.push_back(c3t3.triangulation().segment(e));
    }
    m_segments_aabb_tree.rebuild(m_aabb_segments.begin(), m_aabb_segments.end());
    m_segments_aabb_tree.accelerate_distance_queries();

    // compute epsilon for AABB tree of facets
    const CGAL::Bbox_3& bb = m_triangles_aabb_tree.bbox();
    m_aabb_epsilon = 1e-3 * (std::min)(bb.xmax() - bb.xmin(),
                            (std::min)(bb.ymax() - bb.ymin(),
                                       bb.zmax() - bb.zmin()));
  }

  void collect_incident_cells(const Tr& tr)
  {
    m_inc_cells.clear();
    const std::size_t nbv = tr.number_of_vertices();
    m_inc_cells.resize(nbv, Incident_cells_vector{});
    for (const Cell_handle c : tr.finite_cell_handles())
    {
      for (auto vi : tr.vertices(c))
      {
        const std::size_t idi = vertex_id(vi);
        if(is_free(idi))
          m_inc_cells[idi].push_back(c);
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

  template<typename Gt_>
  Vector_3 compute_normal(const Facet& f,
                          const Vector_3& reference_normal,
                          const Gt_& gt)
  {
    typename Gt_::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();
    typename Gt_::Compute_scalar_product_3
      scalar_product = gt.compute_scalar_product_3_object();

    Vector_3 n = CGAL::Tetrahedral_remeshing::normal(f, gt);
    if (scalar_product(n, reference_normal) < 0.)
      n = opp(n);

    return n;
  }

  template<typename Patch_index>
  std::string debug_to_string(const Patch_index i)
  {
    return std::to_string(i);
  }

  template<typename Patch_index>
  std::string debug_to_string(const std::pair<Patch_index, Patch_index>& pi)
  {
    std::string str = std::to_string(pi.first);
    str.append("_").append(std::to_string(pi.second));
    return str;
  }

  void compute_vertices_normals(const C3t3& c3t3)
  {
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

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::ofstream osf("dump_facet_normals.polylines.txt");
#endif
    for (const auto& [f, n] : fnormals)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      typename Tr::Geom_traits::Point_3 fc = CGAL::centroid(tr.triangle(f));
      osf << "2 " << fc << " " << (fc + n) << std::endl;
#endif
      const Surface_patch_index& surf_i = c3t3.surface_patch_index(f);

      for (const Vertex_handle vi : tr.vertices(f))
      {
        typename Vertices_normals_map::iterator patch_vector_it = m_vertices_normals.find(vi);

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

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    osf.close();
    std::ofstream os("dump_normals.polylines.txt");
    std::unordered_map<Surface_patch_index,
          std::vector<typename Tr::Geom_traits::Segment_3 >, boost::hash<Surface_patch_index> > ons_map;
#endif

    //normalize the computed normals
    for (auto& [v, patch_normals] : m_vertices_normals)
    {
      //value type is map<Surface_patch_index, Vector_3>
      for (auto& [surf_i, n] : patch_normals)
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        auto p = point(v->point());
        os << "2 " << p << " " << (p + n) << std::endl;
#endif

        CGAL::Tetrahedral_remeshing::normalize(n, gt);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        const Surface_patch_index si = surf_i;
        if (ons_map.find(si) == ons_map.end())
          ons_map[si] = std::vector<typename Tr::Geom_traits::Segment_3>();
        ons_map[si].push_back(typename Tr::Geom_traits::Segment_3(p, p + n));
#endif
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    os.close();
    for (auto& kv : ons_map)
    {
      std::ostringstream oss;
      oss << "dump_normals_normalized_["
        << debug_to_string(kv.first) << "].polylines.txt";
      std::ofstream ons(oss.str());
      for (auto s : kv.second)
        ons << "2 " << s.source() << " " << s.target() << std::endl;
      ons.close();
    }
#endif
  }

  void collect_vertices_surface_indices(const C3t3& c3t3)
  {
    m_vertices_surface_indices.clear();
    for (Facet fit : c3t3.facets_in_complex())
    {
      const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);

      for (const Vertex_handle vi : c3t3.triangulation().vertices(fit))
      {
        std::vector<Surface_patch_index>& v_surface_indices = m_vertices_surface_indices[vi];
        if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index) == v_surface_indices.end())
          v_surface_indices.push_back(surface_index);
      }
    }
  }

  void reset_vertex_id_map(const Tr& tr)
  {
    // when flip-smooth steps start,
    // m_vertex_id should already be initialized,
    // done by the last smoothing step.
    // Then, it does not need to be recomputed
    // because no vertices are inserted nor removed anymore
    if(m_flip_smooth_steps)
      return;
    m_vertex_id.clear();
    std::size_t id = 0;
    for (const Vertex_handle v : tr.finite_vertex_handles())
    {
      m_vertex_id[v] = id++;
    }
  }
};


template <typename C3t3, typename SizingFunction, typename CellSelector>
class VertexSmoothOperationBase
{
protected:
  typedef typename C3t3::Triangulation Tr;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Facet Facet;

  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3 Point_3;
  typedef typename Gt::FT FT;

  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  using Context = Vertex_smoothing_context<C3t3, SizingFunction, CellSelector>;

public:
  std::shared_ptr<Context> m_context{nullptr};

  VertexSmoothOperationBase(std::shared_ptr<Context> context)
      : m_context(context) {}

  void set_context(std::shared_ptr<Context> p_context) { m_context = p_context; }

protected:
  Point_3 project_on_tangent_plane(const Point_3& gi, const Point_3& pi, const Vector_3& normal)
  {
    Vector_3 diff(gi, pi);
    return gi + (normal * diff) * normal;
  }

  FT density_along_segment(const Edge& e, const C3t3& c3t3, bool boundary_edge = false) const
  {
    const auto [pt, dim, index] = midpoint_with_info(e, boundary_edge, c3t3);
    const FT s = sizing_at_midpoint(e, pt, dim, index, m_context->m_sizing, c3t3, m_context->m_cell_selector);
    return 1. / s;
  }

  bool is_selected(const Cell_handle c) const { return get(m_context->m_cell_selector, c); }

  template <typename CellRange>
  Dihedral_angle_cosine max_cosine(const Tr& tr, const CellRange& cells) const
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

  // in flip-smooth steps, this function also checks that it improves
  // dihedral angles
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
    const typename Tr::Point backup = v->point();//backup v's position
    const typename Tr::Geom_traits::Point_3 pv = point(backup);

    bool valid_orientation = false;
    bool angles_improved = true;
    double frac = 1.0;
    typename Tr::Geom_traits::Vector_3 move(pv, final_pos);

    const Dihedral_angle_cosine curr_max_cos = m_context->m_flip_smooth_steps
      ? max_cosine(tr, inc_cells)
      : Dihedral_angle_cosine(CGAL::ZERO, 0., 1.);//Dummy unused value

    bool valid_try = true;
    do
    {
      v->set_point(typename Tr::Point(pv + frac * move));

      valid_try = true;
      valid_orientation = true;
      angles_improved = true;

      for (const typename Tr::Cell_handle& ci : inc_cells)
      {
        if (CGAL::POSITIVE != CGAL::orientation(point(ci->vertex(0)->point()),
                                                point(ci->vertex(1)->point()),
                                                point(ci->vertex(2)->point()),
                                                point(ci->vertex(3)->point())))
        {
          frac = 0.5 * frac;
          valid_try = false;
          valid_orientation = false;
          break;
        }
        else if (m_context->m_flip_smooth_steps) //check that dihedral angles get improved
        {
          if(is_selected(ci))
          {
            Dihedral_angle_cosine max_cos_ci = max_cos_dihedral_angle(tr, ci, false);
            if (curr_max_cos < max_cos_ci)
            {
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
    }
    while(!valid_try && frac > 0.1);

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

template <typename C3t3, typename SizingFunction, typename CellSelector>
class ComplexEdgeVertexSmoothOperation
    : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
      public ElementaryOperation<C3t3,
                                 typename C3t3::Triangulation::Vertex_handle,
                                 typename C3t3::Triangulation::Finite_vertex_handles>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using Vertex_handle = typename C3t3::Triangulation::Vertex_handle;
  using Surface_patch_index = typename C3t3::Surface_patch_index;

  using BaseOperation = ElementaryOperation<C3t3,
                                   Vertex_handle,
                                   typename C3t3::Triangulation::Finite_vertex_handles>;
  using ElementType = typename BaseOperation::ElementType;
  static_assert(std::is_same_v<ElementType, Vertex_handle>, "ElementType should be Vertex_handle"); 
  using ElementSource = typename BaseOperation::ElementSource;

  using BaseClass::m_context;

  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;
  using typename BaseClass::Point_3;
  using typename BaseClass::Tr;
  using typename BaseClass::Vector_3;

public:
  ComplexEdgeVertexSmoothOperation(std::shared_ptr<typename BaseClass::Context> context)
      : BaseClass(context) {}

  ElementSource get_element_source(const C3t3& c3t3) const override
  {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override
  {
    auto& tr = c3t3.triangulation();
    const std::size_t vid = m_context->vertex_id(v);
    if (!m_context->is_free(vid) || !is_on_feature(v))
      return false;

    const Point_3 current_pos = point(v->point());
    const auto& moves = m_context->m_moves;

    const std::size_t nb_neighbors = moves[vid].neighbors;
    if (nb_neighbors == 0)
      return false;

    CGAL_assertion(moves[vid].mass > 0);
    const Vector_3 move = (nb_neighbors > 0)
                        ? moves[vid].move / moves[vid].mass
                        : CGAL::NULL_VECTOR;
    const Point_3 smoothed_position = current_pos + move;

#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    Vector_3 sum_projections = CGAL::NULL_VECTOR;
    Point_3 tmp_pos = current_pos;

#ifndef CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION
    const std::vector<Surface_patch_index>& v_surface_indices = m_context->m_vertices_surface_indices.at(v);
    for (const Surface_patch_index& si : v_surface_indices)
    {
      Point_3 normal_projection = BaseClass::project_on_tangent_plane(smoothed_position, current_pos,
                                                                      m_context->m_vertices_normals.at(v).at(si));
      sum_projections += Vector_3(tmp_pos, normal_projection);
      tmp_pos = normal_projection;
    }
#endif

    const Point_3 new_pos = current_pos + sum_projections;
#else
    const Point_3 new_pos = m_context->m_segments_aabb_tree.closest_point(smoothed_position);
#endif

    const auto& inc_cells = m_context->m_inc_cells[vid];
    return BaseClass::check_inversion_and_move(v, new_pos, inc_cells, tr, m_context->m_total_move);
  }

  std::string operation_name() const override { return "Vertex Smooth (Complex Edge Vertices)"; }

  void perform_global_preprocessing(const C3t3& c3t3) const
  {
    auto& tr = c3t3.triangulation();
    auto& moves = m_context->m_moves;

    const std::size_t nbv = tr.number_of_vertices();
    using Move = typename BaseClass::Context::Move;
    const Move default_move{CGAL::NULL_VECTOR, 0 /*neighbors*/, 0. /*mass*/};
    moves.assign(nbv, default_move);

    //collect neighbors
    for (const Edge& e : c3t3.edges_in_complex())
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      CGAL_expensive_assertion(is_on_feature(vh0));
      CGAL_expensive_assertion(is_on_feature(vh1));

      const std::size_t& i0 = m_context->vertex_id(vh0);
      const std::size_t& i1 = m_context->vertex_id(vh1);

      const bool vh0_moving = m_context->is_free(i0);
      const bool vh1_moving = m_context->is_free(i1);

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = BaseClass::density_along_segment(e, c3t3, true);

      if (vh0_moving)
      {
        moves[i0].move += density * Vector_3(p0, p1);
        moves[i0].mass += density;
        ++moves[i0].neighbors;
      }
      if (vh1_moving)
      {
        moves[i1].move += density * Vector_3(p1, p0);
        moves[i1].mass += density;
        ++moves[i1].neighbors;
      }
    }
  }
};

template <typename C3t3, typename SizingFunction, typename CellSelector>
class SurfaceVertexSmoothOperation
    : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
      public ElementaryOperation<C3t3,
                                 typename C3t3::Triangulation::Vertex_handle,
                                 typename C3t3::Triangulation::Finite_vertex_handles>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using BaseOperation = ElementaryOperation<C3t3,
                                   typename C3t3::Triangulation::Vertex_handle,
                                   typename C3t3::Triangulation::Finite_vertex_handles>;
  using ElementType = typename BaseOperation::ElementType;
  static_assert(std::is_same_v<ElementType, typename C3t3::Triangulation::Vertex_handle>,
                               "ElementType should be Vertex_handle");
  using ElementSource = typename BaseOperation::ElementSource;

  using BaseClass::m_context;

  using typename BaseClass::AABB_triangle_tree;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;
  using typename BaseClass::Gt;
  using typename BaseClass::Point_3;
  using typename BaseClass::Surface_patch_index;
  using typename BaseClass::Tr;
  using typename BaseClass::Vector_3;
  using typename BaseClass::Vertex_handle;

private:
  void perform_global_preprocessing(const C3t3& c3t3) const
  {
    auto& tr = c3t3.triangulation();
    auto& moves = m_context->m_moves;
    using Move = typename BaseClass::Context::Move;
    const std::size_t nbv = tr.number_of_vertices();
    const Move default_move{CGAL::NULL_VECTOR, 0/*neighbors*/, 0./*mass*/};
    moves.assign(nbv, default_move);

    for (const Edge& e : tr.finite_edges())
    {
      if (!c3t3.is_in_complex(e) && is_boundary(c3t3, e, m_context->m_cell_selector))
      {
        const Vertex_handle vh0 = e.first->vertex(e.second);
        const Vertex_handle vh1 = e.first->vertex(e.third);

        const std::size_t& i0 = m_context->vertex_id(vh0);
        const std::size_t& i1 = m_context->vertex_id(vh1);

        const bool vh0_moving = !is_on_feature(vh0) && m_context->is_free(i0);
        const bool vh1_moving = !is_on_feature(vh1) && m_context->is_free(i1);

        if (!vh0_moving && !vh1_moving)
          continue;

        const Point_3& p0 = point(vh0->point());
        const Point_3& p1 = point(vh1->point());
        const FT density = BaseClass::density_along_segment(e, c3t3, true);

        if (vh0_moving)
        {
          moves[i0].move += density * Vector_3(p0, p1);
          moves[i0].mass += density;
          ++moves[i0].neighbors;
        }
        if (vh1_moving)
        {
          moves[i1].move += density * Vector_3(p1, p0);
          moves[i1].mass += density;
          ++moves[i1].neighbors;
        }
      }
    }
  }

  std::optional<Point_3> project(const Surface_patch_index& si, const Point_3& gi)
  {
    CGAL_expensive_assertion(m_context->subdomain_FMLS_indices.find(si) != m_context->subdomain_FMLS_indices.end());
    CGAL_assertion(!std::isnan(gi.x()) && !std::isnan(gi.y()) && !std::isnan(gi.z()));

    Vector_3 point_vec(gi.x(), gi.y(), gi.z());
    Vector_3 res_normal = CGAL::NULL_VECTOR;
    Vector_3 result(CGAL::ORIGIN, gi);

    const typename BaseClass::Context::FMLS& fmls = m_context->subdomain_FMLS[m_context->subdomain_FMLS_indices.at(si)];

    int it_nb = 0;
    const int max_it_nb = 5;
    const double epsilon = fmls.getPNScale() / 1000.;
    const double sq_eps = CGAL::square(epsilon);

    do
    {
      point_vec = result;
      fmls.fastProjectionCPU(point_vec, result, res_normal);
      if(std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2]))
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << "MLS error detected si "
                  << "\t(size : " << fmls.getPNSize() << ")"
                  << "\t(point = " << point_vec << " )" << std::endl;
#endif
        return {};
      }
    } while((result - point_vec).squared_length() > sq_eps && ++it_nb < max_it_nb);

    return Point_3(result[0], result[1], result[2]);
  }

public:
  SurfaceVertexSmoothOperation(std::shared_ptr<typename BaseClass::Context> context)
      : BaseClass(context) {}

  ElementSource get_element_source(const C3t3& c3t3) const override
  {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override
  {
    auto& tr = c3t3.triangulation();
    auto& moves = m_context->m_moves;
    const std::size_t vid = m_context->vertex_id(v);
    if (!m_context->is_free(vid) || v->in_dimension() != 2)
      return false;

    const std::size_t nb_neighbors = moves[vid].neighbors;
    const Point_3 current_pos = point(v->point());

    CGAL_assertion(m_context->m_vertices_surface_indices.find(v) != m_context->m_vertices_surface_indices.end());
    const auto& incident_surface_patches = m_context->m_vertices_surface_indices.at(v);

    if (incident_surface_patches.size() > 1)
      return false;

    const Surface_patch_index si = incident_surface_patches[0];
    CGAL_assertion(si != Surface_patch_index());
    CGAL_expensive_assertion_code(auto siv = surface_patch_index(v, c3t3));
    CGAL_expensive_assertion(si == siv);

    Point_3 new_pos;
    bool result = false;

    if (nb_neighbors > 1)
    {
      const Vector_3 move = moves[vid].move / moves[vid].mass;
      const Point_3 smoothed_position = point(v->point()) + move;

#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      Point_3 normal_projection = BaseClass::project_on_tangent_plane(smoothed_position, current_pos,
                                                                      m_context->m_vertices_normals.at(v).at(si));
      std::optional<Point_3> mls_projection = project(si, normal_projection);
      new_pos = (mls_projection != std::nullopt) ? *mls_projection : smoothed_position;
#else
      if(m_context->m_triangles_aabb_tree.squared_distance(smoothed_position) < m_context->m_aabb_epsilon)
      {
        new_pos = m_context->m_triangles_aabb_tree.closest_point(smoothed_position);
      }
      else
      {
        using Ray = typename Tr::Geom_traits::Ray_3;
        using Projection =
            std::optional<typename AABB_triangle_tree::template Intersection_and_primitive_id<Ray>::Type>;

        auto get_intersection_point = [](const Projection& proj) -> std::optional<Point_3> {
          const auto intersection = proj.value().first;
          if(const Point_3* p = std::get_if<Point_3>(&intersection))
            return *p;
          return std::nullopt;
        };

        auto get_intersection_midpoint = [](const Projection& proj) -> std::optional<Point_3> {
          const auto intersection = proj.value().first;
          using Segment = typename Tr::Geom_traits::Segment_3;
          if(const Segment* s = std::get_if<Segment>(&intersection))
            return CGAL::midpoint(s->source(), s->target());
          CGAL_assertion(false);
          return std::nullopt;
        };

        const auto n = m_context->m_vertices_normals.at(v).at(si);
        const Ray ray = tr.geom_traits().construct_ray_3_object()(current_pos, n);
        const Projection proj = m_context->m_triangles_aabb_tree.first_intersection(ray);
        const Projection proj_opp = m_context->m_triangles_aabb_tree.first_intersection(
            tr.geom_traits().construct_opposite_ray_3_object()(ray));

        if(proj != std::nullopt && proj_opp == std::nullopt)
        {
          const auto p = get_intersection_point(proj);
          new_pos = (p != std::nullopt) ? p.value() : get_intersection_midpoint(proj).value();
        }
        else if(proj == std::nullopt && proj_opp != std::nullopt)
        {
          const auto p = get_intersection_point(proj_opp);
          new_pos = (p != std::nullopt)
                    ? p.value()
                    : get_intersection_midpoint(proj_opp).value();
        }
        else if(proj != std::nullopt && proj_opp != std::nullopt)
        {
          const auto op1 = get_intersection_point(proj);
          const auto op2 = get_intersection_point(proj_opp);

          const FT sqd1 = (op1 == std::nullopt) ? 0.
            : CGAL::squared_distance(smoothed_position, op1.value());
          const FT sqd2 = (op2 == std::nullopt) ? 0.
            : CGAL::squared_distance(smoothed_position, op2.value());

          if (sqd1 != 0. && sqd1 < sqd2)
            new_pos = op1.value();
          else if (sqd2 != 0)
            new_pos = op2.value();
          else
            new_pos = smoothed_position;
        }
        else //no valid projection
        {
          new_pos = smoothed_position;
        }
      }
#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      const auto& inc_cells = m_context->m_inc_cells[vid];
      result = BaseClass::check_inversion_and_move(v, new_pos, inc_cells, tr, m_context->m_total_move);
    }
    else if (nb_neighbors > 0)
    {
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      std::optional<Point_3> mls_proj = project(si, current_pos);
      if(mls_proj == std::nullopt)
        continue;

      new_pos = *mls_proj;
#else // AABB_tree projection
      new_pos = m_context->m_segments_aabb_tree.closest_point(current_pos);
#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      const auto& inc_cells = m_context->m_inc_cells[vid];
      result = BaseClass::check_inversion_and_move(v, new_pos, inc_cells, tr, m_context->m_total_move);
    }

    return result;
  }

  std::string operation_name() const override { return "Vertex Smooth (Surface Vertices)"; }
};

template <typename C3t3, typename SizingFunction, typename CellSelector>
class InternalVertexSmoothOperation
    : public VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>,
      public ElementaryOperation<C3t3,
                                 typename C3t3::Triangulation::Vertex_handle,
                                 typename C3t3::Triangulation::Finite_vertex_handles>
{
public:
  using BaseClass = VertexSmoothOperationBase<C3t3, SizingFunction, CellSelector>;
  using BaseOperation = ElementaryOperation<C3t3,
                                   typename C3t3::Triangulation::Vertex_handle,
                                   typename C3t3::Triangulation::Finite_vertex_handles>;
  using ElementType = typename BaseOperation::ElementType;
  static_assert(std::is_same_v<ElementType, typename C3t3::Triangulation::Vertex_handle>,
                               "ElementType should be Vertex_handle");
  using ElementSource = typename BaseOperation::ElementSource;

  using BaseClass::m_context;

  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;
  using typename BaseClass::Point_3;
  using typename BaseClass::Tr;
  using typename BaseClass::Vector_3;
  using typename BaseClass::Vertex_handle;

public:
  InternalVertexSmoothOperation(std::shared_ptr<typename BaseClass::Context> context)
      : BaseClass(context) {}

  void perform_global_preprocessing(const C3t3& c3t3) const
  {
    auto& tr = c3t3.triangulation();
    auto& moves = m_context->m_moves;

    using Move = typename BaseClass::Context::Move;
    const std::size_t nbv = tr.number_of_vertices();
    const Move default_move{CGAL::NULL_VECTOR, 0 /*neighbors*/, 0. /*mass*/};
    moves.assign(nbv, default_move);
    /*for dim 3 vertices, start counting neighbors directly from 0*/

    for (const Edge& e : tr.finite_edges())
    {
      if (is_outside(e, c3t3, m_context->m_cell_selector))
        continue;
      else
      {
        const auto [vh0, vh1] = make_vertex_pair(e);

        const std::size_t& i0 = m_context->vertex_id(vh0);
        const std::size_t& i1 = m_context->vertex_id(vh1);

        const bool vh0_moving = (c3t3.in_dimension(vh0) == 3 && m_context->is_free(i0));
        const bool vh1_moving = (c3t3.in_dimension(vh1) == 3 && m_context->is_free(i1));

        if (!vh0_moving && !vh1_moving)
          continue;

        const Point_3& p0 = point(vh0->point());
        const Point_3& p1 = point(vh1->point());
        const FT density = BaseClass::density_along_segment(e, c3t3);

        if (vh0_moving)
        {
          moves[i0].move += density * Vector_3(p0, p1);
          moves[i0].mass += density;
          ++moves[i0].neighbors;
        }
        if (vh1_moving)
        {
          moves[i1].move += density * Vector_3(p1, p0);
          moves[i1].mass += density;
          ++moves[i1].neighbors;
        }
      }
    }
  }

  ElementSource get_element_source(const C3t3& c3t3) const override
  {
    perform_global_preprocessing(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }

  bool execute_operation(const ElementType& v, C3t3& c3t3) override
  {
    auto& tr = c3t3.triangulation();
    auto& moves = m_context->m_moves;

    const std::size_t vid = m_context->vertex_id(v);
    if (!m_context->is_free(vid))
      return false;

    if (c3t3.in_dimension(v) == 3 && moves[vid].neighbors > 1)
    {
      const Vector_3 move = moves[vid].move / moves[vid].mass;
      const Point_3 new_pos = point(v->point()) + move;
      return BaseClass::check_inversion_and_move(v, new_pos, m_context->incident_cells(vid), tr,
                                                 m_context->m_total_move);
    }
    return false;
  }

  std::string operation_name() const override { return "Vertex Smooth (Internal Vertices)"; }
};

}//namespace internal
}//namespace Tetrahedral_adaptive_remeshing
}//namespace CGAL

#endif //CGAL_INTERNAL_SMOOTH_VERTICES_H
