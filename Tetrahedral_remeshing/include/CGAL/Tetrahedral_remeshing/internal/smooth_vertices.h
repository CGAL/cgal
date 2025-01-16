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
template<typename C3t3, typename SizingFunction, typename CellSelector>
class Tetrahedral_remeshing_smoother
{
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
  typedef  CGAL::Tetrahedral_remeshing::internal::FMLS<Gt> FMLS;
  std::vector<FMLS> subdomain_FMLS;
  std::unordered_map<Surface_patch_index, std::size_t, boost::hash<Surface_patch_index>> subdomain_FMLS_indices;

  Triangle_vec m_aabb_triangles;
  AABB_triangle_tree m_triangles_aabb_tree;
  Segment_vec m_aabb_segments;
  AABB_segment_tree m_segments_aabb_tree;
  FT m_aabb_epsilon;

  const SizingFunction& m_sizing;
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
  Tetrahedral_remeshing_smoother(const SizingFunction& sizing,
                                 const CellSelector& cell_selector,
                                 const bool protect_boundaries,
                                 const bool smooth_constrained_edges)
    : m_sizing(sizing)
    , m_cell_selector(cell_selector)
    , m_protect_boundaries(protect_boundaries)
    , m_smooth_constrained_edges(smooth_constrained_edges)
  {}

  void init(const C3t3& c3t3)
  {
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    //collect a map of vertices surface indices
    std::unordered_map<Vertex_handle, std::vector<Surface_patch_index> > vertices_surface_indices;
    collect_vertices_surface_indices(c3t3, vertices_surface_indices);

    //collect a map of normals at surface vertices
    std::unordered_map<Vertex_handle,
          std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> vertices_normals;
    compute_vertices_normals(c3t3, vertices_normals);

    // Build MLS Surfaces
    createMLSSurfaces(subdomain_FMLS,
                      subdomain_FMLS_indices,
                      vertices_normals,
                      vertices_surface_indices,
                      c3t3);
#else
    // Build AABB tree
    build_aabb_trees(c3t3);
#endif
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

private:

  bool is_selected(const Cell_handle c) const
  {
    return get(m_cell_selector, c);
  }

   std::size_t vertex_id(const Vertex_handle v) const
  {
    CGAL_expensive_assertion(m_vertex_id.find(v) != m_vertex_id.end());
    return m_vertex_id.at(v);
  }

  bool is_free(const Vertex_handle v) const
  {
    return m_free_vertices[vertex_id(v)];
  }
  bool is_free(const std::size_t & vid) const
  {
    return m_free_vertices[vid];
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

  template<typename IncCellsVector>
  void collect_incident_cells(const Tr& tr,
                              IncCellsVector& inc_cells)
  {
    for (const Cell_handle c : tr.finite_cell_handles())
    {
      for (auto vi : tr.vertices(c))
      {
        const std::size_t idi = vertex_id(vi);
        if(is_free(idi))
          inc_cells[idi].push_back(c);
      }
    }
  }

  Point_3 project_on_tangent_plane(const Point_3& gi,
                                    const Point_3& pi,
                                    const Vector_3& normal)
  {
    Vector_3 diff(gi, pi);
    return gi + (normal * diff) * normal;
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

  template<typename VertexNormalsMap>
  void compute_vertices_normals(const C3t3& c3t3,
                                VertexNormalsMap& normals_map)
  {
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
        typename VertexNormalsMap::iterator patch_vector_it = normals_map.find(vi);

        if (patch_vector_it == normals_map.end()
            || patch_vector_it->second.find(surf_i) == patch_vector_it->second.end())
        {
          normals_map[vi][surf_i] = n;
        }
        else
        {
          normals_map[vi][surf_i] += n;
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
    for (typename VertexNormalsMap::iterator vnm_it = normals_map.begin();
         vnm_it != normals_map.end(); ++vnm_it)
    {
      //value type is map<Surface_patch_index, Vector_3>
      for (typename VertexNormalsMap::mapped_type::iterator it = vnm_it->second.begin();
           it != vnm_it->second.end(); ++it)
      {
        Vector_3& n = it->second;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        auto p = point(vnm_it->first->point());
        os << "2 " << p << " " << (p + n) << std::endl;
#endif

        CGAL::Tetrahedral_remeshing::normalize(n, gt);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        const Surface_patch_index si = it->first;
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

  std::optional<Point_3> project(const Surface_patch_index& si,
                                 const Point_3& gi)
  {
    CGAL_expensive_assertion(subdomain_FMLS_indices.find(si) != subdomain_FMLS_indices.end());
    CGAL_assertion(!std::isnan(gi.x()) && !std::isnan(gi.y()) && !std::isnan(gi.z()));

    Vector_3 point(gi.x(), gi.y(), gi.z());
    Vector_3 res_normal = CGAL::NULL_VECTOR;
    Vector_3 result(CGAL::ORIGIN, gi);

    const FMLS& fmls = subdomain_FMLS[subdomain_FMLS_indices.at(si)];

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

  Dihedral_angle_cosine max_cosine(const Tr& tr,
    const boost::container::small_vector<Cell_handle, 40>& cells)
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
                                FT& total_move)
#else
                                FT&)
#endif
  {
    const typename Tr::Point backup = v->point(); //backup v's position
    const typename Tr::Geom_traits::Point_3 pv = point(backup);

    bool valid_orientation = false;
    bool angles_improved = true;
    double frac = 1.0;
    typename Tr::Geom_traits::Vector_3 move(pv, final_pos);

    const Dihedral_angle_cosine curr_max_cos = m_flip_smooth_steps
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
        else if (m_flip_smooth_steps) //check that dihedral angles get improved
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
    bool valid_move =  valid_orientation && angles_improved;

    if(!valid_move)
      v->set_point(backup);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    else
      total_move += CGAL::approximate_sqrt(CGAL::squared_distance(pv, point(v->point())));
#endif

    return valid_move;
  }

  void collect_vertices_surface_indices(
    const C3t3& c3t3,
    std::unordered_map<Vertex_handle,
    std::vector<Surface_patch_index> >& vertices_surface_indices)
  {
    for (Facet fit : c3t3.facets_in_complex())
    {
      const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);

      for (const Vertex_handle vi : c3t3.triangulation().vertices(fit))
      {
        std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[vi];
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


  // boundary_edge is set to false by default because
  // we may not care about this information, for example while collecting
  // weights in the smoothing inside volume step
  FT density_along_segment(const Edge& e,
                           const C3t3& c3t3,
                           const bool boundary_edge = false) const
  {
    const auto mwi = midpoint_with_info(e, boundary_edge, c3t3);
    const FT s = sizing_at_midpoint(e, mwi.dim, mwi.index, m_sizing, c3t3, m_cell_selector);
    const FT density = 1. / s; //density = 1 / size^(dimension)
                 //edge dimension is 1, so density = 1 / size
                 //to have mass = length * density with no dimension
    return density;
  }

  template<typename SurfaceIndices,
           typename IncidentCells, typename NormalsMap>
  std::size_t smooth_edges_in_complex(C3t3& c3t3,
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
                                      const SurfaceIndices& vertices_surface_indices,
#else
                                      const SurfaceIndices&,
#endif
                                      const IncidentCells& inc_cells,
                                      const NormalsMap& vertices_normals,
                                      FT& total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                                    , std::ofstream& os_surf
#endif
                               )
  {
    std::size_t nb_done_1d = 0;
    auto& tr = c3t3.triangulation();

    const std::size_t nbv = tr.number_of_vertices();
    std::vector<Vector_3> moves(nbv, CGAL::NULL_VECTOR);
    std::vector<int> neighbors(nbv, 0);
    std::vector<FT> masses(nbv, 0.);

    //collect neighbors
    for (const Edge& e : c3t3.edges_in_complex())
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      CGAL_expensive_assertion(is_on_feature(vh0));
      CGAL_expensive_assertion(is_on_feature(vh1));

      const std::size_t& i0 = vertex_id(vh0);
      const std::size_t& i1 = vertex_id(vh1);

      const bool vh0_moving = is_free(i0);
      const bool vh1_moving = is_free(i1);

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = density_along_segment(e, c3t3, true);

      if (vh0_moving)
      {
        moves[i0] += density * Vector_3(p0, p1);
        neighbors[i0]++;
        masses[i0] += density;
      }
      if (vh1_moving)
      {
        moves[i1] += density * Vector_3(p1, p0);
        neighbors[i1]++;
        masses[i1] += density;
      }
    }

    // iterate over map of <vertex, id>
    for (auto [v, vid] : m_vertex_id)
    {
      if (!is_free(vid) || !is_on_feature(v))
        continue;

      const Point_3 current_pos = point(v->point());

      const std::size_t nb_neighbors = neighbors[vid];
      if(nb_neighbors == 0)
        continue;

      CGAL_assertion(masses[vid] > 0);
      const Vector_3 move = (nb_neighbors > 0)
                          ? moves[vid] / masses[vid]
                          : CGAL::NULL_VECTOR;

      const Point_3 smoothed_position = current_pos + move;

      CGAL_USE(vertices_normals);
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      Vector_3 sum_projections = CGAL::NULL_VECTOR;
      Point_3 tmp_pos = current_pos;

#ifndef CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION
      const std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices.at(v);
      for (const Surface_patch_index& si : v_surface_indices)
      {
        Point_3 normal_projection = project_on_tangent_plane(smoothed_position,
                                                             current_pos,
                                                             vertices_normals.at(v).at(si));

        sum_projections += Vector_3(tmp_pos, normal_projection);
        tmp_pos = normal_projection;
      }
#endif //CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION

      const Point_3 new_pos = current_pos + sum_projections;

#else // AABB_tree projection

      const Point_3 new_pos = m_segments_aabb_tree.closest_point(smoothed_position);


#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      os_surf << "2 " << current_pos << " " << new_pos << std::endl;
#endif
      // move vertex
      if (check_inversion_and_move(v, new_pos, inc_cells[vid], tr, total_move)){
        nb_done_1d++;
      }
    }
    return nb_done_1d;
  }


template<typename SurfaceIndices,
         typename IncidentCells, typename NormalsMap>
std::size_t smooth_vertices_on_surfaces(C3t3& c3t3,
                                        const SurfaceIndices& vertices_surface_indices,
                                        const IncidentCells& inc_cells,
                                        const NormalsMap& vertices_normals,
                                        FT& total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                                      , std::ofstream& os_surf
                                      , std::ofstream& os_surf0
#endif
                               )
{
  std::size_t nb_done_2d = 0;
  auto& tr = c3t3.triangulation();

  const std::size_t nbv = tr.number_of_vertices();
  std::vector<Vector_3> moves(nbv, CGAL::NULL_VECTOR);
  std::vector<int> neighbors(nbv, 0);
  std::vector<FT> masses(nbv, 0.);

  for (const Edge& e : tr.finite_edges())
  {
    if (!c3t3.is_in_complex(e) && is_boundary(c3t3, e, m_cell_selector))
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const std::size_t& i0 = vertex_id(vh0);
      const std::size_t& i1 = vertex_id(vh1);

      const bool vh0_moving = !is_on_feature(vh0) && is_free(i0);
      const bool vh1_moving = !is_on_feature(vh1) && is_free(i1);

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = density_along_segment(e, c3t3, true);

      if (vh0_moving)
      {
        moves[i0] += density * Vector_3(p0, p1);
        neighbors[i0]++;
        masses[i0] += density;
      }
      if (vh1_moving)
      {
        moves[i1] += density * Vector_3(p1, p0);
        neighbors[i1]++;
        masses[i1] += density;
      }
    }
  }

  // iterate over map of <vertex, id>
  for (auto [v, vid] : m_vertex_id)
  {
    if (!is_free(vid) || v->in_dimension() != 2)
      continue;

    const std::size_t nb_neighbors = neighbors[vid];
    const Point_3 current_pos = point(v->point());

    const auto& incident_surface_patches = vertices_surface_indices.at(v);
    if (incident_surface_patches.size() > 1)
      continue;
    const Surface_patch_index si = incident_surface_patches[0];

    CGAL_assertion(si != Surface_patch_index());
    CGAL_expensive_assertion_code(auto siv = surface_patch_index(v, c3t3));
    CGAL_expensive_assertion(si == siv);

    if (nb_neighbors > 1)
    {
      const Vector_3 move = moves[vid] / masses[vid];
      const Point_3 smoothed_position = point(v->point()) + move;

#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      Point_3 normal_projection = project_on_tangent_plane(smoothed_position,
                                                           current_pos,
                                                           vertices_normals.at(v).at(si));
      std::optional<Point_3> mls_projection = project(si, normal_projection);

      const Point_3 new_pos = (mls_projection != std::nullopt)
                            ? *mls_projection
                            : smoothed_position;

#else // AABB_tree projection
      Point_3 new_pos;
      if (m_triangles_aabb_tree.squared_distance(smoothed_position) < m_aabb_epsilon)
      {
        new_pos = m_triangles_aabb_tree.closest_point(smoothed_position);
      }
      else
      {
        using Ray = typename Tr::Geom_traits::Ray_3;
        using Projection = std::optional<
          typename AABB_triangle_tree::template Intersection_and_primitive_id<Ray>::Type>;

        auto get_intersection_point =
          [](const Projection& proj) -> std::optional<Point_3>
          {
            const auto intersection = proj.value().first;
            if (const Point_3* p = std::get_if<Point_3>(&intersection))
              return *p;
            else
              return std::nullopt;
          };

        // this lambda is called only when we are sure that proj is a Segment
        auto get_intersection_midpoint =
          [](const Projection& proj) -> std::optional<Point_3>
          {
            const auto intersection = proj.value().first;
            using Segment = typename Tr::Geom_traits::Segment_3;
            if (const Segment* s = std::get_if<Segment>(&intersection))
              return CGAL::midpoint(s->source(), s->target());
            else
            {
              CGAL_assertion(false);
              return std::nullopt;
            }
          };

        const auto n = vertices_normals.at(v).at(si);
        const Ray ray = tr.geom_traits().construct_ray_3_object()(current_pos, n);

        const Projection proj = m_triangles_aabb_tree.first_intersection(ray);
        const Projection proj_opp = m_triangles_aabb_tree.first_intersection(
          tr.geom_traits().construct_opposite_ray_3_object()(ray));

        if(proj != std::nullopt && proj_opp == std::nullopt)
        {
          const auto p = get_intersection_point(proj);
          if (p != std::nullopt)
            new_pos = p.value();
          else
            new_pos = get_intersection_midpoint(proj).value();
        }
        else if(proj == std::nullopt && proj_opp != std::nullopt)
        {
          const auto p = get_intersection_point(proj_opp);
          if (p != std::nullopt)
            new_pos = p.value();
          else
            new_pos = get_intersection_midpoint(proj_opp).value();
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
          new_pos = smoothed_position;
      }
#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      if (check_inversion_and_move(v, new_pos, inc_cells[vid], tr, total_move)){
        nb_done_2d++;
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      os_surf << "2 " << current_pos << " " << new_pos << std::endl;
#endif
    }
    else if (nb_neighbors > 0)
    {
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      std::optional<Point_3> mls_proj = project(si, current_pos);
      if (mls_proj == std::nullopt)
        continue;

      const Point_3 new_pos = *mls_proj;
#else // AABB_tree projection
      const Point_3 new_pos = m_segments_aabb_tree.closest_point(current_pos);
#endif // CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      if (check_inversion_and_move(v, new_pos, inc_cells[vid], tr, total_move)){
        nb_done_2d++;
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        os_surf0 << "2 " << current_pos << " " << new_pos << std::endl;
#endif
    }
  }

  return nb_done_2d;
}

template<typename IncidentCells>
std::size_t smooth_internal_vertices(C3t3& c3t3,
                                     const IncidentCells& inc_cells,
                                     FT& total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                                   , std::ofstream& os_vol
#endif
                                     )
{
  std::size_t nb_done_3d = 0;
  auto& tr = c3t3.triangulation();

  const std::size_t nbv = tr.number_of_vertices();
  std::vector<Vector_3> moves(nbv, CGAL::NULL_VECTOR);
  std::vector<int> neighbors(nbv, 0);/*for dim 3 vertices, start counting directly from 0*/
  std::vector<FT> masses(nbv, 0.);

  for (const Edge& e : tr.finite_edges())
  {
    if (is_outside(e, c3t3, m_cell_selector))
      continue;
    else
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const std::size_t& i0 = vertex_id(vh0);
      const std::size_t& i1 = vertex_id(vh1);

      const bool vh0_moving = (c3t3.in_dimension(vh0) == 3 && is_free(i0));
      const bool vh1_moving = (c3t3.in_dimension(vh1) == 3 && is_free(i1));

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = density_along_segment(e, c3t3);

      if (vh0_moving)
      {
        moves[i0] += density * Vector_3(p0, p1);
        neighbors[i0]++;
        masses[i0] += density;
      }
      if (vh1_moving)
      {
        moves[i1] += density * Vector_3(p1, p0);
        neighbors[i1]++;
        masses[i1] += density;
      }
    }
  }

  // iterate over map of <vertex, id>
  for (auto [v, vid] : m_vertex_id)
  {
    if (!is_free(vid))
      continue;

    if (c3t3.in_dimension(v) == 3 && neighbors[vid] > 1)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      os_vol << "2 " << point(v->point());
#endif
      const Vector_3 move = moves[vid] / masses[vid];// static_cast<FT>(neighbors[vid]);
      Point_3 new_pos = point(v->point()) + move;
      if (check_inversion_and_move(v, new_pos, inc_cells[vid], tr, total_move)){
        nb_done_3d++;
      }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      os_vol << " " << point(v->point()) << std::endl;
#endif
    }
  }
  return nb_done_3d;
}

public:
  void smooth_vertices(C3t3& c3t3)
  {
    typedef typename C3t3::Cell_handle            Cell_handle;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::ofstream os_surf("smooth_surfaces.polylines.txt");
    std::ofstream os_surf0("smooth_surfaces0.polylines.txt");
    std::ofstream os_vol("smooth_volume.polylines.txt");
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Smooth vertices...";
    std::cout.flush();

    std::size_t nb_done_3d = 0;
    std::size_t nb_done_2d = 0;
    std::size_t nb_done_1d = 0;
    CGAL::Real_timer timer;
    timer.start();
#endif

    FT total_move = 0.;

    Tr& tr = c3t3.triangulation();

    //collect a map of vertices surface indices
    std::unordered_map<Vertex_handle, std::vector<Surface_patch_index> > vertices_surface_indices;
    if(!m_protect_boundaries)
      collect_vertices_surface_indices(c3t3, vertices_surface_indices);

    //collect a map of normals at surface vertices
    std::unordered_map<Vertex_handle,
          std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> vertices_normals;
    if(!m_protect_boundaries)
      compute_vertices_normals(c3t3, vertices_normals);

    //collect ids
    reset_vertex_id_map(tr);

    //are vertices free to move? indices are in `vertex_id`
    reset_free_vertices(tr);

    //collect incident cells
    using Incident_cells_vector = boost::container::small_vector<Cell_handle, 40>;
    const std::size_t nbv = tr.number_of_vertices();
    std::vector<Incident_cells_vector> inc_cells(nbv, Incident_cells_vector{});
    collect_incident_cells(tr, inc_cells);

    if (!m_protect_boundaries && m_smooth_constrained_edges)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      nb_done_1d =
#endif
      smooth_edges_in_complex(c3t3,
                              vertices_surface_indices, inc_cells, vertices_normals, total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                              , os_surf
#endif
      );
    }

    /////////////// EDGES ON SURFACE, BUT NOT IN COMPLEX //////////////////
    if (!m_protect_boundaries)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      nb_done_2d =
#endif
       smooth_vertices_on_surfaces(c3t3,
                                  vertices_surface_indices, inc_cells, vertices_normals,
                                  total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                                , os_surf, os_surf0
#endif
      );
    }
    CGAL_expensive_assertion(CGAL::Tetrahedral_remeshing::debug::are_cell_orientations_valid(tr));
    ////   end if(!protect_boundaries)

    ////////////// INTERNAL VERTICES ///////////////////////
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    nb_done_3d =
#endif
    smooth_internal_vertices(c3t3, inc_cells,
                             total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                           , os_vol
#endif
    );

    CGAL_expensive_assertion(CGAL::Tetrahedral_remeshing::debug::are_cell_orientations_valid(tr));

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    timer.stop();
    std::size_t nb_done = nb_done_3d + nb_done_2d + nb_done_1d;
    std::cout << " done ("
      << nb_done_1d << "/" << nb_done_2d << "/" << nb_done_3d << " vertices smoothed,"
      << " average move = " << (total_move / nb_done)
      << ", in "<< timer.time() << " seconds)." << std::endl;
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      c3t3.triangulation(), "c3t3_vertices_after_smoothing");
    os_surf.close();
    os_vol.close();
#endif
  }

};//end class Tetrahedral_remeshing_smoother
}//namespace internal
}//namespace Tetrahedral_adaptive_remeshing
}//namespace CGAL

#endif //CGAL_INTERNAL_SMOOTH_VERTICES_H
