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
template<typename C3t3, typename SizingFunction>
class Tetrahedral_remeshing_smoother
{
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename Tr::Edge                  Edge;
  typedef typename Tr::Facet                 Facet;

  typedef typename Tr::Geom_traits           Gt;
  typedef typename Gt::Vector_3              Vector_3;
  typedef typename Gt::Point_3               Point_3;
  typedef typename Gt::FT                    FT;

private:
  typedef  CGAL::Tetrahedral_remeshing::internal::FMLS<Gt> FMLS;
  std::vector<FMLS> subdomain_FMLS;
  std::unordered_map<Surface_patch_index, std::size_t, boost::hash<Surface_patch_index>> subdomain_FMLS_indices;
  bool m_smooth_constrained_edges;
  const SizingFunction& m_sizing;

public:
  Tetrahedral_remeshing_smoother(const SizingFunction& sizing)
    : m_sizing(sizing)
  {}

  template<typename CellSelector>
  void init(const C3t3& c3t3,
            const CellSelector& cell_selector,
            const bool smooth_constrained_edges)
  {
    //collect a map of vertices surface indices
    std::unordered_map<Vertex_handle, std::vector<Surface_patch_index> > vertices_surface_indices;
    collect_vertices_surface_indices(c3t3, vertices_surface_indices);

    //collect a map of normals at surface vertices
    std::unordered_map<Vertex_handle,
          std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> vertices_normals;
    compute_vertices_normals(c3t3, vertices_normals, cell_selector);

    // Build MLS Surfaces
    createMLSSurfaces(subdomain_FMLS,
                      subdomain_FMLS_indices,
                      vertices_normals,
                      vertices_surface_indices,
                      c3t3);

    m_smooth_constrained_edges = smooth_constrained_edges;
  }

private:

  Point_3 project_on_tangent_plane(const Point_3& gi,
                                    const Point_3& pi,
                                    const Vector_3& normal)
  {
    Vector_3 diff(gi, pi);
    return gi + (normal * diff) * normal;
  }

  template<typename CellSelector>
  std::optional<Facet>
  find_adjacent_facet_on_surface(const Facet& f,
                                 const Edge& edge,
                                 const C3t3& c3t3,
                                 const CellSelector& cell_selector)
  {
    CGAL_assertion(is_boundary(c3t3, f, cell_selector));

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
          && is_boundary(c3t3, fi, cell_selector)
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

  template<typename VertexNormalsMap, typename CellSelector>
  void compute_vertices_normals(const C3t3& c3t3,
                                VertexNormalsMap& normals_map,
                                const CellSelector& cell_selector)
  {
    typename Tr::Geom_traits gt = c3t3.triangulation().geom_traits();
    typename Tr::Geom_traits::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();

    const Tr& tr = c3t3.triangulation();

    //collect all facet normals
    std::unordered_map<Facet, Vector_3, boost::hash<Facet>> fnormals;
    for (const Facet& f : tr.finite_facets())
    {
      if (is_boundary(c3t3, f, cell_selector))
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
      CGAL_assertion(is_boundary(c3t3, f, cell_selector));

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

        const typename C3t3::Cell_handle ch = f.first;
        const std::array<std::array<int, 2>, 3> edges
          = {{ {{(ff.second + 1) % 4, (ff.second + 2) % 4}}, //edge 1-2
               {{(ff.second + 2) % 4, (ff.second + 3) % 4}}, //edge 2-3
               {{(ff.second + 3) % 4, (ff.second + 1) % 4}}  //edge 3-1
            }}; //vertex indices in cells

        const Vector_3& ref = fnormals[f];
        for (const std::array<int, 2>& ei : edges)
        {
          Edge edge(ch, ei[0], ei[1]);
          if (std::optional<Facet> neighbor
              = find_adjacent_facet_on_surface(f, edge, c3t3, cell_selector))
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
    for (const auto& fn : fnormals)
    {
      const Facet& f = fn.first;
      const Vector_3& n = fn.second;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      typename Tr::Geom_traits::Point_3 fc
        = CGAL::centroid(point(f.first->vertex(indices(f.second, 0))->point()),
                         point(f.first->vertex(indices(f.second, 1))->point()),
                         point(f.first->vertex(indices(f.second, 2))->point()));
      osf << "2 " << fc << " " << (fc + n) << std::endl;
#endif
      const Surface_patch_index& surf_i = c3t3.surface_patch_index(f);

      for (int i = 0; i < 3; ++i)
      {
        const Vertex_handle vi = f.first->vertex(indices(f.second, i));
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

        CGAL::Tetrahedral_remeshing::normalize(n, c3t3.triangulation().geom_traits());

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
    CGAL_assertion(subdomain_FMLS_indices.find(si) != subdomain_FMLS_indices.end());
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

  template<typename CellRange, typename Tr>
  bool check_inversion_and_move(const typename Tr::Vertex_handle v,
                                const typename Tr::Geom_traits::Point_3& final_pos,
                                const CellRange& inc_cells,
                                const Tr& /* tr */,
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
                                FT& total_move)
#else
                                FT&)
#endif
  {
    const typename Tr::Point backup = v->point(); //backup v's position
    const typename Tr::Geom_traits::Point_3 pv = point(backup);

    bool valid_orientation = false;
    double frac = 1.0;
    typename Tr::Geom_traits::Vector_3 move(pv, final_pos);

    do
    {
      v->set_point(typename Tr::Point(pv + frac * move));

      bool valid_try = true;
      for (const typename Tr::Cell_handle& ci : inc_cells)
      {
        if (CGAL::POSITIVE != CGAL::orientation(point(ci->vertex(0)->point()),
                                                point(ci->vertex(1)->point()),
                                                point(ci->vertex(2)->point()),
                                                point(ci->vertex(3)->point())))
        {
          frac = 0.9 * frac;
          valid_try = false;
          break;
        }
      }
      valid_orientation = valid_try;
    }
    while(!valid_orientation && frac > 0.1);

    if (!valid_orientation) //move failed
      v->set_point(backup);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    else
      total_move += CGAL::approximate_sqrt(CGAL::squared_distance(pv, point(v->point())));
#endif

    return valid_orientation;
  }

  void collect_vertices_surface_indices(
    const C3t3& c3t3,
    std::unordered_map<Vertex_handle,
    std::vector<Surface_patch_index> >& vertices_surface_indices)
  {
    for (Facet fit : c3t3.facets_in_complex())
    {
      const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);

      for (int i = 0; i < 3; i++)
      {
        const Vertex_handle vi = fit.first->vertex(indices(fit.second, i));

        std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[vi];
        if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index) == v_surface_indices.end())
          v_surface_indices.push_back(surface_index);
      }
    }
  }

  auto vertex_id_map(const Tr& tr)
  {
    using Vertex_handle = typename Tr::Vertex_handle;
    std::unordered_map<Vertex_handle, std::size_t> vertex_id;
    std::size_t id = 0;
    for (const Vertex_handle v : tr.finite_vertex_handles())
    {
      vertex_id[v] = id++;
    }
    return vertex_id;
  }

  typename C3t3::Index max_dimension_index(const std::array<Vertex_handle, 2> vs) const
  {
    const int dim0 = vs[0]->in_dimension();
    const int dim1 = vs[1]->in_dimension();

    if(dim0 > dim1)       return vs[0]->index();
    else if(dim1 > dim0)  return vs[1]->index();
    else                  return vs[0]->index(); //arbitrary choice, any of the two should be fine
  }

  FT mass_along_segment(const Edge& e, const Tr& tr) const
  {
    typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    const std::array<Vertex_handle, 2> vs = tr.vertices(e);

    const FT len = CGAL::approximate_sqrt(squared_edge_length(e, tr));

    const int dim = (std::max)(vs[0]->in_dimension(), vs[1]->in_dimension());
    const typename C3t3::Index index = max_dimension_index(vs);

    const FT s = m_sizing(CGAL::midpoint(cp(vs[0]->point()), cp(vs[1]->point())), dim, index);
    const FT density = 1. / (s * s * s); //density = 1 / size^(dimension + 2)
                                         //edge dimension is 1, so density = 1 / size^3
    //const FT mass = len * density;

    return density;// mass;
  }

  template<typename VertexIdMap,
           typename VertexBoolMap, typename SurfaceIndices,
           typename IncidentCells, typename NormalsMap>
  std::size_t smooth_edges_in_complex(C3t3& c3t3,
                                      const VertexIdMap& vertex_id,
                                      const VertexBoolMap& free_vertex,
                                      const SurfaceIndices& vertices_surface_indices,
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
    std::vector<Vector_3> smoothed_positions(nbv, CGAL::NULL_VECTOR);
    std::vector<int> neighbors(nbv, 0);
    std::vector<FT> masses(nbv, 0.);

    //collect neighbors
    for (const Edge& e : c3t3.edges_in_complex())
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const std::size_t& i0 = vertex_id.at(vh0);
      const std::size_t& i1 = vertex_id.at(vh1);

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());

      const auto mass = mass_along_segment(e, tr);

      if (!c3t3.is_in_complex(vh0) && is_on_feature(vh1) && free_vertex[i0])
      {
        smoothed_positions[i0] = smoothed_positions[i0] + mass * Vector_3(p0, p1);
        neighbors[i0]++;
        masses[i0] += mass;
      }
      if (!c3t3.is_in_complex(vh1) && is_on_feature(vh0) && free_vertex[i1])
      {
        smoothed_positions[i1] = smoothed_positions[i1] + mass * Vector_3(p1, p0);
        neighbors[i1]++;
        masses[i1] += mass;
      }
    }

    // iterate over map of <vertex, id>
    for (auto [v, vid] : vertex_id)
    {
      if (!free_vertex[vid])
        continue;
      if(!is_on_feature(v))
        continue;

      const Point_3 current_pos = point(v->point());
      const std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices.at(v);
      const std::size_t count = v_surface_indices.size();

      const std::size_t nb_neighbors = neighbors[vid];
      if(nb_neighbors == 0)
        continue;

      CGAL_assertion(masses[vid] > 0);
      const Vector_3 move = (nb_neighbors > 1)
                          ? smoothed_positions[vid] / masses[vid]//nb_neighbors)
                          : CGAL::NULL_VECTOR;

      const Point_3 smoothed_position = point(v->point()) + move;

      std::size_t nb_projected = 0;
      Vector_3 sum_projections = CGAL::NULL_VECTOR;

      for (const Surface_patch_index& si : v_surface_indices)
      {
        const Point_3 projection = (nb_neighbors > 1)
          ? project_on_tangent_plane(smoothed_position, current_pos, vertices_normals.at(v).at(si))
          : smoothed_position;

        //Check if the mls surface exists to avoid degenerated cases
        std::optional<Point_3> mls_projection = project(si, projection);

        if (mls_projection != std::nullopt)
        {
          sum_projections = sum_projections + Vector_3(current_pos, *mls_projection);
          ++nb_projected;
        }
      }

      Point_3 new_pos = (nb_projected > 0)
                      ? smoothed_position + sum_projections / static_cast<FT>(nb_projected)
                      : smoothed_position;

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


template<typename VertexIdMap,
         typename VertexBoolMap, typename SurfaceIndices,
         typename IncidentCells, typename NormalsMap, typename CellSelector>
std::size_t smooth_vertices_on_surfaces(C3t3& c3t3,
                                        const VertexIdMap& vertex_id,
                                        const VertexBoolMap& free_vertex,
                                        const SurfaceIndices& vertices_surface_indices,
                                        const IncidentCells& inc_cells,
                                        const NormalsMap& vertices_normals,
                                        const CellSelector& cell_selector,
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
  std::vector<Vector_3> smoothed_positions(nbv, CGAL::NULL_VECTOR);
  std::vector<int> neighbors(nbv, 0);
  std::vector<FT> masses(nbv, 0.);

  for (const Edge& e : tr.finite_edges())
  {
    if (!c3t3.is_in_complex(e) && is_boundary(c3t3, e, cell_selector))
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());

      const std::size_t& i0 = vertex_id.at(vh0);
      const std::size_t& i1 = vertex_id.at(vh1);

      const auto mass = mass_along_segment(e, tr);

      if (!is_on_feature(vh0) && free_vertex[i0])
      {
        smoothed_positions[i0] = smoothed_positions[i0] + mass * Vector_3(p0, p1);
        neighbors[i0]++;
        masses[i0] += mass;
      }
      if (!is_on_feature(vh1) && free_vertex[i1])
      {
        smoothed_positions[i1] = smoothed_positions[i1] + mass * Vector_3(p1, p0);
        neighbors[i1]++;
        masses[i1] += mass;
      }
    }
  }

  // iterate over map of <vertex, id>
  for (auto [v, vid] : vertex_id)
  {
    if (!free_vertex[vid] || v->in_dimension() != 2)
      continue;

    const std::size_t nb_neighbors = neighbors[vid];
    const Point_3 current_pos = point(v->point());

    const Surface_patch_index si = surface_patch_index(v, c3t3);
    CGAL_assertion(si != Surface_patch_index());

    if (nb_neighbors > 1)
    {
      const Vector_3 move = smoothed_positions[vid] / masses[vid];
      const Point_3 smoothed_position = point(v->point()) + move;
      Point_3 normal_projection = project_on_tangent_plane(smoothed_position,
                                                           current_pos,
                                                           vertices_normals.at(v).at(si));
      std::optional<Point_3> mls_projection = project(si, normal_projection);

      const Point_3 new_pos = (mls_projection != std::nullopt)
                            ? *mls_projection
                            : smoothed_position;

      if (check_inversion_and_move(v, new_pos, inc_cells[vid], tr, total_move)){
        nb_done_2d++;
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      os_surf << "2 " << current_pos << " " << new_pos << std::endl;
#endif
    }
    else if (nb_neighbors > 0)
    {
      std::optional<Point_3> mls_projection = project(si, current_pos);
      if (mls_projection != std::nullopt)
      {
        const Point_3 new_pos = *mls_projection;
        if (check_inversion_and_move(v, new_pos, inc_cells[vid], tr, total_move)){
          nb_done_2d++;
        }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        os_surf0 << "2 " << current_pos << " " << new_pos << std::endl;
#endif
      }
    }
  }
  return nb_done_2d;
}

template<typename VertexIdMap,
         typename VertexBoolMap, typename SurfaceIndices,
         typename IncidentCells, typename NormalsMap, typename CellSelector>
std::size_t smooth_internal_vertices(C3t3& c3t3,
                                     const VertexIdMap& vertex_id,
                                     const VertexBoolMap& free_vertex,
                                     const SurfaceIndices& vertices_surface_indices,
                                     const IncidentCells& inc_cells,
                                     const NormalsMap& vertices_normals,
                                     const CellSelector& cell_selector,
                                     FT& total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                                   , std::ofstream& os_vol
#endif
                                     )
{
  std::size_t nb_done_3d = 0;
  auto& tr = c3t3.triangulation();

  const std::size_t nbv = tr.number_of_vertices();
  std::vector<Vector_3> smoothed_positions(nbv, CGAL::NULL_VECTOR);
  std::vector<int> neighbors(nbv, 0);/*for dim 3 vertices, start counting directly from 0*/
  std::vector<FT> masses(nbv, 0.);

  for (const Edge& e : tr.finite_edges())
  {
    if (is_outside(e, c3t3, cell_selector))
      continue;
    else
    {
      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const std::size_t& i0 = vertex_id.at(vh0);
      const std::size_t& i1 = vertex_id.at(vh1);

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());

      const auto mass = mass_along_segment(e, tr);

      if (c3t3.in_dimension(vh0) == 3 && free_vertex[i0])
      {
        smoothed_positions[i0] = smoothed_positions[i0] + mass * Vector_3(p0, p1);
        neighbors[i0]++;
        masses[i0] += mass;
      }
      if (c3t3.in_dimension(vh1) == 3 && free_vertex[i1])
      {
        smoothed_positions[i1] = smoothed_positions[i1] + mass * Vector_3(p1, p0);
        neighbors[i1]++;
        masses[i1] += mass;
      }
    }
  }

  // iterate over map of <vertex, id>
  for (auto [v, vid] : vertex_id)
  {
    if (!free_vertex[vid])
      continue;

    if (c3t3.in_dimension(v) == 3 && neighbors[vid] > 1)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      os_vol << "2 " << point(v->point());
#endif
      const Vector_3 move = smoothed_positions[vid] / masses[vid];// static_cast<FT>(neighbors[vid]);
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
  template<typename CellSelector>
  void smooth_vertices(C3t3& c3t3,
                       const bool protect_boundaries,
                       const CellSelector& cell_selector)
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
#endif
    FT total_move = 0.;

    Tr& tr = c3t3.triangulation();

    //collect a map of vertices surface indices
    std::unordered_map<Vertex_handle, std::vector<Surface_patch_index> > vertices_surface_indices;
    if(m_smooth_constrained_edges)
      collect_vertices_surface_indices(c3t3, vertices_surface_indices);

    //collect a map of normals at surface vertices
    std::unordered_map<Vertex_handle,
          std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>> vertices_normals;
    compute_vertices_normals(c3t3, vertices_normals, cell_selector);

    //collect ids
    const std::unordered_map<Vertex_handle, std::size_t>
      vertex_id = vertex_id_map(tr);

    //smooth()
    const std::size_t nbv = tr.number_of_vertices();
    std::vector<bool> free_vertex(nbv, false);//are vertices free to move? indices are in `vertex_id`

    //collect incident cells
    using Incident_cells_vector = boost::container::small_vector<Cell_handle, 40>;
    std::vector<Incident_cells_vector> inc_cells(nbv, Incident_cells_vector());

    for (const Cell_handle c : tr.finite_cell_handles())
    {
      const bool cell_is_selected = get(cell_selector, c);

      for (int i = 0; i < 4; ++i)
      {
        const std::size_t idi = vertex_id.at(c->vertex(i));
        inc_cells[idi].push_back(c);
        if (!cell_is_selected)
          continue;

        const int dim = c3t3.in_dimension(c->vertex(i));
        switch (dim)
        {
        case 3:
          free_vertex[idi] = true;
          break;
        case 2:
          free_vertex[idi] = !protect_boundaries;
          break;
        case 1:
          free_vertex[idi] = !protect_boundaries && m_smooth_constrained_edges;
          break;
        case 0:
          free_vertex[idi] = false;
          break;
        default:
          CGAL_unreachable();
        }
      }
    }

    if (!protect_boundaries && m_smooth_constrained_edges)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      nb_done_1d =
#endif
      smooth_edges_in_complex(c3t3, vertex_id, free_vertex,
                              vertices_surface_indices, inc_cells, vertices_normals, total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                              , os_surf
#endif
      );
    }

    /////////////// EDGES ON SURFACE, BUT NOT IN COMPLEX //////////////////
    if (!protect_boundaries)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      nb_done_2d =
#endif
       smooth_vertices_on_surfaces(c3t3, vertex_id, free_vertex,
                                  vertices_surface_indices, inc_cells, vertices_normals,
                                  cell_selector, total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                                , os_surf, os_surf0
#endif
      );
    }
    CGAL_assertion(CGAL::Tetrahedral_remeshing::debug::are_cell_orientations_valid(tr));
    ////   end if(!protect_boundaries)

    ////////////// INTERNAL VERTICES ///////////////////////
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    nb_done_3d =
#endif
    smooth_internal_vertices(c3t3, vertex_id, free_vertex,
                             vertices_surface_indices, inc_cells, vertices_normals,
                             cell_selector, total_move
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
                           , os_vol
#endif
    );

    CGAL_assertion(CGAL::Tetrahedral_remeshing::debug::are_cell_orientations_valid(tr));

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::size_t nb_done = nb_done_3d + nb_done_2d + nb_done_1d;
    std::cout << " done ("
      << nb_done_1d << "/" << nb_done_2d << "/" << nb_done_3d << " vertices smoothed,"
      << " average move = " << (total_move / nb_done)
      << ")." << std::endl;
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
