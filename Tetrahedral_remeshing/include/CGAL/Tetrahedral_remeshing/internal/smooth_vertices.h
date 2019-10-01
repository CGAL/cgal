// Copyright (c) 2017  GeometryFactory (France).
// All rights reserved.

#ifndef CGAL_INTERNAL_SMOOTH_VERTICES_H
#define CGAL_INTERNAL_SMOOTH_VERTICES_H

#include <CGAL/Vector_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <boost/unordered_map.hpp>  

#include <vector>
#include <cmath>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
  template<typename Gt>
  CGAL::Vector_3<Gt> project_on_tangent_plane(const CGAL::Point_3<Gt>& gi,
                                              const CGAL::Point_3<Gt>& pi,
                                              const CGAL::Vector_3<Gt>& normal)
  {
    typedef typename Gt::Vector_3 Vector_3;
    Vector_3 diff = pi - gi;
    return Vector_3(gi, gi + (normal * diff) * normal);
  }

  template<typename C3t3, typename VertexNormalsMap>
  void compute_vertices_normals(const C3t3& c3t3,
                                VertexNormalsMap& normals_map)
  {
    typedef typename C3t3::Triangulation        Tr;
    typedef typename C3t3::Cell_handle          Cell_handle;
    typedef typename C3t3::Vertex_handle        Vertex_handle;
    typedef typename C3t3::Subdomain_index      Subdomain_index;
    typedef typename C3t3::Surface_patch_index  Surface_patch_index;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
    typedef typename Tr::Gt::Vector_3           Vector_3;

    const Tr& tr = c3t3.triangulation();

    for (Finite_facets_iterator fit = tr.finite_facets_begin();
         fit != tr.finite_facets_end(); ++fit)
    {
      Cell_handle ch = fit->first;
      Cell_handle n_ch = fit->first->neighbor(fit->second);

      Subdomain_index si = ch->subdomain_index();
      Subdomain_index si_mirror = n_ch->subdomain_index();

      if (si != si_mirror || tr.is_infinite(ch) || tr.is_infinite(n_ch))
      {
        Surface_patch_index surf_i = make_surface_patch_index(si, si_mirror);
        for (int i = 0; i < 3; ++i)
        {
          Vertex_handle v_id = fit->first->vertex(indices(fit->second ,i));
          normals_map[v_id][surf_i] = CGAL::NULL_VECTOR;
        }
      }
    }

    for (Finite_facets_iterator fit = tr.finite_facets_begin();
         fit != tr.finite_facets_end(); ++fit)
    {
      Cell_handle ch = fit->first;
      Cell_handle n_ch = fit->first->neighbor(fit->second);

      Subdomain_index si = ch->subdomain_index();
      Subdomain_index si_mirror = n_ch->subdomain_index();

      if (si != si_mirror || tr.is_infinite(ch) || tr.is_infinite(n_ch))
      {
        Surface_patch_index surf_i = make_surface_patch_index(si, si_mirror);

        Vector_3 n = CGAL::normal(*fit, tr.geom_traits());

        if (si < si_mirror || tr.is_infinite(ch))
          n = -1.*n;

        for (int i = 0; i < 3; ++i)
        {
          Vector_3& v_n = normals_map[fit->first->vertex(indices(fit->second, i))][surf_i];
          v_n = v_n + n;
        }
      }
    }

    //normalize the computed normals
    for (typename VertexNormalsMap::iterator vnm_it = normals_map.begin();
         vnm_it != normals_map.end(); ++vnm_it)
    {
      //value type is map<Surface_patch_index, Vector_3>
      for (typename VertexNormalsMap::value_type::iterator it = vnm_it->begin();
           it != vnm_it->end(); ++it)
      {
        Vector_3& n = it->second;
        n = n / CGAL::sqrt(n*n);
      }
    }
  }


  template<typename SurfacePatchIndex, typename Gt>
  bool project(const SurfacePatchIndex& /* si */,
                CGAL::Vector_3<Gt>& gi,
                CGAL::Vector_3<Gt>& projected_point)
  {
//    if (subdomain_FMLS_indices.find(si) == subdomain_FMLS_indices.end())
//      return false;
    typedef typename Gt::Vector_3 Vector_3;
    typedef typename Gt::Point_3  Point_3;

    if (std::isnan(gi.x()) || std::isnan(gi.y()) || isnan(gi.z()))
    {
      std::cout << "Initial point error " << gi << std::endl;
      return false;
    }

    Vector_3 res_normal;
    Point_3 point;
    Point_3 result = CGAL::ORIGIN + gi;

    //FMLS& fmls = subdomain_FMLS[subdomain_FMLS_indices[si]];

    // int it_nb = 0;
    // const int max_it_nb = 5;
    //const float epsilon = fmls.getPNScale() / 1000.;

    //do
    //{
    //  point = result;

    //  //fmls.fastProjectionCPU(point, result, res_normal);

    //  if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2])){
    //    std::cout << "MLS error detected si size " << si.first << " - " << si.second
    //              << " : " << fmls.getPNSize() << std::endl;
    //    return false;
    //  }

    //} while ((result - point).getLength() > epsilon && ++it_nb < max_it_nb);

    projected_point = Vector_3(result.x(), result.y(), result.z());

    return true;
  }

  template<typename Tr, typename K>
  bool check_inversion_and_move(const typename Tr::Vertex_handle v,
                                const CGAL::Vector_3<K>& move,
                                const std::vector<typename Tr::Cell_handle>& cells)
  {
    typedef typename Tr::Cell_handle Cell_handle;

    const typename Tr::Point backup = v->point(); //backup v's position
    typename Tr::Point new_pos(v->point().x() + move.x(),
                               v->point().y() + move.y(),
                               v->point().z() + move.z());
                               //note that weight is lost in case of Regular_triangulation
    v->set_point(new_pos);

    for (Cell_handle ci : cells)
    {
      if (CGAL::POSITIVE != CGAL::orientation(point(ci->vertex(0)->point()),
                                              point(ci->vertex(1)->point()),
                                              point(ci->vertex(2)->point()),
                                              point(ci->vertex(3)->point())))
      {
        v->set_point(backup);
        return false;
      }
    }
    return true;
  }

  template<typename C3T3>
  typename C3T3::Triangulation::Geom_traits::Vector_3
    move_3d(typename C3T3::Vertex_handle v, const C3T3& c3t3)
  {
    typedef typename C3T3::Edge Edge;
    typedef typename C3T3::Vertex_handle Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3 Vector_3;

    Vector_3 move = CGAL::NULL_VECTOR;

    std::vector<Edge> edges;
    c3t3.triangulation().incident_edges(v, std::back_inserter(edges));

    if (edges.empty())
      return move;

    BOOST_FOREACH(Edge e, edges)
    {
      Vertex_handle ve = (e.first->vertex(e.second) != v)
                        ? e.first->vertex(e.second)
                        : e.first->vertex(e.third);
      move = move + Vector_3(CGAL::ORIGIN,  point(ve->point()));
    }

    return 1. / edges.size() * move;
  }

  template<typename C3T3>
  typename C3T3::Triangulation::Geom_traits::Vector_3
    move_2d(typename C3T3::Vertex_handle v,
            const C3T3& c3t3,
            const typename C3T3::Subdomain_index& imaginary_index)
  {
    typedef typename C3T3::Edge Edge;
    typedef typename C3T3::Vertex_handle Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3 Vector_3;

    Vector_3 move = CGAL::NULL_VECTOR;

    std::vector<Edge> edges;
    c3t3.triangulation().incident_edges(v, std::back_inserter(edges));

    if (edges.empty())
      return move;

    std::size_t nbe = 0;
    BOOST_FOREACH(Edge e, edges)
    {
      if (is_on_domain_hull(e, c3t3, imaginary_index))
      {
        Vertex_handle ve = (e.first->vertex(e.second) != v)
                          ? e.first->vertex(e.second)
                          : e.first->vertex(e.third);
        move = move + Vector_3(CGAL::ORIGIN, point(ve->point()));
        ++nbe;
      }
    }

    if (nbe > 0)
      return (1. / nbe) * move;
    else
      return CGAL::NULL_VECTOR;
  }

  template<typename C3T3>
  typename C3T3::Triangulation::Geom_traits::Vector_3
    move_1d(typename C3T3::Vertex_handle v,
      const C3T3& c3t3,
            const typename C3T3::Subdomain_index& /*imaginary_index*/)
  {
    typedef typename C3T3::Edge Edge;
    typedef typename C3T3::Vertex_handle Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits::Vector_3 Vector_3;

    Vector_3 move = CGAL::NULL_VECTOR;

    std::vector<Edge> edges;
    c3t3.triangulation().incident_edges(v, std::back_inserter(edges));

    if (edges.empty())
      return move;

    std::size_t nbe = 0;
    BOOST_FOREACH(Edge e, edges)
    {
      if (!c3t3.is_in_complex(e))
        continue;

      Vertex_handle ve = (e.first->vertex(e.second) != v)
                      ? e.first->vertex(e.second)
                      : e.first->vertex(e.third);

      move = move + Vector_3(CGAL::ORIGIN, point(ve->point()));
      ++nbe;
    }

    if (nbe == 2)
      return 0.5 * move;
    else
      return CGAL::NULL_VECTOR;
  }

  template<typename C3T3, typename CellSelector>
  void smooth_vertices_new(C3T3& c3t3,
                           const typename C3T3::Subdomain_index& imaginary_index,
                           const bool /*protect_boundaries*/,
                           CellSelector cell_selector)
  {
    typedef typename C3T3::Triangulation          Tr;
    typedef typename C3T3::Vertex_handle          Vertex_handle;
    typedef typename C3T3::Cell_handle            Cell_handle;
    typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename C3T3::Point                  Point;
    typedef typename Tr::Geom_traits::Vector_3    Vector_3;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Smooth vertices...";
    std::cout.flush();
#endif

    Tr& tr = c3t3.triangulation();

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      c3t3.triangulation(), "c3t3_vertices_before_smoothing");
#endif

    const std::size_t nbv = tr.number_of_vertices();
    boost::unordered_map<Vertex_handle, std::size_t> vertex_id;
    std::vector<Vector_3> smoothing_vecs(nbv, CGAL::NULL_VECTOR);

    // generate ids for vertices
    std::size_t id = 0;
    for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
      vertex_id[vit] = id++;
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::ofstream ofs_2d("moves_on_surface.polylines.txt");
    std::ofstream ofs_1d("moves_on_features.polylines.txt");
#endif

    // compute move depending on underlying dimension
    for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
      switch (vit->in_dimension())
      {
      case 3:
        if ( is_imaginary(vit, c3t3, imaginary_index)
          || !is_selected(vit, c3t3, cell_selector))
          break;
        else
          smoothing_vecs[vertex_id.at(vit)] = move_3d(vit, c3t3);
        break;

      case 2:
        smoothing_vecs[vertex_id.at(vit)] = move_2d(vit, c3t3, imaginary_index);
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        if (smoothing_vecs[vertex_id.at(vit)] != CGAL::NULL_VECTOR)
          ofs_2d << "2 " << vit->point()
            << " " << (CGAL::ORIGIN + smoothing_vecs[vertex_id.at(vit)]) << std::endl;
#endif
        break;

      case 1:
        smoothing_vecs[vertex_id.at(vit)] = move_1d(vit, c3t3, imaginary_index);
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        if (smoothing_vecs[vertex_id.at(vit)] != CGAL::NULL_VECTOR)
          ofs_1d << "2 " << vit->point()
          << " " << (CGAL::ORIGIN + smoothing_vecs[vertex_id.at(vit)]) << std::endl;
#endif

      default:
        break;
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    ofs_2d.close();
    ofs_1d.close();
#endif

    // apply moves
    for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
      const std::size_t& vid = vertex_id.at(vit);
      const Point new_pos(CGAL::ORIGIN + smoothing_vecs[vid]);
      const Vector_3 move(point(vit->point()), new_pos);

      std::vector<Cell_handle> cells;
      tr.finite_incident_cells(vit, std::back_inserter(cells));

      double frac = 1.;
      while (frac > 0.05   /// 1/16 = 0.0625
          && !check_inversion_and_move<Tr>(vit, frac * move, cells))
      {
        frac = 0.5 * frac;
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      c3t3.triangulation(), "c3t3_vertices_after_smoothing");
#endif
  }

//  template<typename C3T3, typename CellSelector>
//  void smooth_vertices(C3T3& c3t3,
//                       const typename C3T3::Subdomain_index&,
//                       const bool protect_boundaries,
//                       CellSelector cell_selector)
//  {
//    typedef typename C3T3::Surface_patch_index    Surface_patch_index;
//    typedef typename C3T3::Subdomain_index        Subdomain_index;
//    typedef typename C3T3::Triangulation          Tr;
//    typedef typename C3T3::Vertex_handle          Vertex_handle;
//    typedef typename C3T3::Cell_handle            Cell_handle;
//    typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
//    typedef typename Tr::Finite_edges_iterator    Finite_edges_iterator;
//
//    typedef typename Tr::Geom_traits     Gt;
//    typedef typename Gt::Point_3         Point_3;
//    typedef typename Gt::Vector_3        Vector_3;
//    typedef typename Gt::FT              FT;
//
//#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//    std::cout << "Smooth vertices...";
//    std::cout.flush();
//    std::size_t nb_done = 0;
//#endif
//
//    Tr& tr = c3t3.triangulation();
//
//    const std::size_t nbv = tr.number_of_vertices();
//    boost::unordered_map<Vertex_handle, std::size_t> vertex_id;
//    std::vector<Vector_3> smoothing_vecs(nbv, CGAL::NULL_VECTOR);
//    std::vector<int> neighbors(nbv, -1);
//
//    //collect ids
//    std::size_t id = 0;
//    for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
//         vit != tr.finite_vertices_end(); ++vit)
//    {
//      vertex_id[vit] = id++;
//    }
//
//    if (!protect_boundaries)
//    {
//      for (Finite_edges_iterator eit = tr.finite_edges_begin();
//           eit != tr.finite_edges_end(); ++eit)
//      {
//        const Vertex_handle vh0 = eit->first->vertex(eit->second);
//        const Vertex_handle vh1 = eit->first->vertex(eit->third);
//
//        const std::size_t& i0 = vertex_id.at(vh0);
//        const std::size_t& i1 = vertex_id.at(vh1);
//
//        if (/*toRemesh != REMESH_IMAGINARY &&*/ c3t3.is_in_complex(*eit))
//        {
//          if (!is_feature(vh0, c3t3))
//            neighbors[i0] = std::max(0, neighbors[i0]);
//          if (!is_feature(vh1, c3t3))
//            neighbors[i1] = std::max(0, neighbors[i1]);
//
//          bool update_v0 = false, update_v1 = false;
//
//          helpers::get_edge_info(*eit, update_v0, update_v1, c3t3, cell_selector);
//          if (update_v0)
//          {
//            const Point_3& p1 = vh1->point();
//            smoothing_vecs[i0] = smoothing_vecs[i0] + Vector_3(p1.x(), p1.y(), p1.z());
//            neighbors[i0]++;
//          }
//          if (update_v1)
//          {
//            const Point_3& p0 = vh0->point();
//            smoothing_vecs[i1] = smoothing_vecs[i1] + Vector_3(p0.x(), p0.y(), p0.z());
//            neighbors[i1]++;
//          }
//        }
//      }
//
//      //collect a map of vertices subdomain indices
//      boost::unordered_map<Vertex_handle, std::vector<Subdomain_index> > vertices_subdomain_indices;
//      for (typename C3T3::Cell_iterator cit = c3t3.cells_in_complex_begin();
//           cit != c3t3.cells_in_complex_end(); ++cit)
//      {
//        for (int i = 0; i < 4; ++i)
//        {
//          Vertex_handle vi = cit->vertex(i);
//          Subdomain_index si = cit->subdomain_index();
//
//          if (vertices_subdomain_indices.find(vi) == vertices_subdomain_indices.end())
//          {
//            std::vector<Subdomain_index> indices(1);
//            indices[0] = si;
//            vertices_subdomain_indices.insert(std::make_pair(vi, indices));
//          }
//          else
//          {
//            std::vector<Subdomain_index>& v_indices = vertices_subdomain_indices.at(vi);
//            if (std::find(v_indices.begin(), v_indices.end(), si) == v_indices.end())
//              v_indices.push_back(si);
//          }
//        }
//      }
//
//      //collect a map of vertices surface indices
//      boost::unordered_map<Vertex_handle, std::vector<Surface_patch_index> > vertices_surface_indices;
//      for(typename C3T3::Facet_iterator fit = c3t3.facets_in_complex_begin();
//          fit != c3t3.facets_in_complex_end(); ++fit)
//      { 
//        Surface_patch_index surface_index
//          = helpers::make_surface_patch_index(fit->first->subdomain_index(),
//                                              fit->first->neighbor(fit->second)->subdomain_index());
//        for (int i = 0; i < 3; ++i)
//        {
//          Vertex_handle vi = fit->first->vertex(indices(fit->second, i));
//          if (vertices_subdomain_indices.at(vi).size() > 2)
//          {
//            if (vertices_surface_indices.find(vi) == vertices_surface_indices.end())
//            {
//              std::vector<Surface_patch_index> indices(1);
//              indices[0] = surface_index;
//              vertices_surface_indices.insert(std::make_pair(vi, indices));
//            }
//            else
//            {
//              std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices.at(vi);
//              if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index)
//                  == v_surface_indices.end())
//                v_surface_indices.push_back(surface_index);
//            }
//          }
//        }
//      }
//
//      //collect a map of normals at surface vertices
//      boost::unordered_map<Vertex_handle,
//                           boost::unordered_map<Surface_patch_index, Vector_3> > vertices_normals;
//      for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
//           vit != tr.finite_vertices_end(); ++vit)
//      {
//        const std::size_t& vid = vertex_id.at(vit);
//        if (neighbors[vid] > 1)
//        {
//          Point_3 smoothed_position = CGAL::ORIGIN + smoothing_vecs[vid] / neighbors[vid];
//          Vector_3 final_move = CGAL::NULL_VECTOR;
//          Point_3 final_position;
//
//          std::size_t count = 0;
//          Point_3 current_pos = vit->point();
//
//          const std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[vit];
//          for (std::size_t i = 0; i < v_surface_indices.size(); ++i)
//          {
//            const Surface_patch_index& si = v_surface_indices[i];
//
//            Vector_3 normal_projection
//              = project_on_tangent_plane(smoothed_position, current_pos, vertices_normals[vit][si]);
//
//            //Check if the mls surface exists to avoid degenrated cases
//            Vector_3 mls_projection;
//            if (project(si, normal_projection, mls_projection)){
//              final_move = final_move + mls_projection;
//            }
//            else {
//              final_move = final_move + normal_projection;
//            }
//            count++;
//          }
//
//          if (count > 0)
//            final_position = CGAL::ORIGIN + final_move / static_cast<FT>(count);
//          else
//            final_position = smoothed_position;
//
//          // move vertex
//          vit->set_point(final_position);
//
//        }
//        else if (neighbors[vid] > 0)
//        {
//          Vector_3 final_move = CGAL::NULL_VECTOR;
//          Point_3 final_position;
//
//          int count = 0;
//          Vector_3 current_move(CGAL::ORIGIN, vit->point());
//
//          const std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[vit];
//          for (std::size_t i = 0; i < v_surface_indices.size(); ++i)
//          {
//            Surface_patch_index si = v_surface_indices[i];
//            //Check if the mls surface exists to avoid degenrated cases
//
//            Vector_3 mls_projection;
//            if (project(si, current_move, mls_projection)){
//              final_move = final_move + mls_projection;
//            }
//            else {
//              final_move = final_move + current_move;
//            }
//            count++;
//          }
//
//          if (count > 0)
//            final_position = CGAL::ORIGIN + final_move / count;
//          else
//            final_position = CGAL::ORIGIN + current_move;
//
//          // move vertex
//          vit->set_point(final_position);
//        }
//      }
//
//      smoothing_vecs.clear();
//      smoothing_vecs.resize(nbv, CGAL::NULL_VECTOR);
//
//      neighbors.clear();
//      neighbors.resize(nbv, -1);
//
//      for (Finite_edges_iterator eit = tr.finite_edges_begin();
//           eit != tr.finite_edges_end(); ++eit)
//      {
//        const Vertex_handle vh0 = eit->first->vertex(eit->second);
//        const Vertex_handle vh1 = eit->first->vertex(eit->third);
//
//        const std::size_t& i0 = vertex_id.at(vh0);
//        const std::size_t& i1 = vertex_id.at(vh1);
//
//        if ((/*toRemesh != REMESH_IN_COMPLEX &&*/ is_on_hull(*eit, c3t3))
//         || (/*toRemesh != REMESH_IMAGINARY &&*/
//             helpers::is_boundary(c3t3, *eit, cell_selector) && !c3t3.is_in_complex(*eit)))
//        {
//          bool update_v0 = false, update_v1 = false;
//          if (!is_feature(vh0, c3t3))
//            neighbors[i0] = (std::max)(0, neighbors[i0]);
//          if (!is_feature(vh1, c3t3))
//            neighbors[i1] = (std::max)(0, neighbors[i1]);
//
//          helpers::get_edge_info(*eit, update_v0, update_v1, c3t3, cell_selector);
//          if (update_v0)
//          {
//            const Point_3& p1 = vh1->point();
//            smoothing_vecs[i0] = smoothing_vecs[i0] + Vector_3(p1.x(), p1.y(), p1.z());
//            neighbors[i0]++;
//          }
//          if (update_v1)
//          {
//            const Point_3& p0 = vh0->point();
//            smoothing_vecs[i1] = smoothing_vecs[i1] + Vector_3(p0.x(), p0.y(), p0.z());
//            neighbors[i1]++;
//          }
//        }
//      }
//
//      for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
//           vit != tr.finite_vertices_end(); ++vit)
//      {
//        const std::size_t& vid = vertex_id.at(vit);
//
//        if (neighbors[vid] > 1)
//        {
//          Point_3 smoothed_position = CGAL::ORIGIN + smoothing_vecs[vid] / neighbors[vid];
//          Point_3 current_pos = vit->point();
//          Point_3 final_position = CGAL::ORIGIN;
//
//          if (vit->in_dimension() == 3 && is_on_hull(vit, c3t3))
//          {
//            Vector_3 final_move = project_on_tangent_plane(
//              smoothed_position, current_pos, vertices_normals[vit][Surface_patch_index()]);
//            final_position = CGAL::ORIGIN + final_move;
//          }
//          else {
//            // Surface_patch_index si = helpers::make_surface_patch_index(
//            //   vertices_subdomain_indices[vit][0], vertices_subdomain_indices[vit][1]);
//
//            // Vector_3 normal_projection = project_on_tangent_plane(smoothed_position,
//            //                                                       current_pos,
//            //                                                       vertices_normals[vit][si]);
//            //Vector_3 mls_projection;
//            //if (project(si, normal_projection, mls_projection) /*|| project( si, smoothed_position, mls_projection )*/){
//            //  final_position = mls_projection;
//            //  //final_position = smoothed_position;
//            //}
//            //else {
//              final_position = smoothed_position;
//            //}
//            // std::cout << "MLS " << final_position[0] << " - " << final_position[1] << " : " << final_position[2] << std::endl;
//          }
//          /*
//          Normal_iterator it = vertices_normals[vit->info()].end();
//          it--;
//          final_position = final_position + projectOnTangentPlane( smoothed_position, current_pos , it->second );
//          */
//
//          vit->set_point(final_position);
//        }
//        else if (neighbors[vid] > 0)
//        {
//          if (vit->in_dimension() == 2)
//          {
//            // Surface_patch_index si = helpers::make_surface_patch_index(
//            //   vertices_subdomain_indices[vit][0],
//            //   vertices_subdomain_indices[vit][1]);
//
//            Vector_3 current_pos(CGAL::ORIGIN, vit->point());
//            Vector_3 mls_projection;
////            if (project(si, current_pos, mls_projection) /*|| project( si, smoothed_position, mls_projection )*/){
////              vit->set_point(Point_3(mls_projection.x(), mls_projection.y(), mls_projection.z()));
////            }
//          }
//        }
//      }
//    }
//    smoothing_vecs.clear();
//    smoothing_vecs.resize(nbv, CGAL::NULL_VECTOR);
//
//    neighbors.clear();
//    neighbors.resize(nbv, 0);
//
//    for (Finite_edges_iterator eit = tr.finite_edges_begin();
//         eit != tr.finite_edges_end(); ++eit)
//    {
//      //bool in_complex = c3t3.is_in_complex(*eit);
//      //if (  toRemesh == REMESH_ALL
//      //  || (toRemesh == REMESH_IN_COMPLEX && in_complex)
//      //  || (toRemesh == REMESH_IMAGINARY && !in_complex))
//      {
//        const Vertex_handle vh0 = eit->first->vertex(eit->second);
//        const Vertex_handle vh1 = eit->first->vertex(eit->third);
//
//        const std::size_t& i0 = vertex_id.at(vh0);
//        const std::size_t& i1 = vertex_id.at(vh1);
//
//        if (c3t3.in_dimension(vh0) == 3 && !is_on_hull(vh0, c3t3))
//        {
//          const Point_3& p1 = vh1->point();
//          smoothing_vecs[i0] = smoothing_vecs[i0] + Vector_3(CGAL::ORIGIN, p1);
//          neighbors[i0]++;
//        }
//        if (c3t3.in_dimension(vh1) == 3 && !is_on_hull(vh1, c3t3))
//        {
//          const Point_3& p0 = vh0->point();
//          smoothing_vecs[i1] = smoothing_vecs[i1] + Vector_3(CGAL::ORIGIN, p0);
//          neighbors[i1]++;
//        }
//      }
//    }
//
//    for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
//         vit != tr.finite_vertices_end(); ++vit)
//    {
//      const std::size_t& vid = vertex_id.at(vit);
//      if (neighbors[vid] > 1)
//      {
//        if (smoothing_vecs[vid] != CGAL::NULL_VECTOR)
//        {
//#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//          ++nb_done;
//#endif
//          Point_3 new_pos = CGAL::ORIGIN + smoothing_vecs[vid] / neighbors[vid];
//          const Vector_3 move(vit->point(), new_pos);
//
//          std::vector<Cell_handle> cells;
//          tr.finite_incident_cells(vit, std::back_inserter(cells));
//
//          bool selected = true;
//          for (std::size_t i = 0; i < cells.size(); ++i)
//          {
//            if (!cell_selector(cells[i]))
//            {
//              selected = false;
//              break;
//            }
//          }
//          if (!selected)
//            continue;
//
//          double frac = 1.;
//          while (frac > 0.05   /// 1/16 = 0.0625
//                && !check_inversion_and_move(vit, frac * move, cells))
//          {
//            frac = 0.5 * frac;
//          }
//        }
//      }
//    }
//
//#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//    std::cout << " done (" << nb_done << " vertices smoothed)." << std::endl;
//#endif
//  }

}//namespace internal
}//namespace Tetrahedral_adaptive_remeshing
}//namespace CGAL

#endif //CGAL_INTERNAL_SMOOTH_VERTICES_H
