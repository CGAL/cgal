// Copyright (c) 2017  GeometryFactory (France).
// All rights reserved.

#ifndef CGAL_INTERNAL_SMOOTH_VERTICES_H
#define CGAL_INTERNAL_SMOOTH_VERTICES_H

#include <CGAL/Vector_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>

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
        typedef typename Tr::Facet                  Facet;
        typedef typename Tr::Geom_traits::Vector_3  Vector_3;

        typename Tr::Geom_traits::Construct_opposite_vector_3
          opp = c3t3.triangulation().geom_traits().construct_opposite_vector_3_object();
        typename Tr::Geom_traits::Construct_scaled_vector_3
          scale = c3t3.triangulation().geom_traits().construct_scaled_vector_3_object();

        const Tr& tr = c3t3.triangulation();

        for (const Facet& f : tr.finite_facets())
        {
          if (c3t3.is_in_complex(f))
          {
            const Surface_patch_index surf_i = c3t3.surface_patch_index(f);
            for (int i = 0; i < 3; ++i)
            {
              Vertex_handle v_id = f.first->vertex(indices(f.second, i));
              normals_map[v_id][surf_i] = CGAL::NULL_VECTOR;
            }
          }
        }

        for (const Facet& f : tr.finite_facets())
        {
          const Cell_handle ch = f.first;
          const Cell_handle n_ch = f.first->neighbor(f.second);

          const Subdomain_index si = ch->subdomain_index();
          const Subdomain_index si_mirror = n_ch->subdomain_index();

          if (c3t3.is_in_complex(f))
          {
            const Surface_patch_index surf_i = c3t3.surface_patch_index(f);

            Vector_3 n = CGAL::Tetrahedral_remeshing::normal(f, tr.geom_traits());

            if (si < si_mirror || tr.is_infinite(ch)) // todo : fix this condition
              n = opp(n);

            for (int i = 0; i < 3; ++i)
            {
              Vector_3& v_n = normals_map[f.first->vertex(indices(f.second, i))][surf_i];
              v_n = v_n + n;
            }
          }
        }

        //normalize the computed normals
        for (typename VertexNormalsMap::iterator vnm_it = normals_map.begin();
             vnm_it != normals_map.end(); ++vnm_it)
        {
          //value type is map<Surface_patch_index, Vector_3>
          for (typename VertexNormalsMap::mapped_type::iterator it = vnm_it->second.begin();
               it != vnm_it->second.end(); ++it)
          {
            Vector_3& n = it->second;
            n = scale(n, 1. / CGAL::approximate_sqrt(n * n));
          }
        }
      }


      template<typename SurfacePatchIndex,
               typename Gt,
               typename Subdomain__FMLS,
               typename Subdomain__FMLS_indices>
      bool project(const SurfacePatchIndex& si,
                   const CGAL::Vector_3<Gt>& gi,
                   CGAL::Vector_3<Gt>& projected_point,
                   Subdomain__FMLS& subdomain_FMLS,
                   Subdomain__FMLS_indices& subdomain_FMLS_indices)
      {
        if (subdomain_FMLS_indices.find(si) == subdomain_FMLS_indices.end())
          return false;

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

        CGAL::Tetrahedral_remeshing::internal::FMLS&
          fmls = subdomain_FMLS[subdomain_FMLS_indices.at(si)];

        int it_nb = 0;
        const int max_it_nb = 5;
        const float epsilon = fmls.getPNScale() / 1000.;
        const float sq_eps = CGAL::square(epsilon);

        do
        {
          point = result;

          fmls.fastProjectionCPU(point, result, res_normal);

          if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2])){
            std::cout << "MLS error detected si size " << si
                      << " : " << fmls.getPNSize() << std::endl;
            return false;
        }
        } while (CGAL::squared_distance(result, point) > sq_eps && ++it_nb < max_it_nb);

        projected_point = Vector_3(result.x(), result.y(), result.z());

        return true;
      }

      template<typename K, typename CellVector, typename Tr>
      bool check_inversion_and_move(const typename Tr::Vertex_handle v,
                                    const CGAL::Vector_3<K>& move,
                                    const CellVector& cells,
                                    const Tr& tr)
      {
        const typename Tr::Point backup = v->point(); //backup v's position
        const typename Tr::Point new_pos(point(backup) + move);
        v->set_point(new_pos);

        for(const typename CellVector::value_type& ci : cells)
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
      void collect_vertices_subdomain_indices(
        const C3T3& c3t3,
        boost::unordered_map<
          typename C3T3::Vertex_handle,
          std::vector<typename C3T3::Subdomain_index> >& vertices_subdomain_indices)
      {
        typedef typename C3T3::Subdomain_index Subdomain_index;
        typedef typename C3T3::Vertex_handle   Vertex_handle;

        for (typename C3T3::Cell_iterator cit = c3t3.cells_in_complex_begin();
             cit != c3t3.cells_in_complex_end(); ++cit)
        {
          const Subdomain_index si = cit->subdomain_index();
          for (int i = 0; i < 4; ++i)
          {
            const Vertex_handle vi = cit->vertex(i);
            if (vertices_subdomain_indices.find(vi) == vertices_subdomain_indices.end())
            {
              std::vector<Subdomain_index> indices(1);
              indices[0] = si;
              vertices_subdomain_indices.insert(std::make_pair(vi, indices));
            }
            else
            {
              std::vector<Subdomain_index>& v_indices = vertices_subdomain_indices.at(vi);
              if (std::find(v_indices.begin(), v_indices.end(), si) == v_indices.end())
                v_indices.push_back(si);
            }
          }
        }
      }

      template<typename C3T3>
      void collect_vertices_surface_indices(
        const C3T3& c3t3,
        const boost::unordered_map<
            typename C3T3::Vertex_handle,
            std::vector<typename C3T3::Subdomain_index> >& vertices_subdomain_indices,
        boost::unordered_map<
            typename C3T3::Vertex_handle,
            std::vector<typename C3T3::Surface_patch_index> >& vertices_surface_indices)
      {
        typedef typename C3T3::Surface_patch_index Surface_patch_index;
        typedef typename C3T3::Vertex_handle       Vertex_handle;
        typedef typename C3T3::Facet               Facet;

        for (typename C3T3::Facet_iterator fit = c3t3.facets_in_complex_begin();
             fit != c3t3.facets_in_complex_end(); ++fit)
        {
          const Facet& f = *fit;
          const Surface_patch_index surface_index = c3t3.surface_patch_index(f);

          for (int i = 0; i < 3; ++i)
          {
            const Vertex_handle vi = f.first->vertex(indices(f.second, i));
            if (vertices_subdomain_indices.at(vi).size() > 2)
            {
              if (vertices_surface_indices.find(vi) == vertices_surface_indices.end())
              {
                std::vector<Surface_patch_index> indices(1);
                indices[0] = surface_index;
                vertices_surface_indices.insert(std::make_pair(vi, indices));
              }
              else
              {
                std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices.at(vi);
                if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index)
                  == v_surface_indices.end())
                  v_surface_indices.push_back(surface_index);
              }
            }
          }
        }
      }


      template<typename C3T3, typename CellSelector>
      void smooth_vertices(C3T3& c3t3,
                           const bool protect_boundaries,
                           CellSelector cell_selector)
      {
        typedef typename C3T3::Surface_patch_index    Surface_patch_index;
        typedef typename C3T3::Subdomain_index        Subdomain_index;
        typedef typename C3T3::Triangulation          Tr;
        typedef typename C3T3::Vertex_handle          Vertex_handle;
        typedef typename C3T3::Cell_handle            Cell_handle;
        typedef typename C3T3::Facet                  Facet;
        typedef typename Tr::Edge                     Edge;

        typedef typename Tr::Geom_traits     Gt;
        typedef typename Gt::Point_3         Point_3;
        typedef typename Gt::Vector_3        Vector_3;
        typedef typename Gt::FT              FT;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << "Smooth vertices...";
        std::cout.flush();
        std::size_t nb_done = 0;
#endif

        Tr& tr = c3t3.triangulation();

        //collect a map of vertices subdomain indices
        boost::unordered_map<Vertex_handle, std::vector<Subdomain_index> > vertices_subdomain_indices;
        collect_vertices_subdomain_indices(c3t3, vertices_subdomain_indices);

        //collect a map of vertices surface indices
        boost::unordered_map<Vertex_handle, std::vector<Surface_patch_index> > vertices_surface_indices;
        collect_vertices_surface_indices(c3t3, vertices_subdomain_indices, vertices_surface_indices);

        //collect a map of normals at surface vertices
        boost::unordered_map<Vertex_handle,
          boost::unordered_map<Surface_patch_index, Vector_3> > vertices_normals;
        compute_vertices_normals(c3t3, vertices_normals);

        // Build MLS Surfaces
        std::vector < CGAL::Tetrahedral_remeshing::internal::FMLS > subdomain_FMLS;
        boost::unordered_map<Surface_patch_index, unsigned int> subdomain_FMLS_indices;
        createMLSSurfaces(subdomain_FMLS,
                          subdomain_FMLS_indices,
                          vertices_normals,
                          c3t3);

        //smooth()
        const std::size_t nbv = tr.number_of_vertices();
        boost::unordered_map<Vertex_handle, std::size_t> vertex_id;
        std::vector<Vector_3> smoothing_vecs(nbv, CGAL::NULL_VECTOR);
        std::vector<int> neighbors(nbv, -1);

        //collect ids
        std::size_t id = 0;
        for (const Vertex_handle v : tr.finite_vertex_handles())
        {
          vertex_id[v] = id++;
        }

        if (!protect_boundaries)
        {
          for (const Edge& e : tr.finite_edges())
          {
            const Vertex_handle vh0 = e.first->vertex(e.second);
            const Vertex_handle vh1 = e.first->vertex(e.third);

            const std::size_t& i0 = vertex_id.at(vh0);
            const std::size_t& i1 = vertex_id.at(vh1);

            if (c3t3.is_in_complex(e))
            {
              if (!is_feature(vh0, c3t3))
                neighbors[i0] = (std::max)(0, neighbors[i0]);
              if (!is_feature(vh1, c3t3))
                neighbors[i1] = (std::max)(0, neighbors[i1]);

              bool update_v0 = false, update_v1 = false;

              get_edge_info(e, update_v0, update_v1, c3t3, cell_selector);
              if (update_v0)
              {
                const Point_3& p1 = point(vh1->point());
                smoothing_vecs[i0] = smoothing_vecs[i0] + Vector_3(p1.x(), p1.y(), p1.z());
                neighbors[i0]++;
              }
              if (update_v1)
              {
                const Point_3& p0 = point(vh0->point());
                smoothing_vecs[i1] = smoothing_vecs[i1] + Vector_3(p0.x(), p0.y(), p0.z());
                neighbors[i1]++;
              }
            }
          }

          // Smooth
          for (Vertex_handle v : tr.finite_vertex_handles())
          {
            const std::size_t& vid = vertex_id.at(v);
            if (neighbors[vid] > 1)
            {
              Point_3 smoothed_position = CGAL::ORIGIN + smoothing_vecs[vid] / neighbors[vid];
              Vector_3 final_move = CGAL::NULL_VECTOR;
              Point_3 final_position = CGAL::ORIGIN;

              std::size_t count = 0;
              const Point_3 current_pos = point(v->point());

              const std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[v];
              for (std::size_t i = 0; i < v_surface_indices.size(); ++i)
              {
                const Surface_patch_index& si = v_surface_indices[i];

                Vector_3 normal_projection
                  = project_on_tangent_plane(smoothed_position, current_pos, vertices_normals[v][si]);

                //Check if the mls surface exists to avoid degenrated cases
                Vector_3 mls_projection;
                if (project(si, normal_projection, mls_projection, subdomain_FMLS, subdomain_FMLS_indices)) {
                  final_move = final_move + mls_projection;
                }
                else {
                  final_move = final_move + normal_projection;
                }
                count++;
              }

              if (count > 0)
                final_position = CGAL::ORIGIN + final_move / static_cast<FT>(count);
              else
                final_position = smoothed_position;

              // move vertex
              v->set_point(typename Tr::Point(final_position));

            }
            else if (neighbors[vid] > 0)
            {
              Vector_3 final_move = CGAL::NULL_VECTOR;
              Point_3 final_position;

              int count = 0;
              Vector_3 current_move(CGAL::ORIGIN, point(v->point()));

              const std::vector<Surface_patch_index>& v_surface_indices = vertices_surface_indices[v];
              for (std::size_t i = 0; i < v_surface_indices.size(); ++i)
              {
                Surface_patch_index si = v_surface_indices[i];
                //Check if the mls surface exists to avoid degenrated cases

                Vector_3 mls_projection;
                if (project(si, current_move, mls_projection, subdomain_FMLS, subdomain_FMLS_indices)) {
                  final_move = final_move + mls_projection;
                }
                else {
                  final_move = final_move + current_move;
                }
                count++;
              }

              if (count > 0)
                final_position = CGAL::ORIGIN + final_move / count;
              else
                final_position = CGAL::ORIGIN + current_move;

              // move vertex
              v->set_point(typename Tr::Point(final_position));
            }
          }

          smoothing_vecs.clear();
          smoothing_vecs.resize(nbv, CGAL::NULL_VECTOR);

          neighbors.clear();
          neighbors.resize(nbv, -1);

          for (const Edge& e : tr.finite_edges())
          {
            const Vertex_handle vh0 = e.first->vertex(e.second);
            const Vertex_handle vh1 = e.first->vertex(e.third);

            const std::size_t& i0 = vertex_id.at(vh0);
            const std::size_t& i1 = vertex_id.at(vh1);

            if (is_boundary(c3t3, e, cell_selector) && !c3t3.is_in_complex(e))
            {
              bool update_v0 = false, update_v1 = false;
              if (!is_feature(vh0, c3t3))
                neighbors[i0] = (std::max)(0, neighbors[i0]);
              if (!is_feature(vh1, c3t3))
                neighbors[i1] = (std::max)(0, neighbors[i1]);

              get_edge_info(e, update_v0, update_v1, c3t3, cell_selector);
              if (update_v0)
              {
                const Point_3& p1 = point(vh1->point());
                smoothing_vecs[i0] = smoothing_vecs[i0] + Vector_3(p1.x(), p1.y(), p1.z());
                neighbors[i0]++;
              }
              if (update_v1)
              {
                const Point_3& p0 = point(vh0->point());
                smoothing_vecs[i1] = smoothing_vecs[i1] + Vector_3(p0.x(), p0.y(), p0.z());
                neighbors[i1]++;
              }
            }
          }

          for (Vertex_handle v : tr.finite_vertex_handles())
          {
            const std::size_t& vid = vertex_id.at(v);

            if (neighbors[vid] > 1)
            {
              Point_3 smoothed_position = CGAL::ORIGIN + smoothing_vecs[vid] / neighbors[vid];
              const Point_3& current_pos = point(v->point());
              Point_3 final_position = CGAL::ORIGIN;

              if (v->in_dimension() == 3 && is_on_convex_hull(v, c3t3))
              {
                Vector_3 final_move = project_on_tangent_plane(
                  smoothed_position, current_pos, vertices_normals[v][Surface_patch_index()]);
                final_position = CGAL::ORIGIN + final_move;
              }
              else {
                const Surface_patch_index si = surface_patch_index(v, c3t3);

                Vector_3 normal_projection = project_on_tangent_plane(smoothed_position,
                                                                      current_pos,
                                                                      vertices_normals[v][si]);
                Vector_3 mls_projection;
                if (project(si, normal_projection, mls_projection, subdomain_FMLS, subdomain_FMLS_indices)
                    /*|| project( si, smoothed_position, mls_projection )*/){
                  final_position = CGAL::ORIGIN + mls_projection;
                }
                else {
                  final_position = smoothed_position;
                }
                // std::cout << "MLS " << final_position[0] << " - " << final_position[1] << " : " << final_position[2] << std::endl;
              }

              v->set_point(typename Tr::Point(final_position));
            }
            else if (neighbors[vid] > 0)
            {
              if (v->in_dimension() == 2)
              {
                const Surface_patch_index si = surface_patch_index(v, c3t3);

                const Vector_3 current_pos(CGAL::ORIGIN, point(v->point()));
                Vector_3 mls_projection;
                if (project(si, current_pos, mls_projection, subdomain_FMLS, subdomain_FMLS_indices)
                  /*|| project( si, smoothed_position, mls_projection )*/)
                {
                  const typename Tr::Point new_pos(CGAL::ORIGIN + mls_projection);
                  v->set_point(new_pos);
                }
              }
            }
          }
        }
        smoothing_vecs.clear();
        smoothing_vecs.resize(nbv, CGAL::NULL_VECTOR);

        neighbors.clear();
        neighbors.resize(nbv, 0);

        for (const Edge& e : tr.finite_edges())
        {
          if ( !is_outside(e, c3t3, cell_selector))
          {
            const Vertex_handle vh0 = e.first->vertex(e.second);
            const Vertex_handle vh1 = e.first->vertex(e.third);

            const std::size_t& i0 = vertex_id.at(vh0);
            const std::size_t& i1 = vertex_id.at(vh1);

            if (c3t3.in_dimension(vh0) == 3 && !is_on_convex_hull(vh0, c3t3))
            {
              const Point_3& p1 = point(vh1->point());
              smoothing_vecs[i0] = smoothing_vecs[i0] + Vector_3(CGAL::ORIGIN, p1);
              neighbors[i0]++;
            }
            if (c3t3.in_dimension(vh1) == 3 && !is_on_convex_hull(vh1, c3t3))
            {
              const Point_3& p0 = point(vh0->point());
              smoothing_vecs[i1] = smoothing_vecs[i1] + Vector_3(CGAL::ORIGIN, p0);
              neighbors[i1]++;
            }
          }
        }

        for (Vertex_handle v : tr.finite_vertex_handles())
        {
          const std::size_t& vid = vertex_id.at(v);
          if (neighbors[vid] > 1)
          {
            if (smoothing_vecs[vid] != CGAL::NULL_VECTOR)
            {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
              ++nb_done;
#endif
              Point_3 new_pos = CGAL::ORIGIN + smoothing_vecs[vid] / neighbors[vid];
              const Vector_3 move(point(v->point()), new_pos);

              std::vector<Cell_handle> cells;
              tr.finite_incident_cells(v, std::back_inserter(cells));

              bool selected = true;
              for (const Cell_handle ci : cells)
              {
                if (!cell_selector(ci))
                {
                  selected = false;
                  break;
                }
              }
              if (!selected)
                continue;

              double frac = 1.;
              while (frac > 0.05   /// 1/16 = 0.0625
                && !check_inversion_and_move(v, frac * move, cells, tr))
              {
                frac = 0.5 * frac;
              }
            }
          }
        }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << " done (" << nb_done << " vertices smoothed)." << std::endl;
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
          c3t3.triangulation(), "c3t3_vertices_after_smoothing");
#endif
      }

    }//namespace internal
  }//namespace Tetrahedral_adaptive_remeshing
}//namespace CGAL

#endif //CGAL_INTERNAL_SMOOTH_VERTICES_H
