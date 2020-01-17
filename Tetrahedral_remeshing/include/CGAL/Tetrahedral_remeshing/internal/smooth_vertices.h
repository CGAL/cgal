// Copyright (c) 2017  GeometryFactory (France).
// All rights reserved.

#ifndef CGAL_INTERNAL_SMOOTH_VERTICES_H
#define CGAL_INTERNAL_SMOOTH_VERTICES_H

#include <CGAL/Vector_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <CGAL/number_utils.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <limits>
#include <vector>
#include <cmath>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
  template<typename SubdomainIndex>
  std::pair<SubdomainIndex, SubdomainIndex>
    make_surface_index(const SubdomainIndex& s1,
                       const SubdomainIndex& s2)
  {
    if (s1 < s2)
      return std::make_pair(s1, s2);
    else
      return std::make_pair(s2, s1);
  }

  template<typename Gt>
  CGAL::Vector_3<Gt> project_on_tangent_plane(const CGAL::Point_3<Gt>& gi,
                                              const CGAL::Point_3<Gt>& pi,
                                              const CGAL::Vector_3<Gt>& normal)
  {
    typename Gt::Construct_vector_3
      vec = Gt().construct_vector_3_object();
    typename Gt::Construct_scaled_vector_3
      scale = Gt().construct_scaled_vector_3_object();
    return scale(normal, CGAL::scalar_product(normal, vec(gi, pi)));
  }

  template<typename C3t3>
  typename C3t3::Triangulation::Geom_traits::Vector_3
  compute_vertex_normal(const typename C3t3::Vertex_handle v,
                        const C3t3& c3t3)
  {
    typedef typename C3t3::Subdomain_index            Subdomain_index;
    typedef typename C3t3::Triangulation::Facet       Facet;
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;
    typedef typename C3t3::Triangulation::Geom_traits Gt;
    typedef typename Gt::Vector_3                     Vector_3;
    typedef std::pair< Subdomain_index, Subdomain_index> Surface_index;

    typename Gt::Construct_opposite_vector_3
      opp = c3t3.triangulation().geom_traits().construct_opposite_vector_3_object();
    typename Gt::Construct_sum_of_vectors_3
      sum = c3t3.triangulation().geom_traits().construct_sum_of_vectors_3_object();
    typename Gt::Construct_scaled_vector_3
      scale = c3t3.triangulation().geom_traits().construct_scaled_vector_3_object();
    typename Gt::Compute_squared_length_3
      sqlen = c3t3.triangulation().geom_traits().compute_squared_length_3_object();

    std::vector<Facet> facets;
    c3t3.triangulation().incident_facets(v, std::back_inserter(facets));

    Vector_3 normal = CGAL::NULL_VECTOR;

    for (Facet f : facets)
    {
      Cell_handle ch = f.first;
      Cell_handle n_ch = f.first->neighbor(f.second);

      Subdomain_index si = ch->subdomain_index();
      Subdomain_index si_mirror = n_ch->subdomain_index();

      if (si != si_mirror
        || c3t3.triangulation().is_infinite(ch)
        || c3t3.triangulation().is_infinite(n_ch))
      {
        Surface_index surf_i = make_surface_index(si, si_mirror);

        Vector_3 n = facet_normal(c3t3.triangulation(), f);

        if (si < si_mirror || c3t3.triangulation().is_infinite(ch))
          n = opp(n);

        normal = sum(normal, n);
      }
    }

    if (normal != CGAL::NULL_VECTOR)
      return scale(normal, 1. / CGAL::sqrt(sqlen(normal)));
    else
      return CGAL::NULL_VECTOR;
  }

  template<typename C3t3, typename VertexNormalsMap>
  void compute_vertices_normals(const C3t3& c3t3,
                                VertexNormalsMap& normals_map)
  {
    typedef typename C3t3::Triangulation        Tr;
    typedef typename C3t3::Cell_handle          Cell_handle;
    typedef typename C3t3::Vertex_handle        Vertex_handle;
    typedef typename C3t3::Subdomain_index      Subdomain_index;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
    typedef typename Tr::Geom_traits::Vector_3  Vector_3;
    typedef std::pair<Subdomain_index, Subdomain_index> Surface_index;

    const Tr& tr = c3t3.triangulation();

    typename Tr::Geom_traits::Construct_opposite_vector_3
      opp = tr.geom_traits().construct_opposite_vector_3_object();

    for (Finite_facets_iterator fit = tr.finite_facets_begin();
         fit != tr.finite_facets_end(); ++fit)
    {
      Cell_handle ch = fit->first;
      Cell_handle n_ch = fit->first->neighbor(fit->second);

      Subdomain_index si = ch->subdomain_index();
      Subdomain_index si_mirror = n_ch->subdomain_index();

      if (si != si_mirror || tr.is_infinite(ch) || tr.is_infinite(n_ch))
      {
        Surface_index surf_i = make_surface_index(si, si_mirror);
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
        Surface_index surf_i = make_surface_index(si, si_mirror);

        Vector_3 n = CGAL::Tetrahedral_remeshing::facet_normal(tr, *fit);

        if (si < si_mirror || tr.is_infinite(ch))
          n = opp(n);

        for (int i = 0; i < 3; ++i)
        {
          Vertex_handle v_id = fit->first->vertex(indices(fit->second, i));
          Vector_3& v_n = normals_map[v_id][surf_i];
          v_n = v_n + n;
        }
      }
    }

    //normalize the computed normals
    for (typename VertexNormalsMap::iterator vnm_it = normals_map.begin();
         vnm_it != normals_map.end(); ++vnm_it)
    {
      //mapped_type is map<Surface_index, Vector_3>
      for (typename VertexNormalsMap::mapped_type::iterator it = vnm_it->second.begin();
           it != vnm_it->second.end(); ++it)
      {
        Vector_3& n = it->second;
        n = n / CGAL::sqrt(n*n);
      }
    }
  }


  template<typename C3t3>
  std::pair<typename C3t3::Subdomain_index, typename C3t3::Subdomain_index>
  surface_index(const typename C3t3::Vertex_handle v, const C3t3& c3t3)
  {
    typedef typename C3t3::Triangulation::Facet       Facet;
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;
    typedef typename C3t3::Subdomain_index            Subdomain_index;

    std::vector<Facet> facets;
    c3t3.triangulation().incident_facets(v, std::back_inserter(facets));

    for (Facet f : facets)
    {
      Cell_handle ch = f.first;
      Cell_handle n_ch = f.first->neighbor(f.second);

      Subdomain_index si = ch->subdomain_index();
      Subdomain_index si_mirror = n_ch->subdomain_index();

      if (si != si_mirror
        || c3t3.triangulation().is_infinite(ch)
        || c3t3.triangulation().is_infinite(n_ch))
      {
        return make_surface_index(si, si_mirror);
      }
    }
    CGAL_assertion(false);
    return make_surface_index(0, 0);
  }

  template<typename C3t3>
  const boost::unordered_set<typename C3t3::Subdomain_index>
  subdomain_indices(const typename C3t3::Vertex_handle v, const C3t3& c3t3)
  {
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

    std::vector<Cell_handle> cells;
    c3t3.triangulation().incident_cells(v, std::back_inserter(cells));

    boost::unordered_set<typename C3t3::Subdomain_index> res;
    for (Cell_handle c : cells)
    {
      if (c3t3.is_in_complex(c))
        res.insert(c->subdomain_index());
    }
    return res;
  }

  template<typename C3t3, typename FMLSVector, typename SurfaceIndexMap>
  void createMLSSurfaces(const C3t3& c3t3,
                         FMLSVector& subdomain_FMLS,
                         SurfaceIndexMap& subdomain_FMLS_indices)
  {
    typedef typename C3t3::Subdomain_index       Subdomain_index;
    typedef typename C3t3::Triangulation         Tr;
    typedef typename Tr::Geom_traits             Gt;
    typedef typename Tr::Edge                    Edge;
    typedef typename Tr::Vertex_handle           Vertex_handle;
    typedef typename Gt::Point_3                 Point_3;
    typedef typename Gt::Vector_3                Vector_3;

    typedef std::pair<Subdomain_index, Subdomain_index> Surface_index;

    const Tr& tr = c3t3.triangulation();

    SurfaceIndexMap current_subdomain_FMLS_indices;

    SurfaceIndexMap subdomain_sample_numbers;

    //Count the number of vertices for each boundary surface (i.e. one per label)
    for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
      vit != tr.finite_vertices_end(); ++vit)
    {
      if (c3t3.in_dimension(vit) == 2)
      {
        const boost::unordered_set<Subdomain_index>& v_subdomain_indices = subdomain_indices(vit, c3t3);
        if (v_subdomain_indices.size() == 2)
        {
          boost::unordered_set<Subdomain_index>::const_iterator si_it = v_subdomain_indices.cbegin();
          Subdomain_index s1 = *si_it;
          ++si_it;
          Subdomain_index s2 = *si_it;

          subdomain_sample_numbers[make_surface_index(s1, s2)]++;
        }
      }
    }

    std::vector< float* > pns;

    int count = 0;
    //Memory allocation for the point plus normals of the point samples
    for (typename SurfaceIndexMap::iterator it = subdomain_sample_numbers.begin();
      it != subdomain_sample_numbers.end(); ++it)
    {
      current_subdomain_FMLS_indices[it->first] = count;
      pns.push_back(new float[it->second * 6]);
      count++;
    }

    boost::unordered_map<Vertex_handle, 
      boost::unordered_map<Surface_index, Vector_3> > vertices_normals;
    compute_vertices_normals(c3t3, vertices_normals);

    std::vector<int> current_v_count(count, 0);
    std::vector<double> point_spacing(count, 0);
    std::vector<int> point_spacing_count(count, 0);

    //Allocation of the PN
    for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end(); ++vit)
    {
      boost::unordered_set<Subdomain_index> vertices_subdomain_indices
        = subdomain_indices(vit, c3t3);
      if (vertices_subdomain_indices.size() == 2)
      {
        Subdomain_index s1 = *(vertices_subdomain_indices.begin());
        Subdomain_index s2 = *(++vertices_subdomain_indices.begin());

        Surface_index surf_i = make_surface_index(s1, s2);

        int fmls_id = current_subdomain_FMLS_indices[surf_i];

        Point_3& point = vit->point();

        pns[fmls_id][6 * current_v_count[fmls_id]] = point.x();
        pns[fmls_id][6 * current_v_count[fmls_id] + 1] = point.y();
        pns[fmls_id][6 * current_v_count[fmls_id] + 2] = point.z();

        Vector_3& normal = vertices_normals[vit][surf_i];

        pns[fmls_id][6 * current_v_count[fmls_id] + 3] = normal.x();
        pns[fmls_id][6 * current_v_count[fmls_id] + 4] = normal.y();
        pns[fmls_id][6 * current_v_count[fmls_id] + 5] = normal.z();

        current_v_count[fmls_id]++;
      }
    }

    typedef std::pair<Vertex_handle, Vertex_handle> Edge_VV;
    typedef std::map<Edge_VV, unsigned int> EdgeMapIndex;
    EdgeMapIndex edgeMap;

    for (typename C3t3::Facet_iterator fit = c3t3.facets_begin();
         fit != c3t3.facets_end(); ++fit)
    {
      for (int i = 0; i < 2; i++)
      {
        for (int j = i + 1; j < 3; j++)
        {
          Edge edge(fit->first, indices(fit->second, i), indices(fit->second, j));

          Vertex_handle vh0 = edge.first->vertex(edge.second);
          Vertex_handle vh1 = edge.first->vertex(edge.third);
          Edge_VV evv = make_vertex_pair(vh1, vh0);

          if ( subdomain_indices(vh0, c3t3).size() == 2
            && subdomain_indices(vh1, c3t3).size() == 2
            && edgeMap.find(evv) == edgeMap.end())
          {
            edgeMap[evv] = 0;
            Surface_index surf_i = make_surface_index(
                                            fit->first->subdomain_index(),
                                            fit->first->neighbor(fit->second)->subdomain_index());
            int fmls_id = current_subdomain_FMLS_indices[surf_i];

            point_spacing[fmls_id] += CGAL::sqrt(tr.segment(edge).squared_length());
            point_spacing_count[fmls_id] ++;
          }
        }
      }
    }

    int nb_of_mls_to_create = 0;
    double average_point_spacing = 0;

    //Cretaing the actual MLS surfaces
    for (SurfaceIndexMap::iterator it = current_subdomain_FMLS_indices.begin();
      it != current_subdomain_FMLS_indices.end(); ++it)
    {
      if (current_v_count[it->second] > 3)
      {
        nb_of_mls_to_create++;

        double current_point_spacing = point_spacing[it->second] / point_spacing_count[it->second];
        point_spacing[it->second] = current_point_spacing;

        average_point_spacing += current_point_spacing;
      }
    }

    average_point_spacing = average_point_spacing / nb_of_mls_to_create;

    subdomain_FMLS.resize(nb_of_mls_to_create, FMLS());

    count = 0;
    //Cretaing the actual MLS surfaces
    for (SurfaceIndexMap::iterator it = current_subdomain_FMLS_indices.begin();
      it != current_subdomain_FMLS_indices.end(); ++it)
    {
      if (current_v_count[it->second] > 3)
      {
        double current_point_spacing = point_spacing[it->second];

        //subdomain_FMLS[count].toggleHermite(true);
        subdomain_FMLS[count].setPN(pns[it->second], current_v_count[it->second], current_point_spacing);
        //  subdomain_FMLS[count].toggleHermite(true);
        subdomain_FMLS_indices[it->first] = count;

        count++;
      }
      else {
        std::cout << "Problem of number for MLS : " << current_v_count[it->second] << std::endl;
      }
    }
  }

  template<typename C3t3, typename FMLSVector, typename SurfaceIndexMap>
  bool project(const typename C3t3& c3t3,
               const typename C3t3::Vertex_handle& v,
               typename C3t3::Triangulation::Geom_traits::Vector_3& gi,
               typename C3t3::Triangulation::Geom_traits::Vector_3& projected_point,
               FMLSVector& subdomain_FMLS,
               SurfaceIndexMap& subdomain_FMLS_indices)
  {
    typedef typename C3t3::Triangulation::Geom_traits::Point_3  Point_3;
    typedef typename C3t3::Triangulation::Geom_traits::Vector_3 Vector_3;
    typedef typename C3t3::Subdomain_index                      Subdomain_index;

    std::pair<Subdomain_index, Subdomain_index> si = surface_index(v, c3t3);

    if (subdomain_FMLS_indices.find(si) == subdomain_FMLS_indices.end())
      return false;

    Point_3 point(gi.x(), gi.y(), gi.z());

    Vector_3 res_normal;
    Point_3 result(point);

    FMLS& fmls = subdomain_FMLS[subdomain_FMLS_indices[si]];

    int it_nb = 0;

    float epsilon = fmls.getPNScale() / 1000.;
    float sq_eps = epsilon * epsilon;

    do
    {
      point = result;

      fmls.fastProjectionCPU(point, result, res_normal);

      if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2])) {
        std::cout << "MLS error detected si size " << si.first << " - " << si.second
                  << " : " << fmls.getPNSize() << std::endl;
        return false;
      }

      it_nb++;

    } while (CGAL::squared_distance(result, point) > sq_eps && it_nb < 5);

    projected_point = Vector_3(result[0], result[1], result[2]);

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

  template<typename C3t3>
  bool project(const C3t3& c3t3,
    const typename C3t3::Vertex_handle v,
    typename C3t3::Triangulation::Geom_traits::Vector_3& gi,
    typename C3t3::Triangulation::Geom_traits::Vector_3& projected_point )
  {
    typedef typename C3t3::Subdomain_index Subdomain_index;

    const std::pair<Subdomain_index, Subdomain_index> si = surface_index(v);
    if( subdomain_FMLS_indices.find( si ) == subdomain_FMLS_indices.end() )
        return false;

    Vec3Df point( gi.x(), gi.y(), gi.z() );
    if( isnan(point[0]) || isnan(point[1]) || isnan(point[2]) ){
        std::cout << "Initial point error " << point << std::endl;
        return false;
    }

    Vec3Df res_normal;
    Vec3Df result(point);

    FMLS & fmls = subdomain_FMLS[ subdomain_FMLS_indices[ si ] ];

    int it_nb = 0;

    float epsilon = fmls.getPNScale() /1000.;

    do{
        point = result;

        fmls.fastProjectionCPU( point, result, res_normal );

        if( isnan(result[0]) || isnan(result[1]) || isnan(result[2]) ){
            std::cout << "MLS error detected si size " << si.first << " - " << si.second << " : " << fmls.getPNSize() << std::endl;
            return false;
        }

        it_nb++;

    }while ( (result - point).getLength() > epsilon && it_nb < 5 );

    projected_point = K::Vector_3( result[0], result[1], result[2] );

    return true;
}

  template<typename C3T3>
  typename C3T3::Triangulation::Geom_traits::Vector_3
    move_3d(typename C3T3::Vertex_handle v, const C3T3& c3t3)
  {
    typedef typename C3T3::Edge Edge;
    typedef typename C3T3::Vertex_handle Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits Gt;
    typedef typename Gt::Vector_3 Vector_3;

    const Gt& gt = c3t3.triangulation().geom_traits();

    Vector_3 move = CGAL::NULL_VECTOR;

    std::vector<Edge> edges;
    c3t3.triangulation().incident_edges(v, std::back_inserter(edges));

    if (edges.empty())
      return move;

    typename Gt::Construct_vector_3 vec
      = gt.construct_vector_3_object();
    typename Gt::Construct_sum_of_vectors_3 sum
      = gt.construct_sum_of_vectors_3_object();

    BOOST_FOREACH(Edge e, edges)
    {
      Vertex_handle ve = (e.first->vertex(e.second) != v)
                        ? e.first->vertex(e.second)
                        : e.first->vertex(e.third);
      move = sum(move, Vector_3(CGAL::ORIGIN,  point(ve->point())));
    }

    typename Gt::Construct_scaled_vector_3 scale
      = gt.construct_scaled_vector_3_object();
    return scale(move, 1. / edges.size());
  }

  template<typename C3T3, typename CellSelector>
  typename C3T3::Triangulation::Geom_traits::Vector_3
    move_2d(typename C3T3::Vertex_handle v,
            const C3T3& c3t3,
            const typename C3T3::Subdomain_index& imaginary_index,
            const CellSelector cell_selector)
  {
    typedef typename C3T3::Subdomain_index            Subdomain_index;
    typedef typename C3T3::Edge                       Edge;
    typedef typename C3T3::Vertex_handle              Vertex_handle;
    typedef typename C3T3::Triangulation::Geom_traits Gt;
    typedef typename Gt::Vector_3                     Vector_3;
    typedef typename Gt::Point_3                      Point_3;

    const Gt& gt = c3t3.triangulation().geom_traits();

    std::vector<Edge> edges;
    c3t3.triangulation().incident_edges(v, std::back_inserter(edges));

    Vector_3 move = CGAL::NULL_VECTOR;
    if (edges.empty())
      return move;

    typename Gt::Construct_vector_3 vec
      = gt.construct_vector_3_object();
    typename Gt::Construct_sum_of_vectors_3 sum
      = gt.construct_sum_of_vectors_3_object();

    std::size_t nbe = 0;
    for(Edge e : edges)
    {
      if(!c3t3.is_in_complex(e) && is_boundary(c3t3, e, cell_selector))
      {
        Vertex_handle ve = (e.first->vertex(e.second) != v)
                          ? e.first->vertex(e.second)
                          : e.first->vertex(e.third);
        move = sum(move, vec(CGAL::ORIGIN, point(ve->point())));
        ++nbe;
      }
    }

    if (nbe > 0)
    {
      // WIP in this section

      typedef std::pair<Subdomain_index, Subdomain_index> Surface_index;
      typedef std::map<Surface_index, unsigned int/*, Compare*/> SurfaceIndexMap;
      SurfaceIndexMap subdomain_FMLS_indices;
      std::vector< FMLS > subdomain_FMLS;
      //createMLSSurfaces(c3t3, subdomain_FMLS, subdomain_FMLS_indices);

      typename Gt::Construct_scaled_vector_3 scale
        = gt.construct_scaled_vector_3_object();
      move = scale(move, 1. / nbe);

      const Point_3 current_pos = point(v->point());
      const Point_3 smoothed_position = current_pos + move;
      Point_3 final_position = CGAL::ORIGIN;

      Vector_3 normal = compute_vertex_normal(v, c3t3);

      Vector_3 normal_projection = project_on_tangent_plane(
                                     smoothed_position, //smoothed position
                                     current_pos,       //current position
                                     normal);

      Vector_3 mls_projection;
      if (project(c3t3, v, normal_projection, mls_projection,
                  subdomain_FMLS, subdomain_FMLS_indices))
        move = move + mls_projection;
      else
        move = move + normal_projection;

      return move;
    }
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
    typedef typename C3T3::Triangulation::Geom_traits Gt;
    typedef typename Gt::Vector_3 Vector_3;

    const Gt& gt = c3t3.triangulation().geom_traits();

    Vector_3 move = CGAL::NULL_VECTOR;

    std::vector<Edge> edges;
    c3t3.triangulation().incident_edges(v, std::back_inserter(edges));

    if (edges.empty())
      return move;

    typename Gt::Construct_vector_3 vec
      = gt.construct_vector_3_object();
    typename Gt::Construct_sum_of_vectors_3 sum
      = gt.construct_sum_of_vectors_3_object();

    std::size_t nbe = 0;
    BOOST_FOREACH(Edge e, edges)
    {
      if (!c3t3.is_in_complex(e))
        continue;

      Vertex_handle ve = (e.first->vertex(e.second) != v)
                      ? e.first->vertex(e.second)
                      : e.first->vertex(e.third);

      move = sum(move, vec(CGAL::ORIGIN, point(ve->point())));
      ++nbe;
    }

    if (nbe == 2)
    {
      typename Gt::Construct_scaled_vector_3 scale
        = gt.construct_scaled_vector_3_object();
      return scale(move, 0.5);
    }
    else
      return CGAL::NULL_VECTOR;
  }

  template<typename C3T3, typename CellSelector>
  void smooth_vertices_new(C3T3& c3t3,
                           const typename C3T3::Subdomain_index& imaginary_index,
                           const bool protect_boundaries,
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
        if (protect_boundaries)
          break;

        smoothing_vecs[vertex_id.at(vit)] = move_2d(vit, c3t3, imaginary_index, cell_selector);
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        if (smoothing_vecs[vertex_id.at(vit)] != CGAL::NULL_VECTOR)
          ofs_2d << "2 " << vit->point()
            << " " << (CGAL::ORIGIN + smoothing_vecs[vertex_id.at(vit)]) << std::endl;
#endif
        break;

      case 1:
        if (protect_boundaries)
          break;

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


}//namespace internal
}//namespace Tetrahedral_adaptive_remeshing
}//namespace CGAL

#endif //CGAL_INTERNAL_SMOOTH_VERTICES_H
