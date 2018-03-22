// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_BOOST_GRAPH_NEF_POLYHEDRON_TO_POLYGON_MESH_H
#define CGAL_BOOST_GRAPH_NEF_POLYHEDRON_TO_POLYGON_MESH_H

#include <CGAL/license/Nef_3.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/circulator.h>
#include <CGAL/Cartesian_converter.h>
#include <boost/unordered_map.hpp>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Kernel/global_functions_3.h>


namespace CGAL{

namespace nef_to_pm{

// Visitor used to collect and index vertices of a shell
template<class Nef_polyhedron, class Point_3, class Converter>
struct Shell_vertex_index_visitor
{
  typedef boost::unordered_map<typename Nef_polyhedron::Vertex_const_handle, std::size_t> Vertex_index_map;
  std::vector<Point_3>& points;
  Vertex_index_map vertex_indices;
  const Converter& converter;

  Shell_vertex_index_visitor(std::vector<Point_3>& points, const Converter& converter)
    :points(points), converter(converter)
  {}

  void visit(typename Nef_polyhedron::Vertex_const_handle vh)
  {
    std::pair<typename Vertex_index_map::iterator, bool> insert_res =
      vertex_indices.insert( std::make_pair(vh, points.size()) );
    if (insert_res.second)
      points.push_back( converter(vh->point()));
  }
  void visit(typename Nef_polyhedron::Halfedge_const_handle )
  {}
  void visit(typename Nef_polyhedron::Halffacet_const_handle )
  {}
  void visit(typename Nef_polyhedron::SHalfedge_const_handle )
  {}
  void visit(typename Nef_polyhedron::SHalfloop_const_handle )
  {}
  void visit(typename Nef_polyhedron::SFace_const_handle )
  {}
};

struct FaceInfo2
{
  FaceInfo2() : visited(false) {}
  bool visited;
};

//Visitor used to collect polygons
template <class Nef_polyhedron>
struct Shell_polygons_visitor
{
  typedef boost::unordered_map<typename Nef_polyhedron::Vertex_const_handle, std::size_t> Vertex_index_map;
  Vertex_index_map& vertex_indices;
  std::vector< std::vector<std::size_t> >& polygons;
  bool triangulate_all_faces;

  Shell_polygons_visitor(Vertex_index_map& vertex_indices,
                         std::vector< std::vector<std::size_t> >& polygons,
                         bool triangulate_all_faces)
    : vertex_indices( vertex_indices )
    , polygons(polygons)
    , triangulate_all_faces(triangulate_all_faces)
  {}

  std::size_t get_cycle_length( typename Nef_polyhedron::Halffacet_cycle_const_iterator hfc) const
  {
    typename Nef_polyhedron::SHalfedge_around_facet_const_circulator hc_start(hfc), hc_end(hc_start);
    std::size_t i=0;
    CGAL_For_all(hc_start,hc_end)
    {
      ++i;
    }
    return i;
  }

  void visit(typename Nef_polyhedron::Halffacet_const_handle opposite_facet)
  {
    bool is_marked=opposite_facet->incident_volume()->mark();

    CGAL_assertion(Nef_polyhedron::Infi_box::is_standard(opposite_facet->plane()));

    typename Nef_polyhedron::SHalfedge_const_handle se;
    typename Nef_polyhedron::Halffacet_cycle_const_iterator fc;

    typename Nef_polyhedron::Halffacet_const_handle f = opposite_facet->twin();

    if (cpp11::next(f->facet_cycles_begin())==f->facet_cycles_end())
    {
      std::size_t cycle_length = 3;
      if (triangulate_all_faces)
        cycle_length = get_cycle_length(f->facet_cycles_begin());
      switch(cycle_length)
      {
        case 4:
        {
          //collect vertices
          std::vector<typename Nef_polyhedron::Vertex_const_handle> v;
          fc = f->facet_cycles_begin();
          se = typename Nef_polyhedron::SHalfedge_const_handle(fc);
          CGAL_assertion(se!=0);
          typename Nef_polyhedron::SHalfedge_around_facet_const_circulator
            hc_start(se), hc_end(hc_start);
          // insert the vertex indices of the new polygon
          CGAL_For_all(hc_start,hc_end)
            v.push_back(hc_start->source()->center_vertex());

          if( CGAL::cross_product(v[1]->point()-v[0]->point(),v[1]->point()-v[2]->point()) *
              CGAL::cross_product(v[3]->point()-v[2]->point(),v[3]->point()-v[0]->point())
          >
             CGAL::cross_product(v[2]->point()-v[1]->point(),v[3]->point()-v[2]->point()) *
             CGAL::cross_product(v[0]->point()-v[3]->point(),v[1]->point()-v[0]->point()))
          {
            if (is_marked)
            {
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[0]]);
              polygons.back().push_back(vertex_indices[v[1]]);
              polygons.back().push_back(vertex_indices[v[2]]);
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[2]]);
              polygons.back().push_back(vertex_indices[v[0]]);
              polygons.back().push_back(vertex_indices[v[3]]);
            }
            else
            {
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[1]]);
              polygons.back().push_back(vertex_indices[v[0]]);
              polygons.back().push_back(vertex_indices[v[2]]);
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[0]]);
              polygons.back().push_back(vertex_indices[v[2]]);
              polygons.back().push_back(vertex_indices[v[3]]);
            }
          }
          else
          {
            if (is_marked)
            {
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[0]]);
              polygons.back().push_back(vertex_indices[v[1]]);
              polygons.back().push_back(vertex_indices[v[3]]);
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[3]]);
              polygons.back().push_back(vertex_indices[v[1]]);
              polygons.back().push_back(vertex_indices[v[2]]);
            }
            else
            {
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[0]]);
              polygons.back().push_back(vertex_indices[v[3]]);
              polygons.back().push_back(vertex_indices[v[1]]);
              polygons.push_back( std::vector<std::size_t>() );
              polygons.back().push_back(vertex_indices[v[1]]);
              polygons.back().push_back(vertex_indices[v[3]]);
              polygons.back().push_back(vertex_indices[v[2]]);
            }
          }
          return;
        }
        case 3:
        {
          // create a new polygon
          polygons.push_back( std::vector<std::size_t>() );
          fc = f->facet_cycles_begin();
          se = typename Nef_polyhedron::SHalfedge_const_handle(fc);
          CGAL_assertion(se!=0);
          typename Nef_polyhedron::SHalfedge_around_facet_const_circulator
            hc_start(se), hc_end(hc_start);
          // insert the vertex indices of the new polygon
          CGAL_For_all(hc_start,hc_end)
            polygons.back().push_back(vertex_indices[hc_start->source()->center_vertex()]);
          if (!is_marked)
            std::reverse(polygons.back().begin(), polygons.back().end());
          return;
        }
        default:
          break;
      }
    }

    // cases where a cdt is needed
    typedef typename Nef_polyhedron::Kernel Kernel;
    typedef Triangulation_2_projection_traits_3<Kernel>            P_traits;
    typedef Triangulation_vertex_base_with_info_2<std::size_t, P_traits> Vb;
    typedef Triangulation_face_base_with_info_2<FaceInfo2,P_traits>     Fbb;
    typedef Constrained_triangulation_face_base_2<P_traits,Fbb>          Fb;
    typedef Triangulation_data_structure_2<Vb,Fb>                       TDS;
    typedef Exact_predicates_tag                                       Itag;
    typedef Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>   CDT;

    P_traits p_traits(opposite_facet->plane().orthogonal_vector());
    CDT cdt(p_traits);

    // insert each connected component of the boundary of the face as
    // a polygonal constraints
    typename Nef_polyhedron::Halffacet_cycle_const_iterator
      hfc_start=f->facet_cycles_begin(), hfc_end=f->facet_cycles_end();
    CGAL_For_all(hfc_start, hfc_end)
    {
      fc=hfc_start;
      se = typename Nef_polyhedron::SHalfedge_const_handle(fc);
      CGAL_assertion(se!=0);
      typename Nef_polyhedron::SHalfedge_around_facet_const_circulator
        hc_start(se), hc_end(hc_start);
      // collect contour vertices
      std::vector< typename CDT::Vertex_handle > polygon;
      CGAL_For_all(hc_start,hc_end)
      {
        typename CDT::Vertex_handle vh=cdt.insert(hc_start->source()->center_vertex()->point());
        vh->info() = vertex_indices[hc_start->source()->center_vertex()];
        polygon.push_back(vh);
      }
      std::size_t nb_constraints = polygon.size();
      polygon.push_back(polygon.front());
      for(std::size_t i=0; i<nb_constraints; ++i)
        cdt.insert_constraint(polygon[i], polygon[i+1]);
    }
    //look for a triangle inside the domain of the face
    typename CDT::Face_handle fh = cdt.infinite_face();
    fh->info().visited=true;
    std::vector<typename CDT::Edge> queue;
    for (int i=0; i<3; ++i)
      queue.push_back(typename CDT::Edge(fh, i) );
    while(true)
    {
      typename CDT::Edge e = queue.back();
      queue.pop_back();
      e=cdt.mirror_edge(e);
      if (e.first->info().visited) continue;
      if (cdt.is_constrained(e))
      {
        queue.clear();
        queue.push_back(e);
        break;
      }
      else
      {
        for(int i=1; i<3; ++i)
        {
          typename CDT::Edge candidate(e.first, (e.second+i)%3);
          if (!candidate.first->neighbor(candidate.second)->info().visited)
            queue.push_back( candidate );
        }
        e.first->info().visited=true;
      }
    }
    // now extract triangles inside the face
    while(!queue.empty())
    {
      typename CDT::Edge e = queue.back();
      queue.pop_back();
      if (e.first->info().visited) continue;
      e.first->info().visited=true;
      polygons.resize(polygons.size()+1);
      if (is_marked)
        for (int i=2; i>=0; --i)
          polygons.back().push_back(e.first->vertex(i)->info());
      else
        for (int i=0; i<3; ++i)
          polygons.back().push_back(e.first->vertex(i)->info());

      for(int i=1; i<3; ++i)
      {
        typename CDT::Edge candidate(e.first, (e.second+i)%3);
        if (!cdt.is_constrained(candidate) &&
            !candidate.first->neighbor(candidate.second)->info().visited)
        {
          queue.push_back( cdt.mirror_edge(candidate) );
        }
      }
    }
  }

  void visit(typename Nef_polyhedron::SFace_const_handle)
  {}
  void visit(typename Nef_polyhedron::Halfedge_const_handle)
  {}
  void visit(typename Nef_polyhedron::Vertex_const_handle)
  {}
  void visit(typename Nef_polyhedron::SHalfedge_const_handle)
  {}
  void visit(typename Nef_polyhedron::SHalfloop_const_handle)
  {}
};

template <class Point_3, class Nef_polyhedron, class Converter>
void collect_polygon_mesh_info(
  std::vector<Point_3>& points,
  std::vector< std::vector<std::size_t> >& polygons,
  Nef_polyhedron& nef,
  typename Nef_polyhedron::Shell_entry_const_iterator shell,
  const Converter& converter,
  bool triangulate_all_faces)
{
  // collect points and set vertex indices
  Shell_vertex_index_visitor<Nef_polyhedron, Point_3, Converter> vertex_index_visitor(points, converter);
  nef.visit_shell_objects(typename Nef_polyhedron::SFace_const_handle(shell), vertex_index_visitor);

  // collect polygons
  Shell_polygons_visitor<Nef_polyhedron> polygon_visitor(
    vertex_index_visitor.vertex_indices,
    polygons,
    triangulate_all_faces);
  nef.visit_shell_objects(typename Nef_polyhedron::SFace_const_handle(shell), polygon_visitor);
}

} //end of namespace nef_to_pm


template <class Nef_polyhedron, class Polygon_mesh>
void convert_nef_polyhedron_to_polygon_mesh(const Nef_polyhedron& nef, Polygon_mesh& pm, bool triangulate_all_faces = false)
{
  typedef typename Nef_polyhedron::Point_3 Point_3;
  typedef typename boost::property_traits<typename boost::property_map<Polygon_mesh, vertex_point_t>::type>::value_type PM_Point;
  typedef typename Kernel_traits<PM_Point>::Kernel PM_Kernel;
  typedef typename Kernel_traits<Point_3>::Kernel Nef_Kernel;
  typedef Cartesian_converter<Nef_Kernel, PM_Kernel> Converter;

  typename Nef_polyhedron::Volume_const_iterator vol_it = nef.volumes_begin(),
                                                 vol_end = nef.volumes_end();
  if ( Nef_polyhedron::Infi_box::extended_kernel() ) ++vol_it; // skip Infi_box
  CGAL_assertion ( vol_it != vol_end );
  ++vol_it; // skip unbounded volume

  std::vector<PM_Point> points;
  std::vector< std::vector<std::size_t> > polygons;
  Converter to_inexact;
  for (;vol_it!=vol_end;++vol_it)
    nef_to_pm::collect_polygon_mesh_info<PM_Point>(points,
                                                  polygons,
                                                  nef,
                                                  vol_it->shells_begin(),
                                                  to_inexact,
                                                  triangulate_all_faces);

  Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, pm);
}

} //end of namespace CGAL



#endif // CGAL_BOOST_GRAPH_NEF_POLYHEDRON_TO_POLYGON_MESH_H
