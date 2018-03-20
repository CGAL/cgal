// Copyright (c) 2009-2017 GeometryFactory (France).
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
// Author(s)     : Maxime Gimeno

#ifndef CGAL_FACETS_IN_COMPLEX_3_TO_TRIANGLE_MESH_H
#define CGAL_FACETS_IN_COMPLEX_3_TO_TRIANGLE_MESH_H


#include <CGAL/license/Mesh_3.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <map>
#include <boost/unordered_set.hpp>

namespace CGAL {
//! \ingroup PkgMesh_3Functions
//!
//! \brief builds a `TriangleMesh` from the surface facets, with a consistent orientation at the interface of two subdomains.
//!
//! This function exports the surface as a `TriangleMesh` and appends it to `graph`, using
//! `orient_polygon_soup()`.
//!
//! @tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//! @tparam TriangleMesh a model of `MutableFaceGraph` with an internal point property map. The point type should be compatible with the one used in `C3T3`.
//!
//! @param c3t3 an instance of `C3T3`.
//! @param graph an instance of `TriangleMesh`.
template<class C3T3, class TriangleMesh>
void facets_in_complex_3_to_triangle_mesh(const C3T3& c3t3, TriangleMesh& graph) //complexity nlogn(number of facets on surface)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Vertex_handle Vertex_handle;


  typedef typename C3T3::Triangulation Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Weighted_point Weighted_point;
  typedef std::vector<std::size_t> Polygon_3;
  std::vector< Point_3 > points;
  std::vector< Polygon_3 > polygons;

  //add vertices
  //used to set indices of vertices
  std::map<Vertex_handle, int> V;
  int inum = 0;
  boost::unordered_set<Vertex_handle> surface_vertices;
  //only take vertices on surface
  for(typename C3T3::Facets_in_complex_iterator //O(number of facets on surface)
      fit = c3t3.facets_in_complex_begin(),
      end = c3t3.facets_in_complex_end();
      fit != end; ++fit)
  {
    std::size_t s = fit->second;
    for(std::size_t i =1; i<4; ++i)
    {

      typename boost::unordered_set<Vertex_handle>::iterator new_vertex;
      bool is_new;
      boost::tie(new_vertex, is_new) = surface_vertices.insert(fit->first->vertex((s + i) % 4));
      if(is_new)
      {
        const Weighted_point& p = (*new_vertex)->point();
        points.push_back(
              Point_3(CGAL::to_double(p.x()),
                      CGAL::to_double(p.y()),
                      CGAL::to_double(p.z()))
              );
        V.insert(std::make_pair((*new_vertex), inum++));
      }
    }
  }

  //add faces
  for(typename C3T3::Facets_in_complex_iterator //O(number of facets on surface)
      fit = c3t3.facets_in_complex_begin(),
      end = c3t3.facets_in_complex_end();
      fit != end; ++fit)
  {
    std::size_t s = fit->second;
    std::size_t id0(V[(*fit).first->vertex((s + 1) % 4)]),
        id1(V[(*fit).first->vertex((s + 2) % 4)]),
        id2(V[(*fit).first->vertex((s + 3) % 4)]);
    Polygon_3 face;
    face.resize(3);
    //arbitrary choice of orientation. It works better like this on my personnal test data.
    face[0] = id0;
    face[2] = id1;
    face[1] = id2;
    polygons.push_back(face);
  }


  //first skip degenerate polygons
  std::vector<Polygon_3> valid_polygons;
  valid_polygons.reserve(polygons.size());
  BOOST_FOREACH(Polygon_3& polygon, polygons)
  {
    std::set<std::size_t> vids;
    bool to_remove=false;
    BOOST_FOREACH(std::size_t id, polygon)
    {
      if (!vids.insert(id).second){
        to_remove=true;
        break;
      }
    }
    if (!to_remove) valid_polygons.push_back(polygon);
  }
  if (valid_polygons.size()!=polygons.size())
    polygons.swap(valid_polygons);

  CGAL::Polygon_mesh_processing::
      orient_polygon_soup(points, polygons);
  CGAL_assertion(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons));
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, graph);
}//end c3t3_to_face_graph
}//end CGAL
#endif // CGAL_FACETS_IN_COMPLEX_3_TO_TRIANGLE_MESH_H
