// Copyright (c) 2016  GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_IO_FILE_OFF_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_IO_FILE_OFF_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/circulator.h>

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <cstddef>
#include <fstream>
#include <sstream>
#include <vector>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace IO {

template <typename TriangleMesh,
          typename VertexContainer,
          typename FacesContainer,
          typename VertexUVMap,
          typename VertexIndexMap>
void output_uvmap_to_off(const TriangleMesh& mesh,
                         const VertexContainer& vertices,
                         const FacesContainer& faces,
                         const VertexUVMap uvmap,
                         VertexIndexMap vimap,
                         std::ostream& os)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  os << "OFF" << '\n';
  os << vertices.size() << " " << faces.size() << " 0" << '\n';

  std::vector<std::size_t> renumbering_vector(vertices.size());
  std::size_t counter = 0;

  typename VertexContainer::const_iterator vit = vertices.begin();
  typename VertexContainer::const_iterator vend = vertices.end();
  for(; vit!=vend; ++vit){
    vertex_descriptor vd = *vit;
    os << get(uvmap, vd) << " 0" << '\n';

    // in case the vertices in 'vertices' are not in the same order as in vimap
    renumbering_vector[get(vimap, vd)] = counter++;
  }

  typename FacesContainer::const_iterator fit = faces.begin();
  typename FacesContainer::const_iterator fend = faces.end();
  for(; fit!=fend; ++fit){
    face_descriptor fd = *fit;
    halfedge_descriptor hd = halfedge(fd, mesh);

    os << "3";
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, mesh)){
      os << " " << renumbering_vector[get(vimap, vd)];
    }
    os << '\n';
  }
}

template <typename TriangleMesh,
          typename VertexUVMap>
void output_uvmap_to_off(const TriangleMesh& mesh,
                         typename boost::graph_traits<TriangleMesh>::halfedge_descriptor bhd,
                         const VertexUVMap uvmap,
                         std::ostream& os)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef boost::unordered_map<vertex_descriptor, std::size_t> Vertex_index_map;
  Vertex_index_map vium;
  boost::associative_property_map<Vertex_index_map> vimap(vium);

  boost::unordered_set<vertex_descriptor> vertices;
  std::vector<face_descriptor> faces;

  internal::Containers_filler<TriangleMesh> fc(mesh, vertices, &faces);
  Polygon_mesh_processing::connected_component(
                                    face(opposite(bhd, mesh), mesh),
                                    mesh,
                                    boost::make_function_output_iterator(fc));

  std::ostringstream out_vertices, out_faces;
  std::size_t vertices_counter = 0, faces_counter = 0;

  BOOST_FOREACH(vertex_descriptor vd, vertices){
    put(vimap, vd, vertices_counter++);
    out_vertices << get(uvmap, vd) << " 0" << '\n';
  }

  BOOST_FOREACH(face_descriptor fd, faces){
    halfedge_descriptor hd = halfedge(fd, mesh);

    out_faces << "3";
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, mesh)){
      out_faces << " " << get(vimap, vd);
    }
    out_faces << '\n';

    faces_counter++;
  }

  os << "OFF" << '\n';
  os << vertices_counter << " " << faces_counter << " 0" << '\n';
  os << out_vertices.str() << out_faces.str();
}

} // namespace IO

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_IO_FILE_OFF_H
