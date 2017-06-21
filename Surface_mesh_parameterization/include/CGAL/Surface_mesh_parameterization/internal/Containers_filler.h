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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONTAINERS_FILLER_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONTAINERS_FILLER_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <boost/foreach.hpp>
#include "boost/tuple/tuple.hpp"
#include <boost/unordered_set.hpp>

#include <vector>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

// Custom output iterator that fills 'faces' and 'vertices' containers from a mesh
template<typename TriangleMesh,
         typename Vertex_set =
             boost::unordered_set<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>,
         typename Face_vector =
             std::vector<typename boost::graph_traits<TriangleMesh>::face_descriptor> >
class Containers_filler
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  const TriangleMesh& mesh;

  Vertex_set& vertices;
  Face_vector* faces;

public:
  Containers_filler(const TriangleMesh& mesh_,
                    Vertex_set& vertices_,
                    Face_vector* faces_)
    : mesh(mesh_), vertices(vertices_), faces(faces_)
  { }

  Containers_filler(const TriangleMesh& mesh_,
                    Vertex_set& vertices_)
    : mesh(mesh_), vertices(vertices_), faces(NULL)
  { }

  void operator()(face_descriptor fd)
  {
    halfedge_descriptor hd = halfedge(fd, mesh);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, mesh)) {
      vertices.insert(vd);
    }

    if(faces != NULL)
      faces->push_back(fd);
  }
};

template <typename Mesh, typename Map>
struct Index_map_filler
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::face_descriptor    face_descriptor;

  Index_map_filler(const Mesh& mesh, Map& map)
    : mesh(mesh), map(&map), index(0)
  { }

  void operator()(const face_descriptor& fd)
  {
    BOOST_FOREACH(vertex_descriptor vd,
                  vertices_around_face(halfedge(fd, mesh), mesh)) {
      typename Map::iterator it;
      bool new_element;
      boost::tie(it,new_element) = map->insert(std::make_pair(vd,1));
      if(new_element) {
        it->second = index++;
      }
    }
  }

  const Mesh& mesh;
  mutable Map* map;
  int index;
};

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONTAINERS_FILLER_H
