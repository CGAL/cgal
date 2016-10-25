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
//
//
// Author(s)     :

#ifndef CGAL_INTERNAL_SURFACE_MESH_PARAMETERIZATION_CONTAINERS_FILLER_H
#define CGAL_INTERNAL_SURFACE_MESH_PARAMETERIZATION_CONTAINERS_FILLER_H

#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>

#include <vector>

namespace CGAL {

namespace internal {

namespace Surface_mesh_parameterization {

// Custom output iterator that fills 'faces' and 'vertices' containers from a mesh
template<typename TriangleMesh>
class Containers_filler
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor
                                                                face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor
                                                              vertex_descriptor;

  const TriangleMesh& mesh;
  std::vector<face_descriptor>& faces;
  boost::unordered_set<vertex_descriptor>& vertices;

public:
  Containers_filler(const TriangleMesh& mesh,
                    std::vector<face_descriptor>& faces,
                    boost::unordered_set<vertex_descriptor>& vertices)
    : mesh(mesh), faces(faces), vertices(vertices)
  { }

  void operator()(face_descriptor fd)
  {
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd =
                                                              halfedge(fd,mesh);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, mesh)){
      vertices.insert(vd);
    }
    faces.push_back(fd);
  }
};

} // namespace Surface_mesh_parameterization

} // namespace internal

} // namespace CGAL

#endif // CGAL_INTERNAL_SURFACE_MESH_PARAMETERIZATION_CONTAINERS_FILLER_H
