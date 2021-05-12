// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONTAINERS_FILLER_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONTAINERS_FILLER_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/boost/graph/internal/initialized_index_maps_helpers.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <boost/tuple/tuple.hpp>
#include <boost/unordered_set.hpp>
#include <boost/graph/graph_traits.hpp>

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
    : mesh(mesh_), vertices(vertices_), faces(nullptr)
  { }

  void operator()(const face_descriptor fd)
  {
    halfedge_descriptor hd = halfedge(fd, mesh);
    for(vertex_descriptor vd : vertices_around_face(hd, mesh)) {
      vertices.insert(vd);
    }

    if(faces != nullptr)
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

  void operator()(const face_descriptor fd)
  {
    for(vertex_descriptor vd : vertices_around_face(halfedge(fd, mesh), mesh))
    {
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

template <typename TriangleMesh, typename VertexIndexMap>
void fill_index_map_of_cc(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor bhd,
                          const TriangleMesh& mesh,
                          VertexIndexMap vimap)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

  std::vector<face_descriptor> CC_faces;

  // 'reserve' might cause a huge amount of memory to be used for a tiny CC,
  // but if this is a problem as a user, one could simply parameterize a Face_filtered_graph instead.
  CC_faces.reserve(num_faces(mesh));

  PMP::connected_component(face(opposite(bhd, mesh), mesh), mesh, std::back_inserter(CC_faces));

  // If all vertices are involved, avoid walking all the faces
  if(CC_faces.size() == faces(mesh).size())
  {
    BGL::internal::Index_map_initializer<VertexIndexMap, TriangleMesh> id_initializer;
    id_initializer(CGAL::internal_np::vertex_index, vimap, mesh);
  }
  else
  {
    for(vertex_descriptor v : vertices(mesh))
      put(vimap, v, -1);

    int index = 0;
    for(face_descriptor f : CC_faces)
      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh), mesh))
        if(get(vimap, v) == -1)
          put(vimap, v, index++);
  }
}

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_CONTAINERS_FILLER_H
