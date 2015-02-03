// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     :  Sebastien Loriot


#ifndef CGAL_INTERNAL_POLYGON_MESH_SLICER_TRAVERSAL_TRAITS_H
#define CGAL_INTERNAL_POLYGON_MESH_SLICER_TRAVERSAL_TRAITS_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Triangle_3_Ray_3_do_intersect.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/enum.h>

namespace CGAL{
namespace Polygon_mesh_slicer_{

template <typename AL_graph,
          typename TriangleMesh,
          typename VertexPointPmap,
          typename AABBTraits,
          class Traits>
class Traversal_traits
{
/// typedefs
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor     edge_descriptor;
  typedef typename AL_graph::vertex_descriptor                       AL_vertex_descriptor;
  typedef std::pair<const vertex_descriptor, AL_vertex_descriptor>            Vertex_pair;
  typedef std::map<vertex_descriptor, AL_vertex_descriptor>                  Vertices_map;
/// container filled by `intersection()`
  std::set<edge_descriptor>& m_all_coplanar_edges;
  Vertices_map& m_vertices;
  std::vector<edge_descriptor>& m_iedges;
/// data members
  TriangleMesh& m_tmesh;
  const VertexPointPmap& m_vpmap;
  const AABBTraits& m_aabb_traits;
  const Traits& m_traits;
  const typename AL_graph::vertex_descriptor null_vertex;
/// predicates
  typename Traits::Oriented_side_3 oriented_side_3;
  typename Traits::Do_intersect_3 do_intersect_3;

public:

  Traversal_traits( std::set<edge_descriptor>& all_coplanar_edges,
                    std::vector<edge_descriptor>& iedges,
                    Vertices_map& vertices,
                    TriangleMesh& tmesh,
                    const VertexPointPmap& vpmap,
                    const AABBTraits& aabb_traits,
                    const Traits& traits)
    : m_all_coplanar_edges(all_coplanar_edges)
    , m_vertices(vertices)
    , m_iedges(iedges)
    , m_tmesh(tmesh)
    , m_vpmap(vpmap)
    , m_aabb_traits(aabb_traits)
    , m_traits(traits)
    , null_vertex( boost::graph_traits<AL_graph>::null_vertex() )
    , oriented_side_3( m_traits.oriented_side_3_object() )
    , do_intersect_3( m_traits.do_intersect_3_object() )
  {}

  bool go_further() const { return true; }

  void intersection(const typename Traits::Plane_3& plane, const typename AABBTraits::Primitive& primitive)
  {
    typename boost::graph_traits<TriangleMesh>::edge_descriptor ed = primitive.id();

    Oriented_side src = oriented_side_3(plane, get(m_vpmap, source(ed,m_tmesh)) );
    Oriented_side tgt = oriented_side_3(plane, get(m_vpmap, target(ed,m_tmesh)) );

    if (src==ON_ORIENTED_BOUNDARY)
    {
      if (tgt==ON_ORIENTED_BOUNDARY)
        m_all_coplanar_edges.insert(ed);
      else
        m_vertices.insert( Vertex_pair (source(ed,m_tmesh), null_vertex) );
    }
    else{
      if (tgt==ON_ORIENTED_BOUNDARY)
          m_vertices.insert( Vertex_pair (target(ed,m_tmesh), null_vertex) );
      else
        if (src!=tgt)
          return m_iedges.push_back(ed);
    }
  }

  template<class Node>
  bool do_intersect(const typename Traits::Plane_3& plane, const Node& node) const
  {
    return do_intersect_3(plane, node.bbox());
  }
};

} } //end of namespace CGAL::Polygon_mesh_slicer_

#endif // CGAL_INTERNAL_POLYGON_MESH_SLICER_TRAVERSAL_TRAITS_H
