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
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H

namespace CGAL {
namespace Corefinement {

template < class TriangleMesh,
           class VertexPointMap,
           class Node_id,
           class Node_vector,
           class CDT,
           class NewNodeVisitor,
           class NewFaceVisitor  >
void
triangulate_a_face(
  typename boost::graph_traits<TriangleMesh>::face_descriptor current_face,
  TriangleMesh& tm,
  const Node_vector& nodes,
  const std::vector<Node_id>& node_ids,
  typename std::vector<typename boost::graph_traits<TriangleMesh>
                            ::vertex_descriptor>& node_id_to_vertex,
  std::map<std::pair<Node_id,Node_id>,
           typename boost::graph_traits<TriangleMesh>
                ::halfedge_descriptor>& edge_to_hedge,
  const CDT& cdt,
  const VertexPointMap& vpm,
  NewNodeVisitor& new_node_visitor,
  NewFaceVisitor& new_face_visitor)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;

  //insert the intersection point interior to the face inside the polyhedron and
  //save their Polyhedron::vertex_handle
  BOOST_FOREACH(Node_id node_id, node_ids)
  {
    vertex_descriptor v=add_vertex(tm);
    new_node_visitor.new_vertex_added(node_id, v, tm);
    put(vpm, v, nodes[node_id]);
    CGAL_assertion(node_id_to_vertex.size()>node_id);
    node_id_to_vertex[node_id]=v;
  }

  //insert the new halfedge and set their incident vertex

  //grab edges that are not on the convex hull (these have already been created)
  for (typename CDT::Finite_edges_iterator it=cdt.finite_edges_begin();
                                           it!=cdt.finite_edges_end(); ++it)
  {
    typename CDT::Vertex_handle cdt_v0=it->first->vertex( cdt.ccw(it->second) );
    typename CDT::Vertex_handle cdt_v1=it->first->vertex( cdt.cw(it->second) );

    // consider edges not on the convex hull (not on the boundary of the face)
    // and create the corresponding halfedges
    if ( !cdt.is_infinite(it->first->vertex(it->second)) &&
         !cdt.is_infinite(cdt.mirror_vertex(it->first,it->second)) )
    {
      edge_descriptor e=add_edge(tm);
      halfedge_descriptor h=halfedge(e,tm), h_opp=opposite(h,tm);

      Node_id i0=cdt_v0->info(), i1=cdt_v1->info();
      CGAL_assertion( node_id_to_vertex[i0]!=GT::null_vertex());
      CGAL_assertion( node_id_to_vertex[i1]!=GT::null_vertex());

      vertex_descriptor v0=node_id_to_vertex[i0], v1=node_id_to_vertex[i1];

      set_target(h,v0,tm);
      set_target(h_opp,v1,tm);
      set_halfedge(v0,h,tm);
      set_halfedge(v1,h_opp,tm);

      edge_to_hedge[std::make_pair(i0,i1)]=h_opp;
      edge_to_hedge[std::make_pair(i1,i0)]=h;
    }
  }

  //grab triangles.
  new_face_visitor.before_subface_creations(current_face,tm);
  for (typename CDT::Finite_faces_iterator it=cdt.finite_faces_begin(),
                                           it_end=cdt.finite_faces_end();;)
  {
    typename CDT::Vertex_handle cdt_v0=it->vertex(0);
    typename CDT::Vertex_handle cdt_v1=it->vertex(1);
    typename CDT::Vertex_handle cdt_v2=it->vertex(2);

    Node_id i0=cdt_v0->info(), i1=cdt_v1->info(), i2=cdt_v2->info();

    CGAL_assertion(edge_to_hedge.count(std::make_pair(i0,i1)));
    CGAL_assertion(edge_to_hedge.count(std::make_pair(i1,i2)));
    CGAL_assertion(edge_to_hedge.count(std::make_pair(i2,i0)));

    halfedge_descriptor h01=edge_to_hedge[std::make_pair(i0,i1)];
    halfedge_descriptor h12=edge_to_hedge[std::make_pair(i1,i2)];
    halfedge_descriptor h20=edge_to_hedge[std::make_pair(i2,i0)];

    CGAL_assertion(target(h01,tm)==node_id_to_vertex[i1]);
    CGAL_assertion(target(h12,tm)==node_id_to_vertex[i2]);
    CGAL_assertion(target(h20,tm)==node_id_to_vertex[i0]);

    set_next(h01,h12,tm);
    set_next(h12,h20,tm);
    set_next(h20,h01,tm);

    //update face halfedge
    set_halfedge(current_face,h01,tm);

    //update face of halfedges
    set_face(h01,current_face,tm);
    set_face(h12,current_face,tm);
    set_face(h20,current_face,tm);

    if ( ++it!=it_end )
    {
      current_face=add_face(tm);
      new_face_visitor.before_subface_creations(current_face,tm);
    }
    else
      break;
  }
}


} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H
