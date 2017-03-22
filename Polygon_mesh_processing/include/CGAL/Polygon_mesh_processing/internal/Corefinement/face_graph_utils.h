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

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/property_map.h>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <fstream>
#include <sstream>
#include <set>
namespace CGAL {
namespace Corefinement {

template <typename G>
struct No_mark
{
  friend bool get(No_mark<G>,
                  typename boost::graph_traits<G>::edge_descriptor)
  {
    return false;
  }
  friend void put(No_mark<G>,
                  typename boost::graph_traits<G>::edge_descriptor, bool)
  {}
};

template<class G,
         class EdgeMarkMap>
void mark_all_edges(G& tm, EdgeMarkMap& edge_mark_map)
{
  BOOST_FOREACH(typename boost::graph_traits<G>::edge_descriptor ed,
                edges(tm))
  {
    put(edge_mark_map, ed, true);
  }
}

template<class G>
void mark_all_edges(G&, No_mark<G>&)
{} //nothing to do

template<class G,
         class EdgeMarkMap,
         class HalfedgeRange>
void unmark_edges(      G& tm,
                        EdgeMarkMap& edge_mark_map,
                  const HalfedgeRange& hedges)
{
  BOOST_FOREACH(typename boost::graph_traits<G>::halfedge_descriptor hd, hedges)
    put(edge_mark_map, edge(hd, tm), false);
}

template <class G,
          class HalfedgeRange>
void unmark_edges(      G&,
                        No_mark<G>&,
                  const HalfedgeRange&)
{} //nothing to do

template <class G,
          class edge_descriptor,
          class EdgeMarkMapIn,
          class EdgeMarkMapOut>
void copy_edge_mark(edge_descriptor ed_in,
                    edge_descriptor ed_out,
                    const EdgeMarkMapIn& edge_mark_map_in,
                          EdgeMarkMapOut& edge_mark_map_out)
{
  if(get(edge_mark_map_in, ed_in))
    put(edge_mark_map_out, ed_out, true);
}

template <class G,
          class edge_descriptor,
          class EdgeMarkMapOut>
void copy_edge_mark(edge_descriptor,
                    edge_descriptor,
                    No_mark<G>,
                    EdgeMarkMapOut)
{} // nothing to do

template <class G,
          class edge_descriptor,
          class EdgeMarkMapIn>
void copy_edge_mark(edge_descriptor,
                    edge_descriptor,
                    EdgeMarkMapIn,
                    No_mark<G>)
{} // nothing to do

template <class G,
          class edge_descriptor>
void copy_edge_mark(edge_descriptor,
                    edge_descriptor,
                    No_mark<G>,
                    No_mark<G>)
{} // nothing to do

template <class G,
          class EdgeMarkMapIn,
          class EdgeMarkMapOut>
void copy_edge_mark(G& g,
                    const EdgeMarkMapIn & edge_mark_map_in,
                          EdgeMarkMapOut& edge_mark_map_out)
{
  BOOST_FOREACH(typename boost::graph_traits<G>::edge_descriptor ed, edges(g))
    if(get(edge_mark_map_in, ed))
      put(edge_mark_map_out, ed, true);
}

template <class G,
          class EdgeMarkMapOut>
void copy_edge_mark(G&,
                    const No_mark<G> &,
                          EdgeMarkMapOut&)
{} // nothing to do

template <class G,
          class EdgeMarkMapIn>
void copy_edge_mark(G&,
                    const EdgeMarkMapIn&,
                          No_mark<G>&)
{} // nothing to do

template <class G>
void copy_edge_mark(G&,
                    const No_mark<G>&,
                          No_mark<G>&)
{} // nothing to do

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
  Node_vector& nodes,
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
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;

  //insert the intersection point interior to the face inside the mesh and
  //save their vertex_descriptor
  CGAL_assertion( node_ids.size()== std::set<Node_id>(node_ids.begin(), node_ids.end()).size() );
  BOOST_FOREACH(Node_id node_id, node_ids)
  {
    vertex_descriptor v=add_vertex(tm);
    new_node_visitor.new_vertex_added(node_id, v, tm);
    nodes.call_put(vpm, v, node_id, tm);
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

    CGAL_assertion(edge_to_hedge.count(std::make_pair(i0,i1))!= 0);
    CGAL_assertion(edge_to_hedge.count(std::make_pair(i1,i2))!= 0);
    CGAL_assertion(edge_to_hedge.count(std::make_pair(i2,i0))!= 0);

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

template <class PolygonMesh>
class Border_edge_map {
  typedef std::size_t Node_id;
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef boost::unordered_map<edge_descriptor,
                               std::pair<Node_id,Node_id> > Intersection_edge_map;
  const Intersection_edge_map* intersection_edges;
  const PolygonMesh* tm;
public:
  // required by the property forwarder of the filtered_graph in boost
  Border_edge_map()
   : intersection_edges(NULL)
   , tm(NULL)
  {}

  Border_edge_map(const Intersection_edge_map& intersection_edges,
                 const PolygonMesh& tm)
   : intersection_edges(&intersection_edges)
   , tm(&tm)
  {}

  friend
  bool get(const Border_edge_map& map, edge_descriptor e)
  {
    if ( is_border(e,*map.tm) ) return false;
    return map.intersection_edges->count(e)!=0;
  }
};

template<class PolygonMesh>
struct Patch_description{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  std::vector<face_descriptor> faces;
  std::set<vertex_descriptor> interior_vertices;
  std::vector<halfedge_descriptor> interior_edges;
  std::vector<halfedge_descriptor> shared_edges;
  bool is_initialized;

  Patch_description(): is_initialized(false) {};
};

// shared_edges will be filled by halfedges pointing in the patch
// that are inside `is_intersection_edge`, thus mesh boundary halfedges
// are not necessarily inside.
template <class PolygonMesh, class FaceIndexMap, class IsIntersectionEdge>
void extract_patch_simplices(
  std::size_t patch_id,
  PolygonMesh& pm,
  const FaceIndexMap& fids,
  const std::vector<std::size_t>& patch_ids,
  std::vector<typename boost::graph_traits<PolygonMesh>::face_descriptor>& patch_faces,
  std::set<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& interior_vertices,
  std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>& interior_edges,
  std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>& shared_edges,
  const IsIntersectionEdge& is_intersection_edge)
{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  BOOST_FOREACH(face_descriptor f, faces(pm))
  {
    if ( patch_ids[ get(fids, f) ]==patch_id )
    {
      patch_faces.push_back( f );
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_face(halfedge(f, pm),pm))
      {
        if ( !is_intersection_edge.count(edge(h, pm)) )
        {
          if ( h < opposite(h,pm) || is_border(opposite(h,pm),pm) )
            interior_edges.push_back( h );
        }
        else
          shared_edges.push_back(h);
      }
    }
  }

  std::set<vertex_descriptor> border_vertices;
  BOOST_FOREACH(halfedge_descriptor h, shared_edges)
  {
    border_vertices.insert( target(h,pm) );
    // if the model is not closed i.e. patch_border_halfedge is not cycle-only
    border_vertices.insert( source(h,pm) );
  }

  BOOST_FOREACH(halfedge_descriptor h, interior_edges)
  {
    if ( !border_vertices.count( target(h,pm) ) )
      interior_vertices.insert( target(h,pm) );
    if ( !border_vertices.count( source(h,pm) ) )
      interior_vertices.insert( source(h,pm) );
  }
}

template <class PolygonMesh, class FaceIndexMap, class IsIntersectionEdge>
struct Patch_container{
// data members
  std::vector< Patch_description<PolygonMesh> > patches;
// external data members
  PolygonMesh& pm;
  const std::vector<std::size_t>& patch_ids;
  const FaceIndexMap& fids;
  const IsIntersectionEdge& is_intersection_edge;
// constructor
  Patch_container(
    PolygonMesh& pm,
    const std::vector<std::size_t>& patch_ids,
    const FaceIndexMap& fids,
    const IsIntersectionEdge& is_intersection_edge,
    std::size_t nb_patches
  ) : patches(nb_patches)
    , pm(pm)
    , patch_ids(patch_ids)
    , fids(fids)
    , is_intersection_edge(is_intersection_edge)
  {}

  Patch_description<PolygonMesh>& operator[](std::size_t i) {
    if ( !patches[i].is_initialized )
    {
      extract_patch_simplices(
        i, pm,
        fids, patch_ids,
        patches[i].faces, patches[i].interior_vertices,
        patches[i].interior_edges, patches[i].shared_edges,
        is_intersection_edge
      );

      patches[i].is_initialized=true;
    }
    return patches[i];
  }

  /// debug
  std::ostream& dump_patch(std::size_t i, std::ostream& out)
  {
    typedef boost::graph_traits<PolygonMesh> GT;
    typedef typename GT::vertex_descriptor vertex_descriptor;
    typedef typename GT::halfedge_descriptor halfedge_descriptor;
    typedef typename GT::face_descriptor face_descriptor;

    Patch_description<PolygonMesh>& patch=this->operator[](i);
    out << "OFF\n" << patch.interior_vertices.size() +
                      patch.shared_edges.size();
    out << " " << patch.faces.size() << " 0\n";
    std::map<vertex_descriptor, int> vertexid;
    int id=0;
    BOOST_FOREACH(vertex_descriptor vh, patch.interior_vertices)
    {
      vertexid[vh]=id++;
      out << get(vertex_point, pm, vh) << "\n";
    }

    BOOST_FOREACH(halfedge_descriptor hh, patch.shared_edges)
    {
      vertexid[target(hh, pm)]=id++;
      out << get(vertex_point, pm, target(hh, pm)) << "\n";
    }

    BOOST_FOREACH(face_descriptor f, patch.faces)
    {
      out << "3 " << vertexid[source(halfedge(f,pm),pm)] <<
             " "  << vertexid[target(halfedge(f,pm),pm)] <<
             " "  << vertexid[target(next(halfedge(f,pm),pm),pm)] << "\n";
    }

    return out;
  }

  void dump_patches(const boost::dynamic_bitset<>& selection, std::string prefix)
  {
    for (std::size_t i=selection.find_first();
                     i < selection.npos; i = selection.find_next(i))
    {
      std::stringstream ss;
      ss << prefix << "-" << i << ".off";
      std::ofstream output(ss.str().c_str());
      dump_patch(i, output);
    }
  }
};

template <class PolygonMesh,
          class MarkedEdgeSet>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
next_marked_halfedge_around_target_vertex(
   typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
   const PolygonMesh& pm,
   const MarkedEdgeSet& marked_edges)
{
  CGAL_assertion( marked_edges.count(edge(h,pm))!= 0 );
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor nxt =
    next(h, pm);
  while( !marked_edges.count(edge(nxt,pm)) )
  {
    nxt=next(opposite(nxt,pm),pm);
  }
  CGAL_assertion(nxt!=h);
  return nxt;
}

template <class PolygonMesh,
          class EdgeMap,
          class VertexMap,
          class VertexPointMap,
          class IntersectionEdgeMap>
void import_polyline(
  PolygonMesh& output,
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h1,
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h2,
  const PolygonMesh& pm1,
  const PolygonMesh& pm2,
  std::size_t nb_segments,
  EdgeMap& pm1_to_output_edges,
  EdgeMap& pm2_to_output_edges,
  VertexMap& pm1_to_output_vertices,
  const IntersectionEdgeMap& intersection_edges1,
  const IntersectionEdgeMap& intersection_edges2,
  const VertexPointMap& vpm1,
  const VertexPointMap& /*vpm2*/,
  const VertexPointMap& vpm_out,
  std::vector<typename boost::graph_traits<PolygonMesh>
                ::edge_descriptor>& output_shared_edges)
{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  output_shared_edges.push_back(add_edge(output));
  halfedge_descriptor h_out = halfedge(output_shared_edges.back(),output);

  //make sure the first vertex does not already exist
  vertex_descriptor src = GT::null_vertex();
  std::pair< typename VertexMap::iterator, bool > insert_res=
    pm1_to_output_vertices.insert( std::make_pair( source(h1,pm1), src ) );

  if( insert_res.second )
  {
    src = add_vertex(output);
    set_halfedge(src, opposite(h_out,output),output);
    put(vpm_out, src, get(vpm1, source(h1,pm1)));
    insert_res.first->second = src;
  }
  else
    src = insert_res.first->second;

  //make sure the target vertex does not already exist if it is a polyline endpoint
  vertex_descriptor tgt=GT::null_vertex();
  if ( nb_segments==1 )
  {
    insert_res =
      pm1_to_output_vertices.insert( std::make_pair( target(h1,pm1), tgt ) );
    if( insert_res.second )
    {
      tgt = add_vertex(output);
      set_halfedge(tgt, h_out, output);
      put(vpm_out, tgt, get(vpm1, target(h1,pm1)));
      insert_res.first->second = tgt;
    }
    else
      tgt = insert_res.first->second;
  }
  else{
    tgt = add_vertex(output);
    set_halfedge(tgt, h_out, output);
    put(vpm_out, tgt, get(vpm1, target(h1,pm1)));
  }

  //update source and target vertex of the edge created
  set_target(h_out, tgt, output);
  set_target(opposite(h_out,output), src, output);

  halfedge_descriptor prev_out=h_out;
  halfedge_descriptor prev1=h1;
  halfedge_descriptor prev2=h2;

  //set the correspondance
  pm1_to_output_edges.insert(
    std::make_pair(edge(prev1, pm1), edge(prev_out, output)) );
  pm2_to_output_edges.insert(
    std::make_pair(edge(prev2, pm2), edge(prev_out, output)) );

  src=tgt;
  for (std::size_t i=1; i<nb_segments; ++i)
  {
    //create new edge
    output_shared_edges.push_back(add_edge(output));
    h_out = halfedge(output_shared_edges.back(),output);
    //get the new edge
    h1 = next_marked_halfedge_around_target_vertex(prev1, pm1, intersection_edges1);
    h2 = next_marked_halfedge_around_target_vertex(prev2, pm2, intersection_edges2);
    //if this is the final segment, only create a target vertex if it does not exist
    if (i+1!=nb_segments)
    {
      tgt=add_vertex(output);
      set_halfedge(tgt, h_out, output);
      put(vpm_out, tgt, get(vpm1, target(h1,pm1)));
    }
    else{
      std::pair< typename VertexMap::iterator, bool > insert_res =
        pm1_to_output_vertices.insert(std::make_pair(target(h1,pm1), tgt));
      if (insert_res.second)
      {
        tgt=add_vertex(output);
        set_halfedge(tgt, h_out, output);
        put(vpm_out, tgt, get(vpm1, target(h1,pm1)));
        insert_res.first->second = tgt;
      }
      else
        tgt = insert_res.first->second;
    }

    set_target(h_out, tgt, output);
    set_target(opposite(h_out, output), src, output);

    prev_out=h_out;
    prev1 = h1;
    prev2 = h2;
    src = tgt;

    pm1_to_output_edges.insert(
      std::make_pair(edge(prev1, pm1), edge(prev_out, output)) );
    pm2_to_output_edges.insert(
      std::make_pair(edge(prev2, pm2), edge(prev_out, output)) );
  }
}

template <class TriangleMesh, bool reverse_patch_orientation>
struct Triangle_mesh_extension_helper;

template <class TriangleMesh>
struct Triangle_mesh_extension_helper<TriangleMesh, true>
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  typedef boost::unordered_map< edge_descriptor, edge_descriptor> Edge_map;
  Edge_map& tm_to_output_edges;
  const TriangleMesh& tm;
  TriangleMesh& output;

  Triangle_mesh_extension_helper(Edge_map& tm_to_output_edges,
                                 const TriangleMesh& tm,
                                 TriangleMesh& output)
    : tm_to_output_edges(tm_to_output_edges)
    , tm(tm)
    , output(output)
  {}

  halfedge_descriptor get_hedge(halfedge_descriptor h_tm)
  {
    CGAL_assertion( tm_to_output_edges.count(edge(h_tm, tm))!=0 );
    const std::pair<edge_descriptor, edge_descriptor>& key_and_value =
      *tm_to_output_edges.find(edge(h_tm, tm));
    return halfedge(key_and_value.first,tm) != h_tm
           ? halfedge(key_and_value.second, output)
           : opposite(halfedge(key_and_value.second, output), output);
  }

  cpp11::array<halfedge_descriptor,3>
  halfedges(face_descriptor f)
  {
     halfedge_descriptor h=halfedge(f,tm);
     return make_array( get_hedge( h ),
                        get_hedge( prev(h,tm) ),
                        get_hedge( next(h,tm) ) );
  }
};

template <class TriangleMesh>
struct Triangle_mesh_extension_helper<TriangleMesh, false>
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  typedef boost::unordered_map< edge_descriptor, edge_descriptor> Edge_map;
  Edge_map& tm_to_output_edges;
  const TriangleMesh& tm;
  TriangleMesh& output;

  Triangle_mesh_extension_helper(Edge_map& tm_to_output_edges,
                                 const TriangleMesh& tm,
                                 TriangleMesh& output)
    : tm_to_output_edges(tm_to_output_edges)
    , tm(tm)
    , output(output)
  {}

  halfedge_descriptor get_hedge(halfedge_descriptor h_tm)
  {
    CGAL_assertion( tm_to_output_edges.count(edge(h_tm, tm))!=0 );
    const std::pair<edge_descriptor, edge_descriptor>& key_and_value =
      *tm_to_output_edges.find(edge(h_tm, tm));
    return halfedge(key_and_value.first,tm) == h_tm
           ? halfedge(key_and_value.second, output)
           : opposite(halfedge(key_and_value.second, output), output);
  }

  cpp11::array<halfedge_descriptor,3>
  halfedges(face_descriptor f)
  {
     halfedge_descriptor h=halfedge(f,tm);
     return make_array( get_hedge( h ),
                        get_hedge( next(h,tm) ),
                        get_hedge( prev(h,tm) ) );
  }
};


template < bool reverse_patch_orientation,
           class TriangleMesh,
           class PatchContainer,
           class VertexPointMap,
           class EdgeMarkMapOut,
           class EdgeMarkMapIn >
void append_patches_to_triangle_mesh(
  TriangleMesh& output,
  const boost::dynamic_bitset<>& patches_to_append,
  PatchContainer& patches,
  const VertexPointMap& vpm_out,
  const VertexPointMap& vpm_tm,
  EdgeMarkMapOut& edge_mark_map_out,
  const EdgeMarkMapIn& edge_mark_map_in,
  boost::unordered_map<
    typename boost::graph_traits<TriangleMesh>::edge_descriptor,
    typename boost::graph_traits<TriangleMesh>::edge_descriptor
  >& tm_to_output_edges)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  const TriangleMesh& tm = patches.pm;
  Triangle_mesh_extension_helper<TriangleMesh, reverse_patch_orientation>
    helper(tm_to_output_edges, tm, output);

  for (std::size_t i=patches_to_append.find_first();
                   i < patches_to_append.npos;
                   i = patches_to_append.find_next(i))
  {
    #ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning the size of tm_to_output_edges will increase at each step \
             when adding new patches  by the size of internal edges. \
             Maybe the use of a copy would be better since we do not need \
             the internal edges added?
    #endif

    Patch_description<TriangleMesh>& patch = patches[i];

    std::vector<halfedge_descriptor> interior_vertex_halfedges;

    //insert interior halfedges and create interior vertices
    BOOST_FOREACH(halfedge_descriptor h, patch.interior_edges)
    {
      edge_descriptor new_edge = add_edge(output), ed = edge(h,tm);

      // copy the mark on input edge to the output edge
      copy_edge_mark<TriangleMesh>(ed, new_edge,
                                   edge_mark_map_in, edge_mark_map_out);

      halfedge_descriptor new_h = halfedge(new_edge, output);
      tm_to_output_edges[ed] = new_edge;

      set_face(new_h, GT::null_face(), output);
      set_face(opposite(new_h, output), GT::null_face(), output);

      CGAL_assertion(is_border(new_h, output));
      CGAL_assertion(is_border(opposite(new_h,output), output));

      //create a copy of interior vertices only once
      if (  halfedge(target(h,tm),tm)==h &&
            patch.interior_vertices.count(target(h, tm)) )
      {
        vertex_descriptor v = add_vertex(output);
        set_halfedge(v, new_h, output);
        set_target(new_h, v, output);
        put(vpm_out, v, get(vpm_tm, target(h, tm) ) );
        interior_vertex_halfedges.push_back( new_h );
      }
      if (  halfedge(source(h,tm),tm)==opposite(h,tm) &&
            patch.interior_vertices.count(source(h,tm)) )
      {
        vertex_descriptor v = add_vertex(output);
        halfedge_descriptor new_h_opp = opposite(new_h, output);
        set_halfedge(v, new_h_opp, output);
        set_target(new_h_opp, v, output);
        put(vpm_out, v, get(vpm_tm, source(h, tm) ) );
        interior_vertex_halfedges.push_back( new_h_opp );
      }
    }

    //create faces and connect halfedges
    BOOST_FOREACH(face_descriptor f, patch.faces)
    {
      cpp11::array<halfedge_descriptor,3> hedges = helper.halfedges(f);

      face_descriptor new_f = add_face(output);
      set_halfedge(new_f, hedges[0], output);

      for (int i=0;i<3;++i)
      {
        CGAL_assertion(hedges[i]!=GT::null_halfedge());
        set_next(hedges[i], hedges[(i+1)%3], output);
        set_face(hedges[i], new_f, output);
      }
    }

    // handle interior edges that are on the border of the mesh:
    // they do not have a prev/next pointer set since only the pointers
    // of patch interior halfedges part a face have been. In the following
    // (i) we set the next/prev pointer around interior vertices on the mesh
    // boundary and (ii) we collect interior mesh border halfedges incident to
    // a patch border vertex and set their next/prev pointer (possibly of
    // another patch)

    // Containers used for step (ii) for collecting mesh border halfedges
    // with source/target on an intersection polyline that needs it prev/next
    // pointer to be set
    std::vector<halfedge_descriptor> border_halfedges_source_to_link;
    std::vector<halfedge_descriptor> border_halfedges_target_to_link;
    BOOST_FOREACH(halfedge_descriptor h, patch.interior_edges)
      if (is_border_edge(h,tm))
      {
        if (!is_border(h, tm)) h=opposite(h, tm);

        vertex_descriptor src = source(h, tm);
        vertex_descriptor tgt = target(h, tm);
        if (reverse_patch_orientation) std::swap(src, tgt);

        if ( !patch.interior_vertices.count(src) )
          border_halfedges_source_to_link.push_back(helper.get_hedge(h));
        if ( !patch.interior_vertices.count(tgt) ){
          border_halfedges_target_to_link.push_back(helper.get_hedge(h));
          continue; // since the next halfedge should not be in the same patch
        }
        CGAL_assertion( is_border(h, tm) &&
                        is_border(prev(h, tm),tm) &&
                        is_border(next(h, tm),tm));
        // step (i)
        halfedge_descriptor h_out=helper.get_hedge(h);
        halfedge_descriptor h_out_next = reverse_patch_orientation
                                         ? helper.get_hedge(prev(h,tm))
                                         : helper.get_hedge(next(h,tm));
        CGAL_assertion(is_border(h_out,output) && is_border(h_out_next,output));
        set_next(h_out, h_out_next, output);
      }
    // now the step (ii) we look for the candidate halfedge by turning around
    // the vertex in the direction of the interior of the patch
    BOOST_FOREACH(halfedge_descriptor h_out, border_halfedges_target_to_link)
    {
      halfedge_descriptor candidate =
        opposite(prev(opposite(h_out, output), output), output);
      CGAL_assertion_code(halfedge_descriptor start=candidate);
      while (!is_border(candidate, output)){
        candidate=opposite(prev(candidate, output), output);
        CGAL_assertion(candidate!=start);
      }
      set_next(h_out, candidate, output);
    }
    BOOST_FOREACH(halfedge_descriptor h_out, border_halfedges_source_to_link)
    {
      halfedge_descriptor candidate =
      opposite(next(opposite(h_out, output), output), output);
      while (!is_border(candidate, output))
        candidate = opposite(next(candidate, output), output);
      set_next(candidate, h_out, output);
    }

    // For all interior vertices, update the vertex pointer
    // of all but the vertex halfedge
    BOOST_FOREACH(halfedge_descriptor h_out, interior_vertex_halfedges)
    {
      vertex_descriptor v = target(h_out, output);
      halfedge_descriptor next_around_vertex=h_out;
      do{
        CGAL_assertion(next(next_around_vertex, output) != GT::null_halfedge());
        next_around_vertex=opposite(next(next_around_vertex, output), output);
        set_target(next_around_vertex, v, output);
      }while(h_out != next_around_vertex);
    }

    // For all patch boundary vertices, update the vertex pointer
    // of all but the vertex halfedge
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
    {
      //check for a halfedge pointing inside an already imported patch
      halfedge_descriptor h_out = helper.get_hedge(h);
      CGAL_assertion( next(h_out, output)!=GT::null_halfedge() );
      // update the pointers on the target
      halfedge_descriptor next_around_target=h_out;
      vertex_descriptor v=target(h_out, output);
      do{
        next_around_target = opposite(next(next_around_target, output), output);
        set_target(next_around_target, v, output);
      }while(next(next_around_target, output)!=GT::null_halfedge() &&
             next_around_target!=h_out &&
             !is_border(next_around_target, output));
      // update the pointers on the source
      halfedge_descriptor next_around_source=prev(h_out, output);
      CGAL_assertion(next_around_source!=GT::null_halfedge());
      v = source(h_out, output);
      do{
        set_target(next_around_source, v, output);
        next_around_source = prev(opposite(next_around_source, output), output);
      }while( next_around_source!=GT::null_halfedge() &&
              next_around_source!=opposite(h_out, output) &&
              !is_border(next_around_source, output));
    }
  }
}

template < class TriangleMesh,
           class IntersectionEdgeMap,
           class VertexPointMap,
           class EdgeMarkMap1,
           class EdgeMarkMap2,
           class EdgeMarkMapOut,
           class IntersectionPolylines,
           class PatchContainer>
void fill_new_triangle_mesh(
  TriangleMesh& output,
  const boost::dynamic_bitset<>& patches_of_tm1_to_import,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer& patches_of_tm1,
  PatchContainer& patches_of_tm2,
  bool reverse_orientation_of_patches_from_tm1,
  bool reverse_orientation_of_patches_from_tm2,
  const IntersectionPolylines& polylines,
  const IntersectionEdgeMap& intersection_edges1,
  const IntersectionEdgeMap& intersection_edges2,
  const VertexPointMap& vpm1,
  const VertexPointMap& vpm2,
  const VertexPointMap& vpm_out,
  const EdgeMarkMap1& edge_mark_map1,
  const EdgeMarkMap2& edge_mark_map2,
        EdgeMarkMapOut& edge_mark_map_out,
  std::vector< typename boost::graph_traits<TriangleMesh>::edge_descriptor>&
                                                            output_shared_edges)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;

  // this is the miminal number of edges that will be marked (intersection edge).
  // We cannot easily have the total number since some patch interior edges might be marked
  output_shared_edges.reserve(
                              std::accumulate(polylines.lengths.begin(),polylines.lengths.end(),std::size_t(0)) );

  //add a polyline inside O for each intersection polyline
  std::size_t nb_polylines = polylines.lengths.size();
  boost::unordered_map<vertex_descriptor, vertex_descriptor> tm1_to_output_vertices;
  boost::unordered_map<edge_descriptor, edge_descriptor> tm1_to_output_edges,
                                                         tm2_to_output_edges;

  for (std::size_t i=0; i < nb_polylines; ++i)
    if (!polylines.to_skip.test(i))
      import_polyline(output,
                      polylines.tm1[i], polylines.tm2[i],
                      patches_of_tm1.pm,
                      patches_of_tm2.pm,
                      polylines.lengths[i],
                      tm1_to_output_edges, tm2_to_output_edges,
                      tm1_to_output_vertices,
                      intersection_edges1, intersection_edges2,
                      vpm1, vpm2, vpm_out,
                      output_shared_edges);

  //import patches of tm1
  if (reverse_orientation_of_patches_from_tm1)
    append_patches_to_triangle_mesh<true>(output,
                                          patches_of_tm1_to_import,
                                          patches_of_tm1,
                                          vpm_out,
                                          vpm1,
                                          edge_mark_map_out,
                                          edge_mark_map1,
                                          tm1_to_output_edges);
  else
    append_patches_to_triangle_mesh<false>(output,
                                           patches_of_tm1_to_import,
                                           patches_of_tm1,
                                           vpm_out,
                                           vpm1,
                                           edge_mark_map_out,
                                           edge_mark_map1,
                                           tm1_to_output_edges);

  //import patches from tm2
  if (reverse_orientation_of_patches_from_tm2)
    append_patches_to_triangle_mesh<true>(output,
                                          patches_of_tm2_to_import,
                                          patches_of_tm2,
                                          vpm_out,
                                          vpm2,
                                          edge_mark_map_out,
                                          edge_mark_map2,
                                          tm2_to_output_edges);
  else
    append_patches_to_triangle_mesh<false>(output,
                                           patches_of_tm2_to_import,
                                           patches_of_tm2,
                                           vpm_out,
                                           vpm2,
                                           edge_mark_map_out,
                                           edge_mark_map2,
                                           tm2_to_output_edges);
}

template <class TriangleMesh,
          class PatchContainer,
          class EdgeMap>
void disconnect_patches(
  TriangleMesh& tm1,
  const boost::dynamic_bitset<>& patches_to_remove,
  PatchContainer& patches_of_tm1,
  const EdgeMap& tm1_edge_to_tm2_edge, //map intersection edges of tm1 to the equivalent in tm2
        EdgeMap& new_tm1_edge_to_tm2_edge) //map the new intersection edges of tm1 to the equivalent in tm2
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  // disconnect each patch one by one
  for (std::size_t i=patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch=patches_of_tm1[i];
    std::vector<halfedge_descriptor> new_patch_border;

    // start the duplicate step
    std::size_t nb_shared_edges = patch.shared_edges.size();
    new_patch_border.reserve( nb_shared_edges );
    boost::unordered_map<halfedge_descriptor, halfedge_descriptor> old_to_new;

    // save faces inside the patch and set the halfedge
    // to be duplicated on the boundary
    std::vector<face_descriptor> face_backup;
    face_backup.reserve( nb_shared_edges );
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
    {
      face_backup.push_back( face(h, tm1) );
      set_face(h, GT::null_face(), tm1);
    }

    // look for the prev/next halfedges after the patch will be disconnected
    // for the halfedges to be duplicated
    std::vector<halfedge_descriptor> shared_next, shared_prev;
    shared_next.reserve( nb_shared_edges );
    shared_prev.reserve( nb_shared_edges );
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
    {
      halfedge_descriptor nxt=next(h, tm1);
      while(!is_border(nxt, tm1))
        nxt=next(opposite(nxt, tm1), tm1);
      shared_next.push_back(nxt);

      // setting prev is only needed in case tm1 is open
      // and the intersection polyline intersects its boundary
      halfedge_descriptor prv=prev(h, tm1);
      while( !is_border(prv, tm1) )
        prv=prev(opposite(prv, tm1), tm1);
      shared_prev.push_back(prv);

      set_halfedge(target(h, tm1), h, tm1);
      set_halfedge(source(h, tm1), opposite(h, tm1), tm1); //only needed if tm1 is open
    }

    // now duplicate the edge and set its pointers
    for(std::size_t k=0; k<nb_shared_edges; ++k)
    {
      halfedge_descriptor h = patch.shared_edges[k];
      halfedge_descriptor new_hedge = halfedge(add_edge(tm1), tm1);
      set_next(new_hedge, next(h, tm1), tm1);
      set_next(prev(h, tm1), new_hedge, tm1);
      set_face(new_hedge, face_backup[k], tm1);
      set_target(new_hedge, target(h, tm1), tm1);
      set_target(opposite(new_hedge,tm1), source(h, tm1), tm1);
      set_halfedge(face_backup[k], new_hedge, tm1);

      new_patch_border.push_back(new_hedge);
      set_face(opposite(new_hedge, tm1), GT::null_face(), tm1);
      old_to_new.insert( std::make_pair(h, new_hedge) );

      CGAL_assertion( next(opposite(new_hedge, tm1), tm1)==GT::null_halfedge() );
      CGAL_assertion( prev(opposite(new_hedge, tm1), tm1)==GT::null_halfedge() );
      CGAL_assertion( next(prev(new_hedge, tm1), tm1) == new_hedge );
      CGAL_assertion( prev(next(new_hedge, tm1), tm1) == new_hedge );
    }

    // update next/prev pointer of new hedges in case it is one of the new hedges
    BOOST_FOREACH(halfedge_descriptor h, new_patch_border)
      if (is_border(next(h, tm1), tm1))
        set_next(h, old_to_new[next(h,tm1)], tm1);

    // set next/prev pointers on the border of the neighbor patch
    for(std::size_t k=0; k<nb_shared_edges; ++k)
    {
      halfedge_descriptor h = patch.shared_edges[k];
      set_next(h, shared_next[k], tm1);
      set_next(shared_prev[k], h, tm1);
    }

    // update next/prev pointers on the border of the patch
    BOOST_FOREACH(halfedge_descriptor h, new_patch_border)
    {
      halfedge_descriptor h_opp = opposite(h, tm1);
      //set next pointer if not already set
      if ( next(h_opp, tm1)==GT::null_halfedge() )
      {
        // we visit faces inside the patch we consider
        halfedge_descriptor candidate = opposite(prev(h, tm1), tm1);
        while ( !is_border(candidate, tm1) )
          candidate = opposite(prev(candidate, tm1), tm1);
        set_next(h_opp, candidate, tm1);
      }
      CGAL_assertion( prev(next(h_opp, tm1), tm1)==h_opp );

      // set prev pointer if not already set
      if ( prev(h_opp, tm1) == GT::null_halfedge() )
      {
        halfedge_descriptor candidate = opposite(next(h, tm1), tm1);
        while ( !is_border(candidate, tm1) )
          candidate = opposite(next(candidate, tm1), tm1);
        set_next(candidate, h_opp, tm1);
      }

      CGAL_assertion( prev(next(h_opp, tm1), tm1) == h_opp );
      CGAL_assertion( is_border(prev(h_opp, tm1), tm1) );
      CGAL_assertion( is_border(next(h_opp, tm1), tm1) );
    }
    // end of the duplicate step

    CGAL_assertion( new_patch_border.size() == nb_shared_edges );

    for (std::size_t k=0; k < nb_shared_edges; ++k){
      CGAL_assertion( target(patch.shared_edges[k], tm1) == target(new_patch_border[k], tm1) );
      CGAL_assertion( source(patch.shared_edges[k], tm1) == source(new_patch_border[k], tm1) );
      CGAL_assertion( is_border_edge(new_patch_border[k], tm1) );
      CGAL_assertion( !is_border(new_patch_border[k], tm1) );
      CGAL_assertion( is_border(next(opposite(new_patch_border[k], tm1), tm1), tm1) );
      CGAL_assertion( is_border(prev(opposite(new_patch_border[k], tm1), tm1), tm1) );

      typename EdgeMap::const_iterator it_res =
        tm1_edge_to_tm2_edge.find( edge(patch.shared_edges[k], tm1) );
      CGAL_assertion( it_res != tm1_edge_to_tm2_edge.end() );

      new_tm1_edge_to_tm2_edge[
        patch.shared_edges[k]==halfedge(it_res->first, tm1)
        ? edge(new_patch_border[k], tm1)
        : edge(opposite(new_patch_border[k], tm1), tm1) ] = it_res->second;
    }

    patch.shared_edges.swap(new_patch_border);
  }
}

template <class TriangleMesh,
          class PatchContainer,
          class IntersectionPolylines,
          class EdgeMap,
          class VertexPointMap,
          class EdgeMarkMapIn1,
          class EdgeMarkMapIn2,
          class EdgeMarkMapOut>
void compute_inplace_operation_delay_removal_and_insideout(
  TriangleMesh& tm1,
  TriangleMesh& tm2,
  const boost::dynamic_bitset<>& patches_of_tm1_to_keep,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer& patches_of_tm1,
  PatchContainer& patches_of_tm2,
  bool reverse_patch_orientation_tm2,
  const IntersectionPolylines& polylines,
  const VertexPointMap& vpm1,
  const VertexPointMap& vpm2,
        EdgeMarkMapIn1&,
  const EdgeMarkMapIn2& edge_mark_map2,
  const EdgeMarkMapOut& edge_mark_map_out1,
  EdgeMap& disconnected_patches_edge_to_tm2_edge)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef boost::unordered_map<edge_descriptor, edge_descriptor> Edge_map;

  Edge_map tm2_edge_to_tm1_edge, tm1_edge_to_tm2_edge;
  //maps intersection edges from tm2 to tm1
  std::size_t nb_polylines = polylines.lengths.size();
  for(std::size_t i=0; i<nb_polylines; ++i)
  {
    halfedge_descriptor h1 = polylines.tm1[i];
    halfedge_descriptor h2 = polylines.tm2[i];
    std::size_t nb_segments = polylines.lengths[i];

    for (std::size_t k=0;;)
    {
      tm2_edge_to_tm1_edge[edge(h2, tm2)]=edge(h1, tm1);
      tm1_edge_to_tm2_edge[edge(h1, tm1)]=edge(h2, tm2);
      if (++k==nb_segments) break;
      h2 = next_marked_halfedge_around_target_vertex(h2, tm2,
             patches_of_tm2.is_intersection_edge);
      h1=next_marked_halfedge_around_target_vertex(h1, tm1,
           patches_of_tm1.is_intersection_edge);
    }
  }

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
  #warning do not try to disconnect if the patch is isolated? i.e opposite(border_edge_of_patch)->is_border()
#endif
  // disconnect patches inside tm2
  // For the patches scheduled to be removed, their patch descriptions
  // in patches_of_tm1 will be updated so that is_intersection_edge are
  // the newly created halfedges within disconnect_patches.
  // Note that disconnected_patches_edge_to_tm2_edge also refers to those halfedges
  //init the map with the previously filled one (needed when reusing patches in two operations)
  disconnected_patches_edge_to_tm2_edge=tm1_edge_to_tm2_edge;
  disconnect_patches(tm1, ~patches_of_tm1_to_keep, patches_of_tm1,
                     tm1_edge_to_tm2_edge, disconnected_patches_edge_to_tm2_edge);

  //we import patches from tm2
  if (reverse_patch_orientation_tm2)
    append_patches_to_triangle_mesh<true>(tm1,
                                          patches_of_tm2_to_import,
                                          patches_of_tm2,
                                          vpm1,
                                          vpm2,
                                          edge_mark_map_out1,
                                          edge_mark_map2,
                                          tm2_edge_to_tm1_edge);
  else
    append_patches_to_triangle_mesh<false>(tm1,
                                           patches_of_tm2_to_import,
                                           patches_of_tm2,
                                           vpm1,
                                           vpm2,
                                           edge_mark_map_out1,
                                           edge_mark_map2,
                                           tm2_edge_to_tm1_edge);
}

template <class TriangleMesh,
          class PatchContainer,
          class EdgeMarkMap>
void
remove_patches(TriangleMesh& tm,
               const boost::dynamic_bitset<>& patches_to_remove,
               PatchContainer& patches,
               const EdgeMarkMap& edge_mark_map)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  for (std::size_t i=patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch=patches[i];

    // put the halfedges on the boundary of the patch on the boundary of tm
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
      set_face(h, GT::null_face(), tm);

    // set next/prev relationship of border halfedges
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
    {
      halfedge_descriptor nxt=next(h, tm);
      while(!is_border(nxt,tm))
        nxt=next(opposite(nxt, tm), tm);
      set_next(h, nxt, tm);
      set_halfedge(target(h, tm), h, tm);
    }

    // edges removed must be unmarked to avoid issues when adding new elements
    // that could be marked because they retrieve a previously set property
    unmark_edges(tm, edge_mark_map, patch.interior_edges);

    // In case a ccb of the patch is not a cycle (the source and target vertices
    // are border vertices), the first halfedge of that ccb will not have its
    // prev pointer set correctly. To fix that, we consider all interior edges
    // and check for one that is on the border of the patch and that is incident
    // to a border vertex and use it to get the missing prev pointer.
    BOOST_FOREACH(halfedge_descriptor h, patch.interior_edges)
      if(is_border_edge(h, tm))
      {
        if (is_border(h, tm)) h=opposite(h, tm);
        if ( !patch.interior_vertices.count(target(h, tm)) )
        {
          // look for the halfedge belonging to shared_edges
          // having the prev pointer not correctly set
          halfedge_descriptor nxt=next(h, tm);
          while(!is_border(nxt, tm))
            nxt=next(opposite(nxt, tm), tm);
          CGAL_assertion( is_border(nxt, tm) );//we marked it above!
          // now update the prev pointer
          halfedge_descriptor prv=prev(opposite(h, tm), tm);
          set_next(prv, next(h, tm), tm);
          set_halfedge(target(prv, tm), prv, tm);
        }
      }

     //now remove the simplices
    BOOST_FOREACH(halfedge_descriptor h, patch.interior_edges)
      remove_edge(edge(h, tm), tm);
    BOOST_FOREACH(face_descriptor f, patch.faces)
      remove_face(f, tm);
    BOOST_FOREACH(vertex_descriptor v, patch.interior_vertices)
      remove_vertex(v, tm);
  }
}

template <class TriangleMesh,
          class PatchContainer,
          class VertexPointMap,
          class EdgeMarkMapIn1,
          class EdgeMarkMapIn2,
          class EdgeMarkMapOut1>
void compute_inplace_operation(
        TriangleMesh& tm1,
  const TriangleMesh& /*tm2*/,
  const boost::dynamic_bitset<>& patches_of_tm1_to_keep,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer& patches_of_tm1,
  PatchContainer& patches_of_tm2,
  bool reverse_patch_orientation_tm1,
  bool reverse_patch_orientation_tm2,
  const VertexPointMap& vpm1,
  const VertexPointMap& vpm2,
        EdgeMarkMapIn1& edge_mark_map_in1,
  const EdgeMarkMapIn2& edge_mark_map_in2,
        EdgeMarkMapOut1& edge_mark_map_out1,
  boost::unordered_map<
    typename boost::graph_traits<TriangleMesh>::edge_descriptor,
    typename boost::graph_traits<TriangleMesh>::edge_descriptor
  >& tm2_edge_to_tm1_edge)
{
  typedef boost::unordered_map<
      typename boost::graph_traits<TriangleMesh>::edge_descriptor,
      typename boost::graph_traits<TriangleMesh>::edge_descriptor> EdgeMap;
  //clean up patches not kept
  remove_patches(tm1, ~patches_of_tm1_to_keep, patches_of_tm1, edge_mark_map_in1);

  // transfer marks of edges of patches kept to the output edge mark property
  copy_edge_mark<TriangleMesh>(tm1, edge_mark_map_in1, edge_mark_map_out1);

  if (reverse_patch_orientation_tm1){
    Polygon_mesh_processing::
      reverse_face_orientations_of_mesh_with_polylines(tm1);
    // here we need to update the mapping to use the correct border
    // halfedges while appending the patches from tm2
    BOOST_FOREACH(typename EdgeMap::value_type& v, tm2_edge_to_tm1_edge)
      v.second=edge(opposite(halfedge(v.second, tm1), tm1), tm1);
  }

  //we import patches from tm2
  if ( reverse_patch_orientation_tm2 )
    append_patches_to_triangle_mesh<true>(tm1,
                                          patches_of_tm2_to_import,
                                          patches_of_tm2,
                                          vpm1,
                                          vpm2,
                                          edge_mark_map_out1,
                                          edge_mark_map_in2,
                                          tm2_edge_to_tm1_edge);
  else
    append_patches_to_triangle_mesh<false>(tm1,
                                           patches_of_tm2_to_import,
                                           patches_of_tm2,
                                           vpm1,
                                           vpm2,
                                           edge_mark_map_out1,
                                           edge_mark_map_in2,
                                           tm2_edge_to_tm1_edge);
}

template <class TriangleMesh,
          class IntersectionPolylines,
          class PatchContainer,
          class EdgeMap>
void compute_border_edge_map(
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const IntersectionPolylines& polylines,
  PatchContainer& patches_of_tm1,
  PatchContainer& patches_of_tm2,
  EdgeMap& tm2_edge_to_tm1_edge)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  std::size_t nb_polylines = polylines.lengths.size();
  for( std::size_t i=0; i<nb_polylines; ++i)
  {
    if (polylines.to_skip.test(i)) continue;
    halfedge_descriptor h1 = polylines.tm1[i];
    halfedge_descriptor h2 = polylines.tm2[i];
    std::size_t nb_segments = polylines.lengths[i];

    for (std::size_t k=0;;)
    {
      tm2_edge_to_tm1_edge[edge(h2, tm2)]=edge(h1, tm1);
      if (++k==nb_segments) break;
      h2 = next_marked_halfedge_around_target_vertex(
            h2, tm2, patches_of_tm2.is_intersection_edge);
      h1 = next_marked_halfedge_around_target_vertex(
            h1, tm1, patches_of_tm1.is_intersection_edge);
    }
  }
}


template <class TriangleMesh,
          class PatchContainer,
          class IntersectionPolylines,
          class VertexPointMap,
          class EdgeMarkMapIn1,
          class EdgeMarkMapIn2,
          class EdgeMarkMapOut1>
void compute_inplace_operation(
        TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const boost::dynamic_bitset<>& patches_of_tm1_to_keep,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer& patches_of_tm1,
  PatchContainer& patches_of_tm2,
  bool reverse_patch_orientation_tm1,
  bool reverse_patch_orientation_tm2,
  const VertexPointMap& vpm1,
  const VertexPointMap& vpm2,
  const EdgeMarkMapIn1& edge_mark_map_in1,
  const EdgeMarkMapIn2& edge_mark_map_in2,
  const EdgeMarkMapOut1& edge_mark_map_out1,
  const IntersectionPolylines& polylines)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;

  boost::unordered_map<edge_descriptor, edge_descriptor> tm2_edge_to_tm1_edge;

  //maps intersection edges from tm2 to the equivalent in tm1
  compute_border_edge_map(tm1, tm2,
                          polylines,
                          patches_of_tm1, patches_of_tm2,
                          tm2_edge_to_tm1_edge);

  compute_inplace_operation(tm1, tm2,
                            patches_of_tm1_to_keep, patches_of_tm2_to_import,
                            patches_of_tm1, patches_of_tm2,
                            reverse_patch_orientation_tm1,
                            reverse_patch_orientation_tm2,
                            vpm1,
                            vpm2,
                            edge_mark_map_in1,
                            edge_mark_map_in2,
                            edge_mark_map_out1,
                            tm2_edge_to_tm1_edge);
}

// function used to remove polylines imported or kept that are incident only
// to patches not kept for the operation P_ptr is used for storing
// the result. We look for edges with halfedges both on the border of
// the mesh. The vertices incident only to such edges should be removed.
// Here to detect vertices that should be kept, we abuse the fact that
// the halfedge to be removed and incident to a vertex that should not be
// removed will still have its next pointer set to a halfedge part of
// the result.
template <class TriangleMesh, class PatchContainer>
void remove_unused_polylines(
  TriangleMesh& tm,
  const boost::dynamic_bitset<>& patches_to_remove,
  PatchContainer& patches)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;

  std::set<vertex_descriptor> vertices_to_remove;
  std::set<edge_descriptor> edges_to_remove;

  for (std::size_t i = patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch=patches[i];
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
    {
      if (is_border(h, tm) && is_border(opposite(h, tm), tm)){
        vertices_to_remove.insert(target(h, tm));
        vertices_to_remove.insert(source(h, tm));
        edges_to_remove.insert(edge(h, tm));
      }
    }
  }

  BOOST_FOREACH(vertex_descriptor v, vertices_to_remove)
  {
    bool to_remove=true;
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v,tm))
      if (!is_border(h, tm) || !is_border(opposite(h,tm),tm))
      {
        to_remove=false;
        // in case the vertex halfedge was one that is going to remove,
        // update it
        set_halfedge(v, h, tm);
        break;
      }
    if (to_remove)
      remove_vertex(v,tm);
  }
  BOOST_FOREACH(edge_descriptor e, edges_to_remove)
    remove_edge(e,tm);
}

template <class TriangleMesh, class PatchContainer, class EdgeMarkMap>
void remove_disconnected_patches(
  TriangleMesh& tm,
  PatchContainer& patches,
  const boost::dynamic_bitset<>& patches_to_remove,
  EdgeMarkMap& edge_mark_map)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  for (std::size_t i=patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch = patches[i];

    // edges removed must be unmarked to avoid issues when adding new elements
    // that could be marked because they retrieve a previously set property
    unmark_edges(tm, edge_mark_map, patch.interior_edges);

    BOOST_FOREACH(halfedge_descriptor h, patch.interior_edges)
      remove_edge(edge(h, tm), tm);
    // There is no shared halfedge between duplicated patches even
    // if they were before the duplication. Thus the erase that follows is safe.
    // However remember that vertices were not duplicated which is why their
    // removal is not handled here (still in use or to be removed in
    // remove_unused_polylines())
    BOOST_FOREACH(halfedge_descriptor h, patch.shared_edges)
      remove_edge(edge(h, tm), tm);
    BOOST_FOREACH(face_descriptor f, patch.faces)
      remove_face(f, tm);
    BOOST_FOREACH(vertex_descriptor v, patch.interior_vertices)
      remove_vertex(v, tm);
  }
}

} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H
