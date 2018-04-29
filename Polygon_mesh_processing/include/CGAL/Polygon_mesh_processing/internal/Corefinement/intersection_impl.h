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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/Box_with_info_d.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_callbacks.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Intersection_type.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_of_coplanar_triangles_3.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_nodes.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersect_triangle_and_segment_3.h>
#include <CGAL/utility.h>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>

namespace CGAL{
namespace Corefinement {

struct Self_intersection_exception{};
struct Triple_intersection_exception{};
// This functor computes the pairwise intersection of triangle meshes.
// Intersection are given as a set of polylines
// The algorithm works as follow:
// From each triangle mesh, we extract a set of segments from the edges and a
// set of triangles from the faces.
// We use Box_intersection_d to filter intersection between the segments from
// one mesh with the triangles from the other one.
// From this filtered set, for each pair (segment,triangle), we look at the
// intersection type. If not empty, we can have three different cases
//   1)the segment intersect the interior of the triangle:
//        We compute the intersection point and for each face incident to the
//        edge of the segment, we save the fact that the intersection points
//        is common to that face and the face of the intersected triangle
//   2)the segment intersect the triangle on an edge
//        We do the same thing as described above but
//        for all faces incident to the edge intersected
//   3)the segment intersect the triangle at a vertex
//        for each edge incident to the vertex, we do
//        the same operations as in 2)
//
// In case a segment belong to the triangle, we repeat the same procedure
// for each segment incident to this endpoint.
//
// Note that given a pair (segment,triangle)=(S,T), if S belongs
// to the plane of T, we have nothing to do in the following cases:
//   -- no triangle T' contains S such that T and T' are coplanar
//   -- at least one triangle contains S
// Indeed, the intersection points of S and T will be found using segments
// of T or segments adjacent to S.
//
// Coplanar triangles are filtered out and handled separately.
//
template<class TriangleMesh, bool doing_autorefinement = false>
struct Default_surface_intersection_visitor{
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;

  void new_node_added(
    std::size_t,Intersection_type,halfedge_descriptor,halfedge_descriptor,
    const TriangleMesh&,const TriangleMesh&,bool,bool){}
  template<class Graph_node>
  void annotate_graph(std::vector<Graph_node>&){}
  void update_terminal_nodes(std::vector<bool>&){}
  void set_number_of_intersection_points_from_coplanar_faces(std::size_t){};
  void start_new_polyline(std::size_t,std::size_t){}
  void add_node_to_polyline(std::size_t){}
  void input_have_coplanar_faces(){}
  template<class T,class VertexPointMap>
  void finalize(T&,
                const TriangleMesh&, const TriangleMesh&,
                const VertexPointMap&, const VertexPointMap&)
  {}
  void new_node_added_triple_face(std::size_t /* node_id */,
                                  face_descriptor /* f1 */,
                                  face_descriptor /* f2 */,
                                  face_descriptor /* f3 */,
                                  const TriangleMesh& /* tm */)
  {}
  // this is required in autorefinement for the do-intersect of 3 segments.
  // If we implement a predicate only test, we can get rid of it.
  static const bool Predicates_on_constructions_needed = doing_autorefinement;
  static const bool do_need_vertex_graph = false;
};


struct Node_id_set {
  typedef std::size_t Node_id;

  Node_id first;
  Node_id second;
  std::size_t size_;

  Node_id_set()
    : size_(0)
  {}

  void insert(std::size_t v){
    if(size_ == 0){
      first = v;
      ++size_;
    }
    else if((size_ == 1) && (v != first)){
      if(v < first){
        second = first;
        first = v;
      } else {
        second = v;
      }
      ++size_;
    }
    else{
      CGAL_assertion( (size_ == 2 && (v == first || v == second) ) || v==first );
    }
  }

  std::size_t size() const
  {
    return size_;
  }

  std::size_t operator[](std::size_t i) const
  {
    CGAL_assertion( (i < size_) && (i == 0 || i == 1));
    return (i == 0)? first : second;
  }
};

template< class TriangleMesh,
          class VertexPointMap,
          class Node_visitor=Default_surface_intersection_visitor<TriangleMesh>
         >
class Intersection_of_triangle_meshes
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::face_descriptor face_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;

  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor> Box;

  typedef boost::unordered_set<face_descriptor> Face_set;
  typedef boost::unordered_map<edge_descriptor, Face_set> Edge_to_faces;

  static const bool Predicates_on_constructions_needed =
    Node_visitor::Predicates_on_constructions_needed;

  typedef std::pair<face_descriptor, face_descriptor> Face_pair;
  typedef std::set< Face_pair > Coplanar_face_set;

  typedef std::size_t Node_id;

  typedef std::pair<Face_pair,int>                            Face_pair_and_int;
  // we use Face_pair_and_int and not Face_pair to handle coplanar case.
  // Indeed the boundary of the intersection of two coplanar triangles
  // may contain several segments.
  typedef std::map< Face_pair_and_int, Node_id_set >         Faces_to_nodes_map;
  typedef Intersection_nodes<TriangleMesh,
                             VertexPointMap,
                             Predicates_on_constructions_needed>    Node_vector;

// data members
  Edge_to_faces stm_edge_to_ltm_faces; // map edges from the triangle mesh with the smaller address to faces of the triangle mesh with the larger address
  Edge_to_faces ltm_edge_to_stm_faces; // map edges from the triangle mesh with the larger address to faces of the triangle mesh with the smaller address
  // here face descriptor are from tmi and tmj such that &tmi<&tmj
  Coplanar_face_set coplanar_faces;
  Node_vector nodes;
  Node_visitor visitor;
  Faces_to_nodes_map         f_to_node;      //Associate a pair of triangles to their intersection points
  std::vector<Node_id> extra_terminal_nodes; //used only for autorefinement
  CGAL_assertion_code(bool doing_autorefinement;)
// member functions
  void filter_intersections(const TriangleMesh& tm_f,
                            const TriangleMesh& tm_e,
                            const VertexPointMap& vpm_f,
                            const VertexPointMap& vpm_e,
                            bool throw_on_self_intersection)
  {
    std::vector<Box> face_boxes, edge_boxes;
    std::vector<Box*> face_boxes_ptr, edge_boxes_ptr;

    face_boxes.reserve(num_faces(tm_f));
    face_boxes_ptr.reserve(num_faces(tm_f));
    BOOST_FOREACH(face_descriptor fd, faces(tm_f))
    {
      halfedge_descriptor h=halfedge(fd,tm_f);
      face_boxes.push_back( Box(
        get(vpm_f,source(h,tm_f)).bbox() +
        get(vpm_f,target(h,tm_f)).bbox() +
        get(vpm_f,target(next(h,tm_f),tm_f)).bbox(),
        h ) );
      face_boxes_ptr.push_back( &face_boxes.back() );
    }

    edge_boxes.reserve(num_edges(tm_e));
    edge_boxes_ptr.reserve(num_edges(tm_e));
    BOOST_FOREACH(edge_descriptor ed, edges(tm_e))
    {
      halfedge_descriptor h=halfedge(ed,tm_e);
      edge_boxes.push_back( Box(
        get(vpm_e,source(h,tm_e)).bbox() +
        get(vpm_e,target(h,tm_e)).bbox(),
        h ) );
      edge_boxes_ptr.push_back( &edge_boxes.back() );
    }

    /// \todo experiments different cutoff values
    std::ptrdiff_t cutoff = 2 * std::ptrdiff_t(
        std::sqrt(face_boxes.size()+edge_boxes.size()) );

    Edge_to_faces& edge_to_faces = &tm_e < &tm_f
                                 ? stm_edge_to_ltm_faces
                                 : ltm_edge_to_stm_faces;

    #ifdef DO_NOT_HANDLE_COPLANAR_FACES
    typedef Collect_face_bbox_per_edge_bbox<TriangleMesh, Edge_to_faces>
      Callback;
    Callback callback(tm_f, tm_e, edge_to_faces);
    #else
    typedef Collect_face_bbox_per_edge_bbox_with_coplanar_handling<
      TriangleMesh, VertexPointMap, Edge_to_faces, Coplanar_face_set>
     Callback;
    Callback  callback(tm_f, tm_e, vpm_f, vpm_e, edge_to_faces, coplanar_faces);
    #endif
    //using pointers in box_intersection_d is about 10% faster
    if (throw_on_self_intersection){
        Callback_with_self_intersection_report<TriangleMesh, Callback> callback_si(callback);
        CGAL::box_intersection_d( face_boxes_ptr.begin(), face_boxes_ptr.end(),
                                  edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                                  callback_si, cutoff );
        if (callback_si.self_intersections_found())
         throw Self_intersection_exception();
    }
    else
        CGAL::box_intersection_d( face_boxes_ptr.begin(), face_boxes_ptr.end(),
                              edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                              callback, cutoff );
  }

  // for autorefinement
  void filter_intersections(const TriangleMesh& tm,
                            const VertexPointMap& vpm)
  {
    std::vector<Box> face_boxes, edge_boxes;
    std::vector<Box*> face_boxes_ptr, edge_boxes_ptr;

    face_boxes.reserve(num_faces(tm));
    face_boxes_ptr.reserve(num_faces(tm));
    BOOST_FOREACH(face_descriptor fd, faces(tm))
    {
      halfedge_descriptor h=halfedge(fd,tm);
      face_boxes.push_back( Box(
        get(vpm,source(h,tm)).bbox() +
        get(vpm,target(h,tm)).bbox() +
        get(vpm,target(next(h,tm),tm)).bbox(),
        h ) );
      face_boxes_ptr.push_back( &face_boxes.back() );
    }

    edge_boxes.reserve(num_edges(tm));
    edge_boxes_ptr.reserve(num_edges(tm));
    BOOST_FOREACH(edge_descriptor ed, edges(tm))
    {
      halfedge_descriptor h=halfedge(ed,tm);
      edge_boxes.push_back( Box(
        get(vpm,source(h,tm)).bbox() +
        get(vpm,target(h,tm)).bbox(),
        h ) );
      edge_boxes_ptr.push_back( &edge_boxes.back() );
    }

    /// \todo experiments different cutoff values
    std::ptrdiff_t cutoff = 2 * std::ptrdiff_t(
        std::sqrt(face_boxes.size()+edge_boxes.size()) );

    Edge_to_faces& edge_to_faces = stm_edge_to_ltm_faces;

    typedef Collect_face_bbox_per_edge_bbox_with_coplanar_handling_one_mesh<
      TriangleMesh, VertexPointMap, Edge_to_faces, Coplanar_face_set>
     Callback;
    Callback  callback(tm, vpm, edge_to_faces, coplanar_faces);

    //using pointers in box_intersection_d is about 10% faster
    CGAL::box_intersection_d( face_boxes_ptr.begin(), face_boxes_ptr.end(),
                              edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                              callback, cutoff );
  }

  template<class Cpl_inter_pt,class Key>
  std::pair<Node_id,bool>
  get_or_create_node(const Cpl_inter_pt& ipt,
                           Node_id& current_node,
                           std::map<Key,Node_id>& coplanar_node_map,
                     const TriangleMesh& tm1,
                     const TriangleMesh& tm2)
  {
    halfedge_descriptor h1=graph_traits::null_halfedge(),h2=h1;
    switch(ipt.type_1){
      case ON_VERTEX:
        h1=halfedge(target(ipt.info_1,tm1),tm1);
      break;
      case ON_EDGE  :
      {
        h1=opposite(ipt.info_1,tm1);
        if (h1>ipt.info_1) h1=ipt.info_1;
      }
      break;
      case ON_FACE :
        h1=halfedge(face(ipt.info_1,tm1),tm1);
      break;
      default: CGAL_error_msg("Should not get there!");
    }
    switch(ipt.type_2){
      case ON_VERTEX:
        h2=halfedge(target(ipt.info_2,tm2),tm2);
      break;
      case ON_EDGE  :
      {
        h2=opposite(ipt.info_2,tm2);
        if (h2>ipt.info_2) h2=ipt.info_2;
      }
      break;
      case ON_FACE :
        h2=halfedge(face(ipt.info_2,tm2),tm2);
      break;
      default: CGAL_error_msg("Should not get there!");
    }

    Key key(ipt.type_1, ipt.type_2, h1, h2);
    if (&tm1==&tm2 && h1>h2)
      key=Key(ipt.type_2, ipt.type_1, h2, h1);

    std::pair<typename std::map<Key,Node_id>::iterator,bool> res=
      coplanar_node_map.insert(std::make_pair(key,current_node+1));
    if (res.second){ //insert a new node
      nodes.add_new_node(ipt.point);
      return std::pair<Node_id,bool>(++current_node, true);
    }
    return std::pair<Node_id,bool>(res.first->second, false);
  }

  void add_intersection_point_to_face_and_all_edge_incident_faces(face_descriptor f_1,
                                                                  halfedge_descriptor e_2,
                                                                  const TriangleMesh& tm1,
                                                                  const TriangleMesh& tm2,
                                                                  Node_id node_id)
  {
    if (!is_border(e_2, tm2))
    {
      face_descriptor f_2 = face(e_2, tm2);
      if(&tm1!=&tm2 || f_1!=f_2)
      {
        Face_pair face_pair = &tm1==&tm2 ? make_sorted_pair(f_1,f_2):
                                           &tm1<&tm2
                                           ? Face_pair(f_1,f_2)
                                           : Face_pair(f_2,f_1);
        if ( coplanar_faces.count(face_pair)==0 )
          f_to_node[Face_pair_and_int(face_pair,0)].insert(node_id);
      }
    }
    e_2 = opposite(e_2, tm2);
    if (!is_border(e_2, tm2))
    {
      face_descriptor f_2 = face(e_2, tm2);
      if(&tm1!=&tm2 || f_1!=f_2)
      {
        Face_pair face_pair = &tm1==&tm2 ? make_sorted_pair(f_1,f_2):
                                           &tm1<&tm2
                                           ? Face_pair(f_1,f_2)
                                           : Face_pair(f_2,f_1);
        if ( coplanar_faces.count(face_pair) == 0 )
          f_to_node[Face_pair_and_int(face_pair,0)].insert(node_id);
      }
    }
  }

  void cip_handle_case_edge(Node_id node_id,
                            Face_set* fset,
                            halfedge_descriptor e_1,
                            halfedge_descriptor edge_intersected,
                            const TriangleMesh& tm1,
                            const TriangleMesh& tm2)
  {
    //associate the intersection point to all faces incident to the intersected edge using edge
    std::vector<face_descriptor> incident_faces;

    if (!is_border(edge_intersected, tm2))
    {
      face_descriptor f_2 = face(edge_intersected, tm2);
      add_intersection_point_to_face_and_all_edge_incident_faces(f_2,e_1,tm2,tm1,node_id);
      if (fset!=NULL) fset->erase(f_2);
    }
    edge_intersected = opposite(edge_intersected, tm2);
    if (!is_border(edge_intersected, tm2))
    {
      face_descriptor f_2 = face(edge_intersected, tm2);
      add_intersection_point_to_face_and_all_edge_incident_faces(f_2,e_1,tm2,tm1,node_id);
      if (fset!=NULL) fset->erase(f_2);
    }

    //associate the intersection point to all faces incident to edge using the intersected edge
    //at least one pair of faces is already handle above

    Edge_to_faces& tm2_edge_to_tm1_faces = &tm1 < &tm2
                                         ? ltm_edge_to_stm_faces
                                         : stm_edge_to_ltm_faces;

    typename Edge_to_faces::iterator it_fset=tm2_edge_to_tm1_faces.find(edge(edge_intersected,tm2));
    if (it_fset==tm2_edge_to_tm1_faces.end()) return;
    Face_set& fset_bis=it_fset->second;

    if (!is_border(e_1, tm1))
    {
      //the following call is not needed, it is already done in the first loop
      //add_intersection_point_to_face_and_all_edge_incident_faces(f_1,edge_intersected,tm1,tm2,node_id);
      fset_bis.erase(face(e_1, tm1));
    }
    e_1 = opposite(e_1, tm1);
    if (!is_border(e_1, tm1))
    {
      //the following call is not needed, it is already done in the first loop
      //add_intersection_point_to_face_and_all_edge_incident_faces(f_1,edge_intersected,tm1,tm2,node_id);
      fset_bis.erase(face(e_1, tm1));
    }
  }

  void cip_handle_case_vertex(Node_id node_id,
                              Face_set* fset,
                              halfedge_descriptor edge,
                              halfedge_descriptor vertex_intersected,
                              const TriangleMesh& tm1,
                              const TriangleMesh& tm2)
  {
    BOOST_FOREACH(halfedge_descriptor h_2,
                  halfedges_around_target(vertex_intersected,tm2))
    {
      cip_handle_case_edge(node_id,fset,edge,h_2,tm1,tm2);
    }
  }

  void handle_coplanar_case_VERTEX_FACE(halfedge_descriptor v_1,
                                        halfedge_descriptor f_2,
                                        const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        Node_id node_id,
                                        bool is_new_node)
  {
    if(is_new_node)
      visitor.new_node_added(node_id,ON_FACE,v_1,f_2,tm1,tm2,true,false);

    Edge_to_faces& tm1_edge_to_tm2_faces = &tm1 <= &tm2
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    BOOST_FOREACH(halfedge_descriptor h_1,
                  halfedges_around_target(v_1,tm1))
    {
      add_intersection_point_to_face_and_all_edge_incident_faces(face(f_2,tm2),h_1,tm2,tm1,node_id);
      typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(h_1,tm1));
      if (it_ets!=tm1_edge_to_tm2_faces.end()) it_ets->second.erase(face(f_2,tm2));
    }
  }

  void handle_coplanar_case_VERTEX_EDGE(halfedge_descriptor v_1,
                                        halfedge_descriptor e_2,
                                        const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        Node_id node_id,
                                        bool is_new_node)
  {
    if(is_new_node)
      visitor.new_node_added(node_id,ON_VERTEX,e_2,v_1,tm2,tm1,false,false);

    Edge_to_faces& tm1_edge_to_tm2_faces = &tm1 <= &tm2
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    BOOST_FOREACH(halfedge_descriptor h_1,
                  halfedges_around_target(v_1,tm1))
    {
      typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(h_1,tm1));
      Face_set* fset = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):NULL;
      cip_handle_case_edge(node_id,fset,h_1,e_2,tm1,tm2);
    }
  }

  void handle_coplanar_case_VERTEX_VERTEX(halfedge_descriptor v_1,
                                          halfedge_descriptor v_2,
                                          const TriangleMesh& tm1,
                                          const TriangleMesh& tm2,
                                          Node_id node_id,
                                          bool is_new_node)
  {
    if(is_new_node)
      visitor.new_node_added(node_id,ON_VERTEX,v_2,v_1,tm2,tm1,true,false);

    Edge_to_faces& tm1_edge_to_tm2_faces = &tm1 <= &tm2
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    BOOST_FOREACH(halfedge_descriptor h_1,
                  halfedges_around_target(v_1,tm1))
    {
      typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(h_1,tm1));
      Face_set* fset = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):NULL;
      cip_handle_case_vertex(node_id,fset,h_1,v_2,tm1,tm2);
    }
  }

  void compute_intersection_of_coplanar_faces(
    Node_id& current_node,
    const TriangleMesh& tm1,
    const TriangleMesh& tm2,
    const VertexPointMap& vpm1,
    const VertexPointMap& vpm2)
  {
    CGAL_assertion( &tm1 < &tm2 || &tm1==&tm2 );

    typedef cpp11::tuple<Intersection_type,
                         Intersection_type,
                         halfedge_descriptor,
                         halfedge_descriptor> Key;

    typedef std::map<Key,Node_id> Coplanar_node_map;
    Coplanar_node_map coplanar_node_map;

    BOOST_FOREACH(const Face_pair& face_pair, coplanar_faces)
    {
      face_descriptor f1=face_pair.first;
      face_descriptor f2=face_pair.second;

      CGAL_assertion(&tm1!=&tm2 || f1!=f2);

      typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
      typedef Coplanar_intersection<TriangleMesh, EK> Cpl_inter_pt;
      std::list<Cpl_inter_pt> inter_pts;

      // compute the intersection points between the two coplanar faces
      intersection_coplanar_faces(f1, f2, tm1, tm2, vpm1, vpm2, inter_pts);

      std::size_t nb_pts=inter_pts.size();
      std::vector<Node_id> cpln_nodes; cpln_nodes.reserve(nb_pts);

      BOOST_FOREACH(const Cpl_inter_pt& ipt, inter_pts)
      {
        Node_id node_id;
        bool is_new_node;
        cpp11::tie(node_id, is_new_node) =
            get_or_create_node(ipt,current_node,coplanar_node_map,tm1,tm2);
        cpln_nodes.push_back(node_id);

        switch(ipt.type_1){
        case ON_VERTEX:
        {
          switch(ipt.type_2){
          case ON_VERTEX:
            handle_coplanar_case_VERTEX_VERTEX(ipt.info_1,ipt.info_2,tm1,tm2,node_id,is_new_node);
          break;
          case ON_EDGE:
            handle_coplanar_case_VERTEX_EDGE(ipt.info_1,ipt.info_2,tm1,tm2,node_id,is_new_node);
          break;
          case ON_FACE:
            handle_coplanar_case_VERTEX_FACE(ipt.info_1,ipt.info_2,tm1,tm2,node_id,is_new_node);
          break;
          default: CGAL_error_msg("Should not get there!");
          }
        }
        break;
        case ON_EDGE:
        {
          switch(ipt.type_2){
            case ON_VERTEX:
              handle_coplanar_case_VERTEX_EDGE(ipt.info_2,ipt.info_1,tm2,tm1,node_id,is_new_node);
            break;
            case ON_EDGE:
            {
              if(is_new_node)
                visitor.new_node_added(node_id,ON_EDGE,ipt.info_1,ipt.info_2,tm1,tm2,false,false);
              typename Edge_to_faces::iterator it_ets=stm_edge_to_ltm_faces.find(edge(ipt.info_1,tm1));
              Face_set* fset = (it_ets!=stm_edge_to_ltm_faces.end())?&(it_ets->second):NULL;
              cip_handle_case_edge(node_id,fset,ipt.info_1,ipt.info_2,tm1,tm2);
            }
            break;
            default: CGAL_error_msg("Should not get there!");
          }
        }
        break;
        case ON_FACE:
        {
          CGAL_assertion(ipt.type_2==ON_VERTEX);
          handle_coplanar_case_VERTEX_FACE(ipt.info_2,ipt.info_1,tm2,tm1,node_id,is_new_node);
        }
        break;
        default: CGAL_error_msg("Should not get there!");
        }
      }

      switch (nb_pts){
        case 0: break;
        case 1:
        {
            f_to_node[Face_pair_and_int(face_pair,1)].insert(cpln_nodes[0]);
        }
        break;
        default:
        {
          std::size_t stop=nb_pts + (nb_pts<3?-1:0);
          for (std::size_t k=0;k<stop;++k){
            Node_id_set& node_id_set = f_to_node[Face_pair_and_int(face_pair,static_cast<int>(k+1))];
            node_id_set.insert( cpln_nodes[k] );
            node_id_set.insert( cpln_nodes[(k+1)%nb_pts] );
          }
        }
      }
    }
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  void add_new_node(halfedge_descriptor h_1,
                    face_descriptor f_2,
                    const TriangleMesh& tm1,
                    const TriangleMesh& tm2,
                    const VertexPointMap& vpm1,
                    const VertexPointMap& vpm2,
                    cpp11::tuple<Intersection_type,
                                 halfedge_descriptor,
                                 bool,bool> inter_res)
  {
    if ( cpp11::get<3>(inter_res) ) // is edge target in triangle plane
      nodes.add_new_node(get(vpm1, target(h_1,tm1)));
    else{
      if (cpp11::get<3>(inter_res)) // is edge source in triangle plane
        nodes.add_new_node(get(vpm1, source(h_1,tm1)));
      else
        nodes.add_new_node(h_1,f_2,tm1,tm2,vpm1,vpm2);
    }
  }

  void compute_intersection_points(Edge_to_faces& tm1_edge_to_tm2_faces,
                                   const TriangleMesh& tm1,
                                   const TriangleMesh& tm2,
                                   const VertexPointMap& vpm1,
                                   const VertexPointMap& vpm2,
                                   Node_id& current_node)
  {
    typedef cpp11::tuple<Intersection_type, halfedge_descriptor, bool,bool>  Inter_type;


    for(typename Edge_to_faces::iterator it=tm1_edge_to_tm2_faces.begin();
                                         it!=tm1_edge_to_tm2_faces.end();++it)
    {
      edge_descriptor e_1=it->first;
      halfedge_descriptor h_1=halfedge(e_1,tm1);
      Face_set& fset=it->second;
      while (!fset.empty()){
        face_descriptor f_2=*fset.begin();

        Inter_type res=intersection_type(h_1,f_2,tm1,tm2,vpm1,vpm2);
        Intersection_type type=cpp11::get<0>(res);

    //handle degenerate case: one extremity of edge belong to f_2
        std::vector<halfedge_descriptor> all_edges;
        if ( cpp11::get<3>(res) ) // is edge target in triangle plane
          std::copy(halfedges_around_target(h_1,tm1).first,
                    halfedges_around_target(h_1,tm1).second,
                    std::back_inserter(all_edges));
        else{
          if ( cpp11::get<2>(res) ) // is edge source in triangle plane
              std::copy(halfedges_around_source(h_1,tm1).first,
                        halfedges_around_source(h_1,tm1).second,
                        std::back_inserter(all_edges));
          else
            all_edges.push_back(h_1);
        }

        CGAL_precondition(all_edges[0]==h_1 || all_edges[0]==opposite(h_1,tm1));

        // #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
        // check_coplanar_edges(cpp11::next(all_edges.begin()),
        //                      all_edges.end(),CGAL::cpp11::get<1>(res),type);
        // #endif

        typename std::vector<halfedge_descriptor>::iterator it_edge=all_edges.begin();
        switch(type){
          case COPLANAR_TRIANGLES:
            #ifndef DO_NOT_HANDLE_COPLANAR_FACES
            assert(!"COPLANAR_TRIANGLES : this point should never be reached!");
            #else
            //nothing needs to be done, cf. comments at the beginning of the file
            #endif
          break;
          case EMPTY:
            fset.erase(fset.begin());
          break;

          // Case when the edge pierces the face in its interior.
          case ON_FACE:
          {
            CGAL_assertion(f_2==face(cpp11::get<1>(res),tm2));

            Node_id node_id=++current_node;
            add_new_node(h_1,f_2,tm1,tm2,vpm1,vpm2,res);
            visitor.new_node_added(node_id,ON_FACE,h_1,halfedge(f_2,tm2),tm1,tm2,CGAL::cpp11::get<3>(res),CGAL::cpp11::get<2>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              add_intersection_point_to_face_and_all_edge_incident_faces(f_2,*it_edge,tm2,tm1,node_id);
              //erase face from the list to test intersection with it_edge
              if ( it_edge==all_edges.begin() )
                fset.erase(fset.begin());
              else
              {
                typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(*it_edge,tm1));
                if(it_ets!=tm1_edge_to_tm2_faces.end()) it_ets->second.erase(f_2);
              }
            }
          } // end case ON_FACE
          break;

          // Case when the edge intersect one edge of the face.
          case ON_EDGE:
          {
            Node_id node_id=++current_node;
            add_new_node(h_1,f_2,tm1,tm2,vpm1,vpm2,res);
            halfedge_descriptor h_2=cpp11::get<1>(res);
            visitor.new_node_added(node_id,ON_EDGE,h_1,h_2,tm1,tm2,cpp11::get<3>(res),cpp11::get<2>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              if ( it_edge!=all_edges.begin() ){
                typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(*it_edge,tm1));
                Face_set* fset_bis = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):NULL;
                cip_handle_case_edge(node_id,fset_bis,*it_edge,h_2,tm1,tm2);
              }
              else
                cip_handle_case_edge(node_id,&fset,*it_edge,h_2,tm1,tm2);
            }
          } // end case ON_EDGE
          break;

          case ON_VERTEX:
          {
            Node_id node_id=++current_node;
            halfedge_descriptor h_2=cpp11::get<1>(res);
            nodes.add_new_node(get(vpm2, target(h_2,tm2))); //we use the original vertex to create the node
            //before it was ON_FACE but do not remember why, probably a bug...
            visitor.new_node_added(node_id,ON_VERTEX,h_1,h_2,tm1,tm2,cpp11::get<3>(res),cpp11::get<2>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              if ( it_edge!=all_edges.begin() ){
                typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(*it_edge,tm1));
                Face_set* fset_bis = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):NULL;
                cip_handle_case_vertex(node_id,fset_bis,*it_edge,h_2,tm1,tm2);
              }
              else
                cip_handle_case_vertex(node_id,&fset,*it_edge,h_2,tm1,tm2);
            }
          } // end case ON_VERTEX
          break;
        } // end switch on the type of the intersection
      } // end loop on all faces that intersect the edge
    } // end loop on all entries (edges) in 'edge_to_face'
    CGAL_assertion(nodes.size()==unsigned(current_node+1));
  }

  struct Graph_node{
    std::set<Node_id> neighbors;
    unsigned degree;

    Graph_node():degree(0){}

    void insert(Node_id i){
      ++degree;
      CGAL_assertion(!neighbors.count(i));
      neighbors.insert(i);
    }

    void erase(Node_id i){
      CGAL_assertion(neighbors.count(i)!= 0);
      neighbors.erase(i);
    }
    void make_terminal() {if (degree==2) degree=45;}
    bool is_terminal()const {return degree!=2;}
    bool empty() const {return neighbors.empty();}
    Node_id top() const {return *neighbors.begin();}
    void pop() {
      CGAL_assertion(!neighbors.empty());
      neighbors.erase(neighbors.begin());
    }
  };

  /// TODO AUTOREF_TAG replace this by a lexical sort
  struct Less_for_nodes_along_an_edge{
    Node_vector& nodes;
    Node_id ref;
    Less_for_nodes_along_an_edge(Node_vector& nodes, Node_id ref)
      : nodes(nodes), ref(ref)
    {}
    bool operator()(Node_id i, Node_id j) const
    {
      return
        compare_distance_to_point(nodes.exact_node(ref),
                                  nodes.exact_node(i),
                                  nodes.exact_node(j) ) == SMALLER;
    }
  };

  void detect_intersections_in_the_graph(const TriangleMesh& tm,
                                         const VertexPointMap& vpm,
                                         Node_id& current_node)
  {
    boost::unordered_map<face_descriptor,
                         std::vector<face_descriptor> > face_intersections;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();
                                               it!=f_to_node.end();
                                               ++it)
    {
      face_descriptor f1 = it->first.first.first, f2 = it->first.first.second;
      CGAL_assertion( f1 < f2 );
      std::vector<face_descriptor>& inter_faces=face_intersections[f1];
      if (inter_faces.empty() || inter_faces.back()!=f2)
        face_intersections[f1].push_back(f2);
    }

    // TODO AUTOREF_TAG find a better way to avoid too many queries.
    std::map< std::pair< std::pair<face_descriptor, face_descriptor>, int>,
                         std::vector<Node_id> > map_to_process;

    typedef std::pair<const face_descriptor, std::vector<face_descriptor> > Pair_type;
    BOOST_FOREACH(Pair_type& p, face_intersections)
    {
      face_descriptor f1 = p.first;
      std::vector<face_descriptor>& inter_faces=p.second;

      std::size_t nb_faces = inter_faces.size();
      // TODO AUTOREF_TAG handle 4 and more faces intersecting (only 3 right now)
      for(std::size_t i=0; i<nb_faces-1; ++i)
      {
        face_descriptor f2 = p.second[i];
        CGAL_assertion(f1 < f2);

        for(std::size_t j=i+1; j<nb_faces;++j)
        {
          face_descriptor f3 = p.second[j];
          CGAL_assertion(f2 < f3);
          // use lower bound to get the first entry which key is equal or larger
          // the non-coplanar entry case. That way any coplanar entry will be
          // listed afterward
          typename Faces_to_nodes_map::iterator find_it =
            f_to_node.lower_bound(std::make_pair(std::make_pair(f2, f3),0));

          if (find_it!=f_to_node.end() &&
              find_it->first.first.first  == f2 &&
              find_it->first.first.second == f3)
          {
            // look for edges between f1 and f2
            typename Faces_to_nodes_map::iterator it_seg1 =
              f_to_node.lower_bound(std::make_pair(std::make_pair(f1, f2),0));
            CGAL_assertion( it_seg1 != f_to_node.end() &&
                            it_seg1->first.first.first  == f1 &&
                            it_seg1->first.first.second == f2 );

            // look for edges between f1 and f3
            typename Faces_to_nodes_map::iterator it_seg2 =
              f_to_node.lower_bound(std::make_pair(std::make_pair(f1, f3),0));
            CGAL_assertion( it_seg2 != f_to_node.end() &&
                            it_seg2->first.first.first  == f1 &&
                            it_seg2->first.first.second == f3 );

            while(it_seg1 != f_to_node.end() &&
                  it_seg1->first.first.first  == f1 &&
                  it_seg1->first.first.second == f2)
            {
              while(it_seg2 != f_to_node.end() &&
                    it_seg2->first.first.first  == f1 &&
                    it_seg2->first.first.second == f3 )
              {
                const Node_id_set& ns1 = it_seg1->second;
                const Node_id_set& ns2 = it_seg2->second;
                CGAL_assertion (ns1.size()!=0);
                CGAL_assertion (ns2.size()!=0);

                if (ns1.size()==1 || ns2.size()==1){
                  ++it_seg2;
                  continue; /// TODO AUTOREF_TAG shall we ignore tangency points?
                            /// Actually it might be that it is not a tangency point if the third segment was considered!
                            /// so not handling it is a bug
                }

                // TODO AUTOREF_TAG there might be a better test rather than relying on constructions
                typedef typename Node_vector::Exact_kernel::Segment_3 Segment_3;
                Segment_3 s1(nodes.exact_node(ns1[0]), nodes.exact_node(ns1[1])),
                          s2(nodes.exact_node(ns2[0]), nodes.exact_node(ns2[1])),
                          s3;

                if( do_intersect(s1, s2) )
                {
                  /// TODO AUTOREF_TAG it might be the end point of a segment!
                  // we need to find which segment is the third one intersecting
                  // with the point found.
                  typename Faces_to_nodes_map::iterator it_seg3 = find_it;
                  // first check if there is only one such edge (no test is needed then)
                  if (cpp11::next(it_seg3)->first.first == it_seg3->first.first)
                  {
                    while(true)
                    {
                      CGAL_assertion(it_seg3->first.first == find_it->first.first);

                      CGAL_assertion(it_seg3->second.size()!=1 || !"AUTOREF_TAG HANDLE ME");
                      const Node_id_set& ns3 = it_seg3->second;
                      s3=Segment_3(nodes.exact_node(ns3[0]), nodes.exact_node(ns3[1]));
                      if (do_intersect(s1, s3))
                      {
                        CGAL_assertion(do_intersect(s2,s3));
                        break;
                      }
                    }
                  }
                  else
                  {
                    CGAL_assertion(it_seg3->second.size()==2 || !"AUTOREF_TAG TODO: HANDLE ME");
                    s3=Segment_3(nodes.exact_node(it_seg3->second[0]),
                                 nodes.exact_node(it_seg3->second[1]));
                  }
                  CGAL_assertion(it_seg3->first.first == find_it->first.first);

                  // use it_seg3
                  ///TODO AUTOREF_TAG the collinear test could be factorized in the do-intersect
                  /// TODO AUTOREF_TAG if we don't factorise maybe checking for collinearity with
                  ///      cross product would be cheaper?
                  if ( (collinear(s1[0], s1[1], s2[0]) && collinear(s1[0], s1[1], s2[1]) ) ||
                       (collinear(s1[0], s1[1], s3[0]) && collinear(s1[0], s1[1], s3[1]) ) ||
                       (collinear(s2[0], s2[1], s3[0]) && collinear(s2[0], s2[1], s3[1]) ) )
                  {
                    throw Triple_intersection_exception();
                  }

                  // we now considered the refinement of the 3 edges defined
                  // by the intersection of f1, f2, and f3
                  /// TODO AUTOREF_TAG the new intersection point might be the endpoint of an intersection segment!
                  nodes.add_new_node(halfedge(f1, tm),
                                     halfedge(f2, tm),
                                     halfedge(f3, tm),
                                     tm, vpm);
                  Node_id node_id=++current_node;
#ifdef CGAL_DEBUG_AUTOREFINEMENT
                  std::cerr << "New triple node " << node_id << "\n";
#endif
                  visitor.new_node_added_triple_face(node_id, f1, f2, f3, tm);

                  map_to_process[it_seg1->first].push_back(node_id);
                  map_to_process[it_seg2->first].push_back(node_id);
                  map_to_process[it_seg3->first].push_back(node_id);
                }
                ++it_seg2;
              }
              ++it_seg1;
            }
          }
        }
      }
    }

#ifdef CGAL_DEBUG_AUTOREFINEMENT
    std::cout << "\nAt the end of new node creation, current_node is " << current_node << "\n\n";
#endif
    typedef std::pair<const std::pair< std::pair<face_descriptor, face_descriptor>, int>,
                      std::vector<Node_id> > Faces_and_nodes;
    BOOST_FOREACH(Faces_and_nodes& f_n_nids, map_to_process)
    {
      //get the original entry and remove it
      typename Faces_to_nodes_map::iterator find_it =
        f_to_node.find( f_n_nids.first ); // TODO AUTOREF_TAG saving the iterator above would avoid this find
      CGAL_assertion(find_it!=f_to_node.end());
      CGAL_assertion(find_it->second.size()==2); // TODO AUTOREF_TAG handle case size = 1
      Node_id n1 = find_it->second[0];
      Node_id n2 = find_it->second[1];

      // sort node ids along the edge
      std::sort(f_n_nids.second.begin(),
                f_n_nids.second.end(),
                Less_for_nodes_along_an_edge(nodes, n1));

      // insert new segments
      Node_id prev = n1;
      f_n_nids.second.push_back(n2);
      int i=0;
      typename Faces_to_nodes_map::iterator insert_it = find_it;
#ifdef CGAL_DEBUG_AUTOREFINEMENT
      std::cout << n1 << " -> " << n2 << "\n";
#endif
      BOOST_FOREACH(Node_id id, f_n_nids.second)
      {
        insert_it = f_to_node.insert(insert_it, std::make_pair(
          std::make_pair(f_n_nids.first.first,--i), Node_id_set()) ); // I have picked negative int for refined edges
        insert_it->second.insert(prev);
        insert_it->second.insert(id);
        CGAL_assertion(insert_it->second.size()==2);
#ifdef CGAL_DEBUG_AUTOREFINEMENT
        std::cerr <<"  adding " << prev << " " << id << " into "
                  << f_n_nids.first.first.first << " and " << f_n_nids.first.first.second <<  "\n";
#endif
        prev=id;
      }
      f_to_node.erase(find_it);
    }
  }

  template <class Output_iterator>
  void construct_polylines(Output_iterator out){
    typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
    std::size_t nb_nodes=nodes.size();
    std::vector<Graph_node> graph(nb_nodes);
    //counts the number of time each node has been seen
    bool isolated_point_seen=false;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();it!=f_to_node.end();++it){
      const Node_id_set& segment=it->second;
      CGAL_assertion(segment.size()==2 || segment.size()==1);
      if (segment.size()==2){
        Node_id i=segment.first;
        Node_id j=segment.second;
        graph[i].insert(j);
        graph[j].insert(i);
      }
      else{
        CGAL_assertion(segment.size()==1);
        isolated_point_seen=true;
      }
    }

    CGAL_assertion(extra_terminal_nodes.empty() || doing_autorefinement);
    // these nodes are created by pinchements along an edge of the surface.
    // the node ids being the same for the two edges, the degree of the node
    // in the graph is two while it should be 3
    BOOST_FOREACH(Node_id id, extra_terminal_nodes)
      graph[id].make_terminal();

    //visitor call
    visitor.annotate_graph(graph);

    //collect terminal and interior nodes
    boost::dynamic_bitset<> terminal_nodes(nb_nodes), interior_nodes(nb_nodes);
    for (std::size_t i=0;i<nb_nodes;++i)
      if (graph[i].is_terminal())
        terminal_nodes.set(i);
      else
        interior_nodes.set(i);

    //handle isolated points
    if (isolated_point_seen){
      for (std::size_t i=0;i<nb_nodes;++i)
        if (graph[i].degree==0){
          *out++=std::vector<Point_3>(1,nodes[i]);
          visitor.start_new_polyline(i,i);
          terminal_nodes.reset(i);
        }
    }

    //handle polylines
    while(terminal_nodes.any())
    {
      std::size_t i=terminal_nodes.find_first();
      Graph_node& node_i = graph[i];
      std::vector<Point_3> polyline;

      std::size_t j=node_i.top();
      visitor.start_new_polyline(i,j);
      CGAL_assertion(i!=j);
      node_i.pop();
      if (node_i.empty())
        terminal_nodes.reset(i);
      polyline.push_back(nodes[i]);
      while(true){
        Graph_node& node_j=graph[j];
        CGAL_assertion(!node_j.empty());
        node_j.erase(i);
        i=j;
        polyline.push_back(nodes[i]);
        if (node_j.is_terminal())
        {
          if (node_j.empty())
            terminal_nodes.reset(j);
          break;
        }
        else{
          j=node_j.top();
          visitor.add_node_to_polyline(j);
          node_j.pop();
          CGAL_assertion(node_j.empty());
          interior_nodes.reset(i);
        }
      }
      *out++=polyline;
    }

    //handle cycles
    while(interior_nodes.any())
    {
      std::size_t i=interior_nodes.find_first();
      Graph_node& node_i=graph[i];
      std::vector<Point_3> polyline;

      Node_id j=node_i.top();
      visitor.start_new_polyline(i,j);
      interior_nodes.reset(i);
      polyline.push_back(nodes[i]);
      Node_id first=i;
      do{
        Graph_node& node_j=graph[j];
        interior_nodes.reset(j);
        node_j.erase(i);
        i=j;
        polyline.push_back(nodes[i]);
        j=node_j.top();
        visitor.add_node_to_polyline(j);
      }while(j!=first);
      polyline.push_back(nodes[j]);// we duplicate first point for cycles
      *out++=polyline;
    }
  }


  void remove_duplicated_intersecting_edges()
  {
    std::set< std::pair<Node_id,Node_id> > already_seen;
    std::vector<typename Faces_to_nodes_map::iterator> to_erase;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();
          it!=f_to_node.end(); ++it)
    {
      if (it->second.size()==1) continue;
      CGAL_precondition(it->second.size()==2);
      CGAL_assertion((it->second.first)< (it->second.second));
      //it->second is a set so node ids are already sorted
      bool is_new=already_seen.insert(  std::make_pair(
                                       it->second.first,
                                       it->second.second )
                                        ).second;

      if (!is_new)
        to_erase.push_back(it);
    }

    BOOST_FOREACH(typename Faces_to_nodes_map::iterator it, to_erase)
      f_to_node.erase(it);
  }

  void add_common_vertices_for_pairs_of_faces_with_isolated_node(Node_id& current_node)
  {
    const TriangleMesh& tm = nodes.tm1;
    CGAL_assertion(doing_autorefinement);
    std::map<vertex_descriptor, Node_id> vertex_to_node_id;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();
          it!=f_to_node.end(); ++it)
    {
      if (it->second.size()==2) continue;
      CGAL_precondition(it->second.size()==1);
      halfedge_descriptor h1 = halfedge(it->first.first.first, tm);
      halfedge_descriptor h2 = halfedge(it->first.first.second, tm);

      for (int i=0; i<3; ++i)
      {
        for(int j=0; j<3; ++j)
        {
          if ( target(h1, tm)==target(h2,tm) )
          {
            Node_id node_id = current_node+1;
            std::pair< typename std::map<vertex_descriptor, Node_id>::iterator, bool>
              insert_res = vertex_to_node_id.insert(std::make_pair(target(h1,tm), node_id));
            if (insert_res.second)
            {
              ++current_node;
              nodes.add_new_node(get(nodes.vpm1, target(h1,tm)));
              visitor.new_node_added(node_id,ON_VERTEX,h1,h2,tm,tm,true,false);
              extra_terminal_nodes.push_back(node_id);
            }
            else
              node_id = insert_res.first->second;
            it->second.insert(node_id);
            break;
          }
          h2 = next(h2, tm);
        }
        h1 = next(h1, tm);
      }
    }
  }

public:
  Intersection_of_triangle_meshes(const TriangleMesh& tm1,
                                  const TriangleMesh& tm2,
                                  const VertexPointMap& vpm1,
                                  const VertexPointMap& vpm2,
                                  const Node_visitor& v=Node_visitor())
  : nodes(tm1, tm2, vpm1, vpm2)
  , visitor(v)
  {
    CGAL_assertion_code( doing_autorefinement=false; )
  }

  // for autorefinement
  Intersection_of_triangle_meshes(const TriangleMesh& tm,
                                  const VertexPointMap& vpm,
                                  const Node_visitor& v=Node_visitor())
  : nodes(tm, tm, vpm, vpm)
  , visitor(v)
  {
    CGAL_assertion_code( doing_autorefinement=true; )
  }

  template <class OutputIterator>
  OutputIterator operator()(OutputIterator output,
                            bool throw_on_self_intersection,
                            bool build_polylines)
  {
    CGAL_assertion(!doing_autorefinement);

    const TriangleMesh& tm1=nodes.tm1;
    const TriangleMesh& tm2=nodes.tm2;
    const VertexPointMap& vpm1=nodes.vpm1;
    const VertexPointMap& vpm2=nodes.vpm2;

    filter_intersections(tm1, tm2, vpm1, vpm2, throw_on_self_intersection);
    filter_intersections(tm2, tm1, vpm2, vpm1, throw_on_self_intersection);

    Node_id current_node((std::numeric_limits<Node_id>::max)());
    CGAL_assertion(current_node+1==0);

    #ifndef DO_NOT_HANDLE_COPLANAR_FACES
    //first handle coplanar triangles
    if (&tm1<&tm2)
      compute_intersection_of_coplanar_faces(current_node, tm1, tm2, vpm1, vpm2);
    else
      compute_intersection_of_coplanar_faces(current_node, tm2, tm1, vpm2, vpm1);
    visitor.set_number_of_intersection_points_from_coplanar_faces(current_node+1);
    if (!coplanar_faces.empty())
      visitor.input_have_coplanar_faces();
    #endif // not DO_NOT_HANDLE_COPLANAR_FACES

    //compute intersection points of segments and triangles.
    //build the nodes of the graph and connectivity infos
    Edge_to_faces& tm1_edge_to_tm2_faces = (&tm1<&tm2)
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;
    Edge_to_faces& tm2_edge_to_tm1_faces = (&tm1>&tm2)
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    compute_intersection_points(tm1_edge_to_tm2_faces, tm1, tm2, vpm1, vpm2, current_node);
    compute_intersection_points(tm2_edge_to_tm1_faces, tm2, tm1, vpm2, vpm1, current_node);
    if (!build_polylines){
      visitor.finalize(nodes,tm1,tm2,vpm1,vpm2);
      return output;
    }
    //remove duplicated intersecting edges:
    //  In case two faces are incident along an intersection edge coplanar
    //  in a face of another polyhedron (and one extremity inside the face),
    //  the intersection will be reported twice. We kept track
    //  (check_coplanar_edge(s)) of this so that,
    //  we can remove one intersecting edge out of the two
    remove_duplicated_intersecting_edges();

#if 0
    //collect connectivity infos and create polylines
    if ( Node_visitor::do_need_vertex_graph )
#endif
      //using the graph approach (at some point we know all
      // connections between intersection points)
      construct_polylines(output);
#if 0
    else
      construct_polylines_with_info(nodes,out); //direct construction by propagation
#endif

    visitor.finalize(nodes,tm1,tm2,vpm1,vpm2);

    return output;
  }

  //for autorefinement
  template <class OutputIterator>
  OutputIterator operator()(OutputIterator output,
                            bool build_polylines)
  {
    CGAL_assertion(doing_autorefinement);

    const TriangleMesh& tm=nodes.tm1;
    const VertexPointMap& vpm=nodes.vpm1;

    filter_intersections(tm, vpm);

    Node_id current_node((std::numeric_limits<Node_id>::max)());
    CGAL_assertion(current_node+1==0);

    //first handle coplanar triangles
    compute_intersection_of_coplanar_faces(current_node, tm, tm, vpm, vpm);
    if (!coplanar_faces.empty())
      visitor.input_have_coplanar_faces();

    CGAL_assertion(ltm_edge_to_stm_faces.empty());

    //compute intersection points of segments and triangles.
    //build the nodes of the graph and connectivity infos
    compute_intersection_points(stm_edge_to_ltm_faces, tm, tm, vpm, vpm, current_node);

    if (!build_polylines){
      visitor.finalize(nodes,tm,tm,vpm,vpm);
      return output;
    }
    //remove duplicated intersecting edges:
    //  In case two faces are incident along an intersection edge coplanar
    //  in a face of another polyhedron (and one extremity inside the face),
    //  the intersection will be reported twice. We kept track
    //  (check_coplanar_edge(s)) of this so that,
    //  we can remove one intersecting edge out of the two
/// TODO AUTOREF_TAG does this happen in coplanar cases only? + shall we do it have new edge splitting?
    remove_duplicated_intersecting_edges();


    // If a pair of faces defines an isolated node, check if they share a common
    // vertex and create a new node in that case.
    add_common_vertices_for_pairs_of_faces_with_isolated_node(current_node);

    detect_intersections_in_the_graph(tm, vpm, current_node);
#if 0
    //collect connectivity infos and create polylines
    if ( Node_visitor::do_need_vertex_graph )
#endif
      //using the graph approach (at some point we know all
      // connections between intersection points)
      construct_polylines(output);
#if 0
    else
      construct_polylines_with_info(nodes,out); //direct construction by propagation
#endif

    visitor.finalize(nodes,tm,tm,vpm,vpm);

    return output;
  }

};

} } // end of namespace CGAL::Corefinement

#endif //CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_IMPL_H
