// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_callbacks.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Intersection_type.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_of_coplanar_triangles_3.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_nodes.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersect_triangle_and_segment_3.h>
#include <CGAL/Polygon_mesh_processing/Non_manifold_feature_map.h>
#include <CGAL/utility.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/functional/hash.hpp>

#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace CGAL{
namespace Polygon_mesh_processing {
namespace Corefinement {

struct Triple_intersection_exception :
  public std::runtime_error
{
  Triple_intersection_exception()
    : std::runtime_error("Non-handled triple intersection of input triangles")
  {}
};
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
  template<class T>
  void check_no_duplicates(const T&){}
  template<class T,class VPM1,class VPM2>
  void finalize(T&,
                const TriangleMesh&, const TriangleMesh&,
                const VPM1, const VPM2)
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
  void set_non_manifold_feature_map(
    const TriangleMesh&,
    const Non_manifold_feature_map<TriangleMesh>&)
  {}

  // needed for progress tracking
  void start_filtering_intersections() const {}
  void progress_filtering_intersections(double) const{}
  void end_filtering_intersections() const {}

  void start_triangulating_faces(std::size_t) const {}
  void triangulating_faces_step(std::size_t) const {}
  void end_triangulating_faces() const {}

  void start_handling_intersection_of_coplanar_faces(std::size_t) const {}
  void intersection_of_coplanar_faces_step() const {}
  void end_handling_intersection_of_coplanar_faces() const {}

  void start_handling_edge_face_intersections(std::size_t) const {}
  void edge_face_intersections_step() const {}
  void end_handling_edge_face_intersections() const {}

  void start_building_output() const {}
  void end_building_output() const {}
};

struct Node_id_set {
  typedef std::size_t Node_id;

  Node_id first;
  Node_id second;
  std::size_t size_;

  std::vector<std::array<std::size_t, 2> > coplanar_segments;

  Node_id_set()
    : size_(0)
  {}

  void insert(std::size_t i, std::size_t j)
  {
    if (j<i) std::swap(i,j);
    coplanar_segments.push_back( {i,j} );
  }

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

  void get_segments(std::vector<std::array<std::size_t,2> >& ids)
  {
    ids.reserve( (size_==2?1:0) + coplanar_segments.size() );
    if (size_==2)
      ids.push_back({first, second});
    std::copy(coplanar_segments.begin(), coplanar_segments.end(), std::back_inserter(ids));
  }
};

template< class TriangleMesh,
          class VertexPointMap1, class VertexPointMap2,
          class Node_visitor=Default_surface_intersection_visitor<TriangleMesh>
         >
class Intersection_of_triangle_meshes
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::face_descriptor face_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor, Box_policy> Box;

  typedef std::unordered_set<face_descriptor> Face_set;
  typedef std::unordered_map<edge_descriptor, Face_set> Edge_to_faces;

  static const bool Predicates_on_constructions_needed =
    Node_visitor::Predicates_on_constructions_needed;

  typedef std::pair<face_descriptor, face_descriptor> Face_pair;
  typedef std::set< Face_pair > Coplanar_face_set;

  typedef std::size_t Node_id;

  // we use Face_pair_and_int and not Face_pair to handle coplanar case.
  // Indeed the boundary of the intersection of two coplanar triangles
  // may contain several segments.
  typedef std::unordered_map< Face_pair, Node_id_set, boost::hash<Face_pair> >    Faces_to_nodes_map;
  typedef Intersection_nodes<TriangleMesh,
                             VertexPointMap1, VertexPointMap2,
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
  Non_manifold_feature_map<TriangleMesh> non_manifold_feature_map_1,
                                         non_manifold_feature_map_2;
  const TriangleMesh* const_mesh_ptr;
  static const constexpr std::size_t NM_NID = (std::numeric_limits<std::size_t>::max)();
  CGAL_assertion_code(bool doing_autorefinement;)

// member functions
  template <class VPMF, class VPME>
  void filter_intersections(const TriangleMesh& tm_f,
                            const TriangleMesh& tm_e,
                            const VPMF& vpm_f,
                            const VPME& vpm_e,
                            const Non_manifold_feature_map<TriangleMesh>& non_manifold_feature_map,
                            bool throw_on_self_intersection,
                            std::set<face_descriptor>& tm_f_faces,
                            std::set<face_descriptor>& tm_e_faces,
                            bool run_check)
  {
    std::vector<Box> face_boxes, edge_boxes;
    std::vector<Box*> face_boxes_ptr, edge_boxes_ptr;

    face_boxes.reserve(num_faces(tm_f));
    face_boxes_ptr.reserve(num_faces(tm_f));
    for(face_descriptor fd : faces(tm_f))
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
    if (non_manifold_feature_map.non_manifold_edges.empty())
      // general manifold case
      for(edge_descriptor ed : edges(tm_e))
      {
        halfedge_descriptor h=halfedge(ed,tm_e);
        edge_boxes.push_back( Box(
          get(vpm_e,source(h,tm_e)).bbox() +
          get(vpm_e,target(h,tm_e)).bbox(),
          h ) );
        edge_boxes_ptr.push_back( &edge_boxes.back() );
      }
    else
      // non-manifold case
      for(edge_descriptor ed : edges(tm_e))
      {
        std::size_t eid=get(non_manifold_feature_map.e_nm_id, ed);
        halfedge_descriptor h=halfedge(ed,tm_e);
        // insert only one copy of a non-manifold edge
        if (eid!=NM_NID)
        {
          if (non_manifold_feature_map.non_manifold_edges[eid].front()!=ed)
            continue;
          else
            // make sure the halfedge used is consistent with stored one
            h = halfedge(non_manifold_feature_map.non_manifold_edges[eid].front(), tm_e);
        }
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
      TriangleMesh, VPMF, VPME, Edge_to_faces, Coplanar_face_set, Node_visitor>
     Callback;
    Callback  callback(tm_f, tm_e, vpm_f, vpm_e, edge_to_faces, coplanar_faces, visitor);
    #endif
    //using pointers in box_intersection_d is about 10% faster
    if (throw_on_self_intersection){
        Callback_with_self_intersection_report<TriangleMesh, Callback> callback_si(callback, tm_f_faces, tm_e_faces);
        CGAL::box_intersection_d(face_boxes_ptr.begin(), face_boxes_ptr.end(),
                                 edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                                 callback_si, cutoff);
        if (run_check && callback_si.self_intersections_found())
         throw Self_intersection_exception();
    }
    else {
      if (const_mesh_ptr==&tm_e)
      {
        // tm_f might feature degenerate faces
        auto filtered_callback = [&callback](const Box* fb, const Box* eb)
        {
          if (!callback.is_face_degenerated(fb->info()))
            callback(fb, eb);
        };
        CGAL::box_intersection_d( face_boxes_ptr.begin(), face_boxes_ptr.end(),
                                  edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                                  filtered_callback, cutoff );
      }
      else
      {
        if (const_mesh_ptr==&tm_f)
        {
          // tm_e might feature degenerate edges
          auto filtered_callback = [&,this](const Box* fb, const Box* eb)
          {
            if (get(vpm_e, source(eb->info(), tm_e)) != get(vpm_e, target(eb->info(), tm_e)))
              callback(fb, eb);
            else
            {
              halfedge_descriptor hf = fb->info();
              halfedge_descriptor he = eb->info();
              for (int i=0; i<2; ++i)
              {
                if (!is_border(he, tm_e))
                {
                  if ( get(vpm_e, target(next(he, tm_e), tm_e))==get(vpm_e, target(he, tm_e)) &&
                       coplanar(get(vpm_f, source(hf, tm_f)),
                                get(vpm_f, target(hf, tm_f)),
                                get(vpm_f, target(next(hf, tm_f), tm_f)),
                                get(vpm_e, target(he, tm_e))) )
                  {
                    coplanar_faces.insert(
                        &tm_e < &tm_f
                        ? std::make_pair(face(he, tm_e), face(hf, tm_f))
                        : std::make_pair(face(hf, tm_f), face(he, tm_e))
                      );
                  }
                }
                he=opposite(he, tm_e);
              }
            }
          };
          CGAL::box_intersection_d( face_boxes_ptr.begin(), face_boxes_ptr.end(),
                                    edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                                    filtered_callback, cutoff );
        }
        else
          CGAL::box_intersection_d( face_boxes_ptr.begin(), face_boxes_ptr.end(),
                                    edge_boxes_ptr.begin(), edge_boxes_ptr.end(),
                                    callback, cutoff );
      }
    }
  }

  // for autorefinement
  template <class VPM>
  void filter_intersections(const TriangleMesh& tm,
                            const VPM& vpm)
  {
    std::vector<Box> face_boxes, edge_boxes;
    std::vector<Box*> face_boxes_ptr, edge_boxes_ptr;

    face_boxes.reserve(num_faces(tm));
    face_boxes_ptr.reserve(num_faces(tm));
    for(face_descriptor fd : faces(tm))
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
    for(edge_descriptor ed : edges(tm))
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
      TriangleMesh, VPM, Edge_to_faces, Coplanar_face_set>
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
                           const Non_manifold_feature_map<TriangleMesh>& nm_features_map_1,
                           const Non_manifold_feature_map<TriangleMesh>& nm_features_map_2,
                           const TriangleMesh& tm1,
                           const TriangleMesh& tm2)
  {
    halfedge_descriptor h1=graph_traits::null_halfedge(),h2=h1;
    switch(ipt.type_1){
      case ON_VERTEX:
      {
        std::size_t vid1 = nm_features_map_1.non_manifold_vertices.empty()
                         ? NM_NID
                         : get(nm_features_map_1.v_nm_id, target(ipt.info_1,tm1));
        if (vid1==NM_NID)
          h1=halfedge(target(ipt.info_1,tm1),tm1);
        else
          h1=halfedge(nm_features_map_1.non_manifold_vertices[vid1][0],tm1);
      }
      break;
      case ON_EDGE  :
      {
        std::size_t eid1 = nm_features_map_1.non_manifold_edges.empty()
                         ? NM_NID
                         : get(nm_features_map_1.e_nm_id, edge(ipt.info_1,tm1));
        if (eid1==NM_NID)
          h1=ipt.info_1;
        else
          h1=halfedge(nm_features_map_1.non_manifold_edges[eid1][0],tm1);
        h1=(std::max)(h1, opposite(h1, tm1));
      }
      break;
      case ON_FACE :
        h1=halfedge(face(ipt.info_1,tm1),tm1);
      break;
      default: CGAL_error_msg("Should not get there!");
    }
    switch(ipt.type_2){
      case ON_VERTEX:
      {
        std::size_t vid2 = nm_features_map_2.non_manifold_vertices.empty()
                         ? NM_NID
                         : get(nm_features_map_2.v_nm_id, target(ipt.info_2,tm2));
        if (vid2==NM_NID)
          h2=halfedge(target(ipt.info_2,tm2),tm2);
        else
          h2=halfedge(nm_features_map_2.non_manifold_vertices[vid2][0],tm2);
      }
      break;
      case ON_EDGE  :
      {
        std::size_t eid2 = nm_features_map_2.non_manifold_edges.empty()
                         ? NM_NID
                         : get(nm_features_map_2.e_nm_id, edge(ipt.info_2,tm2));
        if (eid2==NM_NID)
          h2=ipt.info_2;
        else
          h2=halfedge(nm_features_map_2.non_manifold_edges[eid2][0],tm2);
        h2=(std::max)(h2, opposite(h2, tm2));
      }
      break;
      case ON_FACE :
        h2=halfedge(face(ipt.info_2,tm2),tm2);
      break;
      default: CGAL_error_msg("Should not get there!");
    }

    Key key(ipt.type_1, ipt.type_2, h1, h2);
    if (&tm1==&tm2 && h2<h1)
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
                                                                  halfedge_descriptor h_2,
                                                                  const TriangleMesh& tm1,
                                                                  const TriangleMesh& tm2,
                                                                  Node_id node_id)
  {
    if (!is_border(h_2, tm2))
    {
      face_descriptor f_2 = face(h_2, tm2);
      if(&tm1!=&tm2 || f_1!=f_2)
      {
        Face_pair face_pair = &tm1==&tm2 ? make_sorted_pair(f_1,f_2):
                                           &tm1<&tm2
                                           ? Face_pair(f_1,f_2)
                                           : Face_pair(f_2,f_1);
        if ( coplanar_faces.count(face_pair)==0 )
          f_to_node[face_pair].insert(node_id);
      }
    }
    h_2 = opposite(h_2, tm2);
    if (!is_border(h_2, tm2))
    {
      face_descriptor f_2 = face(h_2, tm2);
      if(&tm1!=&tm2 || f_1!=f_2)
      {
        Face_pair face_pair = &tm1==&tm2 ? make_sorted_pair(f_1,f_2):
                                           &tm1<&tm2
                                           ? Face_pair(f_1,f_2)
                                           : Face_pair(f_2,f_1);
        if ( coplanar_faces.count(face_pair) == 0 )
          f_to_node[face_pair].insert(node_id);
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
      if (fset!=nullptr) fset->erase(f_2);
    }
    edge_intersected = opposite(edge_intersected, tm2);
    if (!is_border(edge_intersected, tm2))
    {
      face_descriptor f_2 = face(edge_intersected, tm2);
      add_intersection_point_to_face_and_all_edge_incident_faces(f_2,e_1,tm2,tm1,node_id);
      if (fset!=nullptr) fset->erase(f_2);
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
    for(halfedge_descriptor h_2 :
                  halfedges_around_target(vertex_intersected,tm2))
    {
      cip_handle_case_edge(node_id,fset,edge,h_2,tm1,tm2);
    }
  }

  void handle_coplanar_case_VERTEX_FACE(halfedge_descriptor v_1,
                                        halfedge_descriptor f_2,
                                        const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        const Non_manifold_feature_map<TriangleMesh>& nm_features_map_1,
                                        Node_id node_id,
                                        bool is_new_node)
  {
    if(is_new_node)
      visitor.new_node_added(node_id,ON_FACE,v_1,f_2,tm1,tm2,true,false);

    Edge_to_faces& tm1_edge_to_tm2_faces = &tm1 <= &tm2
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    std::vector<vertex_descriptor> tmp_vertices_1(1, target(v_1, tm1));

    std::size_t vid1 = nm_features_map_1.non_manifold_vertices.empty()
      ? NM_NID
      : get(nm_features_map_1.v_nm_id, target(v_1, tm1));

    const std::vector<vertex_descriptor>& vertices_1 = vid1==NM_NID
                                                     ? tmp_vertices_1
                                                     : nm_features_map_1.non_manifold_vertices[vid1];
    for(vertex_descriptor v1 : vertices_1)
      for(halfedge_descriptor h_1 :
                    halfedges_around_target(v1,tm1))
      {
        add_intersection_point_to_face_and_all_edge_incident_faces(face(f_2,tm2),h_1,tm2,tm1,node_id);
        typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(h_1,tm1));
        if (it_ets!=tm1_edge_to_tm2_faces.end()) it_ets->second.erase(face(f_2,tm2));
      }
  }

  void handle_coplanar_case_VERTEX_EDGE(halfedge_descriptor v_1,
                                        halfedge_descriptor h_2,
                                        const TriangleMesh& tm1,
                                        const TriangleMesh& tm2,
                                        const Non_manifold_feature_map<TriangleMesh>& nm_features_map_1,
                                        const Non_manifold_feature_map<TriangleMesh>& nm_features_map_2,
                                        Node_id node_id,
                                        bool is_new_node)
  {
    if(is_new_node)
      visitor.new_node_added(node_id,ON_VERTEX,h_2,v_1,tm2,tm1,false,false);

    Edge_to_faces& tm1_edge_to_tm2_faces = &tm1 <= &tm2
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    std::vector<vertex_descriptor> tmp_vertices_1(1, target(v_1, tm1));

    std::size_t vid1 = nm_features_map_1.non_manifold_vertices.empty()
                     ? NM_NID
                     : get(nm_features_map_1.v_nm_id, target(v_1, tm1));

    const std::vector<vertex_descriptor>& vertices_1 = vid1==NM_NID
                                                      ? tmp_vertices_1
                                                      : nm_features_map_1.non_manifold_vertices[vid1];

    std::vector<edge_descriptor> tmp_edges_2(1, edge(h_2, tm2));

    std::size_t eid2 = nm_features_map_2.non_manifold_edges.empty()
                     ? NM_NID
                     : get(nm_features_map_2.e_nm_id, edge(h_2, tm2));

    const std::vector<edge_descriptor>& edges_2 = eid2==NM_NID
                                                ? tmp_edges_2
                                                : nm_features_map_2.non_manifold_edges[eid2];

    for (vertex_descriptor v1 : vertices_1)
      for(halfedge_descriptor h_1 :
                    halfedges_around_target(v1,tm1))
      {
        typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(h_1,tm1));
        Face_set* fset = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):nullptr;
        for (edge_descriptor e2 : edges_2)
          cip_handle_case_edge(node_id,fset,h_1,halfedge(e2, tm2),tm1,tm2);
      }
  }

  void handle_coplanar_case_VERTEX_VERTEX(halfedge_descriptor v_1,
                                          halfedge_descriptor v_2,
                                          const TriangleMesh& tm1,
                                          const TriangleMesh& tm2,
                                          const Non_manifold_feature_map<TriangleMesh>& nm_features_map_1,
                                          const Non_manifold_feature_map<TriangleMesh>& nm_features_map_2,
                                          Node_id node_id,
                                          bool is_new_node)
  {
    if(is_new_node)
      visitor.new_node_added(node_id,ON_VERTEX,v_2,v_1,tm2,tm1,true,false);

    Edge_to_faces& tm1_edge_to_tm2_faces = &tm1 <= &tm2
                                         ? stm_edge_to_ltm_faces
                                         : ltm_edge_to_stm_faces;

    std::vector<vertex_descriptor> tmp_vertices_1(1, target(v_1, tm1)),
                                   tmp_vertices_2(1, target(v_2, tm2));

    std::size_t vid1 = nm_features_map_1.non_manifold_vertices.empty()
                     ? NM_NID
                     : get(nm_features_map_1.v_nm_id, target(v_1, tm1));

    std::size_t vid2 = nm_features_map_2.non_manifold_vertices.empty()
                     ? NM_NID
                     : get(nm_features_map_2.v_nm_id, target(v_2, tm2));

    const std::vector<vertex_descriptor>& vertices_1 = vid1==NM_NID
                                                     ? tmp_vertices_1
                                                     : nm_features_map_1.non_manifold_vertices[vid1];
    const std::vector<vertex_descriptor>& vertices_2 = vid2==NM_NID
                                                     ? tmp_vertices_2
                                                     : nm_features_map_2.non_manifold_vertices[vid2];

    for (vertex_descriptor v1 : vertices_1)
      for(halfedge_descriptor h_1 :
                    halfedges_around_target(v1,tm1))
      {
        typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(h_1,tm1));
        Face_set* fset = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):nullptr;
        for (vertex_descriptor v2 : vertices_2)
          cip_handle_case_vertex(node_id,fset,h_1,halfedge(v2, tm2),tm1,tm2);
      }
  }

  template <typename VPM1, typename VPM2>
  void compute_intersection_of_coplanar_faces(Node_id& current_node,
                                              const TriangleMesh& tm1,
                                              const TriangleMesh& tm2,
                                              const VPM1& vpm1,
                                              const VPM2& vpm2,
                                              const Non_manifold_feature_map<TriangleMesh>& nm_features_map_1,
                                              const Non_manifold_feature_map<TriangleMesh>& nm_features_map_2)
  {
    CGAL_assertion( &tm1 < &tm2 || &tm1==&tm2 );

    typedef std::tuple<Intersection_type,
                         Intersection_type,
                         halfedge_descriptor,
                         halfedge_descriptor> Key;

    typedef std::map<Key,Node_id> Coplanar_node_map;
    Coplanar_node_map coplanar_node_map;

    visitor.start_handling_intersection_of_coplanar_faces(coplanar_faces.size());
    for(const Face_pair& face_pair : coplanar_faces)
    {
      visitor.intersection_of_coplanar_faces_step();
      face_descriptor f1=face_pair.first;
      face_descriptor f2=face_pair.second;

      typedef typename Node_vector::Exact_kernel EK;
      typedef Coplanar_intersection<TriangleMesh, EK> Cpl_inter_pt;
      std::list<Cpl_inter_pt> inter_pts;

      //handle degenerate faces
      if (const_mesh_ptr)
      {
        halfedge_descriptor h2 = halfedge(f2,tm2);
        if (const_mesh_ptr == &tm1)
        {
          const typename boost::property_traits<VPM2>::reference
            a = get(vpm2, source(h2, tm2)),
            b = get(vpm2, target(h2, tm2)),
            c = get(vpm2, target(next(h2, tm2), tm2));

          if (collinear(a, b, c))
          {
            intersection_coplanar_faces(f2, f1, tm2, tm1, vpm2, vpm1, inter_pts);
            for (Cpl_inter_pt& ipt : inter_pts)
            {
              std::swap(ipt.type_1,ipt.type_2);
              std::swap(ipt.info_1,ipt.info_2);
            }
          }
        }
      }

      // compute the intersection points between the two coplanar faces
      if (inter_pts.empty())
        intersection_coplanar_faces(f1, f2, tm1, tm2, vpm1, vpm2, inter_pts);

      std::size_t nb_pts=inter_pts.size();
      std::vector<Node_id> cpln_nodes; cpln_nodes.reserve(nb_pts);

      for(const Cpl_inter_pt& ipt : inter_pts)
      {
        Node_id node_id;
        bool is_new_node;
        std::tie(node_id, is_new_node) =
            get_or_create_node(ipt,current_node,coplanar_node_map,nm_features_map_1,nm_features_map_2,tm1,tm2);
        cpln_nodes.push_back(node_id);

        switch(ipt.type_1){
        case ON_VERTEX:
        {
          switch(ipt.type_2){
          case ON_VERTEX:
            handle_coplanar_case_VERTEX_VERTEX(ipt.info_1,ipt.info_2,tm1,tm2,nm_features_map_1,nm_features_map_2,node_id,is_new_node);
          break;
          case ON_EDGE:
            handle_coplanar_case_VERTEX_EDGE(ipt.info_1,ipt.info_2,tm1,tm2,nm_features_map_1,nm_features_map_2,node_id,is_new_node);
          break;
          case ON_FACE:
            handle_coplanar_case_VERTEX_FACE(ipt.info_1,ipt.info_2,tm1,tm2,nm_features_map_1,node_id,is_new_node);
          break;
          default: CGAL_error_msg("Should not get there!");
          }
        }
        break;
        case ON_EDGE:
        {
          switch(ipt.type_2){
            case ON_VERTEX:
              handle_coplanar_case_VERTEX_EDGE(ipt.info_2,ipt.info_1,tm2,tm1,nm_features_map_2,nm_features_map_1,node_id,is_new_node);
            break;
            case ON_EDGE:
            {
              std::vector<edge_descriptor> tmp_edges_1(1, edge(ipt.info_1,tm1)),
                                           tmp_edges_2(1, edge(ipt.info_2,tm2));
              std::size_t eid1 = nm_features_map_1.non_manifold_edges.empty()
                 ? NM_NID
                 : get(nm_features_map_1.e_nm_id, edge(ipt.info_1, tm1));
              std::size_t eid2 = nm_features_map_2.non_manifold_edges.empty()
                 ? NM_NID
                 : get(nm_features_map_2.e_nm_id, edge(ipt.info_2, tm2));
              const std::vector<edge_descriptor>& edges_1 = eid1==NM_NID
                                                          ? tmp_edges_1
                                                          : nm_features_map_1.non_manifold_edges[eid1];
              const std::vector<edge_descriptor>& edges_2 = eid2==NM_NID
                                                          ? tmp_edges_2
                                                          : nm_features_map_2.non_manifold_edges[eid2];
              if(is_new_node)
                visitor.new_node_added(node_id,ON_EDGE,halfedge(edges_1.front(), tm1),halfedge(edges_2.front(),tm2),tm1,tm2,false,false);
              for(edge_descriptor e1 : edges_1)
                for(edge_descriptor e2 : edges_2)
                {
                  typename Edge_to_faces::iterator it_ets=stm_edge_to_ltm_faces.find(e1);
                  Face_set* fset = (it_ets!=stm_edge_to_ltm_faces.end())?&(it_ets->second):nullptr;
                  cip_handle_case_edge(node_id,fset,halfedge(e1,tm1),halfedge(e2,tm2),tm1,tm2);
                }
            }
            break;
            default: CGAL_error_msg("Should not get there!");
          }
        }
        break;
        case ON_FACE:
        {
          CGAL_assertion(ipt.type_2==ON_VERTEX);
          handle_coplanar_case_VERTEX_FACE(ipt.info_2,ipt.info_1,tm2,tm1,nm_features_map_2,node_id,is_new_node);
        }
        break;
        default: CGAL_error_msg("Should not get there!");
        }
      }

      switch (nb_pts){
        case 0: break;
        case 1:
        {
            f_to_node[face_pair].insert(cpln_nodes[0]); // TODO: really?
        }
        break;
        default:
        {
          Node_id_set& node_id_set = f_to_node[face_pair];
          std::size_t stop=nb_pts + (nb_pts<3?-1:0);
          for (std::size_t k=0;k<stop;++k)
            node_id_set.insert( cpln_nodes[k], cpln_nodes[(k+1)%nb_pts] );
        }
      }
    }
    visitor.end_handling_intersection_of_coplanar_faces();
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  template <typename VPM1, typename VPM2>
  void add_new_node(halfedge_descriptor h_1,
                    face_descriptor f_2,
                    const TriangleMesh& tm1,
                    const TriangleMesh& tm2,
                    const VPM1& vpm1,
                    const VPM2& vpm2,
                    std::tuple<Intersection_type,
                                 halfedge_descriptor,
                                 bool,bool> inter_res)
  {
    if ( std::get<3>(inter_res) ) // is edge target in triangle plane
      nodes.add_new_node(get(vpm1, target(h_1,tm1)));
    else{
      if (std::get<2>(inter_res)) // is edge source in triangle plane
        nodes.add_new_node(get(vpm1, source(h_1,tm1)));
      else
        nodes.add_new_node(h_1,f_2,tm1,tm2,vpm1,vpm2);
    }
  }

  template <typename VPM1, typename VPM2>
  void compute_intersection_points(Edge_to_faces& tm1_edge_to_tm2_faces,
                                   const TriangleMesh& tm1,
                                   const TriangleMesh& tm2,
                                   const VPM1& vpm1,
                                   const VPM2& vpm2,
                                   const Non_manifold_feature_map<TriangleMesh>& nm_features_map_1,
                                   const Non_manifold_feature_map<TriangleMesh>& nm_features_map_2,
                                   Node_id& current_node)
  {
    typedef std::tuple<Intersection_type, halfedge_descriptor, bool,bool>  Inter_type;

    visitor.start_handling_edge_face_intersections(tm1_edge_to_tm2_faces.size());

    for(typename Edge_to_faces::iterator it=tm1_edge_to_tm2_faces.begin();
                                         it!=tm1_edge_to_tm2_faces.end();++it)
    {
      visitor.edge_face_intersections_step();
      edge_descriptor e_1=it->first;

      halfedge_descriptor h_1=halfedge(e_1,tm1);
      Face_set& fset=it->second;
      while (!fset.empty()){
        face_descriptor f_2=*fset.begin();

        Inter_type res=intersection_type(h_1,f_2,tm1,tm2,vpm1,vpm2);
        Intersection_type type=std::get<0>(res);

    //handle degenerate case: one extremity of edge belong to f_2
        std::vector<halfedge_descriptor> all_edges;
        if ( std::get<3>(res) ) // is edge target in triangle plane
        {
          if (!nm_features_map_1.non_manifold_edges.empty())
          {
            std::size_t vid1 = get(nm_features_map_1.v_nm_id, target(h_1, tm1));
            if (vid1 != NM_NID)
            {
              for (vertex_descriptor vd : nm_features_map_1.non_manifold_vertices[vid1])
              {
                std::copy(halfedges_around_target(vd,tm1).first,
                          halfedges_around_target(vd,tm1).second,
                          std::back_inserter(all_edges));
              }
              if (all_edges.front()!=h_1)
              {
                // restore expected property
                typename std::vector<halfedge_descriptor>::iterator pos =
                  std::find(all_edges.begin(), all_edges.end(), h_1);
                CGAL_assertion(pos!=all_edges.end());
                std::swap(*pos, all_edges.front());
              }
            }
            else
              std::copy(halfedges_around_target(h_1,tm1).first,
                        halfedges_around_target(h_1,tm1).second,
                        std::back_inserter(all_edges));
          }
          else
            std::copy(halfedges_around_target(h_1,tm1).first,
                      halfedges_around_target(h_1,tm1).second,
                      std::back_inserter(all_edges));
        }
        else{
          if ( std::get<2>(res) ) // is edge source in triangle plane
          {
            if (!nm_features_map_1.non_manifold_edges.empty())
            {
              std::size_t vid1 = get(nm_features_map_1.v_nm_id, source(h_1, tm1));
              if (vid1 != NM_NID)
              {
                for (vertex_descriptor vd : nm_features_map_1.non_manifold_vertices[vid1])
                {
                  std::copy(halfedges_around_source(vd,tm1).first,
                            halfedges_around_source(vd,tm1).second,
                            std::back_inserter(all_edges));
                }
                if (all_edges.front()!=h_1)
                {
                  // restore expected property
                  typename std::vector<halfedge_descriptor>::iterator pos =
                    std::find(all_edges.begin(), all_edges.end(), h_1);
                  CGAL_assertion(pos!=all_edges.end());
                  std::swap(*pos, all_edges.front());
                }
              }
              else
                std::copy(halfedges_around_source(h_1,tm1).first,
                          halfedges_around_source(h_1,tm1).second,
                          std::back_inserter(all_edges));
            }
            else
              std::copy(halfedges_around_source(h_1,tm1).first,
                        halfedges_around_source(h_1,tm1).second,
                        std::back_inserter(all_edges));
          }
          else
          {
            all_edges.push_back(h_1);
            edge_descriptor e_1 = edge(h_1, tm1);
            if (!nm_features_map_1.non_manifold_edges.empty())
            {
              std::size_t eid1 = get(nm_features_map_1.e_nm_id, e_1);
              if (eid1 != NM_NID)
              {
                CGAL_assertion( nm_features_map_1.non_manifold_edges[eid1][0]==e_1 );
                for (std::size_t k=1;
                                 k<nm_features_map_1.non_manifold_edges[eid1].size();
                                 ++k)
                {
                  edge_descriptor e_1b = nm_features_map_1.non_manifold_edges[eid1][k];
                  // note that the orientation of the halfedge pushed back is
                  // not relevant for how it is used in the following
                  all_edges.push_back(halfedge(e_1b, tm1));
                }
              }
            }
          }
        }
        CGAL_precondition(all_edges[0]==h_1 || all_edges[0]==opposite(h_1,tm1));

        // #ifdef USE_DETECTION_MULTIPLE_DEFINED_EDGES
        // check_coplanar_edges(std::next(all_edges.begin()),
        //                      all_edges.end(),std::get<1>(res),type);
        // #endif

        typename std::vector<halfedge_descriptor>::iterator it_edge=all_edges.begin();
        switch(type){
          case COPLANAR_TRIANGLES:
            #ifndef DO_NOT_HANDLE_COPLANAR_FACES
            CGAL_error_msg("COPLANAR_TRIANGLES : this point should never be reached!");
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
            CGAL_assertion(f_2==face(std::get<1>(res),tm2));

            Node_id node_id=++current_node;
            add_new_node(h_1,f_2,tm1,tm2,vpm1,vpm2,res);
            visitor.new_node_added(node_id,ON_FACE,h_1,halfedge(f_2,tm2),tm1,tm2,std::get<3>(res),std::get<2>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              add_intersection_point_to_face_and_all_edge_incident_faces(f_2,*it_edge,tm2,tm1,node_id);
              //erase face from the list to test intersection with it_edge
              if ( it_edge==all_edges.begin() )
              {
                fset.erase(fset.begin());
              }
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
            halfedge_descriptor h_2=std::get<1>(res);

            std::size_t eid2 = nm_features_map_2.non_manifold_edges.empty()
                             ? NM_NID
                             : get(nm_features_map_2.e_nm_id, edge(h_2, tm2));

            if (eid2!=NM_NID)
              h_2 = halfedge(nm_features_map_2.non_manifold_edges[eid2].front(), tm2);

            visitor.new_node_added(node_id,ON_EDGE,h_1,h_2,tm1,tm2,std::get<3>(res),std::get<2>(res));
            for (;it_edge!=all_edges.end();++it_edge){
              if ( it_edge!=all_edges.begin() ){
                typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(*it_edge,tm1));
                Face_set* fset_bis = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):nullptr;
                if( eid2 == NM_NID )
                  cip_handle_case_edge(node_id,fset_bis,*it_edge,h_2,tm1,tm2);
                else
                {
                  for (edge_descriptor e2 : nm_features_map_2.non_manifold_edges[eid2])
                    cip_handle_case_edge(node_id,fset_bis,*it_edge,halfedge(e2, tm2),tm1,tm2);
                }
              }
              else
              {
                if( eid2 == NM_NID )
                  cip_handle_case_edge(node_id,&fset,*it_edge,h_2,tm1,tm2);
                else
                  for (edge_descriptor e2 : nm_features_map_2.non_manifold_edges[eid2])
                    cip_handle_case_edge(node_id,&fset,*it_edge,halfedge(e2, tm2),tm1,tm2);
              }
            }
          } // end case ON_EDGE
          break;

          case ON_VERTEX:
          {
            Node_id node_id=++current_node;
            halfedge_descriptor h_2=std::get<1>(res);
            nodes.add_new_node(get(vpm2, target(h_2,tm2))); //we use the original vertex to create the node
            //before it was ON_FACE but do not remember why, probably a bug...
            visitor.new_node_added(node_id,ON_VERTEX,h_1,h_2,tm1,tm2,std::get<3>(res),std::get<2>(res));

            std::size_t vid2 = nm_features_map_2.non_manifold_vertices.empty()
                             ? NM_NID
                             : get(nm_features_map_2.v_nm_id, target(h_2, tm2));

            for (;it_edge!=all_edges.end();++it_edge){
              if ( it_edge!=all_edges.begin() ){
                typename Edge_to_faces::iterator it_ets=tm1_edge_to_tm2_faces.find(edge(*it_edge,tm1));
                Face_set* fset_bis = (it_ets!=tm1_edge_to_tm2_faces.end())?&(it_ets->second):nullptr;
                if( vid2 == NM_NID )
                  cip_handle_case_vertex(node_id,fset_bis,*it_edge,h_2,tm1,tm2);
                else
                  for (vertex_descriptor vd2 : nm_features_map_2.non_manifold_vertices[vid2])
                    cip_handle_case_vertex(node_id,fset_bis,*it_edge,halfedge(vd2, tm2),tm1,tm2);
              }
              else
                if( vid2 == NM_NID )
                  cip_handle_case_vertex(node_id,&fset,*it_edge,h_2,tm1,tm2);
                else
                  for (vertex_descriptor vd2 : nm_features_map_2.non_manifold_vertices[vid2])
                    cip_handle_case_vertex(node_id,&fset,*it_edge,halfedge(vd2, tm2),tm1,tm2);
            }
          } // end case ON_VERTEX
          break;
        } // end switch on the type of the intersection
      } // end loop on all faces that intersect the edge
    } // end loop on all entries (edges) in 'edge_to_face'
    CGAL_assertion(nodes.size()==unsigned(current_node+1));
    visitor.end_handling_edge_face_intersections();
  }

  struct Graph_node{
    boost::container::flat_set<Node_id> neighbors;
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

  template <class VPM>
  void detect_intersections_in_the_graph(const TriangleMesh& tm,
                                         const VPM& vpm,
                                         Node_id& current_node)
  {
    std::unordered_map<face_descriptor,
                       std::vector<face_descriptor> > face_intersections;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();
                                               it!=f_to_node.end();
                                               ++it)
    {
      face_descriptor f1 = it->first.first, f2 = it->first.second;
      CGAL_assertion( f1 < f2 );
      std::vector<face_descriptor>& inter_faces=face_intersections[f1];
      if (inter_faces.empty() || inter_faces.back()!=f2)
        face_intersections[f1].push_back(f2);
    }

    // TODO AUTOREF_TAG find a better way to avoid too many queries.
    std::map< Node_id_set*, std::vector<std::pair<std::array<std::size_t, 2>, Node_id>> > map_to_process;

    typedef std::pair<const face_descriptor, std::vector<face_descriptor> > Pair_type;
    for(Pair_type& p : face_intersections)
    {
      face_descriptor f1 = p.first;
      std::vector<face_descriptor>& inter_faces=p.second;
      std::sort(inter_faces.begin(), inter_faces.end());

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
          typename Faces_to_nodes_map::iterator it_seg23 =
            f_to_node.find(std::make_pair(f2, f3));

          if (it_seg23!=f_to_node.end())
          {
            std::vector<std::array<std::size_t,2> > f2f3_segments;
            it_seg23->second.get_segments(f2f3_segments);

            // look for edges between f1 and f2
            typename Faces_to_nodes_map::iterator it_seg12 =
              f_to_node.find(std::make_pair(f1, f2));
            CGAL_assertion( it_seg12 != f_to_node.end() );
            std::vector<std::array<std::size_t,2> > f1f2_segments;
            it_seg12->second.get_segments(f1f2_segments);

            // look for edges between f1 and f3
            typename Faces_to_nodes_map::iterator it_seg13 =
              f_to_node.find(std::make_pair(f1, f3));
            CGAL_assertion( it_seg13 != f_to_node.end() );
            std::vector<std::array<std::size_t,2> > f1f3_segments;
            it_seg13->second.get_segments(f1f3_segments);

            /// TODO AUTOREF_TAG shall we ignore tangency points?
            /// with the current code, Node_id_set::size()==1 is ignored as we only drop semgents
            /// Actually it might be that it is not a tangency point if the third segment was considered!
            /// so not handling it is a bug

            for (const std::array<std::size_t,2>& ns12 : f1f2_segments)
            {
              for (const std::array<std::size_t,2>& ns13 : f1f3_segments)
              {
                // handle cases of segments sharing an endpoint
                if (ns12[0]==ns13[0] || ns12[0]==ns13[1] ||
                    ns12[1]==ns13[0] || ns12[1]==ns13[1] )
                {
                  Node_id common_nid, nid1, nid2;

                  if (ns12[0]==ns13[0])
                  {
                    common_nid=ns12[0];
                    nid1=ns12[1];
                    nid2=ns13[1];
                  }
                  else
                  {
                    if (ns12[0]==ns13[1])
                    {
                      common_nid=ns12[0];
                      nid1=ns12[1];
                      nid2=ns13[0];
                    }
                    else
                    {
                      if (ns12[1]==ns13[0])
                      {
                        common_nid=ns12[1];
                        nid1=ns12[0];
                        nid2=ns13[1];
                      }
                      else
                      {
                        common_nid=ns12[1];
                        nid1=ns12[0];
                        nid2=ns13[0];
                      }
                    }
                  }

                  if (nid1==nid2)
                  {
                    // shared intersection edge ( TODO what happen for boundary intersection edges)
                    throw Triple_intersection_exception();
                  }

                  typename Node_vector::Exact_kernel::Point_3
                    common_pt = nodes.exact_node(common_nid),
                    pt1 = nodes.exact_node(nid1), pt2 = nodes.exact_node(nid2);
                  if (collinear(common_pt, pt1, pt2))
                  {
                    if ( !collinear_are_ordered_along_line(pt1, common_pt, pt2)  )
                    {
                      throw Triple_intersection_exception();
                    }
                  }
                  continue;
                }

                // TODO AUTOREF_TAG there might be a better test rather than relying on constructions
                typedef typename Node_vector::Exact_kernel::Segment_3 Segment_3;
                Segment_3 s12(nodes.exact_node(ns12[0]), nodes.exact_node(ns12[1])),
                          s13(nodes.exact_node(ns13[0]), nodes.exact_node(ns13[1])),
                          s23;

                if( do_intersect(s12, s13) )
                {
                  /// TODO AUTOREF_TAG it might be the end point of a segment!
                  // we need to find which segment is the third one intersecting
                  // with the point found.

                  // first check if there is only one such edge (no test is needed then)
                  CGAL_assertion(!f2f3_segments.empty() || !"AUTOREF_TAG HANDLE ME");

                  std::array<std::size_t,2> ns23 = f2f3_segments.front();

                  if (f2f3_segments.size()!=1)
                  {
                    std::size_t k=0;
                    for (;k<f2f3_segments.size(); ++k)
                    {
                      ns23=f2f3_segments[k];
                      s23=Segment_3(nodes.exact_node(ns23[0]), nodes.exact_node(ns23[1]));
                      if (do_intersect(s12, s23))
                      {
                        CGAL_assertion(do_intersect(s13,s23));
                        break;
                      }
                    }
                    CGAL_assertion(k!=f2f3_segments.size());
                  }
                  else
                    s23=Segment_3(nodes.exact_node(ns23[0]), nodes.exact_node(ns23[1]));

                  CGAL_assertion(do_intersect(s12, s23));

                  // use s23
                  ///TODO AUTOREF_TAG the collinear test could be factorized in the do-intersect
                  /// TODO AUTOREF_TAG if we don't factorise maybe checking for collinearity with
                  ///      cross product would be cheaper?
                  if ( (collinear(s12[0], s12[1], s13[0]) && collinear(s12[0], s12[1], s13[1]) ) ||
                       (collinear(s12[0], s12[1], s23[0]) && collinear(s12[0], s12[1], s23[1]) ) ||
                       (collinear(s13[0], s13[1], s23[0]) && collinear(s13[0], s13[1], s23[1]) ) )
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

                  map_to_process[&(it_seg12->second)].push_back(std::make_pair(ns12, node_id));
                  map_to_process[&(it_seg13->second)].push_back(std::make_pair(ns13, node_id));
                  map_to_process[&(it_seg23->second)].push_back(std::make_pair(ns23, node_id));
                }
              }
            }
          }
        }
      }
    }

#ifdef CGAL_DEBUG_AUTOREFINEMENT
    std::cout << "\nAt the end of new node creation, current_node is " << current_node << "\n\n";
#endif
    typedef std::pair<Node_id_set* const, std::vector<std::pair<std::array<std::size_t,2>, Node_id> > > Node_id_set_and_nodes;
    for(Node_id_set_and_nodes& n_id_and_new_nodes : map_to_process)
    {
      Node_id_set& nids = *(n_id_and_new_nodes.first);
      CGAL_assertion(nids.size()!=1); // TODO AUTOREF_TAG handle case size = 1

      // let's abuse coplanar_segments to store refinements
      if (nids.size()==2)
      {
        nids.coplanar_segments.push_back({nids[0], nids[1]});
        nids.size_=0;
      }

      std::map<std::array<std::size_t, 2>, std::vector<Node_id> > local_map;
      for (const std::pair<std::array<std::size_t,2>, Node_id>& p : n_id_and_new_nodes.second)
      {
        CGAL_assertion(p.first[0]<p.first[1]);
        local_map[p.first].push_back(p.second);
      }

      for (std::pair<const std::array<std::size_t,2>, std::vector<Node_id> >& p : local_map)
      {
        //get the original entry and remove it
        auto it = std::find(nids.coplanar_segments.begin(), nids.coplanar_segments.end(), p.first);
        CGAL_assertion(it!=nids.coplanar_segments.end());
        nids.coplanar_segments.erase(it);

        std::vector<Node_id>& new_nodes = p.second;
        Node_id n1 = p.first[0];
        Node_id n2 = p.first[1];

        // sort node ids along the edge
        std::sort(new_nodes.begin(),
                  new_nodes.end(),
                  Less_for_nodes_along_an_edge(nodes, n1));

        // insert new segments
        Node_id prev = n1;
        new_nodes.push_back(n2); // add last node id to avoid having another case after the loop
  #ifdef CGAL_DEBUG_AUTOREFINEMENT
        std::cout << n1 << " -> " << n2 << "\n";
  #endif
        for(Node_id id : new_nodes)
        {
          nids.coplanar_segments.push_back( {prev, id} );
  #ifdef CGAL_DEBUG_AUTOREFINEMENT
          std::cerr <<"  adding " << prev << " " << id << " into "
                    << n1 << " and " << n2 <<  "\n";
  #endif
          prev=id;
        }
      }

      // update main entry
      nids.first=nids.coplanar_segments.back()[0];
      nids.second=nids.coplanar_segments.back()[1];
      nids.size_=2;
      nids.coplanar_segments.pop_back();
    }
  }

  template <class Output_iterator>
  void construct_polylines(Output_iterator out){
    typedef typename boost::property_traits<VertexPointMap1>::value_type Point_3;
    std::size_t nb_nodes=nodes.size();
    std::vector<Graph_node> graph(nb_nodes);
    //counts the number of time each node has been seen
    bool isolated_point_seen=false;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();it!=f_to_node.end();++it){
      const Node_id_set& segment=it->second;
      CGAL_assertion(segment.size()<=2);
      if (segment.size()==2){
        Node_id i=segment.first;
        Node_id j=segment.second;
        graph[i].insert(j);
        graph[j].insert(i);
      }
      else{
        if (segment.size()==1)
          isolated_point_seen=true; // NOT TRUE CAN BE END POINT OF POLYLINE FALLING ONTO AN INPUT EDGE
        else
        {
          CGAL_assertion(!it->second.coplanar_segments.empty());
        }
      }

      for (const std::array<std::size_t, 2>& ij : it->second.coplanar_segments)
      {
        graph[ij[0]].insert(ij[1]);
        graph[ij[1]].insert(ij[0]);
      }
    }

    CGAL_assertion(extra_terminal_nodes.empty() || doing_autorefinement);
    // these nodes are created by pinchements along an edge of the surface.
    // the node ids being the same for the two edges, the degree of the node
    // in the graph is two while it should be 3
    for(Node_id id : extra_terminal_nodes)
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
    std::set< std::array<Node_id,2> > already_seen;
    std::vector<typename Faces_to_nodes_map::iterator> to_erase;
    for (typename Faces_to_nodes_map::iterator it=f_to_node.begin();
          it!=f_to_node.end(); ++it)
    {
      if (it->second.size()==2)
      {
        CGAL_assertion(it->second.first < it->second.second);
        if (!already_seen.insert( {it->second.first, it->second.second} ).second)
          it->second.size_=0;
      }

      it->second.coplanar_segments.erase(
        std::remove_if(it->second.coplanar_segments.begin(),
                       it->second.coplanar_segments.end(),
                       [&already_seen](const std::array<std::size_t, 2>& a)
                       {
                         CGAL_assertion(a[0]<a[1]);
                         return !already_seen.insert(a).second;
                       }),
        it->second.coplanar_segments.end());

      if (it->second.size()==0 && it->second.coplanar_segments.empty())
        to_erase.push_back(it);
    }

    for(typename Faces_to_nodes_map::iterator it : to_erase)
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
      if (it->second.size()!=1) continue;

      halfedge_descriptor h1 = halfedge(it->first.first, tm);
      halfedge_descriptor h2 = halfedge(it->first.second, tm);

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
                                  const VertexPointMap1& vpm1,
                                  const VertexPointMap2& vpm2,
                                  const Node_visitor& v=Node_visitor(),
                                  const TriangleMesh* const_mesh_ptr=nullptr)
  : nodes(tm1, tm2, vpm1, vpm2)
  , visitor(v)
  , const_mesh_ptr(const_mesh_ptr)
  {
    CGAL_precondition(is_triangle_mesh(tm1));
    CGAL_precondition(is_triangle_mesh(tm2));
    CGAL_assertion_code( doing_autorefinement=false; )
  }

  // for autorefinement
  Intersection_of_triangle_meshes(const TriangleMesh& tm,
                                  const VertexPointMap1& vpm,
                                  const Node_visitor& v=Node_visitor())
  : nodes(tm, tm, vpm, vpm)
  , visitor(v)
  {
    CGAL_precondition(is_triangle_mesh(tm));
    CGAL_assertion_code( doing_autorefinement=true; )
  }

// setting maps of non manifold features
  void set_non_manifold_feature_map_1(internal_np::Param_not_found){}
  void set_non_manifold_feature_map_2(internal_np::Param_not_found){}
  void set_non_manifold_feature_map_1(const Non_manifold_feature_map<TriangleMesh>& m)
  {
    non_manifold_feature_map_1=m;
    visitor.set_non_manifold_feature_map(nodes.tm1, non_manifold_feature_map_1);
  }
  void set_non_manifold_feature_map_2(const Non_manifold_feature_map<TriangleMesh>& m)
  {
    non_manifold_feature_map_2=m;
    visitor.set_non_manifold_feature_map(nodes.tm2, non_manifold_feature_map_2);
  }

  template <class OutputIterator>
  OutputIterator operator()(OutputIterator output,
                            bool throw_on_self_intersection,
                            bool build_polylines)
  {
    CGAL_assertion(!doing_autorefinement);

    const TriangleMesh& tm1=nodes.tm1;
    const TriangleMesh& tm2=nodes.tm2;
    const VertexPointMap1& vpm1=nodes.vpm1;
    const VertexPointMap2& vpm2=nodes.vpm2;

    // used only if throw_on_self_intersection == true
    std::set<face_descriptor> tm1_faces;
    std::set<face_descriptor> tm2_faces;

    visitor.start_filtering_intersections();
    filter_intersections(tm1, tm2, vpm1, vpm2, non_manifold_feature_map_2, throw_on_self_intersection, tm1_faces, tm2_faces, false);
    filter_intersections(tm2, tm1, vpm2, vpm1, non_manifold_feature_map_1, throw_on_self_intersection, tm2_faces, tm1_faces, true);
    visitor.end_filtering_intersections();

    Node_id current_node((std::numeric_limits<Node_id>::max)());
    CGAL_assertion(current_node+1==0);
// TODO: handle non-manifold edges in coplanar
    #ifndef DO_NOT_HANDLE_COPLANAR_FACES
    //first handle coplanar triangles
    if (&tm1<&tm2)
      compute_intersection_of_coplanar_faces(current_node, tm1, tm2, vpm1, vpm2, non_manifold_feature_map_1, non_manifold_feature_map_2);
    else
      compute_intersection_of_coplanar_faces(current_node, tm2, tm1, vpm2, vpm1, non_manifold_feature_map_2, non_manifold_feature_map_1);

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

    compute_intersection_points(tm1_edge_to_tm2_faces, tm1, tm2, vpm1, vpm2, non_manifold_feature_map_1, non_manifold_feature_map_2, current_node);
    compute_intersection_points(tm2_edge_to_tm1_faces, tm2, tm1, vpm2, vpm1, non_manifold_feature_map_2, non_manifold_feature_map_1, current_node);

    visitor.check_no_duplicates(nodes);

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
    const VertexPointMap1& vpm=nodes.vpm1;

    filter_intersections(tm, vpm);

    Node_id current_node((std::numeric_limits<Node_id>::max)());
    CGAL_assertion(current_node+1==0);

    //first handle coplanar triangles
    compute_intersection_of_coplanar_faces(current_node, tm, tm, vpm, vpm, non_manifold_feature_map_1, non_manifold_feature_map_1);
    if (!coplanar_faces.empty())
      visitor.input_have_coplanar_faces();

    CGAL_assertion(ltm_edge_to_stm_faces.empty());

    //compute intersection points of segments and triangles.
    //build the nodes of the graph and connectivity infos
    compute_intersection_points(stm_edge_to_ltm_faces, tm, tm, vpm, vpm, non_manifold_feature_map_1, non_manifold_feature_map_1, current_node);

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

} } } // CGAL::Polygon_mesh_processing::Corefinement

#endif //CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_IMPL_H
