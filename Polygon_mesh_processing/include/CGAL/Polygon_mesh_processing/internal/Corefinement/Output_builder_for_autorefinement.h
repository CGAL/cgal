// Copyright (c) 2017 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_OUTPUT_BUILDER_FOR_AUTOREFINEMENT_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_OUTPUT_BUILDER_FOR_AUTOREFINEMENT_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/Union_find.h>
#include <CGAL/property_map.h>
#include <CGAL/Default.h>

#include <boost/dynamic_bitset.hpp>

namespace CGAL {
namespace Corefinement {

namespace PMP=Polygon_mesh_processing;
namespace params=PMP::parameters;

template <class TriangleMesh,
          class VertexPointMap,
          class FaceIdMap,
          class Ecm,
          class Kernel_=Default>
class Output_builder_for_autorefinement
{
//Default typedefs
  typedef typename Default::Get<
    Kernel_,
    typename Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type
    >::Kernel >::type                                           Kernel;

// graph_traits typedefs
  typedef TriangleMesh                                              TM;
  typedef boost::graph_traits<TM>                                   GT;
  typedef typename GT::edge_descriptor                 edge_descriptor;
  typedef typename GT::face_descriptor                 face_descriptor;
  typedef typename GT::halfedge_descriptor         halfedge_descriptor;
  typedef typename GT::vertex_descriptor             vertex_descriptor;
// Internal typedefs
  typedef std::size_t                                          Node_id;
  typedef std::pair<Node_id,Node_id>                      Node_id_pair;

  struct Shared_halfedges
  {
    halfedge_descriptor h1, h2;
    Shared_halfedges()
      : h1(boost::graph_traits<TriangleMesh>::null_halfedge())
      , h2(boost::graph_traits<TriangleMesh>::null_halfedge())
    {}

    void add(halfedge_descriptor h)
    {
      if (h1==boost::graph_traits<TriangleMesh>::null_halfedge())
        h1=h;
      else
      {
        if (h1!=h)
        {
          if (h2==boost::graph_traits<TriangleMesh>::null_halfedge())
            h2=h;
          else
          {
            CGAL_assertion(h2==h);
          }
        }
      }
    }
  };

  // to maintain the two halfedges on each polyline
  typedef std::map< Node_id_pair, Shared_halfedges >   An_edge_per_polyline_map;
  typedef boost::unordered_map<vertex_descriptor, Node_id> Node_id_map;
  // typedef boost::unordered_map<edge_descriptor,
  //                              edge_descriptor>               Edge_map;
  typedef boost::unordered_map<Node_id_pair, Shared_halfedges> All_intersection_edges_map;
//Data members
  TriangleMesh &tm;
  // property maps of input mesh
  const VertexPointMap &vpm;
  const FaceIdMap &fids;
  Ecm& ecm;
  // input meshes closed ?
  bool is_tm_closed;
  // orientation of input surface mesh
  bool is_tm_inside_out;
  // constant
  const Node_id NID;
  // boolean indicating if there is an ambiguous or non-manifold situation
  bool all_fixed;
  // for mapping an edge per polyline per triangle mesh
  An_edge_per_polyline_map an_edge_per_polyline;
  // To collect all intersection edges
  All_intersection_edges_map all_intersection_edges_map;

  Node_id get_node_id(vertex_descriptor v,
                      const Node_id_map& node_ids)
  {
    typename Node_id_map::const_iterator it = node_ids.find(v);
    if (it == node_ids.end())
      return NID;
    return it->second;
  }

  bool is_dangling_edge(Node_id src_id, Node_id tgt_id,
                        halfedge_descriptor hedge,
                        TriangleMesh& tm,
                        const boost::dynamic_bitset<>& is_node_of_degree_one) const
  {
    if ( is_node_of_degree_one.test(src_id) )
    {
      bool res=true;
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(hedge, tm))
        if (is_border(h, tm))
        {
          res = false;
          break;
        }
      if (res) return true;
    }
    if ( is_node_of_degree_one.test(tgt_id) )
    {
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(hedge, tm))
        if (is_border(h, tm))
          return false;
      return true;
    }
    return false;
  }

public:

  Output_builder_for_autorefinement(      TriangleMesh& tm,
                                    const VertexPointMap& vpm,
                                          FaceIdMap& fids,
                                          Ecm& ecm)
    : tm(tm)
    , vpm(vpm)
    , fids(fids)
    , ecm(ecm)
    , is_tm_closed( is_closed(tm))
    , is_tm_inside_out( is_tm_closed && !PMP::is_outward_oriented(tm) )
    , NID((std::numeric_limits<Node_id>::max)())
    , all_fixed(true)
  {}

  bool all_self_intersection_fixed() const  { return all_fixed; }

// functions called by the intersection visitor
  void start_new_polyline(Node_id i, Node_id j)
  {
    std::pair<typename An_edge_per_polyline_map::iterator,bool> res=
      an_edge_per_polyline.insert(std::make_pair(make_sorted_pair(i,j),
                                                 Shared_halfedges() ) );
    CGAL_USE(res);
    CGAL_assertion(res.second);
  }

  void add_node_to_polyline(Node_id)
  {}

  void set_edge_per_polyline(TriangleMesh& tm,
                             Node_id_pair indices,
                             halfedge_descriptor hedge)
  {
    if (indices.first>indices.second)
    {
      std::swap(indices.first,indices.second);
      hedge=opposite(hedge,tm);
    }
    typename An_edge_per_polyline_map::iterator it =
      an_edge_per_polyline.find(indices);

    if (it!=an_edge_per_polyline.end())
      it->second.add( hedge );

    //register an intersection halfedge
    all_intersection_edges_map[indices].add(hedge);
  }

  template <class Nodes_vector, class Mesh_to_map_node>
  void operator()(
    const Nodes_vector& nodes,
    bool input_have_coplanar_faces,
    const boost::dynamic_bitset<>& is_node_of_degree_one,
    const Mesh_to_map_node&)
  {
    // this will initialize face indices if the face index map is writable.
    helpers::init_face_indices(tm, fids);

    // first build an unordered_map mapping a vertex to its node id + a set
    // of all intersection edges
    Node_id_map vertex_to_node_id;
    typedef boost::unordered_set<edge_descriptor> Intersection_edge_map;
    Intersection_edge_map intersection_edges;

    typedef std::pair<const Node_id_pair, Shared_halfedges> Pair_type;
    BOOST_FOREACH(const Pair_type& p, all_intersection_edges_map)
    {
      CGAL_assertion(p.second.h1!=boost::graph_traits<TriangleMesh>::null_halfedge());
    // p.second.h2 might be the null halfedge in case two faces sharing an edge
    // intersect (and are obviously coplanar). It is not considered as an intersection
    // and will be discarded later
      if (p.second.h2==boost::graph_traits<TriangleMesh>::null_halfedge())
        continue;
      vertex_to_node_id[source(p.second.h1, tm)] = p.first.first;
      vertex_to_node_id[target(p.second.h1, tm)] = p.first.second;
      vertex_to_node_id[source(p.second.h2, tm)] = p.first.first;
      vertex_to_node_id[target(p.second.h2, tm)] = p.first.second;
      intersection_edges.insert(edge(p.second.h1, tm));
      intersection_edges.insert(edge(p.second.h2, tm));
    }

    // this will initialize face indices if the face index map is writable.
    helpers::init_face_indices(tm, fids);

    // bitset to identify coplanar faces
    boost::dynamic_bitset<> tm_coplanar_faces(num_faces(tm), 0);

    // In the following loop we filter intersection edges that are strictly inside a patch
    // of coplanar facets so that we keep only the edges on the border of the patch.
    // This is not optimal and in an ideal world being able to find the outside edges
    // directly would avoid to compute the intersection of edge/facets inside the patch
    // This loop is done only if the input have some coplanar faces.
    typename An_edge_per_polyline_map::iterator
     epp_it_end=an_edge_per_polyline.end();
    typename An_edge_per_polyline_map::iterator
      epp_it=input_have_coplanar_faces ? an_edge_per_polyline.begin()
                                       : epp_it_end;
    boost::unordered_set<edge_descriptor> inter_edges_to_remove;
    for (;epp_it!=epp_it_end;)
    {
      halfedge_descriptor h1  = epp_it->second.h1;
      halfedge_descriptor h1_opp = opposite(h1, tm);
      halfedge_descriptor h2 = epp_it->second.h2;

      bool to_remove = false;

      if (h2==boost::graph_traits<TriangleMesh>::null_halfedge())
      {
        // we are in the case of two faces sharing an edge and intersecting
        // (coplanar faces)
        CGAL_assertion(get_node_id(target(next(h1, tm), tm), vertex_to_node_id) ==
          get_node_id(target(next(opposite(h1, tm), tm), tm), vertex_to_node_id));
        to_remove=true;
      }
      else{
        halfedge_descriptor h2_opp = opposite(h2, tm);

        if (is_border_edge(h1,tm) || is_border_edge(h2,tm)){
          ++epp_it;
          continue;
        }

        //vertices from tm1
        vertex_descriptor p1 = target(next(h1_opp, tm), tm);
        vertex_descriptor p2 = target(next(h1, tm), tm);
        Node_id index_p1 = get_node_id(p1, vertex_to_node_id);
        Node_id index_p2 = get_node_id(p2, vertex_to_node_id);
        //vertices from tm2
        vertex_descriptor q1 = target(next(h2_opp, tm), tm);
        vertex_descriptor q2 = target(next(h2, tm), tm);
        Node_id index_q1 = get_node_id(q1, vertex_to_node_id);
        Node_id index_q2 = get_node_id(q2, vertex_to_node_id);

        // set boolean for the position of p1 wrt to q1 and q2
        bool p1_eq_q1=false, p1_eq_q2=false;
        if (!is_border(h1_opp, tm) && index_p1!=NID)
        {
          if (!is_border(h2_opp, tm))
            p1_eq_q1 = index_p1 == index_q1;
          if (!is_border(h2, tm))
            p1_eq_q2 = index_p1 == index_q2;
        }

        // set boolean for the position of p2 wrt to q1 and q2
        bool p2_eq_q1=false, p2_eq_q2=false;
        if (!is_border(h1, tm) && index_p2!=NID)
        {
          if (!is_border(h2_opp, tm))
            p2_eq_q1 = index_p2 == index_q1;
          if (!is_border(h2, tm))
            p2_eq_q2 = index_p2 == index_q2;
        }

        //mark coplanar facets if any
        if (p1_eq_q1){
          tm_coplanar_faces.set(get(fids, face(h1_opp, tm)));
          tm_coplanar_faces.set(get(fids, face(h2_opp, tm)));
        }
        if (p1_eq_q2){
          tm_coplanar_faces.set(get(fids, face(h1_opp, tm)));
          tm_coplanar_faces.set(get(fids, face(h2, tm)));
        }
        if (p2_eq_q1){
          tm_coplanar_faces.set(get(fids, face(h1, tm)));
          tm_coplanar_faces.set(get(fids, face(h2_opp, tm)));
        }
        if (p2_eq_q2){
          tm_coplanar_faces.set(get(fids, face(h1, tm)));
          tm_coplanar_faces.set(get(fids, face(h2, tm)));
        }
        if ( (p1_eq_q1 || p1_eq_q2) && (p2_eq_q1 || p2_eq_q2) )
          to_remove = true;
      }
      if (to_remove)
      {
        typename An_edge_per_polyline_map::iterator it_to_rm=epp_it;
        ++epp_it;
        an_edge_per_polyline.erase(it_to_rm);
        inter_edges_to_remove.insert(edge(h1,tm));
        if (h2!=boost::graph_traits<TriangleMesh>::null_halfedge())
          inter_edges_to_remove.insert(edge(h2,tm));
      }
      else
        ++epp_it;
    }

    BOOST_FOREACH(edge_descriptor ed, inter_edges_to_remove)
      intersection_edges.erase(ed);

    // (1) Assign a patch id to each face indicating in which connected
    // component limited by intersection edges of the surface they are.
    // ... for tm
    std::vector<std::size_t> patch_ids( num_faces(tm),NID );
    Boolean_property_map< boost::unordered_set<edge_descriptor> >
      is_intersection(intersection_edges);
    std::size_t nb_patches =
      PMP::connected_components(tm,
                                bind_property_maps(fids,make_property_map(patch_ids)),
                                params::edge_is_constrained_map(
                                    is_intersection)
                                .face_index_map(fids));

    // (2-a) Use the orientation around an edge to classify a patch
    boost::dynamic_bitset<> patches_to_keep(nb_patches);
    boost::dynamic_bitset<> patch_status_not_set(nb_patches);
    boost::dynamic_bitset<> coplanar_patches(nb_patches,false);
    patches_to_keep.set();
    patch_status_not_set.set();
    // use a union-find on patches to track the incidence between patches kept
    typedef Union_find<std::size_t> UF;
    UF uf;
    std::vector<typename UF::handle> patch_handles(nb_patches);
    for (std::size_t p=0; p<nb_patches; ++p)
    {
      patch_handles[p]=uf.make_set(p);
    }

    for (typename An_edge_per_polyline_map::iterator
            it=an_edge_per_polyline.begin(),
            it_end=an_edge_per_polyline.end(); it!=it_end;++it)
    {
      //orientation of faces around the edge (to be sure we can do it)
      const std::pair<Node_id,Node_id>& ids = it->first;
      //const std::pair<bool,int>& polyline_info=it->second.second;

      //get the two halfedges incident to the edge [ids.first,ids.second]
      halfedge_descriptor h1 = it->second.h1;
      halfedge_descriptor h2 = it->second.h2;

      CGAL_assertion(h1!=boost::graph_traits<TriangleMesh>::null_halfedge());
      CGAL_assertion(h2!=boost::graph_traits<TriangleMesh>::null_halfedge());

      CGAL_assertion(ids.first==vertex_to_node_id[source(h1,tm)]);
      CGAL_assertion(ids.second==vertex_to_node_id[target(h1,tm)]);
      CGAL_assertion(ids.first==vertex_to_node_id[source(h2,tm)]);
      CGAL_assertion(ids.second==vertex_to_node_id[target(h2,tm)]);

      // different handling depending on the number of incident
      // triangles to the edge. After sewing there are two, three or
      // four volumes if there are two, three or four incident
      // triangles respectively
      if ( is_border_edge(h1, tm) ){
        if ( is_border_edge(h2,tm) )
        {
          if ( is_border(h1,tm) == is_border(h2,tm) )
          {
            //Orientation issue, nothing done
            all_fixed = false;
          }
          else
          {
            if ( is_border(h1,tm) )
            {
              std::size_t pid1=patch_ids[ get(fids, face(opposite(h1,tm),tm)) ],
                          pid2=patch_ids[ get(fids, face(h2,tm)) ];
              uf.unify_sets(patch_handles[pid1], patch_handles[pid2]);
              patch_status_not_set.reset(pid1);
              patch_status_not_set.reset(pid2);
            }
            else
            {
              std::size_t pid1=patch_ids[ get(fids, face(h1,tm)) ],
                          pid2=patch_ids[ get(fids, face(opposite(h2,tm),tm)) ];
              uf.unify_sets(patch_handles[pid1], patch_handles[pid2]);
              patch_status_not_set.reset(pid1);
              patch_status_not_set.reset(pid2);
            }
          }
        }
        else
        {
          halfedge_descriptor h = is_border(h1, tm) ? opposite(h1, tm) : h1;

          //Sort the three triangle faces around their common edge
          //See the full description in the general case with 4 triangle faces.
          vertex_descriptor p=target(next(h,tm),tm);
          vertex_descriptor q1=target(next(opposite(h2,tm),tm),tm);
          vertex_descriptor q2=target(next(h2,tm),tm);

          Node_id index_p = get_node_id(p, vertex_to_node_id);
          Node_id index_q1 = get_node_id(q1, vertex_to_node_id);
          Node_id index_q2 = get_node_id(q2, vertex_to_node_id);

          std::size_t patch_id_p=patch_ids[ get(fids, face(h,tm)) ];
          std::size_t patch_id_q1=patch_ids[ get(fids, face(opposite(h2,tm),tm)) ];
          std::size_t patch_id_q2=patch_ids[ get(fids, face(h2,tm)) ];

          //indicates that patch status will be updated
          patch_status_not_set.reset(patch_id_p);
          patch_status_not_set.reset(patch_id_q1);
          patch_status_not_set.reset(patch_id_q2);

          bool p_is_between_q1q2 = sorted_around_edge(
            ids.first, ids.second,
            index_q1, index_q2, index_p,
            q1, q2, p,
            vpm, vpm,
            nodes);

          if (p_is_between_q1q2)
          {
            uf.unify_sets(patch_handles[patch_id_q1], patch_handles[patch_id_q2]);
            patches_to_keep.reset(patch_id_p); // even if badly oriented we can
          }                                     // simply discard the patch
          else
          {
            if (h==h1)
            {
              //Orientation issue, nothing done for the incident patches
              all_fixed = false;
            }
            else
            {
              uf.unify_sets(patch_handles[patch_id_p], patch_handles[patch_id_q2]);
              patches_to_keep.reset(patch_id_q1);
            }
          }
        }
      }
      else
        if ( is_border_edge(h2,tm) )
        {
          halfedge_descriptor h = is_border(h2, tm) ? opposite(h2, tm) : h2;

          //Sort the three triangle faces around their common edge
          //See the full description in the general case with 4 triangle faces.
          vertex_descriptor p1=target(next(opposite(h1,tm),tm),tm);
          vertex_descriptor p2=target(next(h1,tm),tm);
          vertex_descriptor q=target(next(h,tm),tm);

          Node_id index_p1 = get_node_id(p1, vertex_to_node_id);
          Node_id index_p2 = get_node_id(p2, vertex_to_node_id);
          Node_id index_q = get_node_id(q, vertex_to_node_id);

          std::size_t patch_id_p1=patch_ids[ get(fids, face(opposite(h1,tm),tm)) ];
          std::size_t patch_id_p2=patch_ids[ get(fids, face(h1,tm)) ];
          std::size_t patch_id_q=patch_ids[ get(fids, face(h,tm)) ];

          //indicates that patch status will be updated
          patch_status_not_set.reset(patch_id_p1);
          patch_status_not_set.reset(patch_id_p2);
          patch_status_not_set.reset(patch_id_q);

          bool q_is_between_p1p2 = sorted_around_edge(
            ids.first, ids.second,
            index_p1, index_p2, index_q,
            p1, p2, q,
            vpm, vpm,
            nodes);

          if (q_is_between_p1p2)
          {
            uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_p2]);
            patches_to_keep.reset(patch_id_q); // even if badly oriented we can
          }                                     // simply discard the patch
          else
          {
            if (h==h2)
            {
              //Orientation issue, nothing done for the incident patches
              all_fixed = false;
            }
            else
            {
              uf.unify_sets(patch_handles[patch_id_q], patch_handles[patch_id_p2]);
              patches_to_keep.reset(patch_id_p1);
            }
          }
        }
        else
        {
          //Sort the four triangle faces around their common edge
          //  we assume that the exterior of the volume is indicated by
          //  counterclockwise oriented faces
          //  (corrected by is_tmi_inside_tmi).
          vertex_descriptor p1=target(next(opposite(h1,tm),tm),tm);
          vertex_descriptor p2=target(next(h1,tm),tm);
          //    when looking from the side of indices.second,
          //    the interior of the first triangle mesh is described
          //    by turning counterclockwise from p1 to p2
          vertex_descriptor q1=target(next(opposite(h2,tm),tm),tm);
          vertex_descriptor q2=target(next(h2,tm),tm);
          //    when looking from the side of indices.second,
          //    the interior of the second volume is described
          //    by turning from q1 to q2

          //check if the third point of each triangular face is an original point (stay NID)
          //or a intersection point (in that case we need the index of the corresponding node to
          //have the exact value of the point)
          Node_id index_p1 = get_node_id(p1, vertex_to_node_id);
          Node_id index_p2 = get_node_id(p2, vertex_to_node_id);
          Node_id index_q1 = get_node_id(q1, vertex_to_node_id);
          Node_id index_q2 = get_node_id(q2, vertex_to_node_id);

          std::size_t patch_id_p1=patch_ids[ get(fids, face(opposite(h1,tm),tm)) ];
          std::size_t patch_id_p2=patch_ids[ get(fids, face(h1,tm)) ];
          std::size_t patch_id_q1=patch_ids[ get(fids, face(opposite(h2,tm),tm)) ];
          std::size_t patch_id_q2=patch_ids[ get(fids, face(h2,tm)) ];

          //indicates that patch status will be updated
          patch_status_not_set.reset(patch_id_p1);
          patch_status_not_set.reset(patch_id_p2);
          patch_status_not_set.reset(patch_id_q1);
          patch_status_not_set.reset(patch_id_q2);

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
          #warning: Factorize the orientation predicates.
#endif //CGAL_COREFINEMENT_POLYHEDRA_DEBUG
          // handle case of coplanar facets
          // We choose that a coplanar patch is classified like the other incident patch since they bound the same volume.
          if ( are_triangles_coplanar_same_side(
                ids.first, ids.second,
                index_p1, index_q1,
                p1, q1,
                vpm, vpm,
                nodes) ) //p1==q1
          {
            coplanar_patches.set(patch_id_p1);
            coplanar_patches.set(patch_id_q1);

            CGAL_assertion(
              !are_triangles_coplanar_same_side(
                ids.first, ids.second,
                index_p2, index_q2,
                p2, q2,
                vpm, vpm,
                nodes) );

            bool q2_is_between_p1p2 = sorted_around_edge(
              ids.first, ids.second,
              index_p1, index_p2, index_q2,
              p1, p2, q2,
              vpm, vpm,
              nodes);
            if ( q2_is_between_p1p2 ){
             //case 1
             patches_to_keep.reset(patch_id_q2);
             if (patch_id_p1<patch_id_q1)
             {
               uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_p2]);
               patches_to_keep.reset(patch_id_q1);
             }
             else
             {
               uf.unify_sets(patch_handles[patch_id_q1], patch_handles[patch_id_p2]);
               patches_to_keep.reset(patch_id_p1);
             }
            }
            else{
              //case 2
              patches_to_keep.reset(patch_id_p2);
              if (patch_id_p1<patch_id_q1)
              {
                uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_q2]);
                patches_to_keep.reset(patch_id_q1);
              }
              else
              {
                uf.unify_sets(patch_handles[patch_id_q1], patch_handles[patch_id_q2]);
                patches_to_keep.reset(patch_id_p1);
              }
            }
            continue;
          }
          else{
            if ( are_triangles_coplanar_same_side(
                   ids.first, ids.second,
                   index_p1, index_q2,
                   p1, q2,
                   vpm, vpm,
                   nodes) ) //p1==q2
            {
              CGAL_assertion( index_p1!=index_p2 || index_p1==Node_id(-1) );
              coplanar_patches.set(patch_id_p1);
              coplanar_patches.set(patch_id_q2);
              bool q1_is_between_p1p2 = sorted_around_edge(
                ids.first, ids.second,
                index_p1, index_p2, index_q1,
                p1, p2, q1,
                vpm, vpm,
                nodes);
              if ( q1_is_between_p1p2 ){
                // case 3
                patches_to_keep.reset(patch_id_p1);
                patches_to_keep.reset(patch_id_p2);
                patches_to_keep.reset(patch_id_q1);
                patches_to_keep.reset(patch_id_q2);
              }
              else{
                // case 4
                uf.unify_sets(patch_handles[patch_id_q1], patch_handles[patch_id_p2]);
                patches_to_keep.reset(patch_id_p1);
                patches_to_keep.reset(patch_id_q2);
              }
              continue;
            }
            else
            {
              if ( are_triangles_coplanar_same_side(
                     ids.first, ids.second,
                     index_p2, index_q1,
                     p2, q1,
                     vpm, vpm,
                     nodes) ) //p2==q1
              {
                coplanar_patches.set(patch_id_p2);
                coplanar_patches.set(patch_id_q1);
                bool q2_is_between_p1p2 = sorted_around_edge(
                  ids.first, ids.second,
                  index_p1, index_p2, index_q2,
                  p1, p2, q2,
                  vpm, vpm,
                  nodes);
                if ( q2_is_between_p1p2 )
                {  //case 5
                  patches_to_keep.reset(patch_id_p1);
                  patches_to_keep.reset(patch_id_p2);
                  patches_to_keep.reset(patch_id_q1);
                  patches_to_keep.reset(patch_id_q2);
                }else{
                  //case 6
                  uf.unify_sets(patch_handles[patch_id_q2], patch_handles[patch_id_p1]);
                  patches_to_keep.reset(patch_id_q1);
                  patches_to_keep.reset(patch_id_p2);
                }
                continue;
              }
              else{
                if ( are_triangles_coplanar_same_side(
                       ids.first, ids.second,
                       index_p2, index_q2,
                       p2, q2,
                       vpm, vpm,
                       nodes) ) //p2==q2
                {
                  coplanar_patches.set(patch_id_p2);
                  coplanar_patches.set(patch_id_q2);
                  bool q1_is_between_p1p2 = sorted_around_edge(
                    ids.first, ids.second,
                    index_p1, index_p2, index_q1,
                    p1, p2, q1,
                    vpm, vpm,
                    nodes);
                  if ( q1_is_between_p1p2 ){
                    //case 7
                    patches_to_keep.reset(patch_id_q1);
                    if(patch_id_p2<patch_id_q2)
                    {
                      uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_p2]);
                      patches_to_keep.reset(patch_id_q2);
                    }
                    else
                    {
                      uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_q2]);
                      patches_to_keep.reset(patch_id_p2);
                    }
                  }
                  else{
                    //case 8
                    patches_to_keep.reset(patch_id_p1);
                    if(patch_id_p2<patch_id_q2)
                    {
                      uf.unify_sets(patch_handles[patch_id_p2], patch_handles[patch_id_q1]);
                      patches_to_keep.reset(patch_id_q2);
                    }
                    else
                    {
                      uf.unify_sets(patch_handles[patch_id_q2], patch_handles[patch_id_q1]);
                      patches_to_keep.reset(patch_id_p2);
                    }
                  }
                  continue;
                }
              }
            }
          }

          CGAL_assertion(
              ( index_p1 == Node_id(-1) ? nodes.to_exact(get(vpm,p1)): nodes.exact_node(index_p1) ) !=
              ( index_q1 == Node_id(-1) ? nodes.to_exact(get(vpm,q1)): nodes.exact_node(index_q1) )
          &&
              ( index_p2 == Node_id(-1) ? nodes.to_exact(get(vpm,p2)): nodes.exact_node(index_p2) ) !=
              ( index_q1 == Node_id(-1) ? nodes.to_exact(get(vpm,q1)): nodes.exact_node(index_q1) )
          &&
              ( index_p1 == Node_id(-1) ? nodes.to_exact(get(vpm,p1)): nodes.exact_node(index_p1) ) !=
              ( index_q2 == Node_id(-1) ? nodes.to_exact(get(vpm,q2)): nodes.exact_node(index_q2) )
          &&
              ( index_p2 == Node_id(-1) ? nodes.to_exact(get(vpm,p2)): nodes.exact_node(index_p2) ) !=
              ( index_q2 == Node_id(-1) ? nodes.to_exact(get(vpm,q2)): nodes.exact_node(index_q2) )
          );

          bool q1_is_between_p1p2 = sorted_around_edge(
            ids.first, ids.second,
            index_p1, index_p2, index_q1,
            p1, p2, q1,
            vpm, vpm,
            nodes);
          bool q2_is_between_p1p2 = sorted_around_edge(
            ids.first, ids.second,
            index_p1, index_p2, index_q2,
            p1, p2, q2,
            vpm, vpm,
            nodes);

          if ( q1_is_between_p1p2 ){
            if( q2_is_between_p1p2 )
            {
              bool p1_is_between_q1q2 = sorted_around_edge(
                ids.first, ids.second,
                index_q1, index_q2, index_p1,
                q1, q2, p1,
                vpm, vpm,
                nodes);
              if (!p1_is_between_q1q2){
                // case (a4)
                uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_p2]);
                patches_to_keep.reset(patch_id_q1);
                patches_to_keep.reset(patch_id_q2);
              }
             else{
               // case (b4)
               patches_to_keep.reset(patch_id_p1);
               patches_to_keep.reset(patch_id_p2);
               patches_to_keep.reset(patch_id_q1);
               patches_to_keep.reset(patch_id_q2);
             }
            }
            else
            {
              //case (c4)
              if ( is_dangling_edge(ids.first, ids.second, h1, tm, is_node_of_degree_one) ||
                   is_dangling_edge(ids.first, ids.second, h2, tm, is_node_of_degree_one) )
              {
                // in case of a surface folding, it might happen that we have a
                // dangling edge that is separating the faces into 2 disjoint
                // patches (think of a sheet twisted at a vertex of the mesh,
                // the degree 1 node being that vertex in the intersection polyline
                // graph).
                // TODO: the condition below is here to avoid removing too many
                //       parts of a mesh (ex: a square in crossing the interior
                //       of a larger). In practice, we could say this is not an
                //       issue if anyway the whole patch would be dropped.
                //       Note that this remark is valid for all cases where
                //       all_fixed is set to false.
                if (patch_id_p1==patch_id_p2 || patch_id_q1==patch_id_q2)
                {
                  all_fixed = false;
                  continue;
                }
              }
              uf.unify_sets(patch_handles[patch_id_p1], patch_handles[patch_id_q2]);
              patches_to_keep.reset(patch_id_p2);
              patches_to_keep.reset(patch_id_q1);
            }
          }
          else
          {
            if( q2_is_between_p1p2 )
            {
              //case (d4)
              if ( is_dangling_edge(ids.first, ids.second, h1, tm, is_node_of_degree_one) ||
                   is_dangling_edge(ids.first, ids.second, h2, tm, is_node_of_degree_one) )
              {
                // same reason as above
                if (patch_id_p1==patch_id_p2 || patch_id_q1==patch_id_q2)
                {
                  all_fixed = false;
                  continue;
                }
              }
              uf.unify_sets(patch_handles[patch_id_p2], patch_handles[patch_id_q1]);
              patches_to_keep.reset(patch_id_q2);
              patches_to_keep.reset(patch_id_p1);
            }
            else
            {
              bool p1_is_between_q1q2 = sorted_around_edge(
                ids.first, ids.second,
                index_q1, index_q2, index_p1,
                q1, q2, p1,
                vpm, vpm,
                nodes);
              if (!p1_is_between_q1q2){
                //case (e4)
                //TODO: This is a "tangency" along an edge here, there is
                //      not much we can do but if one of the two sheets is dropped
                all_fixed = false;
                continue;
              }
              else{
                //case (f4)
                uf.unify_sets(patch_handles[patch_id_q1], patch_handles[patch_id_q2]);
                patches_to_keep.reset(patch_id_p1);
                patches_to_keep.reset(patch_id_p2);
              }
            }
          }
        }
    }

    // the goal here is to remove surface self-intersections. If there are
    // some nested surfaces, they will not be fixed by this function.
    // As a consequence, patch_status_not_set.none() might not be true.

    // use the union-find of incident patches to update the status of patches
    // to keep. If one of the patches in the component is marked as to remove
    // the whole component gets marked.

    // first pass to mark the master patches
    for (typename UF::iterator it=uf.begin(), it_end=uf.end();it!=it_end; ++it)
    {
      if (!patches_to_keep.test(*it))
        patches_to_keep.reset( *uf.find(it) );
    }
    // mark the patch as to be removed if the master is (and not already marked)
    for (typename UF::iterator it=uf.begin(), it_end=uf.end();it!=it_end; ++it)
    {
      if (patches_to_keep.test(*it) && !patches_to_keep.test(*uf.find(it)))
        patches_to_keep.reset( *it );
    }

#ifdef CGAL_COREFINEMENT_DEBUG
    std::cout << "patches_to_keep " <<  patches_to_keep << "\n";
    std::cout << "coplanar_patches " << coplanar_patches << "\n";
    std::cout << "patch_status_not_set " << patch_status_not_set << "\n";
#endif

    //collect edges to stitch before removing patches
    // we could use an_edge_per_polyline but the current version should be cheaper
    // since we don't do any query in intersection_edge_map
    std::vector<edge_descriptor> edges_no_longer_on_intersection;
    std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > hedge_pairs_to_stitch;
    hedge_pairs_to_stitch.reserve(all_intersection_edges_map.size());
    BOOST_FOREACH(const Pair_type& p, all_intersection_edges_map)
    {
      halfedge_descriptor h1 = p.second.h1;
      halfedge_descriptor h2 = p.second.h2;
      halfedge_descriptor h1_opp = opposite(h1, tm);
      halfedge_descriptor h2_opp = opposite(h2, tm);

      if (is_border_edge(h1, tm))
      {
        if (is_border(h1, tm))
        {
          std::size_t patch_id = patch_ids[get(fids,face(opposite(h1,tm),tm))];
          if ( !patches_to_keep.test(patch_id) )
          {
            edges_no_longer_on_intersection.push_back(edge(h1, tm));
            edges_no_longer_on_intersection.push_back(edge(h2, tm));
            continue;
          }
        }
        else
        {
          std::size_t patch_id = patch_ids[get(fids,face(h1,tm))];
          if ( !patches_to_keep.test(patch_id) )
          {
            edges_no_longer_on_intersection.push_back(edge(h1, tm));
            edges_no_longer_on_intersection.push_back(edge(h2, tm));
            continue;
          }
          std::swap(h1, h1_opp);
          std::swap(h2, h2_opp);
        }
      }
      else
      {
        std::size_t patch_id_h1 = patch_ids[get(fids,face(h1,tm))];
        std::size_t patch_id_h1_opp = patch_ids[get(fids,face(h1_opp,tm))];
        if ( patches_to_keep.test(patch_id_h1) ==
             patches_to_keep.test(patch_id_h1_opp))
        {
          edges_no_longer_on_intersection.push_back(edge(h1, tm));
          edges_no_longer_on_intersection.push_back(edge(h2, tm));
          continue;
        }
        if (patches_to_keep.test(patch_id_h1))
        {
          std::swap(h1, h1_opp);
          std::swap(h2, h2_opp);
        }
      }

      if (is_border_edge(h2, tm))
      {
        if (is_border(h2, tm))
        {
          std::size_t patch_id = patch_ids[get(fids,face(opposite(h2,tm),tm))];
          if ( !patches_to_keep.test(patch_id) )
          {
            edges_no_longer_on_intersection.push_back(edge(h1, tm));
            edges_no_longer_on_intersection.push_back(edge(h2, tm));
            continue;
          }
        }
        else
        {
          std::size_t patch_id = patch_ids[get(fids,face(h2,tm))];
          if ( !patches_to_keep.test(patch_id) )
          {
            edges_no_longer_on_intersection.push_back(edge(h1, tm));
            edges_no_longer_on_intersection.push_back(edge(h2, tm));
            continue;
          }
        }
      }
      else
      {
        std::size_t patch_id_h2 = patch_ids[get(fids,face(h2,tm))];
        std::size_t patch_id_h2_opp = patch_ids[get(fids,face(h2_opp,tm))];
        if ( patches_to_keep.test(patch_id_h2) ==
             patches_to_keep.test(patch_id_h2_opp))
        {
          edges_no_longer_on_intersection.push_back(edge(h1, tm));
          edges_no_longer_on_intersection.push_back(edge(h2, tm));
          continue;
        }
      }

      hedge_pairs_to_stitch.push_back( std::make_pair(h1,h2_opp) );
      // unconstrained the halfedge to be removed
      put(ecm, edge(h2_opp,tm), false); // the stitching is always removing the second halfedge
      // force the vertices to be identical for the stitching
      put(vpm, target(h2_opp,tm), get(vpm, source(h1,tm)));
      put(vpm, source(h2_opp,tm), get(vpm, target(h1,tm)));
    }

    // Merge patches to keep only 2: one we keep (1) and one we remove (0)
    const std::size_t PATCH_ID_KEPT = 1;
    BOOST_FOREACH(std::size_t& patch_id, patch_ids)
      patch_id = patches_to_keep.test(patch_id) ? 1 : 0;
    nb_patches=2;
    patches_to_keep=boost::dynamic_bitset<>(2,0);
    patches_to_keep.set(PATCH_ID_KEPT);

    // remove from the set of intersection edges if the patches on both side have
    // the same status.
    BOOST_FOREACH(edge_descriptor e, edges_no_longer_on_intersection)
      intersection_edges.erase(e);

    //store the patch description in a container to avoid recomputing it several times
    typedef Patch_container<TriangleMesh,
                            FaceIdMap,
                            Intersection_edge_map> Patches;
    Patches patches(tm, patch_ids, fids, intersection_edges, nb_patches);
    //remove the extra patch
    remove_patches(tm, ~patches_to_keep,patches, ecm);

    PMP::stitch_borders(tm, hedge_pairs_to_stitch, params::vertex_point_map(vpm));
  }
};


} } // CGAL::Corefinement

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_OUTPUT_BUILDER_FOR_AUTOREFINEMENT_H
