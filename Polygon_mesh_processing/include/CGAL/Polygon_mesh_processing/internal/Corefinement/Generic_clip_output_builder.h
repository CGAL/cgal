// Copyright (c) 2020 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_GENERIC_CLIP_OUTPUT_BUILDER_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_GENERIC_CLIP_OUTPUT_BUILDER_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/property_map.h>
#include <CGAL/Default.h>
#include <CGAL/use.h>

#include <boost/dynamic_bitset.hpp>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace Corefinement {


namespace PMP=Polygon_mesh_processing;
namespace params=PMP::parameters;

template <class TriangleMesh,
          class VertexPointMap1,
          class VertexPointMap2,
          class Ecm1,
          class FaceIdMap1,
          class Kernel_=Default>
class Generic_clip_output_builder
{
//Default typedefs
  typedef typename Default::Get<
    Kernel_,
    typename Kernel_traits<
      typename boost::property_traits<VertexPointMap2>::value_type
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
  // to maintain a halfedge on each polyline per TriangleMesh + pair<bool,size_t>
  // with first = "is the key (pair<Node_id,Node_id>) was reversed?" and
  // second is the number of edges -1 in the polyline
  typedef std::map< Node_id_pair,
                    std::pair< std::map<TriangleMesh*,
                                        halfedge_descriptor>,
                               std::pair<bool,std::size_t> > >
                                              An_edge_per_polyline_map;

  typedef boost::unordered_map<vertex_descriptor, Node_id> Node_id_map;
  typedef boost::unordered_map<edge_descriptor,
                               edge_descriptor>               Edge_map;
//Data members
  TriangleMesh &tm1, &tm2;
  // property maps of input meshes
  const VertexPointMap1 vpm1;
  const VertexPointMap2 vpm2;
  Ecm1 ecm1;
  FaceIdMap1 fids1;
  bool use_compact_clipper;

  // mapping vertex to node id
  Node_id_map vertex_to_node_id1;

  // orientation of input surface meshes
  bool is_tm2_inside_out;
  // constants
  const Node_id NID;

  typename An_edge_per_polyline_map::iterator last_polyline;

  Node_id get_node_id(vertex_descriptor v,
                      const Node_id_map& node_ids)
  {
    typename Node_id_map::const_iterator it = node_ids.find(v);
    if (it == node_ids.end())
      return NID;
    return it->second;
  }

public:

  Generic_clip_output_builder(TriangleMesh& tm1,
                              TriangleMesh& tm2,
                              const VertexPointMap1 vpm1,
                              const VertexPointMap2 vpm2,
                              const Ecm1& ecm1,
                              FaceIdMap1 fids1,
                              bool use_compact_clipper)
    : tm1(tm1), tm2(tm2)
    , vpm1(vpm1), vpm2(vpm2)
    , ecm1(ecm1)
    , fids1(fids1)
    , use_compact_clipper(use_compact_clipper)
    , is_tm2_inside_out( !PMP::is_outward_oriented(tm2, parameters::vertex_point_map(vpm2)) )
    , NID((std::numeric_limits<Node_id>::max)())
  {}

// functions called by the intersection visitor
  void start_new_polyline(Node_id, Node_id) {}
  void add_node_to_polyline(Node_id) {}
  void set_edge_per_polyline(TriangleMesh&, Node_id_pair, halfedge_descriptor){}

  void set_vertex_id(vertex_descriptor v, Node_id node_id, const TriangleMesh& tm)
  {
    CGAL_USE(tm);
    CGAL_assertion(&tm == &tm1);
    vertex_to_node_id1.insert( std::make_pair(v, node_id) );
  }

  template <class Nodes_vector, class Mesh_to_map_node>
  void operator()(
    const Nodes_vector& nodes,
    bool /* input_have_coplanar_faces */,
    const boost::dynamic_bitset<>& /* is_node_of_degree_one */,
    const Mesh_to_map_node&)
  {
    // The property map must be either writable or well-initialized
    if( CGAL::internal::Is_writable_property_map<FaceIdMap1>::value &&
        !BGL::internal::is_index_map_valid(fids1, num_faces(tm1), faces(tm1)) )
    {
      BGL::internal::initialize_face_index_map(fids1, tm1);
    }
    CGAL_assertion(BGL::internal::is_index_map_valid(fids1, num_faces(tm1), faces(tm1)));

    // (1) Assign a patch id to each face indicating in which connected
    // component limited by intersection edges of the surface they are.
    std::vector<std::size_t> tm1_patch_ids( num_faces(tm1),NID );

    std::size_t nb_patches_tm1 =
      PMP::connected_components(tm1,
                                bind_property_maps(fids1,make_property_map(&tm1_patch_ids[0])),
                                params::edge_is_constrained_map(ecm1)
                                       .face_index_map(fids1));

    std::vector <std::size_t> tm1_patch_sizes(nb_patches_tm1, 0);
    for(std::size_t i : tm1_patch_ids)
      if(i!=NID)
        ++tm1_patch_sizes[i];

    // Use the class Side_of_triangle_mesh to classify each patch
    boost::dynamic_bitset<> patch_status_not_set_tm1(nb_patches_tm1);
    patch_status_not_set_tm1.set();

    typedef Side_of_triangle_mesh<TriangleMesh,
                                  Kernel,
                                  VertexPointMap2> Inside_poly_test;

    CGAL::Bounded_side in_tm2 = is_tm2_inside_out
                              ? ON_UNBOUNDED_SIDE : ON_BOUNDED_SIDE;

    Inside_poly_test inside_tm2(tm2, vpm2);
    std::size_t nb_classified=0;
    std::vector<std::size_t> cc_to_keep;
    for(face_descriptor f : faces(tm1))
    {
      const std::size_t f_id = get(fids1, f);
      const std::size_t patch_id = tm1_patch_ids[ f_id ];
      if ( patch_status_not_set_tm1.test( patch_id ) )
      {
        patch_status_not_set_tm1.reset( patch_id );
        halfedge_descriptor h = halfedge(f, tm1);
        Node_id index_p1 = get_node_id(target(h, tm1), vertex_to_node_id1);
        std::array<Node_id, 3> fnids = { index_p1, index_p1, index_p1 };
        if (index_p1 != NID)
        {
          h=next(h, tm1);
          index_p1 = get_node_id(target(h, tm1), vertex_to_node_id1);
          fnids[1]=index_p1;
          if (index_p1 != NID)
          {
            h=next(h, tm1);
            index_p1 = get_node_id(target(h, tm1), vertex_to_node_id1);
            fnids[2]=index_p1;
          }
        }
        if (index_p1 != NID)
        {
          typename Nodes_vector::Exact_kernel ek;
          typedef typename Nodes_vector::Exact_kernel::Point_3 Exact_point_3;
          Exact_point_3 e_centroid = centroid(nodes.exact_node(fnids[0]),
                                              nodes.exact_node(fnids[1]),
                                              nodes.exact_node(fnids[2]));

          Bounded_side position = inside_tm2(e_centroid, ek);
          if ( position==ON_BOUNDARY )
          {
            if (use_compact_clipper)
            {
              cc_to_keep.push_back(patch_id);
            }
          }
          else
            if ( position == in_tm2 )
            {
              cc_to_keep.push_back(patch_id);
            }
        }
        else
        {
          Bounded_side position = inside_tm2( get(vpm1, target(h, tm1)));
          CGAL_assertion( position != ON_BOUNDARY);
          if ( position == in_tm2 )
          {
            cc_to_keep.push_back(patch_id);
          }
        }
        if ( ++nb_classified==nb_patches_tm1) break;
      }
    }

    PMP::keep_connected_components(tm1, cc_to_keep, bind_property_maps(fids1,make_property_map(&tm1_patch_ids[0])));
  }
};


} } } // CGAL::Polygon_mesh_processing::Corefinement

#undef CGAL_COREF_FUNCTION_CALL
#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_GENERIC_CLIP_OUTPUT_BUILDER_H
