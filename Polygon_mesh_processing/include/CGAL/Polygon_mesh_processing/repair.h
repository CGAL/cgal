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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>


#include <set>
#include <vector>
#include <boost/algorithm/minmax_element.hpp>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Union_find.h>
#include <CGAL/algorithm.h>
#include <CGAL/array.h>

// headers for self-intersection removal
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/boost/graph/selection.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/IO/OFF_reader.h>
#include <iostream>
#include <fstream>
#endif

namespace CGAL{
namespace Polygon_mesh_processing {

namespace debug{
  template <class TriangleMesh, class VertexPointMap>
  std::ostream& dump_edge_neighborhood(
    typename boost::graph_traits<TriangleMesh>::edge_descriptor ed,
    TriangleMesh& tmesh,
    const VertexPointMap& vpmap,
    std::ostream& out)
  {
    typedef boost::graph_traits<TriangleMesh> GT;
    typedef typename GT::halfedge_descriptor halfedge_descriptor;
    typedef typename GT::vertex_descriptor vertex_descriptor;
    typedef typename GT::face_descriptor face_descriptor;

    halfedge_descriptor h = halfedge(ed, tmesh);

    std::map<vertex_descriptor, int> vertices;
    std::set<face_descriptor> faces;
    int vindex=0;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(h, tmesh))
    {
      if ( vertices.insert(std::make_pair(source(hd, tmesh), vindex)).second )
        ++vindex;
      if (!is_border(hd, tmesh))
        faces.insert( face(hd, tmesh) );
    }

    h=opposite(h, tmesh);
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(h, tmesh))
    {
      if ( vertices.insert(std::make_pair(source(hd, tmesh), vindex)).second )
        ++vindex;
      if (!is_border(hd, tmesh))
        faces.insert( face(hd, tmesh) );
    }

    std::vector<vertex_descriptor> ordered_vertices(vertices.size());
    typedef std::pair<const vertex_descriptor, int> Pair_type;
    BOOST_FOREACH(const Pair_type& p, vertices)
      ordered_vertices[p.second]=p.first;

    out << "OFF\n" << ordered_vertices.size() << " " << faces.size() << " 0\n";
    BOOST_FOREACH(vertex_descriptor vd, ordered_vertices)
      out << get(vpmap, vd) << "\n";
    BOOST_FOREACH(face_descriptor fd, faces)
    {
      out << "3";
      h=halfedge(fd,tmesh);
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(h, tmesh))
        out << " " << vertices[target(hd, tmesh)];
      out << "\n";
    }
    return out;
  }

  template <class FaceRange, class TriangleMesh>
  void dump_cc_faces(const FaceRange& cc_faces, const TriangleMesh& tm, std::ostream& output)
  {
    typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type Vpm;
    typedef typename boost::property_traits<Vpm>::value_type Point_3;
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

    Vpm vpm = get(boost::vertex_point, tm);

    int id=0;
    std::map<vertex_descriptor, int> vids;
    BOOST_FOREACH(face_descriptor f, cc_faces)
    {
      if ( vids.insert( std::make_pair( target(halfedge(f, tm), tm), id) ).second ) ++id;
      if ( vids.insert( std::make_pair( target(next(halfedge(f, tm), tm), tm), id) ).second ) ++id;
      if ( vids.insert( std::make_pair( target(next(next(halfedge(f, tm), tm), tm), tm), id) ).second ) ++id;
    }
    output << std::setprecision(17);
    output << "OFF\n" << vids.size() << " " << cc_faces.size() << " 0\n";
    std::vector<Point_3> points(vids.size());
    typedef std::pair<const vertex_descriptor, int> Pair_type;
    BOOST_FOREACH(Pair_type p, vids)
      points[p.second]=get(vpm, p.first);
    BOOST_FOREACH(Point_3 p, points)
      output << p << "\n";
    BOOST_FOREACH(face_descriptor f, cc_faces)
    {
      output << "3 "
             << vids[ target(halfedge(f, tm), tm) ] << " "
             << vids[ target(next(halfedge(f, tm), tm), tm) ] << " "
             << vids[ target(next(next(halfedge(f, tm), tm), tm), tm) ] << "\n";
    }
  }
} //end of namespace debug

template <class HalfedgeGraph, class VertexPointMap, class Traits>
struct Less_vertex_point{
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor vertex_descriptor;
  const Traits& m_traits;
  const VertexPointMap& m_vpmap;
  Less_vertex_point(const Traits& traits, const VertexPointMap& vpmap)
    : m_traits(traits)
    , m_vpmap(vpmap) {}
  bool operator()(vertex_descriptor v1, vertex_descriptor v2) const{
    return m_traits.less_xyz_3_object()(get(m_vpmap, v1), get(m_vpmap, v2));
  }
};

///\cond SKIP_IN_MANUAL

template <class Traits, class TriangleMesh, class VertexPointMap, class OutputIterator>
OutputIterator
degenerate_faces(const TriangleMesh& tm,
                 const VertexPointMap& vpmap,
                 const Traits& traits,
                 OutputIterator out)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  BOOST_FOREACH(face_descriptor fd, faces(tm))
  {
    if ( is_degenerate_triangle_face(fd, tm, vpmap, traits) )
      *out++=fd;
  }
  return out;
}

template <class TriangleMesh, class OutputIterator>
OutputIterator
degenerate_faces(const TriangleMesh& tm, OutputIterator out)
{
  typedef typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type Vpm;
  typedef typename boost::property_traits<Vpm>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;

  return degenerate_faces(tm, get(vertex_point, tm), Kernel(), out);
}

// this function remove a border edge even if it does not satisfy the link condition.
// The only limitation is that the length connected component of the boundary this edge
// is strictly greater than 3
template <class TriangleMesh>
typename boost::graph_traits<TriangleMesh>::vertex_descriptor
remove_a_border_edge(typename boost::graph_traits<TriangleMesh>::edge_descriptor ed,
                     TriangleMesh& tm)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  halfedge_descriptor h=halfedge(ed,tm);

  if ( is_border(h,tm) )
    h=opposite(h,tm);
  halfedge_descriptor opp_h = opposite(h,tm);
  CGAL_assertion(is_border(opp_h,tm));
  CGAL_assertion(!is_border(h,tm));

  CGAL_assertion(next(next(opp_h, tm), tm) !=opp_h); // not working for a hole made of 2 edges
  CGAL_assertion(next(next(next(opp_h, tm), tm), tm) !=opp_h); // not working for a hole make of 3 edges

  if (CGAL::Euler::does_satisfy_link_condition(edge(h,tm),tm))
    return CGAL::Euler::collapse_edge(ed, tm);

  // collect edges that have one vertex in the link of
  // the vertices of h and one of the vertex of h as other vertex
  std::set<edge_descriptor> common_incident_edges;
  BOOST_FOREACH(halfedge_descriptor hos, halfedges_around_source(h, tm))
    BOOST_FOREACH(halfedge_descriptor hot, halfedges_around_target(h, tm))
    {
      if( target(hos, tm) == source(hot, tm) )
      {
        common_incident_edges.insert( edge(hot, tm) );
        common_incident_edges.insert( edge(hos, tm) );
      }
    }

  // in the following loop, we visit define a connected component of
  // faces bounded by edges in common_incident_edges and h. We look
  // for the maximal one. This set of faces is the one that will
  // disappear while collapsing ed
  std::set<face_descriptor> marked_faces;

  std::vector<halfedge_descriptor> queue;
  queue.push_back( opposite(next(h,tm), tm) );
  queue.push_back( opposite(prev(h,tm), tm) );
  marked_faces.insert( face(h, tm) );

  do{
    std::vector<halfedge_descriptor> boundary;
    while(!queue.empty())
    {
      halfedge_descriptor back=queue.back();
      queue.pop_back();
      face_descriptor fback=face(back,tm);
      if (common_incident_edges.count(edge(back,tm)))
      {
        boundary.push_back(back);
        continue;
      }
      if ( !marked_faces.insert(fback).second )
        continue;
      queue.push_back( opposite(next(back,tm), tm) );
      queue.push_back( opposite(prev(back,tm), tm) );
    }
    CGAL_assertion( boundary.size() == 2 );
    common_incident_edges.erase( edge(boundary[0], tm) );
    common_incident_edges.erase( edge(boundary[1], tm) );

    queue.push_back(boundary[0]);
    queue.push_back(boundary[1]);
  }
  while(!common_incident_edges.empty());

  // hk1 and hk2 are bounding the region that will be removed.
  // The edge of hk2 will be removed and hk2 will be replaced
  // by the opposite edge of hk1
  halfedge_descriptor hk1=queue.front();
  halfedge_descriptor hk2=queue.back();
  if ( target(hk1,tm)!=source(hk2,tm) )
    std::swap(hk1, hk2);

  CGAL_assertion( target(hk1,tm)==source(hk2,tm) );
  CGAL_assertion( source(hk1,tm)==source(h,tm) );
  CGAL_assertion( target(hk2,tm)==target(h,tm) );


  // collect vertices and edges to remove and do remove faces
  std::set<edge_descriptor> edges_to_remove;
  std::set<vertex_descriptor> vertices_to_remove;
  BOOST_FOREACH(face_descriptor fd, marked_faces)
  {
    halfedge_descriptor hd=halfedge(fd, tm);
    for(int i=0; i<3; ++i)
    {
      edges_to_remove.insert( edge(hd, tm) );
      vertices_to_remove.insert( target(hd,tm) );
      hd=next(hd, tm);
    }
  }

  vertex_descriptor vkept=source(hk1,tm);

  //back-up next, prev halfedge to be restore pointers after removal
  halfedge_descriptor hp=prev(opp_h, tm);
  halfedge_descriptor hn=next(opp_h, tm);
  halfedge_descriptor hk1_opp_next = next(hk2, tm);
  halfedge_descriptor hk1_opp_prev = prev(hk2, tm);
  face_descriptor hk1_opp_face = face(hk2,tm);

  // we will remove the target of hk2, update vertex pointers
  BOOST_FOREACH(halfedge_descriptor hot,
                halfedges_around_target(hk2, tm))
  {
    set_target(hot, vkept, tm);
  }

  // update halfedge pointers since hk2 will be removed
  set_halfedge(vkept, opposite(hk1, tm), tm);
  set_halfedge(target(hk1,tm), hk1, tm);

  // do not remove hk1 and its vertices
  vertices_to_remove.erase( vkept );
  vertices_to_remove.erase( target(hk1, tm) );
  edges_to_remove.erase( edge(hk1,tm) );

  bool hk2_equals_hp = hk2==hp;
  CGAL_assertion( is_border(hk2, tm) == hk2_equals_hp );

  /*
  - case hk2!=hp

         /\      /
     hk1/  \hk2 /
       /    \  /
  ____/______\/____
  hn   h_opp   hp

  - case hk2==hp

         /\
     hk1/  \hk2 == hp
       /    \
  ____/______\
  hn   h_opp
  */

  // remove vertices
  BOOST_FOREACH(vertex_descriptor vd, vertices_to_remove)
    remove_vertex(vd, tm);
  // remove edges
  BOOST_FOREACH(edge_descriptor ed, edges_to_remove)
    remove_edge(ed, tm);
  // remove faces
  BOOST_FOREACH(face_descriptor fd, marked_faces)
    remove_face(fd, tm);

  // now update pointers
  set_face(opposite(hk1, tm), hk1_opp_face, tm);
  if (!hk2_equals_hp)
  {
    set_next(hp, hn, tm);
    set_next(opposite(hk1, tm), hk1_opp_next, tm);
    set_next(hk1_opp_prev, opposite(hk1, tm), tm);
    set_halfedge(hk1_opp_face, opposite(hk1, tm), tm);
  }
  else
  {
    set_next(hk1_opp_prev, opposite(hk1, tm), tm);
    set_next(opposite(hk1, tm), hn, tm);
  }
  return vkept;
}

template <class EdgeRange, class TriangleMesh, class NamedParameters>
std::size_t remove_null_edges(
                       const EdgeRange& edge_range,
                       TriangleMesh& tmesh,
                       const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tmesh));

  using boost::get_param;
  using boost::choose_param;

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  typedef typename GetVertexPointMap<TM, NamedParameters>::type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_property_map(vertex_point, tmesh));
  typedef typename GetGeomTraits<TM, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  std::size_t nb_deg_faces = 0;

  // collect edges of length 0
  std::set<edge_descriptor> null_edges_to_remove;
  BOOST_FOREACH(edge_descriptor ed, edge_range)
  {
    if ( traits.equal_3_object()(get(vpmap, target(ed, tmesh)), get(vpmap, source(ed, tmesh))) )
      null_edges_to_remove.insert(ed);
  }

  #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
  std::cout << "Found " << null_edges_to_remove.size() << " null edges.\n";
  #endif

  while (!null_edges_to_remove.empty())
  {
    edge_descriptor ed = *null_edges_to_remove.begin();
    null_edges_to_remove.erase(null_edges_to_remove.begin());

    halfedge_descriptor h = halfedge(ed, tmesh);

    if (CGAL::Euler::does_satisfy_link_condition(ed,tmesh))
    {
      // remove edges that could also be set for removal
      if ( face(h, tmesh)!=GT::null_face() )
      {
        ++nb_deg_faces;
        null_edges_to_remove.erase(edge(prev(h, tmesh), tmesh));
      }
      if (face(opposite(h, tmesh), tmesh)!=GT::null_face())
      {
        ++nb_deg_faces;
        null_edges_to_remove.erase(edge(prev(opposite(h, tmesh), tmesh), tmesh));
      }
      //now remove the edge
      CGAL::Euler::collapse_edge(ed, tmesh);
    }
    else{
      //handle the case when the edge is incident to a triangle hole
      //we first fill the hole and try again
      if ( is_border(ed, tmesh) )
      {
        halfedge_descriptor hd = halfedge(ed,tmesh);
        if (!is_border(hd,tmesh)) hd=opposite(hd,tmesh);
        if (is_triangle(hd, tmesh))
        {
          Euler::fill_hole(hd, tmesh);
          null_edges_to_remove.insert(ed);
          continue;
        }
      }

      // When the edge does not satisfy the link condition, it means that it cannot be
      // collapsed as is. In the following we assume that there is no topological issue
      // with contracting the edge (no volume will disappear).
      // We start by marking the faces that are incident to an edge endpoint.
      // If the set of marked faces is a topologically disk, then we simply remove all the simplicies
      // inside the disk and star the hole with the edge vertex kept.
      // If the set of marked faces is not a topological disk, it has some non-manifold vertices
      // on its boundary. We need to mark additional faces to make it a topological disk.
      // We can then apply the star hole procedure.
      // Right now we additionally mark the smallest connected components of non-marked faces
      // (using the numnber of faces)

      //backup central point
      typename Traits::Point_3 pt = get(vpmap, source(ed, tmesh));

      // mark faces of the link of each endpoints of the edge which collapse is not topologically valid
      std::set<face_descriptor> marked_faces;
      //   first endpoint
      BOOST_FOREACH( halfedge_descriptor hd, CGAL::halfedges_around_target(halfedge(ed,tmesh), tmesh) )
        if (!is_border(hd,tmesh)) marked_faces.insert( face(hd, tmesh) );
      //   second endpoint
      BOOST_FOREACH( halfedge_descriptor hd, CGAL::halfedges_around_target(opposite(halfedge(ed, tmesh), tmesh), tmesh) )
        if (!is_border(hd,tmesh)) marked_faces.insert( face(hd, tmesh) );

      // extract the halfedges on the boundary of the marked region
      std::vector<halfedge_descriptor> border;
      BOOST_FOREACH(face_descriptor fd, marked_faces)
        BOOST_FOREACH(halfedge_descriptor hd, CGAL::halfedges_around_face(halfedge(fd,tmesh), tmesh))
        {
          halfedge_descriptor hd_opp = opposite(hd, tmesh);
          if ( is_border(hd_opp, tmesh) ||
               marked_faces.count( face(hd, tmesh) )!=
               marked_faces.count( face(hd_opp, tmesh) ) )
          {
            border.push_back( hd );
          }
        }
      CGAL_assertion( !border.empty() ); // a whole connected component got selected and will disappear (not handled for now)
      // define cc of border halfedges: two halfedges are in the same cc
      // if they are on the border of the cc of non-marked faces.
      typedef CGAL::Union_find<halfedge_descriptor> UF_ds;
      UF_ds uf;
      std::map<halfedge_descriptor, typename UF_ds::handle> handles;
      // one cc per border halfedge
      BOOST_FOREACH(halfedge_descriptor hd, border)
        handles.insert( std::make_pair(hd, uf.make_set(hd)) );

      // join cc's
      BOOST_FOREACH(halfedge_descriptor hd, border)
      {
        CGAL_assertion( marked_faces.count( face( hd, tmesh) ) > 0);
        CGAL_assertion( marked_faces.count( face( opposite(hd, tmesh), tmesh) ) == 0 );
        halfedge_descriptor candidate = hd;

        do{
          candidate = prev( opposite(candidate, tmesh), tmesh );
        } while( !marked_faces.count( face( opposite(candidate, tmesh), tmesh) ) );
        uf.unify_sets( handles[hd], handles[opposite(candidate, tmesh)] );
      }

      std::size_t nb_cc = uf.number_of_sets();
      if ( nb_cc != 1 )
      {
        // if more than one connected component is found then the patch
        // made of marked faces contains "non-manifold" vertices.
        // The smallest components need to be marked so that the patch
        // made of marked faces is a topological disk

        // we will explore in parallel the connected components and will stop
        // when all but one connected component have been entirely explored.
        // We add one face at a time for each cc in order to not explore a
        // potentially very large cc.
        std::vector< std::vector<halfedge_descriptor> > stacks_per_cc(nb_cc);
        std::vector< std::set<face_descriptor> > faces_per_cc(nb_cc);
        std::vector< bool > exploration_finished(nb_cc, false);


        // init the stacks of halfedges using the cc of the boundary
        std::size_t index=0;
        std::map< halfedge_descriptor, std::size_t > ccs;
        typedef std::pair<const halfedge_descriptor, typename UF_ds::handle> Pair_type;
        BOOST_FOREACH(Pair_type p, handles)
        {
          halfedge_descriptor opp_hedge = opposite(p.first, tmesh);
          if (is_border(opp_hedge, tmesh)) continue; // nothing to do on the boundary

          typedef typename std::map< halfedge_descriptor, std::size_t >::iterator Map_it;
          std::pair<Map_it, bool> insert_res=
            ccs.insert( std::make_pair(*uf.find( p.second ), index) );
          if (insert_res.second) ++index;

          stacks_per_cc[ insert_res.first->second ].push_back( prev(opp_hedge, tmesh) );
          stacks_per_cc[ insert_res.first->second ].push_back( next(opp_hedge, tmesh) );
          faces_per_cc[ insert_res.first->second ].insert( face(opp_hedge, tmesh) );
        }

        std::size_t nb_ccs_to_be_explored = nb_cc;
        index=0;
        //explore the cc's
        do{
          // try to extract one more face for a given cc
          do{
            CGAL_assertion( !exploration_finished[index] );
            halfedge_descriptor hd = stacks_per_cc[index].back();
            stacks_per_cc[index].pop_back();
            hd = opposite(hd, tmesh);
            if ( !is_border(hd,tmesh) && !marked_faces.count(face(hd, tmesh) ) )
            {
              if ( faces_per_cc[index].insert( face(hd, tmesh) ).second )
              {
                stacks_per_cc[index].push_back( next(hd, tmesh) );
                stacks_per_cc[index].push_back( prev(hd, tmesh) );
                break;
              }
            }
            if (stacks_per_cc[index].empty()) break;
          }
          while(true);
          // the exploration of a cc is finished when its stack is empty
          exploration_finished[index]=stacks_per_cc[index].empty();
          if ( exploration_finished[index] ) --nb_ccs_to_be_explored;
          if ( nb_ccs_to_be_explored==1 ) break;
          while ( exploration_finished[(++index)%nb_cc] );
          index=index%nb_cc;
        }while(true);

        /// \todo use the area criteria? this means maybe continue exploration of larger cc
        // mark faces of completetly explored cc
        for (index=0; index< nb_cc; ++index)
          if( exploration_finished[index] )
          {
            BOOST_FOREACH(face_descriptor fd, faces_per_cc[index])
              marked_faces.insert(fd);
          }
      }

      // collect simplices to be removed
      std::set<vertex_descriptor> vertices_to_keep;
      std::set<halfedge_descriptor> halfedges_to_keep;
      BOOST_FOREACH(halfedge_descriptor hd, border)
        if (  !marked_faces.count(face(opposite(hd, tmesh), tmesh)) )
        {
          halfedges_to_keep.insert( hd );
          vertices_to_keep.insert( target(hd, tmesh) );
        }

      // backup next,prev relationships to set after patch removal
      std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > next_prev_halfedge_pairs;
      halfedge_descriptor first_border_hd=*( halfedges_to_keep.begin() );
      halfedge_descriptor current_border_hd=first_border_hd;
      do{
        halfedge_descriptor prev_border_hd=current_border_hd;
        current_border_hd=next(current_border_hd, tmesh);
        while( marked_faces.count( face( opposite(current_border_hd, tmesh), tmesh) ) )
          current_border_hd=next(opposite(current_border_hd, tmesh), tmesh);
        next_prev_halfedge_pairs.push_back( std::make_pair(prev_border_hd, current_border_hd) );
      }while(current_border_hd!=first_border_hd);

      // collect vertices and edges to remove and do remove faces
      std::set<edge_descriptor> edges_to_remove;
      std::set<vertex_descriptor> vertices_to_remove;
      BOOST_FOREACH(face_descriptor fd, marked_faces)
      {
        halfedge_descriptor hd=halfedge(fd, tmesh);
        for(int i=0; i<3; ++i)
        {
          if ( !halfedges_to_keep.count(hd) )
            edges_to_remove.insert( edge(hd, tmesh) );
          if ( !vertices_to_keep.count(target(hd,tmesh)) )
            vertices_to_remove.insert( target(hd,tmesh) );
          hd=next(hd, tmesh);
        }
        remove_face(fd, tmesh);
      }

      // remove vertices
      BOOST_FOREACH(vertex_descriptor vd, vertices_to_remove)
        remove_vertex(vd, tmesh);
      // remove edges
      BOOST_FOREACH(edge_descriptor ed, edges_to_remove)
      {
        null_edges_to_remove.erase(ed);
        remove_edge(ed, tmesh);
      }

      // add a new face, set all border edges pointing to it
      // and update halfedge vertex of patch boundary vertices
      face_descriptor new_face = add_face(tmesh);
      typedef std::pair<halfedge_descriptor, halfedge_descriptor> Pair_type;
      BOOST_FOREACH(const Pair_type& p, next_prev_halfedge_pairs)
      {
        set_face(p.first, new_face, tmesh);
        set_next(p.first, p.second, tmesh);
        set_halfedge(target(p.first, tmesh), p.first, tmesh);
      }
      set_halfedge(new_face, first_border_hd, tmesh);
      // triangulate the new face and update the coordinate of the central vertex
      halfedge_descriptor new_hd=Euler::add_center_vertex(first_border_hd, tmesh);
      put(vpmap, target(new_hd, tmesh), pt);

      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(new_hd, tmesh))
        if ( traits.equal_3_object()(get(vpmap, target(hd, tmesh)), get(vpmap, source(hd, tmesh))) )
          null_edges_to_remove.insert(edge(hd, tmesh));

      CGAL_assertion( is_valid(tmesh) );
    }
  }

  return nb_deg_faces;
}

template <class EdgeRange, class TriangleMesh>
std::size_t remove_null_edges(
                       const EdgeRange& edge_range,
                       TriangleMesh& tmesh)
{
  return remove_null_edges(edge_range, tmesh,
                           parameters::all_default());
}

/// \ingroup PMP_repairing_grp
/// removes the degenerate faces from a triangulated surface mesh.
/// A face is considered degenerate if two of its vertices share the same location,
/// or more generally if all its vertices are collinear.
///
/// @pre `CGAL::is_triangle_mesh(tmesh)`
///
/// @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param tmesh the  triangulated surface mesh to be repaired
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`. The type of this map is model of `ReadWritePropertyMap`.
/// If this parameter is omitted, an internal property map for
/// `CGAL::vertex_point_t` should be available in `TriangleMesh`
/// \cgalParamEnd
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///       The traits class must provide the nested type `Point_3`,
///       and the nested functors :
///         - `Compare_distance_3` to compute the distance between 2 points
///         - `Collinear_are_ordered_along_line_3` to check whether 3 collinear points are ordered
///         - `Collinear_3` to check whether 3 points are collinear
///         - `Less_xyz_3` to compare lexicographically two points
///         - `Equal_3` to check whether 2 points are identical
///         -  for each functor Foo, a function `Foo foo_object()`
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \todo the function might not be able to remove all degenerate faces.
///       We should probably do something with the return type.
/// \return number of removed degenerate faces
template <class TriangleMesh, class NamedParameters>
std::size_t remove_degenerate_faces(TriangleMesh& tmesh,
                                    const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tmesh));

  using boost::get_param;
  using boost::choose_param;

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  typedef typename GetVertexPointMap<TM, NamedParameters>::type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_property_map(vertex_point, tmesh));
  typedef typename GetGeomTraits<TM, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref;
// First remove edges of length 0
  std::size_t nb_deg_faces = remove_null_edges(edges(tmesh), tmesh, np);

  #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
  {
  std::cout <<"Done with null edges.\n";
  std::ofstream output("/tmp/no_null_edges.off");
  output << std::setprecision(17) << tmesh << "\n";
  output.close();
  }
  #endif

// Then, remove triangles made of 3 collinear points
  std::set<face_descriptor> degenerate_face_set;
  BOOST_FOREACH(face_descriptor fd, faces(tmesh))
    if ( is_degenerate_triangle_face(fd, tmesh, vpmap, traits) )
      degenerate_face_set.insert(fd);
  nb_deg_faces+=degenerate_face_set.size();

  // first remove degree 3 vertices that are part of a cap
  // (only the vertex in the middle of the opposite edge)
  // This removal does not change the shape of the mesh.
  while (!degenerate_face_set.empty())
  {
    std::set<vertex_descriptor> vertices_to_remove;
    BOOST_FOREACH(face_descriptor fd, degenerate_face_set)
    {
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
        vertex_descriptor vd = target(hd, tmesh);
        if (degree(vd, tmesh) == 3)
        {
          vertices_to_remove.insert(vd);
          break;
        }
      }
    }

    BOOST_FOREACH(vertex_descriptor vd, vertices_to_remove)
    {
      halfedge_descriptor hd=halfedge(vd, tmesh);
      BOOST_FOREACH(halfedge_descriptor hd2, halfedges_around_target(hd, tmesh))
        if (!is_border(hd2, tmesh))
          degenerate_face_set.erase( face(hd2, tmesh) );
      // remove the central vertex and check if the new face is degenerated
      hd=CGAL::Euler::remove_center_vertex(hd, tmesh);
      if (is_degenerate_triangle_face(face(hd, tmesh), tmesh, vpmap, traits))
      {
        degenerate_face_set.insert( face(hd, tmesh) );
      }
    }
    if (vertices_to_remove.empty()) break;
  }

  while (!degenerate_face_set.empty())
  {
    #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
    std::cout << "Loop on removing deg faces\n";
    // ensure the mesh is not broken
    {
      std::ofstream out("/tmp/out.off");
      out << tmesh;
      out.close();

      std::vector<typename Traits::Point_3> points;
      std::vector<std::vector<std::size_t> > triangles;
      std::ifstream in("/tmp/out.off");
      CGAL::read_OFF(in, points, triangles);
      if (!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles))
      {
        std::cerr << "ERROR: got a polygon soup!\n";
        exit(EXIT_FAILURE);
      }
    }
    #endif

    face_descriptor fd = *degenerate_face_set.begin();

    // look whether an incident triangle is also degenerated
    bool detect_cc_of_degenerate_triangles = false;
    BOOST_FOREACH(halfedge_descriptor hd,
                  halfedges_around_face(halfedge(fd, tmesh), tmesh) )
    {
      face_descriptor adjacent_face = face( opposite(hd, tmesh), tmesh );
      if ( adjacent_face!=GT::null_face() &&
           degenerate_face_set.count(adjacent_face) )
      {
        detect_cc_of_degenerate_triangles = true;
        break;
      }
    }

    if (!detect_cc_of_degenerate_triangles)
    {
      #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      std::cout << "  no degenerate neighbors, using a flip.\n";
      #endif
      degenerate_face_set.erase(degenerate_face_set.begin());
    // flip the longest edge of the triangle
      Point_ref p1 = get(vpmap, target( halfedge(fd, tmesh), tmesh) );
      Point_ref p2 = get(vpmap, target(next(halfedge(fd, tmesh), tmesh), tmesh) );
      Point_ref p3 = get(vpmap, source( halfedge(fd, tmesh), tmesh) );

      CGAL_assertion(p1!=p2 && p1!=p3 && p2!=p3);

      typename Traits::Compare_distance_3 compare_distance = traits.compare_distance_3_object();

      halfedge_descriptor edge_to_flip;
      if (compare_distance(p1,p2, p1,p3) != CGAL::SMALLER) // p1p2 > p1p3
      {
        if (compare_distance(p1,p2, p2,p3) != CGAL::SMALLER) // p1p2 > p2p3
          // flip p1p2
          edge_to_flip = next( halfedge(fd, tmesh), tmesh );
        else
          // flip p2p3
          edge_to_flip = prev( halfedge(fd, tmesh), tmesh );
      }
      else
        if (compare_distance(p1,p3, p2,p3) != CGAL::SMALLER) // p1p3>p2p3
          //flip p3p1
          edge_to_flip = halfedge(fd, tmesh);
        else
          //flip p2p3
          edge_to_flip = prev( halfedge(fd, tmesh), tmesh );

      face_descriptor opposite_face=face( opposite(edge_to_flip, tmesh), tmesh);
      if ( opposite_face == GT::null_face() )
        // simply remove the face
        Euler::remove_face(edge_to_flip, tmesh);
      else
      {
        // condition for the flip to be valid (the edge to be created do not already exists)
        if ( !halfedge(target(next(edge_to_flip, tmesh), tmesh),
                       target(next(opposite(edge_to_flip, tmesh), tmesh), tmesh),
                       tmesh).second )
        {
          Euler::flip_edge(edge_to_flip, tmesh);
        }
        #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        else{
          std::cout << "  WARNING: flip is not possible\n";
          // \todo Let p and q be the vertices opposite to `edge_to_flip`, and let
          //       r be the vertex of `edge_to_flip` that is the furthest away from
          //       the edge `pq`. In that case I think we should remove all the triangles
          //       so that the triangle pqr is in the mesh.
        }
        #endif
      }
    }
    else
    {
    // Process a connected component of degenerate faces
      // get all the faces from the connected component
      // and the boundary edges
      std::set<face_descriptor> cc_faces;
      std::vector<face_descriptor> queue;
      std::vector<halfedge_descriptor> boundary_hedges;
      std::vector<halfedge_descriptor> inside_hedges;
      queue.push_back(fd);
      cc_faces.insert(fd);

      while(!queue.empty())
      {
        face_descriptor top=queue.back();
        queue.pop_back();
        BOOST_FOREACH(halfedge_descriptor hd,
                      halfedges_around_face(halfedge(top, tmesh), tmesh) )
        {
          face_descriptor adjacent_face = face( opposite(hd, tmesh), tmesh );
          if ( adjacent_face==GT::null_face() ||
               degenerate_face_set.count(adjacent_face)==0 )
            boundary_hedges.push_back(hd);
          else
          {
            if (cc_faces.insert(adjacent_face).second)
              queue.push_back(adjacent_face);
            if ( hd < opposite(hd, tmesh) )
              inside_hedges.push_back(hd);
          }
        }
      }

      #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      std::cout << "  Deal with a cc of " << cc_faces.size() << " degenerate faces.\n";
      /// dump cc_faces
      {
      int id=0;
      std::map<vertex_descriptor, int> vids;
      BOOST_FOREACH(face_descriptor f, cc_faces)
      {
        if ( vids.insert( std::make_pair( target(halfedge(f, tmesh), tmesh), id) ).second ) ++id;
        if ( vids.insert( std::make_pair( target(next(halfedge(f, tmesh), tmesh), tmesh), id) ).second ) ++id;
        if ( vids.insert( std::make_pair( target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh), id) ).second ) ++id;
      }
      std::ofstream output("/tmp/cc_faces.off");
      output << std::setprecision(17);
      output << "OFF\n" << vids.size() << " " << cc_faces.size() << " 0\n";
      std::vector<typename Traits::Point_3> points(vids.size());
      typedef std::pair<const vertex_descriptor, int> Pair_type;
      BOOST_FOREACH(Pair_type p, vids)
        points[p.second]=get(vpmap, p.first);
      BOOST_FOREACH(typename Traits::Point_3 p, points)
        output << p << "\n";
      BOOST_FOREACH(face_descriptor f, cc_faces)
      {
        output << "3 "
               << vids[ target(halfedge(f, tmesh), tmesh) ] << " "
               << vids[ target(next(halfedge(f, tmesh), tmesh), tmesh) ] << " "
               << vids[ target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh) ] << "\n";
      }

      for (std::size_t pid=2; pid!=points.size(); ++pid)
      {
        CGAL_assertion(collinear(points[0], points[1], points[pid]));
      }
      }
      #endif

      // find vertices strictly inside the cc
      std::set<vertex_descriptor> boundary_vertices;
      BOOST_FOREACH(halfedge_descriptor hd, boundary_hedges)
        boundary_vertices.insert( target(hd, tmesh) );
      std::set<vertex_descriptor> inside_vertices;
      BOOST_FOREACH(halfedge_descriptor hd, inside_hedges)
      {
        if (!boundary_vertices.count( target(hd, tmesh) ))
          inside_vertices.insert( target(hd, tmesh) );
        if (!boundary_vertices.count( source(hd, tmesh) ))
          inside_vertices.insert( source(hd, tmesh) );
      }

      // v-e+f = 1 for a topological disk and e = (3f+#boundary_edges)/2
      if (boundary_vertices.size()+inside_vertices.size() -
          (cc_faces.size()+boundary_hedges.size())/2 != 1)
      {
        //cc_faces does not define a topological disk
        /// \todo Find to way to handle that case
        #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "  WARNING: Cannot remove the component of degenerate faces: not a topological disk.\n";
        #endif
        BOOST_FOREACH(face_descriptor f, cc_faces)
          degenerate_face_set.erase(f);
        continue;
      }
    // preliminary step to check if the operation is possible
      // sort the boundary points along the common supporting line
      //    we first need a reference point
      typedef Less_vertex_point<TriangleMesh, VertexPointMap, Traits> Less_vertex;
      std::pair<
        typename std::set<vertex_descriptor>::iterator,
        typename std::set<vertex_descriptor>::iterator > ref_vertices =
        boost::minmax_element( boundary_vertices.begin(),
                               boundary_vertices.end(),
                               Less_vertex(traits, vpmap) );

      //    and then we sort the vertices using this reference point
      typedef std::set<typename Traits::Point_3> Sorted_point_set;
      Sorted_point_set sorted_points;
      BOOST_FOREACH(vertex_descriptor v, boundary_vertices)
        sorted_points.insert( get(vpmap,v) );

      CGAL_assertion(sorted_points.size()==
                     std::set<typename Traits::Point_3>(sorted_points.begin(),
                                                        sorted_points.end()).size());

      CGAL_assertion( get( vpmap, *ref_vertices.first)==*sorted_points.begin() );
      CGAL_assertion( get( vpmap, *ref_vertices.second)==*cpp11::prev(sorted_points.end()) );

      const typename Traits::Point_3& xtrm1 = *sorted_points.begin();
      const typename Traits::Point_3& xtrm2 = *cpp11::prev(sorted_points.end());

      // recover halfedges on the hole, bounded by the reference vertices
      std::vector<halfedge_descriptor> side_one, side_two;

      // look for the outgoing border halfedge of the first extreme point
      BOOST_FOREACH(halfedge_descriptor hd, boundary_hedges)
        if ( get(vpmap, source(hd, tmesh)) == xtrm1 )
        {
          side_one.push_back(hd);
          break;
        }
      CGAL_assertion(side_one.size()==1);

      while( get(vpmap, target(side_one.back(), tmesh)) != xtrm2 )
      {
        vertex_descriptor prev_vertex = target(side_one.back(), tmesh);
        BOOST_FOREACH(halfedge_descriptor hd, boundary_hedges)
          if ( source(hd, tmesh) == prev_vertex )
          {
            side_one.push_back(hd);
            break;
          }
      }

      // look for the outgoing border halfedge of second extreme vertex
      BOOST_FOREACH(halfedge_descriptor hd, boundary_hedges)
        if ( source(hd, tmesh) == target(side_one.back(), tmesh) )
        {
          side_two.push_back(hd);
          break;
        }
      CGAL_assertion(side_two.size()==1);

      while( target(side_two.back(), tmesh) != source(side_one.front(), tmesh) )
      {
        vertex_descriptor prev_vertex = target(side_two.back(), tmesh);
        BOOST_FOREACH(halfedge_descriptor hd, boundary_hedges)
          if ( source(hd, tmesh) == prev_vertex )
          {
            side_two.push_back(hd);
            break;
          }
      }

      CGAL_assertion( side_one.size()+side_two.size()==boundary_hedges.size() );

      // reverse the order of the second side so as to follow
      // the same order than side one
      std::reverse(side_two.begin(), side_two.end());
      BOOST_FOREACH(halfedge_descriptor& h, side_two)
        h=opposite(h, tmesh);

      //make sure the points of the vertices along side_one are correctly sorted
      std::vector<Point_3> side_points;
      side_points.reserve(side_one.size()+1);
      side_points.push_back(get(vpmap,source(side_one.front(), tmesh)));
      BOOST_FOREACH(halfedge_descriptor h, side_one)
        side_points.push_back(get(vpmap,target(h, tmesh)));
      CGAL_assertion(get(vpmap,source(side_one.front(), tmesh))==side_points.front());
      CGAL_assertion(get(vpmap,target(side_one.back(), tmesh))==side_points.back());
      std::sort(side_points.begin(), side_points.end());
      //\todo the reordering could lead to the apparition of null edges.
      CGAL_assertion(std::unique(side_points.begin(), side_points.end())==side_points.end());
      for(std::size_t i=0;i<side_one.size()-1;++i)
        put(vpmap,target(side_one[i], tmesh), side_points[i+1]);

      //same thing for side_two
      side_points.clear();
      side_points.reserve(side_two.size()+1);
      side_points.push_back(get(vpmap,source(side_two.front(), tmesh)));
      BOOST_FOREACH(halfedge_descriptor h, side_two)
        side_points.push_back(get(vpmap,target(h, tmesh)));
      CGAL_assertion(get(vpmap,source(side_two.front(), tmesh))==side_points.front());
      CGAL_assertion(get(vpmap,target(side_two.back(), tmesh))==side_points.back());
      std::sort(side_points.begin(), side_points.end());
      //\todo the reordering could lead to the apparition of null edges.
      CGAL_assertion(std::unique(side_points.begin(), side_points.end())==side_points.end());
      for(std::size_t i=0;i<side_two.size()-1;++i)
        put(vpmap,target(side_two[i], tmesh), side_points[i+1]);


      CGAL_assertion( source(side_one.front(), tmesh) == *ref_vertices.first );
      CGAL_assertion( source(side_two.front(), tmesh) == *ref_vertices.first );
      CGAL_assertion( target(side_one.back(), tmesh) == *ref_vertices.second );
      CGAL_assertion( target(side_two.back(), tmesh) == *ref_vertices.second );

      typename Sorted_point_set::iterator it_pt = cpp11::next(sorted_points.begin()),
                                          it_pt_end = cpp11::prev(sorted_points.end());

      bool non_collapsable = false;
      typename std::vector<halfedge_descriptor>::iterator side_one_it = side_one.begin();
      typename std::vector<halfedge_descriptor>::iterator side_two_it = side_two.begin();
      for(;it_pt!=it_pt_end;++it_pt)
      {
        // check if it_pt is the point of the target of one or two halfedges
        bool target_of_side_one = get(vpmap, target(*side_one_it, tmesh))==*it_pt;
        bool target_of_side_two = get(vpmap, target(*side_two_it, tmesh))==*it_pt;

        if (target_of_side_one && target_of_side_two)
        {
          BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(*side_one_it, tmesh))
          {
            if (source(h, tmesh)==target(*side_two_it, tmesh))
            {
              non_collapsable=true;
              break;
            }
          }
        }
        else{
          CGAL_assertion(target_of_side_one || target_of_side_two);
          vertex_descriptor v1 = target_of_side_one ? target(*side_one_it, tmesh)
                                                    : target(*side_two_it, tmesh);
          vertex_descriptor v2 = target_of_side_two ? target(next(opposite(*side_one_it, tmesh), tmesh), tmesh)
                                                    : target(next(*side_two_it, tmesh), tmesh);
          BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v1, tmesh))
            if (source(h, tmesh)==v2)
            {
              non_collapsable=true;
              break;
            }
        }
        if(non_collapsable) break;
        if (target_of_side_one) ++side_one_it;
        if (target_of_side_two) ++side_two_it;
      }
      if (non_collapsable){
        BOOST_FOREACH(face_descriptor f, cc_faces)
          degenerate_face_set.erase(f);
        #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "  WARNING: cannot remove a connected components of degenerate faces.\n";
        #endif
        continue;
      }

    // now proceed to the fix
      // update the face and halfedge vertex pointers on the boundary
      BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
      {
        set_face(h, GT::null_face(), tmesh);
        set_halfedge(target(h,tmesh), h, tmesh);
      }
      // update next/prev pointers of boundary_hedges
      BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
      {
        halfedge_descriptor next_candidate = next( h, tmesh);
        while (face(next_candidate, tmesh)!=GT::null_face())
          next_candidate = next( opposite( next_candidate, tmesh), tmesh);
        set_next(h, next_candidate, tmesh);
      }
      // remove degenerate faces
      BOOST_FOREACH(face_descriptor f, cc_faces)
      {
        degenerate_face_set.erase(f);
        remove_face(f, tmesh);
      }
      // remove interior edges
      BOOST_FOREACH(halfedge_descriptor h, inside_hedges)
        remove_edge(edge(h, tmesh), tmesh);
      // remove interior vertices
      BOOST_FOREACH(vertex_descriptor v, inside_vertices)
        remove_vertex(v, tmesh);

      #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      std::cout << "  side_one.size() " << side_one.size() << "\n";
      std::cout << "  side_two.size() " << side_two.size() << "\n";
      #endif

      CGAL_assertion( source(side_one.front(), tmesh) == *ref_vertices.first );
      CGAL_assertion( source(side_two.front(), tmesh) == *ref_vertices.first );
      CGAL_assertion( target(side_one.back(), tmesh) == *ref_vertices.second );
      CGAL_assertion( target(side_two.back(), tmesh) == *ref_vertices.second );

      // now split each side to contains the same sequence of points
      //    first side
      int hi=0;
      for (typename Sorted_point_set::iterator it=cpp11::next(sorted_points.begin()),
                                               it_end=sorted_points.end(); it!=it_end; ++it)
      {
        CGAL_assertion( *cpp11::prev(it) == get(vpmap, source(side_one[hi], tmesh) ) );
        if( *it != get(vpmap, target(side_one[hi], tmesh) ) ){
          // split the edge and update the point
          halfedge_descriptor h1 = next(opposite(side_one[hi], tmesh), tmesh);
          put(vpmap,
              target(Euler::split_edge(side_one[hi], tmesh), tmesh),
              *it);
          // split_edge updates the halfedge of the source vertex of h,
          // since we reuse later the halfedge of the first refernce vertex
          // we must set it as we need.
          if ( source(h1,tmesh) == *ref_vertices.first)
            set_halfedge(*ref_vertices.first, prev( prev(side_one[hi], tmesh), tmesh), tmesh );
          // retriangulate the opposite face
          if ( face(h1, tmesh) != GT::null_face())
            Euler::split_face(h1, opposite(side_one[hi], tmesh), tmesh);
        }
        else
          ++hi;
      }
      //    second side
      hi=0;
      for (typename Sorted_point_set::iterator it=cpp11::next(sorted_points.begin()),
                                               it_end=sorted_points.end(); it!=it_end; ++it)
      {
        CGAL_assertion( *cpp11::prev(it) == get(vpmap, source(side_two[hi], tmesh) ) );
        if( *it != get(vpmap, target(side_two[hi], tmesh) ) ){
          // split the edge and update the point
          halfedge_descriptor h2 = Euler::split_edge(side_two[hi], tmesh);
          put(vpmap, target(h2, tmesh), *it);
          // split_edge updates the halfedge of the source vertex of h,
          // since we reuse later the halfedge of the first refernce vertex
          // we must set it as we need.
          if ( source(h2,tmesh) == *ref_vertices.first)
            set_halfedge(*ref_vertices.first, opposite( h2, tmesh), tmesh );
          // retriangulate the face
          if ( face(h2, tmesh) != GT::null_face())
            Euler::split_face(h2, next(side_two[hi], tmesh), tmesh);
        }
        else
          ++hi;
      }

      CGAL_assertion( target(halfedge(*ref_vertices.first, tmesh), tmesh) == *ref_vertices.first );
      CGAL_assertion( face(halfedge(*ref_vertices.first, tmesh), tmesh) == GT::null_face() );

      #ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      {
        halfedge_descriptor h_side2 = halfedge(*ref_vertices.first, tmesh);
        halfedge_descriptor h_side1 = next(h_side2, tmesh);

        do
        {
          CGAL_assertion( get(vpmap, source(h_side1, tmesh)) == get(vpmap, target(h_side2, tmesh)) );
          CGAL_assertion( get(vpmap, target(h_side1, tmesh)) == get(vpmap, source(h_side2, tmesh)) );

          if ( target(next(opposite(h_side1, tmesh), tmesh), tmesh) ==
               target(next(opposite(h_side2, tmesh), tmesh), tmesh) )
          {
            CGAL_assertion(!"Forbidden simplification");
          }

          h_side2 = prev(h_side2, tmesh);
          h_side1 = next(h_side1, tmesh);
        }
        while( target(h_side1, tmesh) != *ref_vertices.second );
      }
      #endif

      // remove side1 and replace its opposite hedges by those of side2
      halfedge_descriptor h_side2 = halfedge(*ref_vertices.first, tmesh);
      halfedge_descriptor h_side1 = next(h_side2, tmesh);
      while(true)
      {
        CGAL_assertion( get(vpmap, source(h_side1, tmesh)) == get(vpmap, target(h_side2, tmesh)) );
        CGAL_assertion( get(vpmap, target(h_side1, tmesh)) == get(vpmap, source(h_side2, tmesh)) );
        // backup target vertex
        vertex_descriptor vertex_to_remove = target(h_side1, tmesh);
        if (vertex_to_remove!=*ref_vertices.second){
          vertex_descriptor replacement_vertex = source(h_side2, tmesh);
          // replace the incident vertex
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(h_side1, tmesh))
            set_target(hd, replacement_vertex, tmesh);
        }
        // prev side2 hedge for next loop
        halfedge_descriptor h_side2_for_next_turn = prev(h_side2, tmesh);
        // replace the opposite of h_side1 by h_side2
        halfedge_descriptor opposite_h_side1 = opposite( h_side1, tmesh);
        face_descriptor the_face = face(opposite_h_side1, tmesh);
        set_face(h_side2, the_face, tmesh);
        if (the_face!=GT::null_face()) set_halfedge(the_face, h_side2, tmesh);
        set_next(h_side2, next(opposite_h_side1, tmesh), tmesh);
        set_next(prev(opposite_h_side1, tmesh), h_side2, tmesh);
        // take the next hedges
        edge_descriptor edge_to_remove = edge(h_side1, tmesh);
        h_side1 = next(h_side1, tmesh);
        // now remove the extra edge
        remove_edge(edge_to_remove, tmesh);
        // ... and the extra vertex if it's not the second reference
        if (vertex_to_remove==*ref_vertices.second)
        {
          // update the halfedge pointer of the last vertex (others were already from side 2)
          CGAL_assertion( target(opposite(h_side2, tmesh), tmesh) == vertex_to_remove );
          set_halfedge(vertex_to_remove, opposite(h_side2, tmesh), tmesh);
          break;
        }
        else
          remove_vertex(vertex_to_remove , tmesh);
        h_side2 = h_side2_for_next_turn;
      }
    }
  }

  return nb_deg_faces;
}


template<class TriangleMesh>
std::size_t remove_degenerate_faces(TriangleMesh& tmesh)
{
  return remove_degenerate_faces(tmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <class TriangleMesh, class Vpm>
std::size_t duplicate_non_manifold_vertices(TriangleMesh& tm, Vpm vpm)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  boost::unordered_set<vertex_descriptor> vertices_handled;
  boost::unordered_set<halfedge_descriptor> halfedges_handled;

  std::size_t nb_new_vertices=0;

  std::vector<halfedge_descriptor> non_manifold_cones;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(tm))
  {
    if (halfedges_handled.insert(h).second)
    {
      vertex_descriptor vd = target(h, tm);
      if ( !vertices_handled.insert(vd).second )
      {
        non_manifold_cones.push_back(h);
      }
      else
        set_halfedge(vd, h, tm);
      halfedge_descriptor start=opposite(next(h, tm), tm);
      h=start;
      do{
        halfedges_handled.insert(h);
        h=opposite(next(h, tm), tm);
      }while(h!=start);
    }
  }

  if (!non_manifold_cones.empty())  {
    BOOST_FOREACH(halfedge_descriptor h, non_manifold_cones)
    {
      halfedge_descriptor start = h;
      vertex_descriptor new_vd = add_vertex(tm);
      ++nb_new_vertices;
      put(vpm, new_vd, get(vpm, target(h, tm)));
      set_halfedge(new_vd, h, tm);
      do{
        set_target(h, new_vd, tm);
        h=opposite(next(h, tm), tm);
      } while(h!=start);
    }
  }

  return nb_new_vertices;
}

template <class TriangleMesh>
std::size_t duplicate_non_manifold_vertices(TriangleMesh& tm)
{
  return duplicate_non_manifold_vertices(tm, get(vertex_point, tm));
}

/// \endcond


/// \ingroup PMP_repairing_grp
/// removes the isolated vertices from any polygon mesh.
/// A vertex is considered isolated if it is not incident to any simplex
/// of higher dimension.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
///
/// @param pmesh the polygon mesh to be repaired
///
/// @return number of removed isolated vertices
///
template <class PolygonMesh>
std::size_t remove_isolated_vertices(PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> to_be_removed;

  BOOST_FOREACH(vertex_descriptor v, vertices(pmesh))
  {
    if (CGAL::halfedges_around_target(v, pmesh).first
      == CGAL::halfedges_around_target(v, pmesh).second)
      to_be_removed.push_back(v);
  }
  std::size_t nb_removed = to_be_removed.size();
  BOOST_FOREACH(vertex_descriptor v, to_be_removed)
  {
    remove_vertex(v, pmesh);
  }
  return nb_removed;
}

/// \cond SKIP_IN_MANUAL
template <class TriangleMesh,  class face_descriptor, class VertexPointMap>
std::pair< bool, bool >
remove_self_intersections_one_step(TriangleMesh& tm,
                                   std::set<face_descriptor>& faces_to_remove,
                                   VertexPointMap& vpmap,
                                   int step,
                                   bool preserve_genus,
                                   bool verbose)
{
  std::set<face_descriptor> faces_to_remove_copy = faces_to_remove;

  if (verbose)
    std::cout << "DEBUG: running remove_self_intersections_one_step, step " << step
              << " with " << faces_to_remove.size() << " intersecting faces\n";

  CGAL_assertion(tm.is_valid());

  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;

  bool something_was_done = false; // indicates if a region was successfully remeshed
  bool all_fixed = true; // indicates if all removal went well
  // indicates if a removal was not possible because the region handle has
  // some boundary cycle of halfedges
  bool topology_issue = false;
  if (verbose)
  {
    std::cout << "  DEBUG: is_valid in one_step(tm)? ";
    std::cout.flush();
    std::cout << is_valid(tm) << "\n";
  }

  if(!faces_to_remove.empty()){

    while(!faces_to_remove.empty())
    {
      // Process a connected component of faces to remove.
      // collect all the faces from the connected component
      std::set<face_descriptor> cc_faces;
      std::vector<face_descriptor> queue(1, *faces_to_remove.begin()); // temporary queue
      cc_faces.insert(queue.back());
      while(!queue.empty())
      {
        face_descriptor top=queue.back();
        queue.pop_back();
        halfedge_descriptor h = halfedge(top,tm);
        for (int i=0;i<3; ++i)
        {
          face_descriptor adjacent_face = face( opposite(h, tm), tm );
          if ( adjacent_face!=boost::graph_traits<TriangleMesh>::null_face())
          {
            if (faces_to_remove.count(adjacent_face) != 0 &&
                cc_faces.insert(adjacent_face).second)
              queue.push_back(adjacent_face);
          }
          h = next(h, tm);
        }
      }

      // expand the region to be filled
      if (step > 0)
        expand_face_selection(cc_faces, tm, step,
                              make_boolean_property_map(cc_faces),
                              Emptyset_iterator());

      // try to compactify the selection region by also selecting all the faces included
      // in the bounding box of the initial selection
      std::vector<halfedge_descriptor> stack_for_expension;
      Bbox_3 bb;
      BOOST_FOREACH(face_descriptor fd, cc_faces)
      {
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(fd, tm), tm))
        {
          bb += get(vpmap, target(h, tm)).bbox();
          face_descriptor nf = face(opposite(h, tm), tm);
          if (nf != boost::graph_traits<TriangleMesh>::null_face() &&
              cc_faces.count(nf)==0)
          {
            stack_for_expension.push_back(opposite(h, tm));
          }
        }
      }

      while(!stack_for_expension.empty())
      {
        halfedge_descriptor h=stack_for_expension.back();
        stack_for_expension.pop_back();
        if ( cc_faces.count(face(h,tm))==1) continue;
        if ( do_overlap(bb, get(vpmap, target(next(h, tm), tm)).bbox()) )
        {
          cc_faces.insert(face(h,tm));
          halfedge_descriptor candidate = opposite(next(h, tm), tm);
          if ( face(candidate, tm) != boost::graph_traits<TriangleMesh>::null_face() )
            stack_for_expension.push_back( candidate );
          candidate = opposite(prev(h, tm), tm);
          if ( face(candidate, tm) != boost::graph_traits<TriangleMesh>::null_face() )
            stack_for_expension.push_back( candidate );
        }
      }

      // remove faces from the set to process
      BOOST_FOREACH(face_descriptor f, cc_faces)
        faces_to_remove.erase(f);

      if (cc_faces.size()==1) continue; // it is a triangle nothing better can be done

      //Check for non-manifold vertices in the selection and remove them by selecting all incident faces:
      //  extract the set of halfedges that is on the boundary of the holes to be
      //  made. In addition, we make sure no hole to be created contains a vertex
      //  visited more than once along a hole border (pinched surface)
      //  We save the size of boundary_hedges to make sur halfedges added
      //  from non_filled_hole are not removed.
      do{
        bool non_manifold_vertex_removed = false; //here non-manifold is for the 1D polyline
        std::vector<halfedge_descriptor> boundary_hedges;
        BOOST_FOREACH(face_descriptor fh, cc_faces)
        {
          halfedge_descriptor h = halfedge(fh,tm);
          for (int i=0;i<3; ++i)
          {
            if ( is_border( opposite(h, tm), tm) ||
                 cc_faces.count( face( opposite(h, tm), tm) ) == 0)
            {
              boundary_hedges.push_back(h);
            }
            h=next(h, tm);
          }
        }

        // detect vertices visited more than once along
        // a hole border. We then remove all faces incident
        // to such a vertex to force the removal of the vertex.
        // Actually even if two holes are sharing a vertex, this
        // vertex will be removed. It is not needed but since
        // we do not yet have one halfedge per hole it is simpler
        // and does not harm
        std::set<vertex_descriptor> border_vertices;
        BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
        {
          if (!border_vertices.insert(target(h,tm)).second){
            BOOST_FOREACH(halfedge_descriptor hh, halfedges_around_target(h,tm)){
              if (!is_border(hh, tm))
              {
                cc_faces.insert(face(hh, tm)); // add the face to the current selection
                faces_to_remove.erase(face(hh, tm));
              }
            }
            non_manifold_vertex_removed=true;
          }
        }

        if (!non_manifold_vertex_removed)
        {
          break;
        }
      }
      while(true);

      // Collect halfedges on the boundary of the region to be selected
      // (pointing inside the domain to be remeshed)
      std::vector<halfedge_descriptor> cc_border_hedges;
      BOOST_FOREACH(face_descriptor fd, cc_faces)
      {
        halfedge_descriptor h = halfedge(fd, tm);
        for (int i=0; i<3;++i)
        {
          if ( is_border(opposite(h, tm), tm) ||
               cc_faces.count( face(opposite(h, tm), tm) )== 0)
          {
            cc_border_hedges.push_back(h);
          }
          h=next(h, tm);
        }
      }

      if(!is_selection_a_topological_disk(cc_faces, tm))
      {
        // check if the selection contains cycles of border halfedges
        bool only_border_edges = true;
        std::set<halfedge_descriptor> mesh_border_hedge;

        BOOST_FOREACH(halfedge_descriptor h, cc_border_hedges)
        {
          if ( !is_border(opposite(h, tm), tm) )
            only_border_edges = false;
          else
            mesh_border_hedge.insert( opposite(h, tm) );
        }
        int nb_cycles=0;
        while(!mesh_border_hedge.empty())
        {
          // we must count the number of cycle of boundary edges
          halfedge_descriptor h_b = *mesh_border_hedge.begin(), h=h_b;
          mesh_border_hedge.erase( mesh_border_hedge.begin() );
          do{
            h=next(h, tm);
            if (h==h_b)
            {
              // found a cycle
              ++nb_cycles;
              break;
            }
            else
            {
              typename std::set<halfedge_descriptor>::iterator it =
                mesh_border_hedge.find(h);
              if ( it == mesh_border_hedge.end() )
                break; // not a cycle
              mesh_border_hedge.erase(it);
            }
          }while(true);
        }

        if(nb_cycles > (only_border_edges ? 1 : 0) )
        {
          if(verbose)
            std::cout << "  DEBUG: CC not handled due to the presence of  "
                      << nb_cycles << " of boundary edges\n";
          topology_issue = true;
          continue;
        }
        else
        {
          if (preserve_genus)
          {
            if(verbose)
              std::cout << "  DEBUG: CC not handled because it is not a topological disk (preserve_genus=true)\n";
            all_fixed = false;
            continue;
          }
          // count the number of cycles of halfedges of the boundary
          std::map<vertex_descriptor, vertex_descriptor> bhs;
          BOOST_FOREACH(halfedge_descriptor h, cc_border_hedges)
          {
            bhs[source(h, tm)]=target(h, tm);
          }
          int nbc=0;
          while(!bhs.empty())
          {
            ++nbc;
            std::pair<vertex_descriptor, vertex_descriptor > top=*bhs.begin();
            bhs.erase(bhs.begin());
            do
            {
              typename std::map<vertex_descriptor, vertex_descriptor>::iterator
                it_find = bhs.find(top.second);
              if (it_find == bhs.end()) break;
              top = *it_find;
              bhs.erase(it_find);
            }
            while(true);
          }
          if (nbc!=1){
            if(verbose)
              std::cout << "  DEBUG: CC not handled because it is not a topological disk("
                        << nbc << " boundary cycles)<<\n";
            all_fixed = false;
            continue;
          }
          else
          {
            if(verbose)
              std::cout << "  DEBUG: CC that is not a topological disk but has only one boundary cycle(preserve_genus=false)\n";
          }
        }
      }

      // sort halfedges so that they describe the sequence
      // of halfedges of the hole to be made
      CGAL_assertion( cc_border_hedges.size() > 2 );
      for(std::size_t i=0; i < cc_border_hedges.size()-2; ++i)
      {
        vertex_descriptor tgt = target(cc_border_hedges[i], tm);
        for(std::size_t j=i+1; j<cc_border_hedges.size(); ++j)
        {
          if(tgt == source(cc_border_hedges[j], tm))
          {
            std::swap(cc_border_hedges[i+1], cc_border_hedges[j]);
            break;
          }
          CGAL_assertion(j!=cc_border_hedges.size()-1);
        }
      }
      CGAL_assertion( source(cc_border_hedges.front(), tm) ==
                      target(cc_border_hedges.back(), tm) );

      // collect vertices and edges inside the current selection cc
      std::set<vertex_descriptor> cc_interior_vertices;
      std::set<edge_descriptor>  cc_interior_edges;

      // first collect all vertices and edges incident to the faces to remove
      BOOST_FOREACH(face_descriptor fh, cc_faces)
      {
        BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(fh,tm),tm))
        {
          if (halfedge(target(h, tm), tm)==h) // limit the number of insertions
            cc_interior_vertices.insert(target(h, tm));
          cc_interior_edges.insert(edge(h,tm));
        }
      }
      // and then remove those on the boundary
      BOOST_FOREACH(halfedge_descriptor h, cc_border_hedges)
      {
        cc_interior_vertices.erase(target(h, tm));
        cc_interior_edges.erase(edge(h,tm));
      }

      if (verbose)
      {
        std::cout << "  DEBUG: is_valid(tm) in one_step, before mesh changes? ";
        std::cout << is_valid(tm) << std::endl;
      }

      //try hole_filling.
      typedef CGAL::Triple<int, int, int> Face_indices;
      typedef typename boost::property_traits<VertexPointMap>::value_type Point;
      std::vector<Point> hole_points, third_points;
      hole_points.reserve(cc_border_hedges.size());
      third_points.reserve(cc_border_hedges.size());
      std::vector<vertex_descriptor> border_vertices;
      BOOST_FOREACH(halfedge_descriptor h, cc_border_hedges)
      {
        vertex_descriptor v = source(h, tm);
        hole_points.push_back( get(vpmap, v) );
        border_vertices.push_back(v);
        third_points.push_back(get(vpmap, target(next(opposite(h, tm), tm), tm))); // TODO fix me for mesh border edges
      }
      CGAL_assertion(hole_points.size() >= 3);

      // try to triangulate the hole using default parameters
      //(using Delaunay search space if CGAL_HOLE_FILLING_DO_NOT_USE_DT3 is not defined)
      std::vector<Face_indices> patch;
      if (hole_points.size()>3)
        triangulate_hole_polyline(hole_points,
                                  third_points,
                                  std::back_inserter(patch));
      else
        patch.push_back(Face_indices(0,1,2)); // trivial hole filling

      if(patch.empty())
      {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
        if (verbose)
          std::cout << "  DEBUG: Failed to fill a hole using Delaunay search space.\n";
        triangulate_hole_polyline(hole_points,
                                  third_points,
                                  std::back_inserter(patch),
                                  parameters::use_delaunay_triangulation(false));
#endif
        if (patch.empty())
        {
          if (verbose)
            std::cout << "  DEBUG: Failed to fill a hole using the whole search space.\n";
          all_fixed = false;
          continue;
        }
      }

      // make sure that the hole filling is valid, we check that no
      // edge already in the mesh is present in patch.
      bool non_manifold_edge_found = false;
      BOOST_FOREACH(const Face_indices& triangle, patch)
      {
        cpp11::array<int, 6> edges =
          make_array(triangle.first, triangle.second,
                     triangle.second, triangle.third,
                     triangle.third, triangle.first);
        for (int k=0; k<3; ++k)
        {
          int vi=edges[2*k], vj=edges[2*k+1];
          // ignore boundary edges
          if (vi+1==vj || (vj==0 && static_cast<std::size_t>(vi)==border_vertices.size()-1) )
            continue;
          halfedge_descriptor h = halfedge(border_vertices[vi], border_vertices[vj], tm).first;
          if (h!=boost::graph_traits<TriangleMesh>::null_halfedge() &&
              cc_interior_edges.count(edge(h, tm))==0)
          {
            non_manifold_edge_found=true;
            break;
          }
        }
        if (non_manifold_edge_found) break;
      }
      if (non_manifold_edge_found)
      {
        if (verbose)
          std::cout << "  DEBUG: Triangulation produced is non-manifold when plugged into the mesh.\n";
        all_fixed = false;
        continue;
      }

      // plug the new triangles in the mesh, reusing previous edges and faces
      std::vector<edge_descriptor> edge_stack(cc_interior_edges.begin(), cc_interior_edges.end());
      std::vector<face_descriptor> face_stack(cc_faces.begin(), cc_faces.end());

      std::map< std::pair<int, int>, halfedge_descriptor > halfedge_map;
      int i=0;
      // register border halfedges
      BOOST_FOREACH(halfedge_descriptor h, cc_border_hedges)
      {
        int j = static_cast<int>( std::size_t(i+1)%cc_border_hedges.size() );
        halfedge_map.insert(std::make_pair( std::make_pair(i, j), h) );
        set_halfedge(target(h, tm), h, tm); // update vertex halfedge pointer
        CGAL_assertion( border_vertices[i] == source(h, tm) &&
                        border_vertices[j] == target(h, tm) );
        ++i;
      }

      std::vector<halfedge_descriptor> hedges;
      hedges.reserve(4);
      face_descriptor f = boost::graph_traits<TriangleMesh>::null_face();
      BOOST_FOREACH(const Face_indices& triangle, patch)
      {
        // get the new face
        if (face_stack.empty())
          f=add_face(tm);
        else
        {
          f=face_stack.back();
          face_stack.pop_back();
        }

        cpp11::array<int, 4> indices =
          make_array( triangle.first,
                      triangle.second,
                      triangle.third,
                      triangle.first );
        for (int i=0; i<3; ++i)
        {
          // get the corresponding halfedge (either a new one or an already created)
          typename std::map< std::pair<int, int> , halfedge_descriptor >::iterator insert_res =
            halfedge_map.insert(
              std::make_pair( std::make_pair(indices[i], indices[i+1]),
                              boost::graph_traits<TriangleMesh>::null_halfedge() ) ).first;
          if (insert_res->second == boost::graph_traits<TriangleMesh>::null_halfedge())
          {
            if (edge_stack.empty())
              insert_res->second = halfedge(add_edge(tm), tm);
            else
            {
              insert_res->second = halfedge(edge_stack.back(), tm);
              edge_stack.pop_back();
            }

            halfedge_map[std::make_pair(indices[i+1], indices[i])] =
              opposite(insert_res->second, tm);
          }
          hedges.push_back(insert_res->second);
        }
        hedges.push_back(hedges.front());
        // update halfedge connections + face pointers
        for(int i=0; i<3;++i)
        {
          set_next(hedges[i], hedges[i+1], tm);
          set_face(hedges[i], f, tm);
          set_target(hedges[i], border_vertices[indices[i+1]], tm);
        }
        set_halfedge(f, hedges[0], tm);
        hedges.clear();
      }

      // now remove remaining edges,
      BOOST_FOREACH(edge_descriptor e, edge_stack)
        remove_edge(e, tm);
      // vertices,
      BOOST_FOREACH(vertex_descriptor vh, cc_interior_vertices)
        remove_vertex(vh, tm);
      // and remaning faces
      BOOST_FOREACH(face_descriptor f, face_stack)
        remove_face(f, tm);

      if (verbose)
        std::cout << "  DEBUG: " << cc_faces.size() << " triangles removed, "
                  << patch.size() << " created\n";

      CGAL_assertion(is_valid(tm));

      something_was_done = true;
    }
  }
  if (!something_was_done)
  {
    faces_to_remove.swap(faces_to_remove_copy);
    if (verbose)
      std::cout<<"  DEBUG: Nothing was changed during this step, self-intersections won't be recomputed."<<std::endl;
  }
  return std::make_pair(all_fixed, topology_issue);
}

template <class TriangleMesh, class NamedParameters>
bool remove_self_intersections(TriangleMesh& tm, const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::face_descriptor face_descriptor;

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                           get_property_map(vertex_point, tm));

  const int max_steps = boost::choose_param(boost::get_param(np, internal_np::number_of_iterations), 7);
  bool verbose = boost::choose_param(boost::get_param(np, internal_np::verbosity_level), 0) > 0;
  bool preserve_genus = boost::choose_param(boost::get_param(np, internal_np::preserve_genus), true);

  if (verbose)
    std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tm)? " << is_valid(tm) << "\n";

  // first handle the removal of degenerate faces
  remove_degenerate_faces(tm, np);

  if (verbose)
    std::cout << "DEBUG: After degenerate faces removal, is_valid(tm)? " << is_valid(tm) << "\n";

  // Look for self-intersections in the polyhedron and remove them
  int step=-1;
  bool all_fixed = true; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> faces_to_remove;
  while( ++step<max_steps )
  {
    if (faces_to_remove.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;
      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(tm, std::back_inserter(self_inter));

      BOOST_FOREACH(Face_pair fp, self_inter)
      {
        faces_to_remove.insert(fp.first);
        faces_to_remove.insert(fp.second);
      }
    }

    if ( faces_to_remove.empty() && all_fixed){
      if (verbose)
        std::cout<<"DEBUG: There is no more face to remove."<<std::endl;
      break;
    }

    cpp11::tie(all_fixed, topology_issue) =
      remove_self_intersections_one_step(tm, faces_to_remove, vpm, step, preserve_genus, verbose);
    if (all_fixed && topology_issue)
    {
      if (verbose)
        std::cout<< "DEBUG: Process stopped because of boundary cycles"
                    " of boundary edges involved in self-intersections.\n";
      return false;
    }
  }

  return step<max_steps;
}

template <class TriangleMesh>
bool remove_self_intersections(TriangleMesh& tm)
{
  return remove_self_intersections(tm, parameters::all_default());
}
/// \endcond

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
