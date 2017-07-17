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

// headers for self-intersection removal
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/boost/graph/selection.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

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

template <class Traits>
struct Less_along_ray{
  const Traits& m_traits;
  typename Traits::Point_3 m_source;
  Less_along_ray(const Traits& traits,
                 const typename Traits::Point_3& s)
    : m_traits(traits)
    , m_source(s)
  {};
  bool operator()( const typename Traits::Point_3& p1,
                   const typename Traits::Point_3& p2) const
  {
    return m_traits.collinear_are_ordered_along_line_3_object()(m_source, p1, p2);
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
/// @tparam NamedParameters a sequence of \ref namedparameters
///
/// @param tmesh the  triangulated surface mesh to be repaired
/// @param np optional \ref namedparameters described below
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

// First remove edges of length 0
  std::size_t nb_deg_faces = remove_null_edges(edges(tmesh), tmesh, np);

// Then, remove triangles made of 3 collinear points
  std::set<face_descriptor> degenerate_face_set;
  BOOST_FOREACH(face_descriptor fd, faces(tmesh))
    if ( is_degenerate_triangle_face(fd, tmesh, vpmap, traits) )
      degenerate_face_set.insert(fd);
  nb_deg_faces+=degenerate_face_set.size();

  while (!degenerate_face_set.empty())
  {
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
      degenerate_face_set.erase(degenerate_face_set.begin());
    // flip the longest edge of the triangle
      const typename Traits::Point_3& p1 = get(vpmap, target( halfedge(fd, tmesh), tmesh) );
      const typename Traits::Point_3& p2 = get(vpmap, target(next(halfedge(fd, tmesh), tmesh), tmesh) );
      const typename Traits::Point_3& p3 = get(vpmap, source( halfedge(fd, tmesh), tmesh) );

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
        Euler::flip_edge(edge_to_flip, tmesh);
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

      #if 0
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
      output << std::setprecision(44);
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
      typedef Less_along_ray<Traits> Less_point;
      typedef std::set<typename Traits::Point_3, Less_point> Sorted_point_set;
      Sorted_point_set sorted_points( Less_point( traits, get(vpmap, *ref_vertices.first) ) );
      BOOST_FOREACH(vertex_descriptor v, boundary_vertices)
        sorted_points.insert( get(vpmap,v) );

      CGAL_assertion( get( vpmap, *ref_vertices.first)==*sorted_points.begin() );
      CGAL_assertion( get( vpmap, *ref_vertices.second)==*cpp11::prev(sorted_points.end()) );

      // recover halfedges on the hole, bounded by the reference vertices
      std::vector<halfedge_descriptor> side_one, side_two;
      side_one.push_back( next( halfedge(*ref_vertices.first, tmesh), tmesh) );
      while( target(side_one.back(), tmesh)!=*ref_vertices.second)
        side_one.push_back( next(side_one.back(), tmesh) );
      side_two.push_back( next(side_one.back(), tmesh) );
      while( target(side_two.back(), tmesh)!=*ref_vertices.first )
        side_two.push_back( next(side_two.back(), tmesh) );
      // reverse the order of the second side so as to follow
      // the same order than side one
      std::reverse(side_two.begin(), side_two.end());
      BOOST_FOREACH(halfedge_descriptor& h, side_two)
        h=opposite(h, tmesh);

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
namespace internal{
template <class Descriptor>
struct Is_selected{
  std::set<Descriptor>& selection;

  Is_selected(std::set<Descriptor>& sel)
    :selection(sel)
  {}

  friend bool get(Is_selected is, Descriptor d){
    return is.selection.count(d);
  }

  friend void put(Is_selected is, Descriptor d, bool b){
    if (b)
      is.selection.insert(d);
    else
      is.selection.erase(d);
  }
};
} // end of namespace internal

template <class TriangleMesh>
bool remove_self_intersections(TriangleMesh& tm, const int max_steps = 7, bool verbose=false)
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;

// Look for self-intersections in the polyhedron and remove them
  int step=-1;
  std::vector<halfedge_descriptor> non_filled_hole;
  bool no_hole_was_filled=false; // indicates if the filling of all previously
                                 // created holes failed. If true then no new
                                 // self-intersection have been created and
                                 // checking for it is useless.
  while( ++step<max_steps )
  {
    if (verbose)
      std::cout << "DEBUG: is_valid(tm)? " << is_valid(tm) << "\n";

    typedef std::pair<face_descriptor, face_descriptor> Face_pair;
    std::vector<Face_pair> self_inter;
    std::vector<halfedge_descriptor> one_halfedge_per_border;
    if (!no_hole_was_filled)
      self_intersections(tm, std::back_inserter(self_inter));
    no_hole_was_filled=true;

    if(!self_inter.empty() || !non_filled_hole.empty()){
      if (verbose)
        std::cout << "DEBUG: Iterative self-intersection removal step " << step
                  << " - non_filled_hole.size() = " << non_filled_hole.size() << std::endl;
      std::set<face_descriptor> faces_to_remove;
      BOOST_FOREACH(Face_pair fp, self_inter)
      {
        faces_to_remove.insert(fp.first);
        faces_to_remove.insert(fp.second);
      }

      // expand the region to be filled
      internal::Is_selected<face_descriptor> is_selected(faces_to_remove);
      expand_face_selection(faces_to_remove, tm, step+1, is_selected, Emptyset_iterator());
      // try to avoid non-manifold vertices (morpho-math)
      reduce_face_selection(faces_to_remove, tm, 1, is_selected, Emptyset_iterator());

      // now expand holes than were not filled
      std::vector<halfedge_descriptor> boundary_hedges; // this container will contain the halfedges of
                                                       // all created holes
      std::set<halfedge_descriptor> border_created; // track border halfedges that were previously created
                                                    // to avoid considering them as original mesh border
                                                    // edges that should be kept
      BOOST_FOREACH(halfedge_descriptor h, non_filled_hole)
      {
        select_incident_faces(halfedges_around_face(h,tm), tm,
          std::inserter(faces_to_remove, faces_to_remove.begin()) );
        BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_face(h,tm))
        {
          CGAL_assertion(is_border(h2, tm));
          border_created.insert(h2);
          if ( is_border(opposite(h2, tm), tm) )
            boundary_hedges.push_back(h2); // mesh border halfedges should be re-added

        }
      }
      non_filled_hole.clear();

      /// save the halfedges that will get on the boundary
      bool border_edges_found=false; // indicates at least one face incident to
                                     // the mesh border will be removed. In that
                                     // case, border halfedges should not be removed.

      // extract the set of halfedges that is on the boundary of the holes to be
      // made. In addition, we make sure no hole to be created contains a vertex
      // visited more than once along a hole border (pinched surface)
      //  We save the size of boundary_hedges to make sur halfedges added
      // from non_filled_hole are not removed.
      bool non_manifold_vertex_removed; //here non-manifold is for the 1D polyline
      std::size_t boundary_hedges_initial_size=boundary_hedges.size();
      do{
        non_manifold_vertex_removed=false;
        boundary_hedges.resize(boundary_hedges_initial_size);
        BOOST_FOREACH(face_descriptor fh, faces_to_remove)
        {
          halfedge_descriptor h = halfedge(fh,tm);
          for (int i=0;i<3; ++i)
          {
            if ( is_border( opposite(h, tm), tm) ){
              if (!border_created.count(opposite(h, tm))){
                // only border halfedges that were not created by a previous face
                // removal should be considered as hole boundary
                boundary_hedges.push_back(h);
                border_edges_found=true;
              }
            }
            else
              if ( !faces_to_remove.count( face( opposite(h, tm), tm) ) )
                boundary_hedges.push_back(h);
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
                faces_to_remove.insert(face(hh, tm));
            }
            non_manifold_vertex_removed=true;
          }
        }
      }
      while(non_manifold_vertex_removed);

    /// remove the selection
      if (border_edges_found){
        // When at least one face incident to the mesh border is set
        // to be removed, we should pay attention not to remove
        // the border halfedge that is part of the mesh.

        // first collect all vertices and edges incident to the faces to remove
        std::set<vertex_descriptor> vertices_to_remove;
        std::set<edge_descriptor>  edges_to_remove;
        BOOST_FOREACH(face_descriptor fh, faces_to_remove)
        {
          BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(halfedge(fh,tm),tm))
          {
            if (halfedge(target(h, tm), tm)==h) // limit the number of insertions
              vertices_to_remove.insert(target(h, tm));
            edges_to_remove.insert(edge(h,tm));
          }
        }

        // detect cycles of input border halfedges:
        // if such a cycle is found input border halfedges are removed to
        // remove small border cycle that are imposing a self-intersection
        // that could not be fixed if kept. (This also prevent small island
        // if a small hole is incident to the faces to be removed).
        if (border_edges_found){
          std::set<halfedge_descriptor> cycles;
          std::set<halfedge_descriptor>  boundary_set;
          BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
          {
            if ( is_border(opposite(h, tm), tm) )
              boundary_set.insert(opposite(h, tm));
          }

          BOOST_FOREACH(halfedge_descriptor hd, boundary_set)
          {
            CGAL_assertion(is_border(hd,tm));
            if(cycles.count(hd)) continue;
            halfedge_descriptor nhd=next(hd,tm);
            bool remove_it=true;
            do{
              if (!boundary_set.count(nhd))
              {
                remove_it=false;
                break;
              }
              nhd=next(nhd,tm);
            }
            while(nhd!=hd);
            if (remove_it){
              BOOST_FOREACH(halfedge_descriptor h, halfedges_around_face(hd,tm))
                cycles.insert(h);
            }
          }

          //remove cycle edges
          if (!cycles.empty())
          {
            std::vector<halfedge_descriptor> tmp;
            tmp.reserve(boundary_hedges.size()-cycles.size());
            BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
            {
              if (!cycles.count(opposite(h, tm)))
                tmp.push_back(h);
            }
            tmp.swap(boundary_hedges);
          }
        }

        // do not remove edges on the boundary of the selection of faces,
        // nor its vertices
        BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
        {
          vertices_to_remove.erase(target(h, tm));
          edges_to_remove.erase(edge(h,tm));
        }
        // now remove edges,
        BOOST_FOREACH(edge_descriptor e, edges_to_remove)
          remove_edge(e, tm);
        // and vertices,
        BOOST_FOREACH(vertex_descriptor vh, vertices_to_remove)
          remove_vertex(vh, tm);
        // and finally facets
        BOOST_FOREACH(face_descriptor f, faces_to_remove)
          remove_face(f, tm);
        // set new border_vertices to the boundary and update
        // the halfedge pointer of the border vertices
        BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
        {
          set_face(h,graph_traits::null_face(), tm);
          set_halfedge(target(h, tm), h, tm);
          set_next(h,h,tm); // set himself as next to track edges of the holes
        }
        // update next/prev relationships of the hole
        BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
        {
          halfedge_descriptor nh=next(opposite(h, tm), tm);
          while( !is_border(opposite(nh, tm), tm) ||
                  next(opposite(nh, tm), tm) != opposite(nh, tm)) //this part makes sure we consider halfedges of the hole
          {
            nh=next(opposite(nh, tm), tm);
            CGAL_assertion(nh!=h);
          }
          CGAL_assertion(next(opposite(nh, tm), tm)==opposite(nh,tm));
          set_next(opposite(nh, tm), h, tm);
        }
      }
      else
        /// \todo check whether this is more expensive than the previous code above
        BOOST_FOREACH(face_descriptor f, faces_to_remove)
          Euler::remove_face(halfedge(f, tm), tm);

      if (verbose)
        std::cout << "  DEBUG: " << faces_to_remove.size() << " triangles removed" << std::endl;

      /// now get one halfedge per hole
      std::set<halfedge_descriptor> visited;
      BOOST_FOREACH(halfedge_descriptor h, boundary_hedges)
      {
        if (visited.insert(h).second)
        {
          one_halfedge_per_border.push_back(h);
          BOOST_FOREACH(halfedge_descriptor hh, halfedges_around_face(h, tm))
          {
            CGAL_assertion_code(bool insert_ok =)
            visited.insert(hh)
            CGAL_assertion_code(.second);
            CGAL_assertion(insert_ok || h==hh);
          }
        }
      }
    }

    if (!one_halfedge_per_border.empty()){
      BOOST_FOREACH(halfedge_descriptor h, one_halfedge_per_border)
      {
        std::size_t nb_new_triangles = 0;
        Counting_output_iterator out(&nb_new_triangles);
        triangulate_hole(tm, h, out);
        if (!nb_new_triangles)
        {
          if (verbose)
            std::cout << "  DEBUG: Failed to fill a hole!!!" << std::endl;
          non_filled_hole.push_back(h);
        }
        else
          no_hole_was_filled=false;
      }
      if (verbose)
        std::cout << "  DEBUG: Number of holes " << one_halfedge_per_border.size() << std::endl;
    }
    else{
      if (verbose)
        std::cout << "INFO: All self-intersections were corrected\n";
      break;
    }
  }

  return step<max_steps;
}
/// \endcond

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
