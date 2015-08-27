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

#include <set>
#include <vector>
#include <boost/algorithm/minmax_element.hpp>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Union_find.h>

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

template <class Traits, class TriangleMesh, class VertexPointMap>
bool is_degenerated(
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
  TriangleMesh& tmesh,
  const VertexPointMap& vpmap,
  const Traits& traits)
{
  const typename Traits::Point_3& p1 = get(vpmap, target( hd, tmesh) );
  const typename Traits::Point_3& p2 = get(vpmap, target(next(hd, tmesh), tmesh) );
  const typename Traits::Point_3& p3 = get(vpmap, source( hd, tmesh) );
  return traits.collinear_3_object()(p1, p2, p3);
}

template <class Traits, class TriangleMesh, class VertexPointMap>
bool is_degenerated(
  typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
  TriangleMesh& tmesh,
  const VertexPointMap& vpmap,
  const Traits& traits)
{
  return is_degenerated(halfedge(fd,tmesh), tmesh, vpmap, traits);
}

///\cond SKIP_IN_MANUAL
template <class EdgeRange, class TriangleMesh, class NamedParameters>
std::size_t remove_null_edges(
                       const EdgeRange& edge_range,
                       TriangleMesh& tmesh,
                       const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tmesh));

  using boost::choose_const_pmap;
  using boost::get_param;
  using boost::choose_param;

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  typedef typename GetVertexPointMap<TM, NamedParameters>::type VertexPointMap;
  VertexPointMap vpmap = choose_pmap(get_param(np, boost::vertex_point),
                                     tmesh,
                                     boost::vertex_point);
  typedef typename GetGeomTraits<TM, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, geom_traits), Traits());

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

/// \ingroup PkgPolygonMeshProcessing
/// removes the degenerate faces from a triangulated surface mesh.
/// A face is considered degenerate if two of its vertices share the same location,
/// or more generally if all its vertices are collinear.
///
/// @pre `CGAL::is_triangle_mesh(tmesh)`
///
/// @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
///        that has an internal property map for `boost::vertex_point_t`
/// @tparam NamedParameters a sequence of \ref namedparameters
///
/// @param tmesh the  triangulated surface mesh to be repaired
/// @param np optional \ref namedparameters described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`. The type of this map is model of `ReadWritePropertyMap` \cgalParamEnd
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
/// \endcond
template <class TriangleMesh, class NamedParameters>
std::size_t remove_degenerate_faces(TriangleMesh& tmesh,
                                    const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tmesh));

  using boost::choose_const_pmap;
  using boost::get_param;
  using boost::choose_param;

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  typedef typename GetVertexPointMap<TM, NamedParameters>::type VertexPointMap;
  VertexPointMap vpmap = choose_pmap(get_param(np, boost::vertex_point),
                                     tmesh,
                                     boost::vertex_point);
  typedef typename GetGeomTraits<TM, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, geom_traits), Traits());

// First remove edges of length 0
  std::size_t nb_deg_faces = remove_null_edges(edges(tmesh), tmesh, np);

// Then, remove triangles made of 3 collinear points
  std::set<face_descriptor> degenerate_face_set;
  BOOST_FOREACH(face_descriptor fd, faces(tmesh))
    if ( is_degenerated(fd, tmesh, vpmap, traits) )
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
               !degenerate_face_set.count(adjacent_face) )
            boundary_hedges.push_back(hd);
          else
            if (cc_faces.insert(adjacent_face).second)
            {
              inside_hedges.push_back(hd);
              queue.push_back(adjacent_face);
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
        degenerate_face_set.erase(f);
      BOOST_FOREACH(face_descriptor f, cc_faces)
        remove_face(f, tmesh);
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
            set_halfedge(*ref_vertices.first, opposite( prev(side_two[hi], tmesh), tmesh), tmesh );
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

/// \cond SKIP_IN_MANUAL
template<class TriangleMesh>
std::size_t remove_degenerate_faces(TriangleMesh& tmesh)
{
  return remove_degenerate_faces(tmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}
/// \endcond

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
