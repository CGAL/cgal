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

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL{
namespace Polygon_mesh_processing {

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

namespace internal {

  // hbase is an edge of length 0 we cannot contract because
  // the link condition is not satisfied.
  // In this function we look whether the condition is not
  // satisfied because of some faces on the "side" of hbase.
  // Here we look for two halfedges that together with hbase
  // enclose a set of degenerate faces so as to replace that
  // set with only one triangle
  template <class Traits, class TriangleMesh, class VertexPointMap>
  boost::optional<
    std::pair<
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor,
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor >
  >
  find_larger_triangle(
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hbase,
    TriangleMesh& tmesh,
    const VertexPointMap& vpmap,
    const Traits& traits)
  {
    typedef typename boost::graph_traits<TriangleMesh> GT;
    typedef typename GT::halfedge_descriptor halfedge_descriptor;

    bool no_pb=false;
    // hbase = v0 -> v1
    // consider hedges with v0 as target
    halfedge_descriptor left_hedge = prev(hbase, tmesh), o_right_hedge=left_hedge;
    do
    {
      left_hedge = prev( opposite(left_hedge, tmesh), tmesh );
      if ( face(left_hedge, tmesh) == GT::null_face() ||
           !is_degenerated(left_hedge, tmesh, vpmap, traits) )
      {
        no_pb=true;
        break;
      }
      // look for a halfedge with v1 as target
      BOOST_FOREACH(o_right_hedge, halfedges_around_source(left_hedge, tmesh) )
      {
        if ( target(o_right_hedge, tmesh) == target(hbase, tmesh) )
        {
          if ( !is_degenerated( opposite(o_right_hedge, tmesh), tmesh, vpmap, traits) )
            no_pb=true;
          break;
        }
      }
      if (  target(o_right_hedge, tmesh) == target(hbase, tmesh) )
        break;
    }
    while(true);

    if (!no_pb)
      return std::make_pair(left_hedge, opposite(o_right_hedge, tmesh));
    return boost::none;
  }

  // h2 is a halfedge of length 0 but collapsing the edge
  // is not directly possible because the link condition is not
  // satisfied. We remove all simplices bounded by h1, h2 and h3
  // so as to keep only that triangle
  template <class TriangleMesh, class EdgeMap>
  void
  remove_faces_inside_triangle(
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h1,
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h2,
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h3,
    TriangleMesh& tmesh,
    EdgeMap& edge_map1,
    EdgeMap& edge_map2)
  {
    typedef typename boost::graph_traits<TriangleMesh> GT;
    typedef typename GT::edge_descriptor edge_descriptor;
    typedef typename GT::halfedge_descriptor halfedge_descriptor;
    typedef typename GT::face_descriptor face_descriptor;
    typedef typename GT::vertex_descriptor vertex_descriptor;

    CGAL_assertion( source(h1,tmesh) == target(h3, tmesh) );
    CGAL_assertion( source(h2,tmesh) == target(h1, tmesh) );
    CGAL_assertion( source(h3,tmesh) == target(h2, tmesh) );

    std::set<face_descriptor> all_faces;
    all_faces.insert(face(h1, tmesh));
    all_faces.insert(face(h2, tmesh));
    all_faces.insert(face(h3, tmesh));

    std::vector<halfedge_descriptor> queue;
    queue.push_back( opposite(prev(h1,tmesh), tmesh) );
    queue.push_back( opposite(next(h1,tmesh), tmesh) );
    queue.push_back( opposite(prev(h2,tmesh), tmesh) );
    queue.push_back( opposite(next(h2,tmesh), tmesh) );
    queue.push_back( opposite(prev(h3,tmesh), tmesh) );
    queue.push_back( opposite(next(h3,tmesh), tmesh) );

    std::set<edge_descriptor> all_edges;
    while(!queue.empty())
    {
      halfedge_descriptor back=queue.back();
      queue.pop_back();
      all_edges.insert( edge(back, tmesh) );

      if ( all_faces.insert( face(back, tmesh) ).second )
      {
        queue.push_back( opposite(prev(back,tmesh), tmesh) );
        queue.push_back( opposite(next(back,tmesh), tmesh) );
      }
    }

    std::set<vertex_descriptor> all_vertices;
    BOOST_FOREACH(edge_descriptor ed, all_edges)
    {
      all_vertices.insert( source(ed, tmesh) );
      all_vertices.insert( target(ed, tmesh) );
    }
    all_vertices.erase( target(h1, tmesh) );
    all_vertices.erase( target(h2, tmesh) );
    all_vertices.erase( target(h3, tmesh) );

    // create the triangle
    face_descriptor remaining_face = *all_faces.begin();
    all_faces.erase(all_faces.begin());
    //   update next-prev pointers
    set_next(h1, h2, tmesh);
    set_next(h2, h3, tmesh);
    set_next(h3, h1, tmesh);
    //   update vertex pointers
    set_halfedge(target(h1,tmesh), h1, tmesh);
    set_halfedge(target(h2,tmesh), h2, tmesh);
    set_halfedge(target(h3,tmesh), h3, tmesh);
    //  update face-halfedge pointers
    set_face(h1, remaining_face, tmesh);
    set_face(h2, remaining_face, tmesh);
    set_face(h3, remaining_face, tmesh);
    set_halfedge(remaining_face, h2, tmesh);

    // remove interior simplices
    BOOST_FOREACH(edge_descriptor ed, all_edges){
      edge_map1.erase(ed);
      edge_map2.erase(ed);
      remove_edge(ed, tmesh);
    }
    BOOST_FOREACH(vertex_descriptor vd, all_vertices)
      remove_vertex(vd, tmesh);
    BOOST_FOREACH(face_descriptor fd, all_faces)
      remove_face(fd, tmesh);
  }
} // end of namespace internal

/// \ingroup PkgPolygonMeshProcessing
/// removes the degenerate faces from a triangle mesh.
///
/// @pre `CGAL::is_pure_triangle(tmesh)`
///
/// @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
///        that has a property map for `boost::vertex_point_t`
/// @tparam NamedParameters a sequence of \ref namedparameters
///
/// @param tmesh the triangle mesh to be repaired
/// @param np optional \ref namedparameters described below
///
/// \b Named \b parameters
/// <ul>
/// <li>\b vertex_point_map the property map with the points associated to the vertices of `pmesh`. The type of this mad is model of `ReadWritePropertyMap`.
/// <li>\b kernel a geometric traits class instance. 
/// The traits class must provide the nested types :
///     - `Point_3`,
///     - `Compare_distance_3` to compute the distance between 2 points
///     - `Collinear_are_ordered_along_line_3` to check whether 3 collinear points are ordered
///     - `Collinear_3` to check whether 3 points are collinear
///     - `Less_xyz_3` to compare lexicographically two points
///     - `Equal_3` to check whether 2 points are identical
///     -  for each functor Foo, a function `Foo foo_object()`
/// </ul>

/// \return number of degenerate faces found
///
template <class TriangleMesh, class NamedParameters>
std::size_t remove_degenerate_faces(TriangleMesh& tmesh,
                                    const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_pure_triangle(tmesh));

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
  VertexPointMap vpmap = choose_const_pmap(get_param(np, boost::vertex_point),
                                           tmesh,
                                           boost::vertex_point);
  typedef typename GetKernel<TM, NamedParameters>::type Traits;
  Traits traits = choose_param(get_param(np, geom_traits), Traits());

  std::size_t nb_deg_faces = 0;

// First remove edges of length 0
  std::set<edge_descriptor> edges_to_remove;

  BOOST_FOREACH(edge_descriptor ed, edges(tmesh))
  {
    if ( traits.equal_3_object()(get(vpmap, target(ed, tmesh)),
                                 get(vpmap, source(ed, tmesh))) )
      edges_to_remove.insert(ed);
  }

  std::size_t nb_edges_previously_skipped = 0;
  do{
    std::set<edge_descriptor> edges_to_remove_skipped;
    while (!edges_to_remove.empty())
    {
      edge_descriptor ed = *edges_to_remove.begin();
      edges_to_remove.erase(edges_to_remove.begin());

      halfedge_descriptor h = halfedge(ed, tmesh);

      if (CGAL::Euler::satisfies_link_condition(ed,tmesh))
      {
        // remove edges that could also be set for removal
        if ( face(h, tmesh)!=GT::null_face() )
        {
          ++nb_deg_faces;
          edges_to_remove.erase(edge(prev(h, tmesh), tmesh));
          edges_to_remove_skipped.erase(edge(prev(h, tmesh), tmesh));
        }
        if (face(opposite(h, tmesh), tmesh)!=GT::null_face())
        {
          ++nb_deg_faces;
          edges_to_remove.erase(edge(prev(opposite(h, tmesh), tmesh), tmesh));
          edges_to_remove_skipped.erase(edge(prev(opposite(h, tmesh), tmesh), tmesh));
        }
        //now remove the edge
        CGAL::Euler::collapse_edge(ed, tmesh);
      }
      else{
        // the following test filters in particular cases where the
        // detection of the vertex that breaks the link condition
        // is not correctly detected by find_larger_triangle
        // (i.e. when the hedges detected should be used with the opposite edge)
        if ( traits.collinear_3_object()(
              get(vpmap, target(h, tmesh)),
              get(vpmap, target(next(h, tmesh), tmesh)),
              get(vpmap, target(next(opposite(h, tmesh), tmesh), tmesh)) ) )
        {
          // we handle this case in the next loop removing set of degenerated
          // triangles which points are collinear
          edges_to_remove_skipped.insert(ed);
          continue;
        }

        do{
          // improve the link condition on h's side
          boost::optional< std::pair<halfedge_descriptor, halfedge_descriptor> >
            res = internal::find_larger_triangle(h, tmesh, vpmap, traits);
          if (res)
            internal::remove_faces_inside_triangle(res->first, h, res->second, tmesh, edges_to_remove, edges_to_remove_skipped);

          // improve the link condition on opposite(h)'s side
          h=opposite(h, tmesh);
          res = internal::find_larger_triangle(h, tmesh, vpmap, traits);
          if (res)
            internal::remove_faces_inside_triangle(res->first, h, res->second, tmesh, edges_to_remove, edges_to_remove_skipped);
          h=opposite(h, tmesh);
          ed=edge(h, tmesh);

        } while(!CGAL::Euler::satisfies_link_condition(ed,tmesh));

        edges_to_remove.insert( ed );
      }
    }

    // check if some edges were skipped due to link condition not satisfied
    // that could now be satisfied
    if (edges_to_remove_skipped.empty() ||
        edges_to_remove_skipped.size()==nb_edges_previously_skipped) break;
    edges_to_remove.swap(edges_to_remove_skipped);
    nb_edges_previously_skipped = edges_to_remove.size();
  } while(true);

// remove triangles made of 3 collinear points
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
      if (!compare_distance(p1,p2, p1,p3)) // p1p2 > p1p3
      {
        if (!compare_distance(p1,p2, p2,p3)) // p1p2 > p2p3
          // flip p1p2
          edge_to_flip = next( halfedge(fd, tmesh), tmesh );
        else
          // flip p2p3
          edge_to_flip = prev( halfedge(fd, tmesh), tmesh );
      }
      else
        if (!compare_distance(p1,p3, p2,p3)) // p1p3>p2p3
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
