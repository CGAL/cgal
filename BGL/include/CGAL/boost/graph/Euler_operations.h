// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Philipp Moeller

#ifndef CGAL_EULER_OPERATIONS_H
#define CGAL_EULER_OPERATIONS_H

#include <stdexcept>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/internal/helpers.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

namespace EulerImpl {

template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
join_face(typename boost::graph_traits<Graph>::halfedge_descriptor h,
          Graph& g)
{
 typedef typename boost::graph_traits<Graph> Traits;
  typedef typename Traits::halfedge_descriptor           halfedge_descriptor;

  typedef typename Traits::face_descriptor               face_descriptor;


  halfedge_descriptor hop = opposite(h,g);
  halfedge_descriptor hprev = prev(h, g), gprev = prev(hop, g);
  face_descriptor f = face(h, g), f2 = face(hop, g);

  internal::remove_tip(hprev, g);

  internal::remove_tip(gprev, g);

  if(! is_border(hop,g)){
    remove_face(f2, g);
  }
  bool fnull = is_border(h,g);

  
  halfedge_descriptor hprev2 = hprev;
  while(hprev2 != gprev) {
    hprev2 = next(hprev2, g);
    set_face(hprev2, f, g);
  }
  
  if (! fnull)
    set_halfedge(f, hprev, g);
  set_halfedge(target(hprev,g), hprev, g);
  set_halfedge(target(gprev,g), gprev, g);
  //  internal::set_constant_vertex_is_border(g, target(h, g));
  //  internal::set_constant_vertex_is_border(g, target(opposite(h, g), g));

  remove_edge(edge(h, g), g);
  return hprev;

}
} // namespace EulerImpl

/// \endcond

  namespace Euler {
/// \ingroup PkgBGLEulerOperations
/// @{


/**  
 * joins the two vertices incident to `h`, (that is `source(h, g)` and
 * `target(h, g)`) and removes `source(h,g)`. Returns the predecessor
 * of `h` around the vertex, i.e., `prev(opposite(h,g))`.  The
 * invariant `join_vertex(split_vertex(h,g),g)` returns `h`.  The
 * time complexity is linear in the degree of the vertex removed.
 *
 * \image html join_vertex.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \param g the graph
 * \param h the halfedge which incident vertices are joint
 *
 * \returns `prev(opposite(h,g))`
 *
 * \pre The size of the faces incident to `h` and `opposite(h,g)` is at least 4.
 *
 * \post `source(h, g)` is invalidated
 * \post `h` is invalidated 
 * 
 * \sa `split_vertex()`
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
join_vertex(typename boost::graph_traits<Graph>::halfedge_descriptor h,
            Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;
  typedef Halfedge_around_target_iterator<Graph>           halfedge_around_vertex_iterator;

  halfedge_descriptor hop = opposite(h, g)
    , hprev = prev(hop, g)
    , gprev = prev(h, g)
    , hnext = next(hop, g)
    , gnext = next(h, g);
  vertex_descriptor v_to_remove = target(hop, g)
    , v = target(h, g);

  // this assertion fires needlessly
  // CGAL_precondition(std::distance(
  //                     halfedges_around_face(e, g).first,
  //                     halfedges_around_face(e, g).second) >= 4);

  CGAL_assertion( halfedge(v_to_remove, v, g).first == h );

  halfedge_around_vertex_iterator ieb, iee;
  for(boost::tie(ieb, iee) = halfedges_around_target(hop, g); ieb != iee; ++ieb) {
    CGAL_assertion( target(*ieb,g) == v_to_remove);
    set_target(*ieb ,v , g);
  }

  set_next(hprev, hnext, g);
  set_next(gprev, gnext, g);
  set_halfedge(v, gprev, g);
  // internal::set_constant_vertex_is_border(g, v);

  remove_edge(edge(h, g), g);
  remove_vertex(v_to_remove, g);

  return hprev;
}



/** 
 * splits the target vertex `v` of `h1` and `h2`, and connects the new vertex
 * and `v` with a new edge. Let `hnew` be `opposite(next(h1, g), g)` after the
 * split. The split regroups the halfedges around the two vertices. The
 * edge sequence `hnew`, `opposite(next(h2, g), g)`, ..., `h1`
 * remains around the old vertex, while the halfedge sequence
 * `opposite(hnew, g)`, `opposite(next(h1, g), g)` (before the
 * split), ..., `h2` is regrouped around the new vertex. The split
 * returns `hnew`, i.e., the new edge incident to vertex `v`. The
 * time is proportional to the distance from `h1` to `h2` around the
 * vertex.
 *
 * \image html split_vertex.svg
 *
 * \tparam Graph must be a model of  `MutableFaceGraph`
 *
 * \param g the graph
 * \param h1 halfedge descriptor
 * \param h2 halfedge descriptor
 *
 * \returns `hnew`
 *
 * \pre `target(h1, g) == target(h2, g)`, that is  `h1` and `h2` are incident to the same vertex
 * \pre `h1 != h2`, that is no antennas
 *
 * \sa `join_vertex()`
 *
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
split_vertex(typename boost::graph_traits<Graph>::halfedge_descriptor h1,
             typename boost::graph_traits<Graph>::halfedge_descriptor h2,
             Graph& g)
{
  CGAL_assertion(h1 != h2);
  CGAL_assertion(target(h1, g) == target(h2, g));

  typename boost::graph_traits<Graph>::halfedge_descriptor 
    hnew = halfedge(add_edge(g), g),
    hnewopp = opposite(hnew, g);
  typename boost::graph_traits<Graph>::vertex_descriptor
    vnew = add_vertex(g);
  internal::insert_halfedge(hnew, h2, g);
  internal::insert_halfedge(hnewopp, h1, g);
  set_target(hnew, target(h1, g), g);
  
  typename boost::graph_traits<Graph>::halfedge_descriptor 
    end = hnewopp;
  do
  {
    set_target(hnewopp, vnew, g);
    hnewopp = opposite(next(hnewopp, g), g);
  } while (hnewopp != end);

  internal::set_vertex_halfedge(hnew, g);
  // internal::set_constant_vertex_is_border(g, target(hnew, g));
  internal::set_vertex_halfedge(hnewopp, g);
  // internal::set_constant_vertex_is_border(g, target(hnewopp, g));
  return hnew;
}

/**
 * splits the halfedge `h` into two halfedges inserting a new vertex that is a copy of `vertex(opposite(h,g),g)`.
 * Is equivalent to `opposite(split_vertex( prev(h,g), opposite(h,g),g), g)`.
 * \returns the new halfedge `hnew` pointing to the inserted vertex. The new halfedge is followed by the old halfedge, i.e., `next(hnew,g) == h`.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
split_edge(typename boost::graph_traits<Graph>::halfedge_descriptor h, Graph& g)
{ return opposite(split_vertex(prev(h,g), opposite(h,g),g), g); }


/**
 * joins the two faces incident to `h` and `opposite(h,g)`. 
 * The faces may be holes.
 *
 * If `Graph` is a model of `MutableFaceGraph`
 * the face incident to `opposite(h,g)` is removed.
 *
 * `join_face()` and `split_face()` are inverse operations, that is
 * `join_face(split_face(h,g),g)` returns `h`.
 *
 * \image html join_face.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`.
 * \param g the graph
 * \param h the halfedge incident to one of the faces to be joined.
 *
 * \returns `prev(h,g)`
 *
 * \pre `out_degree(source(h,g)), g)) >= 3`
 * \pre `out_degree(target(h,g)) >= 3`
 *
 * \sa `split_face()`
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
join_face(typename boost::graph_traits<Graph>::halfedge_descriptor h,
          Graph& g)
{
  return EulerImpl::join_face(h,g);
}



/**
 * splits the face incident to `h1` and `h2`.  Creates the opposite
 * halfedges `h3` and `h4`, such that `next(h1,g) == h3` and `next(h2,g) == h4`.
 * Performs the inverse operation to `join_face()`. 
 *
 * If `Graph` is a model of `MutableFaceGraph` and if the update of faces is not disabled
 * a new face incident to `h4` is added. 
 *
 * \image html split_face.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \param g the graph
 * \param h1
 * \param h2
 *
 * \returns `h3`
 *
 * \pre `h1` and `h2` are incident to the same face
 * \pre `h1 != h2`
 * \pre `next(h1,g) != h2` and `next(h2,g) != h1` (no loop)
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
split_face(typename boost::graph_traits<Graph>::halfedge_descriptor h1,
           typename boost::graph_traits<Graph>::halfedge_descriptor h2,
           Graph& g)
{ 
  typedef typename boost::graph_traits<Graph> Traits;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Traits::face_descriptor face_descriptor;
  halfedge_descriptor hnew = halfedge(add_edge(g), g);
  face_descriptor fnew = add_face(g);
  internal::insert_tip( hnew, h2, g);
  internal::insert_tip( opposite(hnew, g), h1, g);
  set_face( hnew, face(h1,g), g);
  internal::set_face_in_face_loop(opposite(hnew,g), fnew, g);
  set_halfedge(face(hnew,g), hnew, g);
  set_halfedge(face(opposite(hnew,g),g), opposite(hnew,g), g);
  return hnew; 
}

  
/**
 * glues the cycle of halfedges of `h1` and `h2` together.
 * The vertices in the cycle of `h2` get removed.
 * If `h1` or `h2` are not border halfedges their faces get removed.
 * The vertices on the face cycle of `h1` get removed.
 * The invariant `join_loop(h1, split_loop(h1,h2,h3,g), g)` returns `h1` and keeps 
 * the graph unchanged.
 * 
 * \image html join_loop.svg
 *
 * \tparam Graph must be a `MutableFaceGraph`
 *
 * \returns `h1`. 
 *
 * \pre The faces incident to `h` and `g` are different and have equal number of edges.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
join_loop(typename boost::graph_traits<Graph>::halfedge_descriptor h1,
          typename boost::graph_traits<Graph>::halfedge_descriptor h2,
          Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;

  CGAL_precondition( is_border(h1,g) || face(h1, g) != face(h2, g));
  if (! is_border(h1,g))
    remove_face(face(h1, g), g);
  if (! is_border(h2,g))
    remove_face(face(h2,g), g);
  halfedge_descriptor hi = h1;
  halfedge_descriptor gi = h2;
  CGAL_assertion_code( std::size_t termination_count = 0;)
  do {
    CGAL_assertion( ++termination_count != 0);
    halfedge_descriptor hii = hi;
    halfedge_descriptor gii = gi;
    hi = next(hi, g);
    // gi = find_prev(gi); // Replaced by search around vertex.
    set_face( hii, face( opposite(gii, g), g), g);
    set_halfedge(face(hii, g), hii, g);
    remove_vertex(target(opposite(gii, g), g), g);
    if ( next(opposite(next(opposite(gii,g), g), g), g) == gii) {
      gi = opposite(next(opposite(gii,g),g), g);
    } else {
      set_next(hii, next(opposite(gii,g), g), g);
      gii = opposite(next(opposite(gii, g), g), g);
      set_target( gii, target(hii, g), g);
      while ( next(opposite(next(gii, g), g), g) != gi) {
        CGAL_assertion( ++termination_count != 0);
        gii = opposite(next(gii,g), g);
        set_target( gii, target(hii, g), g);
      }
      gi = opposite(next(gii,g), g);
      set_next(gii, hi, g);
    }
  } while ( hi != h1);
  CGAL_assertion( gi == h2);
  do {
    halfedge_descriptor gii = gi;
    gi = next(gi, g);
    remove_edge(edge(gii,g), g);
  } while ( gi != h2);
  return h1;
}


/**
 * cuts the graph along the cycle `(h1,h2,h3)` changing the genus 
 * (halfedge `h3` runs on the backside of the three dimensional figure below).
 * Three new vertices, three new pairs of halfedges,
 * and two new triangular faces are created.
 *
 * `h1`, `h2`, and `h3` will be incident to the first new face. 
 *
 * Note that `split_loop()` does not deal with properties of new vertices, halfedges, and faces.
 *
 * \image html split_loop.svg
 * 
 * \tparam Graph must be a `MutableFaceGraph`
 *
 * \returns the halfedge incident to the second new face.
 *
 * \pre `h1`, `h2`, and `h3` denote distinct, consecutive halfedges of the graph 
 * and form a cycle: i.e., `target(h1) == target(opposite(h2,g),g)`, … , 
 * `target(h3,g) == target(opposite(h1,g),g)`. 
 * \pre The six faces incident to `h1`, `h2`, and `h3` are all distinct.
 */
  template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
split_loop(typename boost::graph_traits<Graph>::halfedge_descriptor h1,
           typename boost::graph_traits<Graph>::halfedge_descriptor h2,
           typename boost::graph_traits<Graph>::halfedge_descriptor h3,
           Graph& g)
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::face_descriptor                 face_descriptor;

  halfedge_descriptor h = h1, i = h2, j = h3;
   CGAL_precondition( h != i);
        CGAL_precondition( h != j);
        CGAL_precondition( i != j);
        CGAL_precondition( target(h,g) == target(opposite(i,g),g));
        CGAL_precondition( target(i,g) == target(opposite(j,g),g));
        CGAL_precondition( target(j,g) == target(opposite(h,g),g));
        // Create a copy of the triangle.
        halfedge_descriptor hnew = internal::copy(h,g);
        halfedge_descriptor inew = internal::copy(i,g);
        halfedge_descriptor jnew = internal::copy(j,g);
        internal::close_tip( hnew, add_vertex(g), g);
        internal::close_tip( inew, add_vertex(g), g);
        internal::close_tip( jnew, add_vertex(g), g);
        internal::insert_tip( opposite(inew, g), hnew, g);
        internal::insert_tip( opposite(jnew, g), inew, g);
        internal::insert_tip( opposite(hnew, g), jnew, g);
        // Make the new incidences with the old stucture.
        CGAL_assertion_code( std::size_t termination_count = 0;)
          if ( next(h,g) != i) {
            halfedge_descriptor nh = next(h, g);
            set_next(h, i, g);
            set_next(hnew, nh, g);
            nh = opposite(nh, g);
            while ( next(nh, g) != i) {
                CGAL_assertion( ++termination_count != 0);
                set_target( nh, target(hnew,g), g);
                nh = opposite(next(nh, g), g);
            }
            set_target( nh, target(hnew,g), g);
            set_next(nh, inew, g);
        }
        if ( next(i, g) != j) {
          halfedge_descriptor nh = next(i, g);
          set_next(i, j, g);
          set_next(inew, nh, g);
          nh = opposite(nh,g);
          while ( next(nh,g) != j) {
                CGAL_assertion( ++termination_count != 0);
                set_target( nh, target(inew, g), g);
                nh = opposite(next(nh, g), g);
            }
          set_target( nh, target(inew, g), g);
          set_next(nh, jnew, g);

        }
        if ( next(j,g) != h) {
          halfedge_descriptor nh = next(j, g);
          set_next(j, h, g);
          set_next(jnew, nh, g);
          nh = opposite(nh, g);
          while ( next(nh,g) != h) {
                CGAL_assertion( ++termination_count != 0);
                set_target( nh, target(jnew, g), g);
                nh = opposite(next(nh, g), g);
            }
          set_target(nh, target(jnew, g), g);
          set_next(nh, hnew, g);
        }
        // Fill the holes with two new faces.
        face_descriptor f = add_face(g);
        set_face( h, f, g);
        set_face( i, f, g);
        set_face( j, f, g);
        set_halfedge(face(h,g), h, g);
        f = add_face(g);
        set_face( opposite(hnew, g), f, g);
        set_face( opposite(inew, g), f, g);
        set_face( opposite(jnew, g), f, g);
        set_halfedge(face(opposite(hnew,g),g), opposite(hnew,g), g);
        // Take care of maybe changed halfedge pointers.
        set_halfedge(face(hnew, g), hnew, g);
        set_halfedge(face(inew, g), inew, g);
        set_halfedge(face(jnew, g), jnew, g);
        set_halfedge(target(hnew, g), hnew, g);
        set_halfedge(target(inew, g), inew, g);
        set_halfedge(target(jnew, g), jnew, g);
        return opposite(hnew, g);
}


/**
 * removes the incident face of `h` and changes all halfedges incident to the face into border halfedges 
 * or removes them from the graph if they were already border halfedges.
 *
 * If this creates isolated vertices they get removed as well. 
 *
 * \image html remove_face.svg
 * \image html remove_face_and_vertex.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \pre `h` is not a border halfedge
 *
 * \sa `make_hole()` for a more specialized variant.
 */
template< typename Graph >
void remove_face(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                 Graph& g)
{
  typedef typename boost::graph_traits<Graph>            Traits;
  typedef typename Traits::halfedge_descriptor           halfedge_descriptor;
  typedef typename Traits::face_descriptor               face_descriptor;

  CGAL_precondition(! is_border(h,g));
  face_descriptor f = face(h, g);

  halfedge_descriptor end = h;
  do {
    internal::set_border(h,g);
    halfedge_descriptor nh = next(h, g);
    bool h_border = is_border(opposite(h, g),g);
    bool nh_bborder = is_border(opposite(nh, g),g);

    if(h_border && nh_bborder && next(opposite(nh, g), g) == opposite(h, g)) {
      remove_vertex(target(h, g), g);
      if(h != end)
        remove_edge(edge(h, g), g);
    } else {
      if(nh_bborder) {
        internal::set_vertex_halfedge(opposite(next(opposite(nh, g), g), g), g);
        internal::remove_tip(h, g);
        //internal::set_constant_vertex_is_border(g, target(h, g));
      }
      if(h_border) {
        internal::set_vertex_halfedge(opposite(next(h, g), g), g);
        internal::remove_tip(prev(opposite(h, g), g), g);
        //internal::set_constant_vertex_is_border(g, target(prev(opposite(h, g), g), g));
        if(h != end)
          remove_edge(edge(h, g), g);
      }
    }
    h = nh;
  } while(h != end);
  remove_face(f, g);

  if(is_border(opposite(h, g),g))
    remove_edge(edge(h, g), g);
}

/**
* adds and returns the edge `e` connecting `s` and `t`
* halfedge(e, g) has s as source and t as target
*/
template<typename Graph>
typename boost::graph_traits<Graph>::edge_descriptor
add_edge(typename boost::graph_traits<Graph>::vertex_descriptor s,
         typename boost::graph_traits<Graph>::vertex_descriptor t,
         Graph& g)
{
  typename boost::graph_traits<Graph>::edge_descriptor e = add_edge(g);
  set_target(halfedge(e, g), t, g);
  set_target(opposite(halfedge(e, g), g), s, g);
  return e;
}

/**
* adds a new face defined by a range of vertices (identified by their descriptors,
* `boost::graph_traits<Graph>::%vertex_descriptor`).
* For each pair of consecutive vertices, the corresponding halfedge
* is added in `g` if new, and its connectivity is updated otherwise.
* The face can be added only at the boundary of `g`, or as a new connected component.
*
* @pre `vr` contains at least 3 vertices
* @returns the added face descriptor, or `boost::graph_traits<Graph>::%null_face()` if the face could not be added.
*/
template< typename Graph, typename VertexRange >
typename boost::graph_traits<Graph>::face_descriptor
add_face(const VertexRange& vr, Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor     face_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor     edge_descriptor;

  std::vector<vertex_descriptor> vertices(vr.begin(), vr.end()); // quick and dirty copy
  unsigned int n = (unsigned int)vertices.size();
  // don't allow degenerated faces
  CGAL_assertion(n > 2);

  std::vector<halfedge_descriptor> halfedges(n);
  std::vector<bool>                is_new(n);

  for (unsigned int i = 0, ii = 1; i<n; ++i, ++ii, ii %= n)
  {
    if ( ! internal::is_isolated(vertices[i], g)
      && ! is_border(vertices[i], g))
      return boost::graph_traits<Graph>::null_face();

    std::pair<halfedge_descriptor, bool> he
      = halfedge(vertices[i], vertices[ii], g);
    halfedges[i] = he.first;//collect if exists
    is_new[i] = !(he.second/*true if exists*/);

    if (!is_new[i] && !is_border(halfedges[i], g))
      return boost::graph_traits<Graph>::null_face();
  }

  halfedge_descriptor inner_next, inner_prev,
                      outer_next, outer_prev,
                      border_next, border_prev,
                      patch_start, patch_end;
  // cache for set_next and vertex' set_halfedge
  typedef std::pair<halfedge_descriptor, halfedge_descriptor> NextCacheEntry;
  typedef std::vector<NextCacheEntry>    NextCache;
  NextCache next_cache;
  next_cache.reserve(3 * n);

  // re-link patches if necessary
  for (unsigned int i = 0, ii = 1; i<n; ++i, ++ii, ii %= n)
  {
    if (!is_new[i] && !is_new[ii])
    {
      inner_prev = halfedges[i];
      inner_next = halfedges[ii];

      if (next(inner_prev, g) != inner_next)
      {
        // here comes the ugly part... we have to relink a whole patch

        // search a free gap
        // free gap will be between border_prev and border_next
        outer_prev = opposite(inner_next, g);
        outer_next = opposite(inner_prev, g);
        border_prev = outer_prev;
        do{
          border_prev = opposite(next(border_prev, g), g);
        }while (!is_border(border_prev, g) || border_prev == inner_prev);
        border_next = next(border_prev, g);
        CGAL_assertion(is_border(border_prev, g));
        CGAL_assertion(is_border(border_next, g));

        if (border_next == inner_next)
          return boost::graph_traits<Graph>::null_face();

        // other halfedges' indices
        patch_start = next(inner_prev, g);
        patch_end   = prev(inner_next, g);

        // relink
        next_cache.push_back(NextCacheEntry(border_prev, patch_start));
        next_cache.push_back(NextCacheEntry(patch_end, border_next));
        next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
      }
    }
  }
  // create missing edges
  for (unsigned int i = 0, ii = 1; i<n; ++i, ++ii, ii %= n)
  {
    if (is_new[i])
    {
      edge_descriptor ne = add_edge(vertices[i], vertices[ii], g);
      halfedges[i] = halfedge(ne, g);
      CGAL_assertion(halfedges[i] != boost::graph_traits<Graph>::null_halfedge());

      set_face(opposite(halfedges[i], g), boost::graph_traits<Graph>::null_face(), g); // as it may be recycled we have to reset it  
      CGAL_assertion(source(halfedges[i], g) == vertices[i]);
    }
  }
  // create the face
  face_descriptor f = add_face(g);
  set_halfedge(f, halfedges[n - 1], g);

  // setup halfedges
  for (unsigned int i = 0, ii = 1; i<n; ++i, ++ii, ii %= n)
  {
    vertex_descriptor v = vertices[ii];
    inner_prev = halfedges[i];
    inner_next = halfedges[ii];

    unsigned int id = 0;
    if (is_new[i])  id |= 1;
    if (is_new[ii]) id |= 2;

    if (id)
    {
      outer_prev = opposite(inner_next, g);
      outer_next = opposite(inner_prev, g);

      // set outer links
      switch (id)
      {
      case 1: // prev is new, next is old
        border_prev = prev(inner_next, g);
        next_cache.push_back(NextCacheEntry(border_prev, outer_next));
        set_halfedge(v, border_prev, g);
        break;

      case 2: // next is new, prev is old
        border_next = next(inner_prev, g);
        next_cache.push_back(NextCacheEntry(outer_prev, border_next));
        set_halfedge(v, outer_prev, g);
        break;

      case 3: // both are new
        if (halfedge(v, g) == boost::graph_traits<Graph>::null_halfedge())
        {
          set_halfedge(v, outer_prev, g);
          next_cache.push_back(NextCacheEntry(outer_prev, outer_next));
        }
        else
        {
          border_prev = halfedge(v, g);
          border_next = next(border_prev, g);
          next_cache.push_back(NextCacheEntry(border_prev, outer_next));
          next_cache.push_back(NextCacheEntry(outer_prev, border_next));
        }
        break;
      }

      // set inner link
      next_cache.push_back(NextCacheEntry(inner_prev, inner_next));
    }

    // set face index
    set_face(halfedges[i], f, g);
  }

  // process next halfedge cache
  typename NextCache::const_iterator ncIt(next_cache.begin()), ncEnd(next_cache.end());
  for (; ncIt != ncEnd; ++ncIt)
    set_next(ncIt->first, ncIt->second, g);

  // adjust vertices' halfedge index
  for (unsigned int i = 0; i<n; ++i)
    internal::adjust_incoming_halfedge(vertices[i], g);

  return f;
}


  /**
   * removes the incident face of `h` and changes all halfedges incident to the face into border halfedges. See `remove_face(g,h)` for a more generalized variant.
   *
   * \pre None of the incident edges of the face is a border edge.
   */
template< typename Graph>
void make_hole(typename boost::graph_traits<Graph>::halfedge_descriptor h,
               Graph& g)
{
  typedef typename boost::graph_traits<Graph>            Traits;
  typedef typename Traits::face_descriptor               face_descriptor;
  typedef Halfedge_around_face_iterator<Graph>           halfedge_around_face_iterator;

  CGAL_precondition(! is_border(h,g));
  face_descriptor fd = face(h, g);
  halfedge_around_face_iterator hafib, hafie;
  for(boost::tie(hafib, hafie) = halfedges_around_face(h, g); 
      hafib != hafie; 
      ++hafib){
    CGAL_assertion(! is_border(opposite(*hafib,g),g));
    internal::set_border(*hafib, g);
  }
  remove_face(fd,g);  
}


    /** fills the hole incident to `h`.
     * \pre `h` must be a border halfedge
     */
template< typename Graph>
void fill_hole(typename boost::graph_traits<Graph>::halfedge_descriptor h,
               Graph& g)
{
  typedef typename boost::graph_traits<Graph>  Traits;
  typedef typename Traits::face_descriptor     face_descriptor;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;

  face_descriptor f = add_face(g);
  BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(h,g)){
    set_face(hd, f,g);
  }
  set_halfedge(f,h,g);
}


/** 
 * creates a barycentric triangulation of the face incident to `h`. Creates a new
 * vertex and connects it to each vertex incident to `h` and splits `face(h, g)` 
 * into triangular faces.
 * `h` remains incident to
 * the original face. The time complexity is linear in the size of the face.
 *
 * \image html add_center_vertex.svg
 *
 * \returns the halfedge `next(h, g)` after the
 * operation, i.e., a halfedge pointing to the new vertex. 
 *
 * Note that `add_center_vertex()` does not deal with properties of new vertices, 
 * halfedges, and faces. 
 *  \pre `h` is not a border halfedge.
 *
 * \param g the graph
 * \param h halfedge descriptor
 * \tparam Graph must be a model of `MutableFaceGraph`
 * \sa `remove_center_vertex()`
 *
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
add_center_vertex(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                  Graph& g) 
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::vertex_descriptor               vertex_descriptor;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;
  typedef typename Traits::face_descriptor                 face_descriptor;

  halfedge_descriptor hnew = halfedge(add_edge(g),g);
  vertex_descriptor vnew = add_vertex(g);
  internal::close_tip(hnew, vnew, g);
  internal::insert_tip(opposite(hnew, g), h, g);
  set_face(hnew, face(h, g), g);
  set_halfedge(face(h,g), h, g);
  halfedge_descriptor h2 = next(opposite(hnew, g), g);
  while ( next(h2, g) != hnew) {
    halfedge_descriptor gnew = halfedge(add_edge(g),g);
    internal::insert_tip( gnew, hnew, g);
    internal::insert_tip( opposite(gnew,g), h2, g);
    face_descriptor fnew = add_face(g);
    set_face( h2, fnew, g);
    set_face( gnew, fnew, g);
    set_face( next(gnew,g), fnew, g);
    set_halfedge(face(h2, g), h2, g);
    h2 = next(opposite(gnew, g), g);
  }
  set_face(next(hnew,g), face(hnew,g), g);
  internal::set_vertex_halfedge(hnew, g);
  return hnew;
}

/**
 * removes the vertex `target(h, g)` and all incident halfedges thereby merging all
 * incident faces.   The resulting face may not be triangulated. 
 * This function is the inverse operation of `add_center_vertex()`. 
 * The invariant `h == remove_center_vertex(add_center_vertex(h,g),g)` 
 * holds, if `h` is not a border halfedge. 
 *
 * \image html remove_center_vertex.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \param g the graph
 * \param h halfedge descriptor
 *
 * \returns `prev(h, g)`
 *
 * \pre None of the incident faces of `target(h,g)` is a
 * hole. There are at least two distinct faces incident to the faces
 * that are incident to `target(h,g)`. (This prevents the
 * operation from collapsing a volume into two faces glued together
 * with opposite orientations, such as would happen with any vertex of
 * a tetrahedron.)
 *
 * \sa `add_center_vertex()`
 *
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
remove_center_vertex(typename boost::graph_traits<Graph>::halfedge_descriptor h,
                     Graph& g) 
{
  typedef typename boost::graph_traits<Graph>              Traits;
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;

  // h points to the vertex that gets removed
  halfedge_descriptor h2    = opposite(next(h, g), g);
  halfedge_descriptor hret = prev(h, g);
  while (h2 != h) {
    halfedge_descriptor gprev = prev(h2, g);
    internal::set_vertex_halfedge(gprev, g);
    internal::remove_tip(gprev, g);

    remove_face(face(h2, g), g);

    halfedge_descriptor gnext = opposite(next(h2, g), g);
    remove_edge(edge(h2,g), g);
    h2 = gnext;
  }
  internal::set_vertex_halfedge(hret, g);
  internal::remove_tip(hret, g);
  remove_vertex(target(h, g), g);
  remove_edge(edge(h, g), g);
  internal::set_face_in_face_loop(hret, face(hret, g), g);
  set_halfedge(face(hret, g), hret, g);
  return hret;
}

/**
 * appends a new face to the border halfedge `h2` by connecting 
 * the tip of `h2` with the tip of `h1` with two new halfedges and a new vertex 
 * and creating a new face that is incident to `h2`. 
 * Note that `add_vertex_and_face_to_border()` does not deal with properties of new 
 * vertices, halfedges, and faces.
 * 
 * \image html add_vertex_and_face_to_border.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \returns the halfedge of the new edge that is incident to the new face
 * and the new vertex.
 *
 * \pre `h1` and `h2` are border halfedges
 * \pre `h1 != h2`,
 * \pre `h1` and `h2` are on the same border.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
add_vertex_and_face_to_border(typename boost::graph_traits<Graph>::halfedge_descriptor h1,
                              typename boost::graph_traits<Graph>::halfedge_descriptor h2,
                              Graph& g)
{
  typename boost::graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
  typename boost::graph_traits<Graph>::face_descriptor f = add_face(g);
  typename boost::graph_traits<Graph>::edge_descriptor e1 = add_edge(g);
  typename boost::graph_traits<Graph>::edge_descriptor e2 = add_edge(g);
  typename boost::graph_traits<Graph>::halfedge_descriptor he1= halfedge(e1, g);
  typename boost::graph_traits<Graph>::halfedge_descriptor he2= halfedge(e2, g);
  typename boost::graph_traits<Graph>::halfedge_descriptor ohe1= opposite(he1, g);
  typename boost::graph_traits<Graph>::halfedge_descriptor ohe2= opposite(he2, g);

  set_next(he1, next(h1,g),g);
  set_next(h1,ohe1,g);
  set_target(he1,target(h1,g),g);
  set_target(ohe1,v,g);

  set_next(he2,he1,g);
  set_next(ohe1,ohe2,g);
  set_target(he2,v,g);
  set_halfedge(v,ohe1,g);
  set_next(ohe2,next(h2,g),g);
  set_target(ohe2,target(h2,g),g);
  set_next(h2,he2,g);
  internal::set_border(ohe1,g);
  internal::set_border(ohe2,g);

  CGAL::Halfedge_around_face_iterator<Graph> hafib,hafie;
  for(boost::tie(hafib, hafie) = halfedges_around_face(he1, g);
      hafib != hafie;
      ++hafib){
    set_face(*hafib, f, g);
  }
  set_halfedge(f, he1, g);
  return ohe2;
}


/**
 * appends a new face incident to the border halfedge `h1` and `h2` by connecting the vertex `target(h2,g)` 
 * and the vertex `target(h1,g)` with a new halfedge, and filling this separated part of the hole 
 * with a new face, such that the new face is incident to `h2`. 
 *
 * \image html add_face_to_border.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \returns the halfedge of the new edge that is incident to the new face.
 *
 * \pre  `h1` and `h2` are border halfedges, 
 * \pre `h1 != h2`, 
 * \pre `next(h1,g) != h2`,
 * \pre `h1` and `h2` are on the same border.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
add_face_to_border(typename boost::graph_traits<Graph>::halfedge_descriptor h1,
                   typename boost::graph_traits<Graph>::halfedge_descriptor h2,
                   Graph& g)
{
  CGAL_precondition(is_border(h1,g) == true);
  CGAL_precondition(is_border(h2,g) == true);
  CGAL_precondition(h1 != h2);
  CGAL_precondition(next(h1, g) != h2);

  typename boost::graph_traits<Graph>::face_descriptor f = add_face(g);
  typename boost::graph_traits<Graph>::edge_descriptor e = add_edge(g);
  typename boost::graph_traits<Graph>::halfedge_descriptor 
      newh= halfedge(e, g)
    , newhop = opposite(newh, g);

  set_next(newhop, next(h2, g), g);

  set_next(h2, newh, g);
  
  set_next(newh, next(h1, g), g);

  set_next(h1, newhop, g);
  
  set_target(newh, target(h1, g), g);
  set_target(newhop, target(h2, g), g);

  // make the vertices point to the border halfedge
  set_halfedge(target(h2,g), newhop, g);
  internal::set_border(newhop, g);

  CGAL::Halfedge_around_face_iterator<Graph> hafib,hafie;
  for(boost::tie(hafib, hafie) = halfedges_around_face(newh, g); 
      hafib != hafie; 
      ++hafib){
    set_face(*hafib, f, g);
  }

  set_halfedge(f, newh, g);

  return newh;
}


/**
 * collapses an edge in a graph.
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 * Let `v0` and `v1` be the source and target vertices, and let `e` and `e'` be the halfedges of edge `v0v1`.
 *
 * For `e`, let `en` and `ep` be the next and previous
 * halfedges, that is `en = next(e, g)`, `ep = prev(e, g)`, and let
 * `eno` and `epo` be their opposite halfedges, that is
 * `eno = opposite(en, g)` and `epo = opposite(ep, g)`.
 * Analoguously, for `e'` define  `en'`, `ep'`, `eno'`, and  `epo'`.
 *
 * Then, after the collapse of edge `v0v1` the following holds for `e` (and analoguously for `e'`)
 *
 * <UL>
 *   <LI>The edge `v0v1` is no longer in `g`.
 *   <LI>The faces incident to edge `v0v1` are no longer in `g`.
 *   <LI>Either `v0`, or `v1` is no longer in `g` while the other remains.
 *       Let `vgone` be the removed vertex and `vkept` be the remaining vertex.
 *   <LI>If `e` was a border halfedge, that is `is_border(e, g) == true`, then `next(ep,g) == en`, and `prev(en,g) == ep`.
 *   <LI>If `e` was not a border halfedge, that is `is_border(e, g) == false`, then `ep` and `epo` are no longer in `g` while `en` and `eno` are kept in `g`.
 *   <LI>For all halfedges `hv` in `halfedges_around_target(vgone, g)`, `target(*hv, g) == vkept` and `source(opposite(*hv, g), g) == vkept`.
 *   <LI>No other incidence information has changed in `g`.
 * </UL>
 * \returns vertex `vkept` (which can be either `v0` or `v1`).
 * \pre g must be a triangulated graph
 * \pre `does_satisfy_link_condition(v0v1,g) == true`.
 */
template<typename Graph>
typename boost::graph_traits<Graph>::vertex_descriptor
collapse_edge(typename boost::graph_traits<Graph>::edge_descriptor v0v1,
              Graph& g)
{
  typedef boost::graph_traits< Graph > Traits;
  typedef typename Traits::vertex_descriptor          vertex_descriptor;
  typedef typename Traits::halfedge_descriptor            halfedge_descriptor;

  halfedge_descriptor pq = halfedge(v0v1,g);
  halfedge_descriptor qp = opposite(pq, g);
  halfedge_descriptor pt = opposite(prev(pq, g), g);
  halfedge_descriptor qb = opposite(prev(qp, g), g);
  
  bool lTopFaceExists         = ! is_border(pq,g);
  bool lBottomFaceExists      = ! is_border(qp,g);
  bool lTopLeftFaceExists     = lTopFaceExists    && ! is_border(pt,g);
  bool lBottomRightFaceExists = lBottomFaceExists && ! is_border(qb,g);

  CGAL_precondition( !lTopFaceExists    || (lTopFaceExists    && ( degree(target(pt, g), g) > 2 ) ) ) ;
  CGAL_precondition( !lBottomFaceExists || (lBottomFaceExists && ( degree(target(qb, g), g) > 2 ) ) ) ;

  vertex_descriptor q = target(pq, g);
  vertex_descriptor p = source(pq, g);
#if 0
  if(lTopLeftFaceExists && lBottomRightFaceExists){
    std::cerr <<    " // do it low level" << std::endl;
    halfedge_descriptor qt = next(pq,g);
    halfedge_descriptor pb = next(qp,g);
    halfedge_descriptor ppt = prev(pt,g);
    halfedge_descriptor pqb = prev(qb,g);
    if(halfedge(q,g) == pq){
      set_halfedge(q, pqb,g);
    }
    vertex_descriptor t = target(qt,g);
    if(halfedge(t,g) == pt){
      set_halfedge(t, qt,g);
    } 
    vertex_descriptor b = target(pb,g);
    if(halfedge(b,g) == qb){
      set_halfedge(t, pb,g);
    }
    set_face(qt, face(pt,g),g);
    set_halfedge(face(qt,g),qt,g);
    set_face(pb, face(qb,g),g);
    set_halfedge(face(pb,g),pb,g);
    set_next(qt, next(pt,g),g);
    set_next(pb, next(qb,g),g);
    set_next(ppt, qt,g);
    set_next(pqb,pb,g);
    remove_face(face(pq,g),g);
    remove_face(face(qp,g),g);
    remove_edge(v0v1,g);
    remove_edge(edge(pt,g),g);
    remove_edge(edge(qb,g),g);
    remove_vertex(p,g);
    Halfedge_around_target_circulator<Graph> beg(ppt,g), end(pqb,g);
    while(beg != end){
      CGAL_assertion(target(*beg,g) == p);
      set_target(*beg,q,g);
      --beg;
    }

    return q;
    // return the vertex kept
  }
#endif

  bool lP_Erased = false, lQ_Erased = false ;

  if ( lTopFaceExists )
  { 
    CGAL_precondition( ! is_border(opposite(pt, g),g) ) ; // p-q-t is a face of the mesh
    if ( lTopLeftFaceExists )
    {
      //CGAL_ECMS_TRACE(3, "Removing p-t E" << pt.idx() << " (V" 
      //                << p.idx() << "->V" << target(pt, g).idx() 
      //                << ") by joining top-left face" ) ;

      join_face(pt,g);
    }
    else
    {
      //CGAL_ECMS_TRACE(3, "Removing p-t E" << pt.idx() << " (V" << p.idx() 
      //                << "->V" << target(pt, g).idx() << ") by erasing top face" ) ;

      remove_face(opposite(pt, g),g);

      if ( !lBottomFaceExists )
      {
        //CGAL_ECMS_TRACE(3, "Bottom face doesn't exist so vertex P already removed" ) ;

        lP_Erased = true ;
      }  
    } 
  }

  if ( lBottomFaceExists )
  {   
    CGAL_precondition( ! is_border(opposite(qb, g),g) ) ; // p-q-b is a face of the mesh
    if ( lBottomRightFaceExists )
    {
      //CGAL_ECMS_TRACE(3, "Removing q-b E" << qb.idx() << " (V" 
      //                << q.idx() << "->V" << target(qb, g).idx() 
      //                << ") by joining bottom-right face" ) ;

      join_face(qb,g);
    }
    else
    {
      //CGAL_ECMS_TRACE(3, "Removing q-b E" << qb.idx() << " (V" 
      //                << q.idx() << "->V" << target(qb, g).idx() 
      //                << ") by erasing bottom face" ) ;

      remove_face(opposite(qb, g),g);

      if ( !lTopFaceExists )
      {
        //CGAL_ECMS_TRACE(3, "Top face doesn't exist so vertex Q already removed" ) ;
        lQ_Erased = true ;
      }  
    }
  }

  CGAL_assertion( !lP_Erased || !lQ_Erased ) ;

  if ( !lP_Erased && !lQ_Erased )
  {
    //CGAL_ECMS_TRACE(3, "Removing vertex P by joining pQ" ) ;

    join_vertex(pq,g);
    lP_Erased = true ;
  }    

  CGAL_expensive_assertion(is_valid(g));

  return lP_Erased ? q : p ;
}

/**
 * Collapses the edge `v0v1` replacing it with v0 or v1, as described in the paragraph above
 * and guarantees that an edge `e2`, for which `get(edge_is_constrained_map, e2)==true`, 
 * is not removed after the collapse.
 * 
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 * \tparam EdgeIsConstrainedMap mut be a model of `ReadablePropertyMap` with the edge descriptor of `Graph` 
 *       as key type and a Boolean as value type. It indicates if an edge is constrained or not. 
 *
 * \pre This function requires `g` to be an oriented 2-manifold with or without boundaries. 
 *       Furthermore, the edge `v0v1` must satisfy the link condition, which guarantees that the surface mesh is also 2-manifold after the edge collapse. 
 * \pre `get(edge_is_constrained_map, v0v1)==false`. 
 * \pre  `v0` and `v1` are not both incident to a constrained edge. 
 */

template<typename Graph, typename EdgeIsConstrainedMap>
typename boost::graph_traits<Graph>::vertex_descriptor
collapse_edge(typename boost::graph_traits<Graph>::edge_descriptor v0v1,
              Graph& g
              , EdgeIsConstrainedMap Edge_is_constrained_map)
{
  typedef boost::graph_traits< Graph > Traits;
  typedef typename Traits::vertex_descriptor          vertex_descriptor;
  typedef typename Traits::halfedge_descriptor            halfedge_descriptor;

  halfedge_descriptor pq = halfedge(v0v1,g);
  CGAL_assertion( !get(Edge_is_constrained_map,v0v1) );

  halfedge_descriptor qp = opposite(pq,g);
  halfedge_descriptor pt = opposite(prev(pq,g),g);
  halfedge_descriptor qb = opposite(prev(qp,g),g);
  halfedge_descriptor tq = opposite(next(pq,g),g);
  halfedge_descriptor bp = opposite(next(qp,g),g);

  bool lTopFaceExists         = ! is_border(pq,g) ;
  bool lBottomFaceExists      = ! is_border(qp,g) ;

  vertex_descriptor q = target(pq,g);
  vertex_descriptor p = source(pq,g);

  //used to collect edges to remove from the surface
  halfedge_descriptor edges_to_erase[2];
  halfedge_descriptor* edges_to_erase_ptr=edges_to_erase;

  // If the top facet exists, we need to choose one out of the two edges which one disappears:
  //   p-t if it is not constrained and t-q otherwise
  if ( lTopFaceExists )
  {
    if ( !get(Edge_is_constrained_map,edge(pt,g)) )
    {
      *edges_to_erase_ptr++=pt;
    }
    else
    {
      *edges_to_erase_ptr++=tq;
    }
  }

  // If the bottom facet exists, we need to choose one out of the two edges which one disappears:
  //   q-b if it is not constrained and b-p otherwise
  if ( lBottomFaceExists )
  {
    if ( !get(Edge_is_constrained_map,edge(qb,g)) )
    {
      *edges_to_erase_ptr++=qb;
    }
    else{
      *edges_to_erase_ptr++=bp;
    }
  }

  if (lTopFaceExists && lBottomFaceExists)
  {
    if ( face(edges_to_erase[0],g) == face(edges_to_erase[1],g)
         && (! is_border(edges_to_erase[0],g)) )
    {
      // the vertex is of valence 3 and we simply need to remove the vertex
      // and its indicent edges
      bool lP_Erased = false;
      halfedge_descriptor edge =
        next(edges_to_erase[0],g) == edges_to_erase[1]?
          edges_to_erase[0]:edges_to_erase[1];
      if (target(edge,g) == p)
        lP_Erased = true;
      remove_center_vertex(edge,g);
      return lP_Erased? q : p;
    }
    else
    {
      if (!(is_border(edges_to_erase[0],g)))
        join_face(edges_to_erase[0],g);
      else
        remove_face(opposite(edges_to_erase[0],g),g);
      if (!is_border(edges_to_erase[1],g))
        join_face(edges_to_erase[1],g);
      else
        remove_face(opposite(edges_to_erase[1],g),g);
      join_vertex(pq,g);
      return q;
    }
  }
  else
  {
      if (lTopFaceExists)
      {
        if (!(is_border(edges_to_erase[0],g))){
          join_face(edges_to_erase[0],g);
          join_vertex(pq,g);
          return q;
        }
        bool lQ_Erased = is_border(opposite(next(pq,g),g),g);
        remove_face(opposite(edges_to_erase[0],g),g);
        return lQ_Erased?p:q;
      }

      if (! (is_border(edges_to_erase[0],g))){
        join_face(edges_to_erase[0],g);
        join_vertex(qp,g);
        return p;
      }
      bool lP_Erased= is_border(opposite(next(qp,g),g),g);
      remove_face(opposite(edges_to_erase[0],g),g);
      return lP_Erased?q:p;
  };
}

/// performs an edge flip, rotating the edge pointed by
/// `h` by one vertex in the direction of the face orientation.
/// \pre Both faces incident to `h` are triangles.
template<typename Graph>
void
flip_edge(typename boost::graph_traits<Graph>::halfedge_descriptor h,
          Graph& g)
{
  typedef boost::graph_traits<Graph> Traits;
  typedef typename Traits::vertex_descriptor   vertex_descriptor;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Traits::face_descriptor     face_descriptor;

  vertex_descriptor s = source(h,g); 
  vertex_descriptor t = target(h,g); 
  halfedge_descriptor nh = next(h,g), nnh = next(nh,g), oh = opposite(h,g), noh = next(oh,g), nnoh = next(noh,g);
  vertex_descriptor s2 = target(nh,g), t2 = target(noh,g);
  face_descriptor fh = face(h,g), foh = face(oh,g);

  CGAL_assertion(fh != Traits::null_face() && foh != Traits::null_face());

  if(halfedge(s,g) == oh){
    set_halfedge(s,nnh,g);
  }
  if(halfedge(t,g) == h){
    set_halfedge(t,nnoh,g);
  }
  set_next(h,nnoh,g);
  set_next(oh,nnh,g);
  set_target(h,t2,g);
  set_target(oh,s2,g);
  set_next(nh,h,g);
  set_next(noh,oh,g);
  set_next(nnoh,nh,g);
  set_next(nnh,noh,g);
  set_face(nnoh,fh,g);
  set_face(nnh,foh,g);
  set_halfedge(fh,h,g);
  set_halfedge(foh,oh,g);
}

/**
 *  \returns `true` if `e` satisfies the *link condition* \cgalCite{degn-tpec-98}, which guarantees that the surface is also 2-manifold after the edge collapse.
 */
  template<typename Graph>
bool
  does_satisfy_link_condition(typename boost::graph_traits<Graph>::edge_descriptor e,
                           Graph& g)
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
    typedef CGAL::Halfedge_around_source_iterator<Graph> out_edge_iterator;

    halfedge_descriptor v0_v1 = halfedge(e,g);
    halfedge_descriptor v1_v0 = opposite(v0_v1,g);
    
    vertex_descriptor v0 = target(v1_v0,g), v1 = target(v0_v1,g);

    vertex_descriptor vL = target(next(v0_v1,g),g);
    vertex_descriptor vR = target(next(v1_v0,g),g);

    out_edge_iterator eb1, ee1 ; 
    out_edge_iterator eb2, ee2 ; 


  // The following loop checks the link condition for v0_v1.
  // Specifically, that for every vertex 'k' adjacent to both 'p and 'q', 'pkq' is a face of the mesh.
  // 
  for ( boost::tie(eb1,ee1) = halfedges_around_source(v0,g) ;  eb1 != ee1 ; ++ eb1 )
  {
    halfedge_descriptor v0_k = *eb1;
    
    if ( v0_k != v0_v1 )
    {
      vertex_descriptor k = target(v0_k,g);
      
      for ( boost::tie(eb2,ee2) =  halfedges_around_source(k,g) ; eb2 != ee2 ; ++ eb2 )
      {
        halfedge_descriptor k_v1 = *eb2;

        if ( target(k_v1,g) == v1 )
        {
          // At this point we know p-q-k are connected and we need to determine if this triangle is a face of the mesh.
          //
          // Since the mesh is known to be triangular there are at most two faces sharing the edge p-q.
          //
          // If p->q is NOT a border edge, the top face is p->q->t where t is target(next(p->q))
          // If q->p is NOT a border edge, the bottom face is q->p->b where b is target(next(q->p))
          //
          // If k is either t or b then p-q-k *might* be a face of the mesh. It won't be if k==t but p->q is border
          // or k==b but q->b is a border (because in that case even though there exists triangles p->q->t (or q->p->b)
          // they are holes, not faces)
          // 
     
          bool lIsFace =   ( vL == k && (! is_border(v0_v1,g)) )
            || ( vR == k && (! is_border(v1_v0,g)) ) ;
                        
         

          if ( !lIsFace )
          {
            // CGAL_ECMS_TRACE(3,"  k=V" << get(Vertex_index_map,k) << " IS NOT in a face with p-q. NON-COLLAPSABLE edge." ) ;
            return false ;
          }  
          else 
          {
            //CGAL_ECMS_TRACE(4,"  k=V" << get(Vertex_index_map,k) << " is in a face with p-q") ;
          }
        }
      }  
    }
  }   
     
  
    if ( is_border(v0_v1,g) )
    {
      if ( next(next(next(v0_v1,g),g),g) == v0_v1 )
      {
        //CGAL_ECMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge." ) ;
        return false ;
      }
    }
    else if ( is_border(v1_v0,g) )
    {
      if ( next(next(next(v1_v0,g),g),g) == v1_v0 )
      {
        //CGAL_ECMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge." ) ;
        return false ;
      }
    }
    else
    {
      if ( is_border(v0,g) && is_border(v1,g) )
      {
        //CGAL_ECMS_TRACE(3,"  both p and q are boundary vertices but p-q is not. NON-COLLAPSABLE edge." ) ;
        return false ;
      }  
      else
      {
        if ( is_tetrahedron(v0_v1,g) )
        {
          //CGAL_ECMS_TRACE(3,"  p-q belongs to a tetrahedron. NON-COLLAPSABLE edge." ) ;
          return false ;
        }
      }
    }

  
  return true ;
}

#ifndef CGAL_NO_DEPRECATED_CODE
/// \cond SKIP_IN_MANUAL
template<typename Graph>
bool
  satisfies_link_condition(typename boost::graph_traits<Graph>::edge_descriptor e,
                           Graph& g)
{
  return does_satisfy_link_condition(e, g);
}
/// \endcond
#endif
/// @}

} // CGAL

} // CGAL


#endif /* CGAL_EULER_OPERATIONS_H */
