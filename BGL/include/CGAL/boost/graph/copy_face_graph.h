// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_H
#define CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_H

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <boost/unordered_map.hpp>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/iterator.h>

namespace CGAL {

/*!
  \ingroup PkgBGLHelperFct

  copies a source FaceGraph into another FaceGraph of different
  type. OutputIterators can be provided to produce a mapping between
  source and target elements.

  \tparam SourceMesh a `FaceListGraph`
  \tparam TargetMesh a `FaceListGraph`
  \tparam V2V an `OutputIterator` accepting `std::pair<sm_vertex_descriptor, tm_vertex_descriptor>`
  \tparam H2H an `OutputIterator` accepting `std::pair<sm_halfedge_descriptor, tm_halfedge_descriptor>`
  \tparam F2F an `OutputIterator` accepting `std::pair<sm_face_descriptor, tm_face_descriptor>`

  where the prefixx `sm_` and `tm_` mean belonging to the source or
  target mesh respectively.

  \param sm the source mesh of the copy operation
  \param tm the target mesh of the copy operation
  \param v2v pairs of `vertex_descriptorS` from `sm` and corresponding `vertex_descriptorS` in `tm` are added to `v2v`
  \param h2h pairs of `halfedge_descriptorS` from `sm` and corresponding `halfedge_descriptorS` in `tm` are added to `h2h`
  \param f2f pairs of `face_descriptorS` from `sm` and corresponding `face_descriptorS` in `tm` are added to `f2f`


  This function assumes that both graphs have an internal property
  `vertex_point` and that the `vertex_point` values of `tm` are
  constructible from the ones of `sm`.
*/
#if defined(CGAL_CXX11) || defined(DOXYGEN_RUNNING) // Use template default arguments
template <typename SourceMesh, typename TargetMesh,
          typename V2V = Emptyset_iterator,
          typename H2H = Emptyset_iterator,
          typename F2F = Emptyset_iterator>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm,
                     V2V v2v = V2V(), H2H h2h = H2H(), F2F f2f = F2F())
#else // use the overloads
template <typename SourceMesh, typename TargetMesh,
          typename V2V, typename H2H, typename F2F>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm,
                     V2V v2v, H2H h2h, F2F f2f)
#endif
{
  typedef typename boost::graph_traits<SourceMesh>::vertex_descriptor sm_vertex_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::vertex_descriptor tm_vertex_descriptor;

  typedef typename boost::graph_traits<SourceMesh>::face_descriptor sm_face_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::face_descriptor tm_face_descriptor;

  typedef typename boost::graph_traits<SourceMesh>::halfedge_descriptor sm_halfedge_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::halfedge_descriptor tm_halfedge_descriptor;

  typedef typename boost::property_map<SourceMesh, vertex_point_t>::const_type sm_PMap;
  typedef typename boost::property_map<TargetMesh, vertex_point_t>::type tm_PMap;

  sm_PMap sm_pmap = get(vertex_point, sm);
  tm_PMap tm_pmap = get(vertex_point, tm);


  BOOST_FOREACH(sm_vertex_descriptor svd, vertices(sm)){
    tm_vertex_descriptor tvd = add_vertex(tm);
    *v2v++ = std::make_pair(svd, tvd);
    put(tm_pmap, tvd, get(sm_pmap, svd));
  }

  // internal f2f
  boost::unordered_map<sm_face_descriptor, tm_face_descriptor> f2f_;
  BOOST_FOREACH(sm_face_descriptor sfd, faces(sm)){
    std::vector<tm_vertex_descriptor> tv;
    BOOST_FOREACH(sm_vertex_descriptor svd, vertices_around_face(halfedge(sfd,sm),sm)){
      tv.push_back(v2v.at(svd));
    }
    tm_face_descriptor new_face = Euler::add_face(tv,tm);
    f2f_[sfd] = new_face;
    *f2f++ = std::make_pair(sfd, new_face);
  }
  
  BOOST_FOREACH(sm_face_descriptor sfd, faces(sm)){
    sm_halfedge_descriptor shd = halfedge(sfd,sm), done(shd);
    tm_halfedge_descriptor thd = halfedge(f2f_[sfd],tm);
    tm_vertex_descriptor tvd = v2v.at(target(shd,sm));
    while(target(thd,tm) != tvd){
      thd = next(thd,tm);
    }
    do {
      *h2h++ = std::make_pair(shd, thd);

      if (face(opposite(shd, sm), sm) == boost::graph_traits<SourceMesh>::null_face())
        h2h.insert(std::make_pair(opposite(shd, sm), opposite(thd, tm)));

      shd = next(shd,sm);
      thd = next(thd,tm);
    }while(shd != done);
  }
  
}

#if !defined(CGAL_CXX11)
template <typename SourceMesh, typename TargetMesh>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm)
{ copy_face_graph(sm, tm, Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator()); }

template <typename SourceMesh, typename TargetMesh, typename V2V>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm, V2V v2v)
{ copy_face_graph(sm, tm, v2v, Emptyset_iterator(), Emptyset_iterator()); }

template <typename SourceMesh, typename TargetMesh, typename V2V, typename H2H>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm, V2V v2v, H2H h2h)
{ copy_face_graph(sm, tm, v2v, h2h, Emptyset_iterator()); }
#endif

} // namespace CGAL

#endif //  CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_H
