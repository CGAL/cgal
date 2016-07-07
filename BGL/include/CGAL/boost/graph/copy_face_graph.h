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

#include <CGAL/config.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>

#include <boost/unordered_map.hpp>

namespace CGAL {

/*!
  \ingroup PkgBGLHelperFct

  copies a source model of `FaceListGraph` into a target model of a
  `FaceListGraph`. `OutputIterators` can be provided to produce a
  mapping between source and target elements. The target graph is not
  cleared.

  \tparam SourceMesh a model of `FaceListGraph`. 
          The descriptor types `boost::graph_traits<SourceMesh>::%vertex_descriptor`
          and `boost::graph_traits<SourceMesh>::%face_descriptor` must be
          models of `Hashable`.
  \tparam TargetMesh a model of `FaceListGraph`
  \tparam V2V a model of `OutputIterator` accepting `std::pair<sm_vertex_descriptor, tm_vertex_descriptor>`
  \tparam H2H a model of `OutputIterator` accepting `std::pair<sm_halfedge_descriptor, tm_halfedge_descriptor>`
  \tparam F2F a model of `OutputIterator` accepting `std::pair<sm_face_descriptor, tm_face_descriptor>`

  where the prefix `sm_` and `tm_` mean belonging to the source and
  target mesh respectively.

  The types `sm_vertex_descriptor` and `sm_face_descriptor` must be models of the concept `Hashable`.

  \param sm the source mesh
  \param tm the target mesh
  \param v2v pairs of `vertex_descriptors` from `sm` and corresponding `vertex_descriptors` in `tm` are added to `v2v`
  \param h2h pairs of `halfedge_descriptors` from `sm` and corresponding `halfedge_descriptors` in `tm` are added to `h2h`
  \param f2f pairs of `face_descriptors` from `sm` and corresponding `face_descriptors` in `tm` are added to `f2f`

  This function assumes that both graphs have an internal property
  `vertex_point`. Values of that property are converted using
  `CGAL::Cartesian_converter<SourceKernel,
  TargetKernel>`. `SourceKernel` and `TargetKernel` are deduced using
  `CGAL::Kernel_traits`.

  Other properties are not copied.
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

  Cartesian_converter<typename Kernel_traits<typename boost::property_traits<sm_PMap>::value_type>::type, 
                      typename Kernel_traits<typename boost::property_traits<tm_PMap>::value_type>::type > 
    conv;

  sm_PMap sm_pmap = get(vertex_point, sm);
  tm_PMap tm_pmap = get(vertex_point, tm);

  // internal f2f and v2v
  boost::unordered_map<sm_vertex_descriptor, tm_vertex_descriptor> v2v_;
  boost::unordered_map<sm_face_descriptor, tm_face_descriptor> f2f_;

  BOOST_FOREACH(sm_vertex_descriptor svd, vertices(sm)){
    tm_vertex_descriptor tvd = add_vertex(tm);
    v2v_[svd] = tvd;
    *v2v++ = std::make_pair(svd, tvd);
    put(tm_pmap, tvd, conv(get(sm_pmap, svd)));
  }

  BOOST_FOREACH(sm_face_descriptor sfd, faces(sm)){
    std::vector<tm_vertex_descriptor> tv;
    BOOST_FOREACH(sm_vertex_descriptor svd, vertices_around_face(halfedge(sfd,sm),sm)){
      tv.push_back(v2v_.at(svd));
    }
    tm_face_descriptor new_face = Euler::add_face(tv,tm);
    f2f_[sfd] = new_face;
    *f2f++ = std::make_pair(sfd, new_face);
  }
  
  BOOST_FOREACH(sm_face_descriptor sfd, faces(sm)){
    sm_halfedge_descriptor shd = halfedge(sfd,sm), done(shd);
    tm_halfedge_descriptor thd = halfedge(f2f_[sfd],tm);
    tm_vertex_descriptor tvd = v2v_.at(target(shd,sm));
    while(target(thd,tm) != tvd){
      thd = next(thd,tm);
    }
    do {
      *h2h++ = std::make_pair(shd, thd);
      if (face(opposite(shd, sm), sm) == boost::graph_traits<SourceMesh>::null_face()){
        *h2h++  = std::make_pair(opposite(shd, sm), opposite(thd, tm));
      }
      shd = next(shd,sm);
      thd = next(thd,tm);
    }while(shd != done);
  }
  
}

#if !defined(CGAL_CXX11)  && !defined(DOXYGEN_RUNNING)
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
