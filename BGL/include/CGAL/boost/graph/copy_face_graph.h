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
#include <CGAL/property_map.h>
#include <boost/unordered_map.hpp>

namespace CGAL {

namespace internal {

template <typename SourceMesh, typename TargetMesh,
          typename Hmap,
          typename V2V, typename H2H, typename F2F,
          typename Src_vpm, typename Tgt_vpm>
void copy_face_graph_impl(const SourceMesh& sm, TargetMesh& tm,
                          Hmap hmap,
                          V2V v2v, H2H h2h, F2F f2f,
                          Src_vpm sm_vpm, Tgt_vpm tm_vpm )
{
  typedef typename boost::graph_traits<SourceMesh>::vertex_descriptor sm_vertex_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::vertex_descriptor tm_vertex_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::vertex_iterator   tm_vertex_iterator;

  typedef typename boost::graph_traits<SourceMesh>::face_descriptor sm_face_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::face_descriptor tm_face_descriptor;

  typedef typename boost::graph_traits<SourceMesh>::halfedge_descriptor sm_halfedge_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::halfedge_descriptor tm_halfedge_descriptor;

  typedef typename boost::graph_traits<SourceMesh>::edge_descriptor sm_edge_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::edge_descriptor tm_edge_descriptor;

  Cartesian_converter<typename Kernel_traits<typename boost::property_traits<Src_vpm>::value_type>::type,
                      typename Kernel_traits<typename boost::property_traits<Tgt_vpm>::value_type>::type >
    conv;

  std::vector<tm_halfedge_descriptor> tm_border_halfedges;
  std::vector<sm_halfedge_descriptor> sm_border_halfedges;

  tm_face_descriptor tm_null_face = boost::graph_traits<TargetMesh>::null_face();

  //insert halfedges and create each vertex when encountering its halfedge
  BOOST_FOREACH(sm_edge_descriptor sm_e, edges(sm))
  {
    tm_edge_descriptor tm_e = add_edge(tm);
    sm_halfedge_descriptor sm_h = halfedge(sm_e, sm), sm_h_opp = opposite(sm_h, sm);
    tm_halfedge_descriptor tm_h = halfedge(tm_e, tm), tm_h_opp = opposite(tm_h, tm);

    put(hmap, sm_h, tm_h);
    put(hmap, sm_h_opp, tm_h_opp);
    *h2h++=std::make_pair(sm_h, tm_h);
    *h2h++=std::make_pair(sm_h_opp, tm_h_opp);

    if ( is_border(sm_h, sm) ){
      tm_border_halfedges.push_back( tm_h );
      sm_border_halfedges.push_back( sm_h );
      set_face(tm_h, tm_null_face, tm);
      CGAL_assertion(next(tm_h, tm) == boost::graph_traits<TargetMesh>::null_halfedge() );
    }

    if( is_border(sm_h_opp, sm) ){
      tm_border_halfedges.push_back( tm_h_opp );
      sm_border_halfedges.push_back( sm_h_opp );
      set_face(tm_h_opp, tm_null_face, tm);
      CGAL_assertion(next(tm_h_opp, tm) == boost::graph_traits<TargetMesh>::null_halfedge() );
    }

    //create a copy of interior vertices only once
    sm_vertex_descriptor sm_h_src = source(sm_h,sm), sm_h_tgt = target(sm_h,sm);
    if ( halfedge(sm_h_tgt,sm)==sm_h )
    {
      tm_vertex_descriptor tm_h_tgt = add_vertex(tm);
      *v2v++=std::make_pair(sm_h_tgt, tm_h_tgt);
      set_halfedge(tm_h_tgt, tm_h, tm);
      set_target(tm_h, tm_h_tgt, tm);
      put(tm_vpm, tm_h_tgt, conv(get(sm_vpm, sm_h_tgt)));
    }
    if ( halfedge(sm_h_src,sm)==sm_h_opp )
    {
      tm_vertex_descriptor tm_h_src = add_vertex(tm);
      *v2v++=std::make_pair(sm_h_src, tm_h_src);
      set_halfedge(tm_h_src, tm_h_opp, tm);
      set_target(tm_h_opp, tm_h_src, tm);
      put(tm_vpm, tm_h_src, conv(get(sm_vpm, sm_h_src)));
    }
  }
  //create faces and connect halfedges
  BOOST_FOREACH(sm_face_descriptor sm_f, faces(sm))
  {
    tm_face_descriptor tm_f = add_face(tm);
    *f2f++=std::make_pair(sm_f, tm_f);

    sm_halfedge_descriptor sm_h_i=halfedge(sm_f, sm);
    tm_halfedge_descriptor tm_h_prev = get(hmap, prev(sm_h_i, sm));
    set_halfedge(tm_f, tm_h_prev, tm);

    CGAL_precondition(*halfedges_around_face(sm_h_i, sm).first == sm_h_i);
    BOOST_FOREACH(sm_halfedge_descriptor sm_h, halfedges_around_face(sm_h_i, sm))
    {
      tm_halfedge_descriptor tm_h = get(hmap, sm_h);
      set_next(tm_h_prev, tm_h, tm);
      set_face(tm_h, tm_f, tm);
      tm_h_prev=tm_h;
    }
  }

  // update next/prev of tm border halfedges
  std::size_t nb_border_hedges = tm_border_halfedges.size();
  for (std::size_t i=0; i< nb_border_hedges; ++i)
  {
    tm_halfedge_descriptor tm_h = tm_border_halfedges[i];

    if ( next(tm_h, tm) != boost::graph_traits<TargetMesh>::null_halfedge() )
      continue; //already set

    tm_halfedge_descriptor tm_h_prev = tm_h;
    CGAL_precondition(*halfedges_around_face(sm_border_halfedges[i], sm).first == sm_border_halfedges[i]);
    BOOST_FOREACH(sm_halfedge_descriptor sm_h,
                  halfedges_around_face(next(sm_border_halfedges[i], sm), sm))
    {
      CGAL_assertion(next(tm_h_prev, tm) == boost::graph_traits<TargetMesh>::null_halfedge());
      tm_h = get(hmap, sm_h);
      set_next(tm_h_prev, tm_h, tm);
      tm_h_prev=tm_h;
    }
  }
  // update halfedge vertex of all but the vertex halfedge
  for(tm_vertex_iterator vit = vertices(tm).first;
      vit != vertices(tm).second; ++vit)
  {
    tm_vertex_descriptor v = *vit;
    tm_halfedge_descriptor h = halfedge(v, tm);
    tm_halfedge_descriptor next_around_vertex=h;
    do{
      next_around_vertex=opposite(next(next_around_vertex, tm), tm);
      set_target(next_around_vertex, v, tm);
    }while(h != next_around_vertex);
  }
}

template <typename SourceMesh, typename TargetMesh,
          typename V2V, typename H2H, typename F2F,
          typename Src_vpm, typename Tgt_vpm>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm,
                     Tag_false,
                     V2V v2v, H2H h2h, F2F f2f,
                     Src_vpm sm_vpm, Tgt_vpm tm_vpm )
{
  typedef typename boost::graph_traits<SourceMesh>::halfedge_descriptor sm_halfedge_descriptor;
  typedef typename boost::graph_traits<TargetMesh>::halfedge_descriptor tm_halfedge_descriptor;

  boost::unordered_map<sm_halfedge_descriptor,
                       tm_halfedge_descriptor> hash_map(num_halfedges(sm));
  copy_face_graph_impl(sm, tm,
                       boost::make_assoc_property_map(hash_map),
                       v2v, h2h, f2f,
                       sm_vpm, tm_vpm);
}

template <typename SourceMesh, typename TargetMesh,
          typename V2V, typename H2H, typename F2F,
          typename Src_vpm, typename Tgt_vpm>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm,
                     Tag_true,
                     V2V v2v, H2H h2h, F2F f2f,
                     Src_vpm sm_vpm, Tgt_vpm tm_vpm )
{
  typedef typename boost::graph_traits<TargetMesh>::halfedge_descriptor tm_halfedge_descriptor;
  std::vector<tm_halfedge_descriptor> hedges(num_halfedges(sm));

  // init halfedge index map
  /// \TODO shall we keep that?
  helpers::init_halfedge_indices(const_cast<SourceMesh&>(sm),
                                 get(boost::halfedge_index, sm));

  copy_face_graph_impl(sm, tm,
                       bind_property_maps(get(boost::halfedge_index, sm),
                                          make_property_map(hedges)),
                       v2v, h2h, f2f,
                       sm_vpm, tm_vpm);
}

} // end of namespace internal

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
  \tparam Src_vpm a model of `ReadablePropertyMap` with `sm_vertex_descriptor` as key
  \tparam Tgt_vpm a model of `WritablePropertyMap` with `tm_vertex_descriptor` as key
  where the prefix `sm_` and `tm_` mean belonging to the source and
  target mesh respectively.

  The types `sm_vertex_descriptor` and `sm_face_descriptor` must be models of the concept `Hashable`.

  \param sm the source mesh
  \param tm the target mesh
  \param v2v pairs of `vertex_descriptors` from `sm` and corresponding `vertex_descriptors` in `tm` are added to `v2v`
  \param h2h pairs of `halfedge_descriptors` from `sm` and corresponding `halfedge_descriptors` in `tm` are added to `h2h`
  \param f2f pairs of `face_descriptors` from `sm` and corresponding `face_descriptors` in `tm` are added to `f2f`
  \param sm_vpm vertex point map for `sm`
  \param tm_vpm vertex point map for `tm`

  The points from `sm` to `tm` are converted using
  `CGAL::Cartesian_converter<SourceKernel, TargetKernel>`.
  `SourceKernel` and `TargetKernel` are deduced using `CGAL::Kernel_traits`
  from the value types of `Src_vpm` and `Tgt_vpm`.

  Other properties are not copied.
*/
#if defined(DOXYGEN_RUNNING) // Use template default arguments
template <typename SourceMesh, typename TargetMesh,
          typename V2V = Emptyset_iterator,
          typename H2H = Emptyset_iterator,
          typename F2F = Emptyset_iterator,
          typename Src_vpm = typename boost::property_map<SourceMesh, vertex_point_t>::const_type,
          typename Tgt_vpm = typename boost::property_map<TargetMesh, vertex_point_t>::type>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm,
                     V2V v2v = V2V(), H2H h2h = H2H(), F2F f2f = F2F(),
                     Src_vpm sm_vpm = get(vertex_point, sm),
                     Tgt_vpm tm_vpm = get(vertex_point, tm) )
#else // use the overloads
template <typename SourceMesh, typename TargetMesh,
          typename V2V, typename H2H, typename F2F,
          typename Src_vpm, typename Tgt_vpm>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm,
                     V2V v2v, H2H h2h, F2F f2f,
                     Src_vpm sm_vpm, Tgt_vpm tm_vpm )
#endif
{
  internal::copy_face_graph(sm, tm,
                            boost::graph_has_property<SourceMesh,boost::halfedge_index_t>(),
                            v2v, h2h, f2f,
                            sm_vpm, tm_vpm);
}

#if !defined(DOXYGEN_RUNNING)
template <typename SourceMesh, typename TargetMesh>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm)
{ copy_face_graph(sm, tm, Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename TargetMesh, typename V2V>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm, V2V v2v)
{ copy_face_graph(sm, tm, v2v, Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename TargetMesh, typename V2V, typename H2H>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm, V2V v2v, H2H h2h)
{ copy_face_graph(sm, tm, v2v, h2h, Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename TargetMesh, typename V2V, typename H2H, typename F2F>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f)
{ copy_face_graph(sm, tm, v2v, h2h, f2f,
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename TargetMesh, typename V2V, typename H2H, typename F2F, typename Src_vpm>
void copy_face_graph(const SourceMesh& sm, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f, Src_vpm sm_vpm)
{ copy_face_graph(sm, tm, v2v, h2h, f2f,
                  sm_vpm, get(vertex_point, tm)); }
#endif

} // namespace CGAL

#endif //  CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_H
