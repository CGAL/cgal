// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Maxime Gimeno

#ifndef  CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_PATCH_H
#define  CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_PATCH_H

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Connected_component_graph.h>
#include <CGAL/boost/graph/halfedge_graph_traits.h>
#include <CGAL/boost/graph/Connected_component_graph.h>

namespace CGAL {

/*!
  \ingroup PkgBGLHelperFct

  copies a connected component of a source model of `FaceListGraph` into a target model of a
  `FaceListGraph`. `OutputIterators` can be provided to produce a
  mapping between source and target elements. The target graph is not
  cleared.

  \tparam SourceMesh a model of `FaceListGraph`.
          The descriptor types `boost::graph_traits<SourceMesh>::%vertex_descriptor`
          and `boost::graph_traits<SourceMesh>::%face_descriptor` must be
          models of `Hashable`.
  \tparam TargetMesh a model of `FaceListGraph`
  \tparam FaceComponentMap a model of `WritablePropertyMap` with
          `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
          `graph_traits<PolygonMesh>::faces_size_type` as value type.
  \tparam V2V a model of `OutputIterator` accepting `std::pair<sm_vertex_descriptor, tm_vertex_descriptor>`
  \tparam H2H a model of `OutputIterator` accepting `std::pair<sm_halfedge_descriptor, tm_halfedge_descriptor>`
  \tparam F2F a model of `OutputIterator` accepting `std::pair<sm_face_descriptor, tm_face_descriptor>`
  \tparam Src_vpm a model of `ReadablePropertyMap` with `sm_vertex_descriptor` as key
  \tparam Tgt_vpm a model of `WritablePropertyMap` with `tm_vertex_descriptor` as key
  where the prefix `sm_` and `tm_` mean belonging to the source and
  target mesh respectively.

  The types `sm_vertex_descriptor` and `sm_face_descriptor` must be models of the concept `Hashable`.

  \param sm the source mesh
  \param fccmap the property map containing the faces of `sm` and the index of their connected component.
  \param pid the index in `fccmap` of the connected component to copy
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
          typename FaceComponentMap,
          typename V2V = Emptyset_iterator,
          typename H2H = Emptyset_iterator,
          typename F2F = Emptyset_iterator,
          typename Src_vpm = typename boost::property_map<SourceMesh, vertex_point_t>::const_type,
          typename Tgt_vpm = typename boost::property_map<TargetMesh, vertex_point_t>::type>
void copy_face_graph_patch(const SourceMesh& sm,
                     FaceComponentMap fccmap,
                     typename boost::graph_traits<SourceMesh>::faces_size_type pid,
                     TargetMesh& tm,
                     V2V v2v = V2V(), H2H h2h = H2H(), F2F f2f = F2F(),
                     Src_vpm sm_vpm = get(vertex_point, sm),
                     Tgt_vpm tm_vpm = get(vertex_point, tm) )
#else // use the overloads
template <typename SourceMesh, typename TargetMesh,
          typename FaceComponentMap, typename V2V,
          typename H2H, typename F2F,
          typename Src_vpm, typename Tgt_vpm >
void copy_face_graph_patch(const SourceMesh& sm, FaceComponentMap fccmap, typename boost::graph_traits<SourceMesh>::faces_size_type pid, TargetMesh& tm,
                           V2V v2v, H2H h2h, F2F f2f,
                           Src_vpm sm_vpm, Tgt_vpm tm_vpm )
#endif
{
    typedef CGAL::Connected_component_graph<SourceMesh, FaceComponentMap> Adapter;

    copy_face_graph(Adapter(sm, fccmap, pid), tm, v2v, h2h, f2f, sm_vpm, tm_vpm);
}


/*!
  \ingroup PkgBGLHelperFct

  copies the connected component containg `seed_face` of a source model of `FaceListGraph` into a target model of a
  `FaceListGraph`. `OutputIterators` can be provided to produce a
  mapping between source and target elements. The target graph is not
  cleared.

  \tparam SourceMesh a model of `FaceListGraph`.
          The descriptor types `boost::graph_traits<SourceMesh>::%vertex_descriptor`
          and `boost::graph_traits<SourceMesh>::%face_descriptor` must be
          models of `Hashable`.
  \tparam TargetMesh a model of `FaceListGraph`
  \tparam NamedParameters a sequence of \ref namedparameters
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
  \param seed_face a face of `sm` belonging to the connected component to copy
  \param np optional \ref namedparameters described below
*
*  \cgalNamedParamsBegin
*     \cgalParamBegin{edge_is_constrained_map}  a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
*  \cgalNamedParamsEnd

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
          typename NamedParameters,
          typename V2V = Emptyset_iterator,
          typename H2H = Emptyset_iterator,
          typename F2F = Emptyset_iterator,
          typename Src_vpm = typename boost::property_map<SourceMesh, vertex_point_t>::const_type,
          typename Tgt_vpm = typename boost::property_map<TargetMesh, vertex_point_t>::type>
void copy_face_graph_patch(const SourceMesh& sm,
                     const typename boost::graph_traits<SourceMesh>::
                                           face_descriptor& seed_face,
                     const NamedParameters& np = boost::parameter::all_default(),
                     TargetMesh& tm,
                     V2V v2v = V2V(), H2H h2h = H2H(), F2F f2f = F2F(),
                     Src_vpm sm_vpm = get(vertex_point, sm),
                     Tgt_vpm tm_vpm = get(vertex_point, tm) )
#else // use the overloads

template <typename SourceMesh, typename TargetMesh,
          typename NamedParameters,
          typename V2V,typename H2H, typename F2F,
          typename Src_vpm, typename Tgt_vpm >
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face,
                           const NamedParameters& np, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f,
                           Src_vpm sm_vpm, Tgt_vpm tm_vpm )
#endif
{

    typedef typename boost::graph_traits<SourceMesh>::face_descriptor g_face_descriptor;
    typedef typename boost::graph_traits<SourceMesh>::faces_size_type faces_size_t;
    typedef boost::associative_property_map< boost::unordered_map< g_face_descriptor, faces_size_t>  >FCMap;
    typedef CGAL::Connected_component_graph<SourceMesh, FCMap> Adapter;

    boost::unordered_map<g_face_descriptor, faces_size_t> map(CGAL::num_faces(sm));
    Polygon_mesh_processing::connected_components(sm, boost::make_assoc_property_map(map), np);
    copy_face_graph(Adapter(sm, boost::make_assoc_property_map(map), map[seed_face]), tm, v2v, h2h, f2f, sm_vpm, tm_vpm);
}

#if !defined(DOXYGEN_RUNNING)
template <typename SourceMesh, typename FaceComponentMap, typename TargetMesh>
void copy_face_graph_patch(const SourceMesh& sm, FaceComponentMap fccmap, typename boost::graph_traits<SourceMesh>::faces_size_type pid, TargetMesh& tm)
{ copy_face_graph_patch(sm, fccmap, pid, tm, Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename FaceComponentMap, typename TargetMesh, typename V2V>
void copy_face_graph_patch(const SourceMesh& sm, FaceComponentMap fccmap, typename boost::graph_traits<SourceMesh>::faces_size_type pid, TargetMesh& tm, V2V v2v)
{ copy_face_graph_patch(sm, fccmap, pid, tm, v2v, Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename FaceComponentMap, typename TargetMesh, typename V2V, typename H2H>
void copy_face_graph_patch(const SourceMesh& sm, FaceComponentMap fccmap, typename boost::graph_traits<SourceMesh>::faces_size_type pid, TargetMesh& tm, V2V v2v, H2H h2h)
{ copy_face_graph_patch(sm, fccmap, pid, tm, v2v, h2h, Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename FaceComponentMap, typename TargetMesh, typename V2V, typename H2H, typename F2F>
void copy_face_graph_patch(const SourceMesh& sm, FaceComponentMap fccmap, typename boost::graph_traits<SourceMesh>::faces_size_type pid, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f)
{ copy_face_graph_patch(sm, fccmap, pid, tm, v2v, h2h, f2f,
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename FaceComponentMap, typename TargetMesh, typename V2V, typename H2H, typename F2F, typename Src_vpm>
void copy_face_graph_patch(const SourceMesh& sm, FaceComponentMap fccmap, typename boost::graph_traits<SourceMesh>::faces_size_type pid, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f, Src_vpm sm_vpm)
{ copy_face_graph_patch(sm, fccmap, pid, tm, v2v, h2h, f2f,
                  sm_vpm, get(vertex_point, tm)); }


template <typename SourceMesh, typename TargetMesh>
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face, TargetMesh& tm)
{ copy_face_graph_patch(sm, seed_face, Polygon_mesh_processing::parameters::all_default(), tm, Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename NamedParameters, typename TargetMesh>
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face,
                           const NamedParameters& np, TargetMesh& tm)
{ copy_face_graph_patch(sm, seed_face, np, tm, Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename NamedParameters, typename TargetMesh, typename V2V>
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face,
                           const NamedParameters& np, TargetMesh& tm, V2V v2v)
{ copy_face_graph_patch(sm, seed_face, np, tm, v2v, Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename NamedParameters, typename TargetMesh, typename V2V, typename H2H>
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face,
                           const NamedParameters& np, TargetMesh& tm, V2V v2v, H2H h2h)
{ copy_face_graph_patch(sm, seed_face, np, tm, v2v, h2h, Emptyset_iterator(),
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename NamedParameters, typename TargetMesh, typename V2V, typename H2H, typename F2F>
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face,
                           const NamedParameters& np, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f)
{ copy_face_graph_patch(sm, seed_face, np, tm, v2v, h2h, f2f,
                  get(vertex_point, sm), get(vertex_point, tm)); }

template <typename SourceMesh, typename NamedParameters, typename TargetMesh, typename V2V, typename H2H, typename F2F, typename Src_vpm>
void copy_face_graph_patch(const SourceMesh& sm, const typename boost::graph_traits<SourceMesh>::
                           face_descriptor& seed_face,
                           const NamedParameters& np, TargetMesh& tm, V2V v2v, H2H h2h, F2F f2f, Src_vpm sm_vpm)
{ copy_face_graph_patch(sm, seed_face, np, tm, v2v, h2h, f2f,
                  sm_vpm, get(vertex_point, tm)); }
#endif
}//end CGAL
#endif //  CGAL_BOOST_GRAPH_COPY_FACE_GRAPH_PATCH_H
