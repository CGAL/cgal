// Copyright (c) 2024 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BGL_IO_OM_H
#define CGAL_BGL_IO_OM_H

#if defined(CGAL_USE_OPENMESH) || defined(DOXYGEN_RUNNING)

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/property_map.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <fstream>
#include <map>

namespace CGAL {
namespace IO {

namespace internal {
template <typename Graph, typename VPM, typename VFeaturePM, typename EFeaturePM>
bool read_OM(const std::string& fname, Graph& g, VPM vpm, VFeaturePM vfpm, EFeaturePM efpm)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> OMesh;
  typedef typename boost::graph_traits<OMesh>::vertex_descriptor om_vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<OMesh>::halfedge_descriptor om_halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  OMesh omesh;
  OpenMesh::IO::Options options = OpenMesh::IO::Options::Status;
  bool ok = OpenMesh::IO::read_mesh(omesh, fname, options);
  if(! ok){
    return false;
  }

  std::map<om_vertex_descriptor, vertex_descriptor> v2v;
  auto v2vpmap = boost::make_assoc_property_map(v2v);

  std::map<om_halfedge_descriptor, halfedge_descriptor> h2h;
  auto h2hpmap = boost::make_assoc_property_map(h2h);

  CGAL::copy_face_graph<OMesh,Graph>(omesh, g,
                                     CGAL::parameters::vertex_to_vertex_map(v2vpmap)
                                                      .halfedge_to_halfedge_map(h2hpmap),
                                     CGAL::parameters::vertex_point_map(vpm));
  if(options.vertex_has_status()){
    for(auto v : vertices(omesh)){
      put(vfpm, v2v[v], omesh.status(v).feature());
    }
  }

  if(options.edge_has_status()){
    for(auto e : edges(omesh)){
      auto sme = edge(h2h[halfedge(e,omesh)], g);
      put(efpm, sme , omesh.status(OpenMesh::EdgeHandle(e.idx())).feature());
    }
  }
  return true;
}

template <typename Graph, typename VPM, typename VFeaturePM, typename EFeaturePM>
bool write_OM(std::string fname, Graph& g, VPM vpm, VFeaturePM vfpm, EFeaturePM efpm)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> OMesh;
  typedef typename boost::graph_traits<OMesh>::vertex_descriptor om_vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<OMesh>::halfedge_descriptor om_halfedge_descriptor;
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

  std::map<vertex_descriptor, om_vertex_descriptor> v2v;
  auto v2vpmap = boost::make_assoc_property_map(v2v);

  std::map<halfedge_descriptor, om_halfedge_descriptor> h2h;
  auto h2hpmap = boost::make_assoc_property_map(h2h);

  OMesh omesh;
  CGAL::copy_face_graph<Graph, OMesh>(g, omesh,
                                      CGAL::parameters::vertex_point_map(vpm)
                                                       .vertex_to_vertex_map(v2vpmap)
                                                       .halfedge_to_halfedge_map(h2hpmap));
  omesh.request_edge_status();
  omesh.request_vertex_status();

  for (auto h : halfedges(g))
  {
    om_halfedge_descriptor omh = h2h.at(h);
    const bool isfeature = get(efpm, edge(h, g));
    omesh.status(omesh.edge_handle(omh)).set_feature(isfeature);
  }
  for (auto v : vertices(g))
  {
    auto omv = v2v.at(v);
    const bool isfeature = get(vfpm, v);
    omesh.status(omv).set_feature(isfeature);
  }

  return OpenMesh::IO::write_mesh(omesh, fname, OpenMesh::IO::Options::Status);
}
} // end of internal namespace

/*!
  \ingroup PkgBGLIoFuncsOM

  \brief reads the graph `g` from the file `fname`, using the \ref IOStreamOM.

  The data is expected to represent a 2-manifold (possibly with borders).

  \attention The graph `g` is not cleared, and the data from the file are appended.

  \note This function is only available if OpenMesh is available (`CGAL_USE_OPENMESH` is defined or CMake target is linked with `CGAL::OpenMesh_support`).

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the input file
  \param g the graph to be built from the input data
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd
    \cgalParamNBegin{edge_is_constrained_map}
      \cgalParamDescription{a property map containing the feature-or-not status of each edge of `g` to be filled by the reader}
      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%edge_descriptor`
                     as key type and `bool` as value type.}
      \cgalParamDefault{a default property map where no edge is marked as feature}
    \cgalParamNEnd
    \cgalParamNBegin{vertex_is_constrained_map}
      \cgalParamDescription{a property map containing the feature-or-not status of each vertex of `g` to be filled by the reader}
      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `bool` as value type.}
      \cgalParamDefault{a default property map where no vertex is marked as feature}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if reading was successful and the resulting mesh is valid, `false` otherwise.
*/
template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OM(const std::string& fname,
             Graph& g,
             const CGAL_NP_CLASS& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;

  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;
  using Default_vfmap = Static_boolean_property_map<vertex_descriptor, false>;
  using Default_efmap = Static_boolean_property_map<edge_descriptor, false>;
  auto vfpm = choose_parameter<Default_vfmap>(get_parameter(np, internal_np::vertex_is_constrained));
  auto efpm = choose_parameter<Default_efmap>(get_parameter(np, internal_np::edge_is_constrained));
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_property_map(vertex_point, g));
  return internal::read_OM(fname, g, vpm, vfpm, efpm);
}


/*!
\ingroup PkgBGLIoFuncsOM

  \brief writes the graph `g` into a file named `fname`, using the \ref IOStreamOM.

  \note This function is only available if OpenMesh is available (`CGAL_USE_OPENMESH` is defined or CMake target is linked with `CGAL::OpenMesh_support`).

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the output file
  \param g the graph to be written
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd
    \cgalParamNBegin{edge_is_constrained_map}
      \cgalParamDescription{a property map containing the feature-or-not status of each edge of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%edge_descriptor`
                     as key type and `bool` as value type.}
    \cgalParamNEnd
    \cgalParamNBegin{vertex_is_constrained_map}
      \cgalParamDescription{a property map containing the feature-or-not status of each vertex of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `bool` as value type.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful, `false` otherwise.
*/
template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OM(const std::string& fname,
              const Graph& g,
              const CGAL_NP_CLASS& np = parameters::default_values())
{
  using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;

  using CGAL::parameters::get_parameter;
  using CGAL::parameters::choose_parameter;
  auto vfpm = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                               CGAL::Constant_property_map<vertex_descriptor, bool>(false));
  auto efpm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                               CGAL::Constant_property_map<edge_descriptor, bool>(false));
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, g));
  return internal::write_OM(fname, g, vpm, vfpm, efpm);
}


} // namespace IO
} // namespace CGAL

#endif // CGAL_USE_OPENMESH || DOXYGEN_RUNNING
#endif // CGAL_IO_OM
