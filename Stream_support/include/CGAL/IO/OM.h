// Copyright (c) 2024 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_IO_OM_H
#define CGAL_IO_OM_H

#if defined(CGAL_USE_OPENMESH) || defined(DOXYGEN_RUNNING)

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/property_map.h>
#include <iostream>
#include <fstream>
#include <map>

namespace CGAL {
namespace IO {

template <typename SM, typename VFeaturePM, typename EFeaturePM>
bool read_OM(std::string fname, SM& sm, VFeaturePM vfpm, EFeaturePM efpm)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> OMesh;
  typedef typename boost::graph_traits<OMesh>::vertex_descriptor om_vertex_descriptor;
  typedef typename boost::graph_traits<SM>::vertex_descriptor sm_vertex_descriptor;
  typedef typename boost::graph_traits<OMesh>::halfedge_descriptor om_halfedge_descriptor;
  typedef typename boost::graph_traits<SM>::halfedge_descriptor sm_halfedge_descriptor;

  OMesh omesh;
  OpenMesh::IO::Options options = OpenMesh::IO::Options::Status;
  bool ok = OpenMesh::IO::read_mesh(omesh, fname, options);
  if(! ok){
    return false;
  }

  std::map<om_vertex_descriptor,sm_vertex_descriptor> v2v;
  auto v2vpmap = boost::make_assoc_property_map(v2v);

  std::map<om_halfedge_descriptor,sm_halfedge_descriptor> h2h;
  auto h2hpmap = boost::make_assoc_property_map(h2h);

  CGAL::copy_face_graph<OMesh,SM>(omesh, sm, CGAL::parameters::vertex_to_vertex_map(v2vpmap).halfedge_to_halfedge_map(h2hpmap));

  if(options.vertex_has_status()){
    for(auto v : vertices(omesh)){
        put(vfpm, v2v[v], omesh.status(v).feature());
    }
  }

  if(options.edge_has_status()){
    for(auto e : edges(omesh)){
        auto sme = edge(h2h[halfedge(e,omesh)], sm);
        put(efpm, sme , omesh.status(OpenMesh::EdgeHandle(e.idx())).feature());
    }
  }
  return true;
}

template <typename SM, typename VFeaturePM, typename EFeaturePM>
bool write_OM(std::string fname, SM& sm, VFeaturePM vfpm, EFeaturePM efpm)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> OMesh;
  typedef typename boost::graph_traits<OMesh>::vertex_descriptor om_vertex_descriptor;
  typedef typename boost::graph_traits<SM>::vertex_descriptor sm_vertex_descriptor;
  typedef typename boost::graph_traits<OMesh>::halfedge_descriptor om_halfedge_descriptor;
  typedef typename boost::graph_traits<SM>::halfedge_descriptor sm_halfedge_descriptor;

  std::map<sm_vertex_descriptor, om_vertex_descriptor> v2v;
  auto v2vpmap = boost::make_assoc_property_map(v2v);

  std::map<sm_halfedge_descriptor, om_halfedge_descriptor> h2h;
  auto h2hpmap = boost::make_assoc_property_map(h2h);

  OMesh omesh;
  CGAL::copy_face_graph<SM, OMesh>(sm, omesh,
    CGAL::parameters::vertex_to_vertex_map(v2vpmap).halfedge_to_halfedge_map(h2hpmap));

  omesh.request_edge_status();
  omesh.request_vertex_status();

  for (auto h : halfedges(sm))
  {
    om_halfedge_descriptor omh = h2h.at(h);
    const bool isfeature = get(efpm, edge(h, sm));
    omesh.status(omesh.edge_handle(omh)).set_feature(isfeature);
  }
  for (auto v : vertices(sm))
  {
    auto omv = v2v.at(v);
    const bool isfeature = get(vfpm, v);
    omesh.status(omv).set_feature(isfeature);
  }

  return OpenMesh::IO::write_mesh(omesh, fname, OpenMesh::IO::Options::Status);
}

} // namespace IO
} // namespace CGAL

#endif // CGAL_USE_OPENMESH || DOXYGEN_RUNNING
#endif // CGAL_IO_OM
