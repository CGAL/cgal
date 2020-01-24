// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_WRL_H
#define CGAL_BGL_IO_WRL_H

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/container/flat_map.hpp>

#include <fstream>

namespace CGAL {

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the wrl format (VRML 2.0).

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \see \ref IOStreamWRL
*/
template <typename FaceGraph, typename NamedParameters>
bool write_WRL(std::ostream& os,
               const FaceGraph& g,
               const NamedParameters& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor    face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typename CGAL::GetVertexPointMap<FaceGraph, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  boost::container::flat_map<vertex_descriptor,vertices_size_type> reindex;
  int n = 0;

  os << "#VRML V2.0 utf8\n"
        "Group {\n"
        "children [\n"
        "Shape {\n"
        "appearance DEF A1 Appearance {\n"
        "material Material {\n"
        "diffuseColor .6 .5 .9\n"
        "}\n"
        "}\n"
        "appearance\n"
        "Appearance {\n"
        "material DEF Material Material {}\n"
        "}\n"
        "}\n"
        "Group {\n"
        "children [\n"
        "Shape {\n"
        "appearance Appearance { material USE Material }\n"
        "geometry IndexedFaceSet {\n"
        "convex FALSE\n"
        "solid  FALSE\n"
        "coord  Coordinate {\n"
        "point [\n";

  for(vertex_descriptor v : vertices(g))
  {
    os << get(vpm,v) << ",\n";
    reindex[v] = n++;
  }

  os << "] #point\n"
        "} #coord Coordinate\n"
        "coordIndex  [\n";

  for(face_descriptor f : faces(g))
  {
    for(vertex_descriptor v : vertices_around_face(halfedge(f, g), g))
      os << reindex[v] << ",";
    os << "-1,\n";
  }

  os << "] #coordIndex\n"
        "} #geometry\n"
        "} #Shape\n"
        "] #children\n"
        "} #group\n"
        "]\n"
        "}\n";

  return os.good();
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_WRL(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_WRL(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool write_WRL(std::ostream& os, const FaceGraph& g) { return write_WRL(os, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_WRL(const char* fname, const FaceGraph& g) { return write_WRL(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_WRL(const std::string& fname, const FaceGraph& g) { return write_WRL(fname, g, parameters::all_default()); }


} // namespace CGAL

#endif // CGAL_BGL_IO_WRL_H
