// Copyright (c) 2015-2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_INP_H
#define CGAL_BGL_IO_INP_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/container/flat_map.hpp>

#include <fstream>
#include <string>

namespace CGAL {

namespace IO {

/// \cond SKIP_IN_MANUAL

template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_INP(std::ostream& os,
               const std::string& name,
               const std::string& type,
               const Graph& g,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor                  vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::face_descriptor                    face_descriptor;
  typedef typename boost::graph_traits<Graph>::vertices_size_type                 vertices_size_type;

  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_NP_CLASS>::const_type  VPM;
  typedef typename boost::property_traits<VPM>::reference                         Point_ref;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  if(!os.good())
    return false;

  os << "*Part, name=" << name << "\n*Node\n";
  boost::container::flat_map<vertex_descriptor,vertices_size_type> reindex;
  int n = 1;
  for(const vertex_descriptor v : vertices(g))
  {
    Point_ref p = get(vpm,v);
    os << n << ", " << p.x() << ", " << p.y() << ", " << p.z() << '\n';
    reindex[v] = n++;
  }

  n = 1;
  os << "*Element, type=" << type << std::endl;
  for(const face_descriptor f : faces(g))
  {
    os << n++;
    for(const vertex_descriptor v : CGAL::vertices_around_face(halfedge(f, g), g))
      os << ", " << reindex[v];

    os << '\n';
  }

  os << "*End Part"<< std::endl;

  return os.good();
}

template <typename Graph, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_INP(const std::string& fname,
               const std::string& type,
               const Graph& g,
               const CGAL_NP_CLASS& np = parameters::default_values() )
{
  std::ofstream os(fname);
  return write_INP(os, fname, type, g, np);
}

#ifndef CGAL_NO_DEPRECATED_CODE
template <typename FaceGraph, typename NamedParameters = parameters::Default_named_parameters>
CGAL_DEPRECATED bool write_inp(std::ostream& os,
               const FaceGraph& g,
               std::string name,
               std::string type,
               const NamedParameters& np = parameters::default_values())
{
  return write_INP(os, name, type, g, np);
}
#endif
/// \endcond

}} // namespace CGAL::IO

#endif // CGAL_BGL_IO_INP_H
