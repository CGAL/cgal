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

#ifndef CGAL_BGL_IO_INP_H
#define CGAL_BGL_IO_INP_H

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/container/flat_map.hpp>

namespace CGAL {

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_INP(std::ostream& os,
               const std::string name,
               const std::string type,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor                  vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor                    face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type                 vertices_size_type;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::const_type  VPM;
  typedef typename boost::property_traits<VPM>::reference                             Point_ref;

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

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_INP(const char* fname,
               const std::string type,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream out(fname);
  return write_INP(out, fname, type, g, np);
}

template <typename FaceGraph>
bool write_INP(std::ostream& os, const std::string name, const std::string type, const FaceGraph& g)
{
  return write_INP(os, name, type, g, parameters::all_default());
}

template <typename FaceGraph>
bool write_INP(const char* fname, const std::string type, const FaceGraph& g)
{
  return write_INP(fname, type, g, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_BGL_IO_INP_H
