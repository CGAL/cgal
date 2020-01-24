// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BGL_IO_OFF_H
#define CGAL_BGL_IO_OFF_H

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/container/flat_map.hpp>

#include <fstream>
#include <iostream>
#include <string>

// @todo reintroduce deprecated versions of the functions using lower case file formats

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {
namespace read_off_tools {

inline bool is_whitespace(const std::string& s)
{
  for(unsigned int i=0; i<s.size(); ++i)
  {
    if(s[i] != ' ' && s[i] != '\t')
      return false;
  }

  return true;
}

inline std::string next_non_comment(std::istream& is)
{
  std::string line;
  do
  {
    std::getline(is, line);
  }
  while(line[0] == '#' || is_whitespace(line));

  return line;
}

} // namespace read_off_tools
} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash, and lines with whitespace.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
    \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream are added.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
              FaceGraph& g,
              const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type          vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor             face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type             faces_size_type;

  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type                     Point;

  using namespace IO::internal::read_off_tools;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, g));

  vertices_size_type nv, nvf;
  faces_size_type nf;
  int ignore;

  std::string line = next_non_comment(is);
  {
    std::istringstream iss(line);
    std::string off;
    iss >> off;
    CGAL_assertion( off == "OFF" || off == "COFF");
  }

  line = next_non_comment(is);
  {
    std::istringstream iss(line);
    iss >> nv >> nf >> ignore;
  }

  std::vector<vertex_descriptor> vertices(nv);
  Point p;

  for(vertices_size_type i=0; i<nv; ++i)
  {
    line = next_non_comment(is);
    std::istringstream iss(line);
    iss >> p;
    vertices[i] = add_vertex(g);
    put(vpm, vertices[i], p);
  }

  for(faces_size_type i=0; i<nf; ++i)
  {
    line = next_non_comment(is);
    std::istringstream iss(line);
    iss >> nvf;
    std::vector<vertex_descriptor> face(nvf);
    for(vertices_size_type j=0; j<nvf; ++j)
    {
      faces_size_type fvi;
      iss >> fvi;
      face[j] = vertices[fvi];
    }

    face_descriptor f = CGAL::Euler::add_face(face, g);
    if(f == boost::graph_traits<FaceGraph>::null_face())
      return false;
  }

  return true;
}

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream are added.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const char* fname,
              FaceGraph& g,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  return read_OFF(in, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return read_OFF(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool read_OFF(std::istream& is, FaceGraph& g) { return read_OFF(is, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_OFF(const char* fname, FaceGraph& g) { return read_OFF(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_OFF(const std::string& fname, FaceGraph& g) { return read_OFF(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the OFF format.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor    face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<FaceGraph>::faces_size_type    faces_size_type;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, g));

  vertices_size_type nv = static_cast<vertices_size_type>(std::distance(vertices(g).first, vertices(g).second));
  faces_size_type nf = static_cast<faces_size_type>(std::distance(faces(g).first, faces(g).second));

  if(!os.good())
    return false;

  os << "OFF\n"
     << nv << " " << nf << " 0\n";

  boost::container::flat_map<vertex_descriptor, vertices_size_type> reindex;
  int n = 0;
  for(vertex_descriptor v : vertices(g))
  {
    os << get(vpm, v) << '\n';
    reindex[v] = n++;
  }

  for(face_descriptor f : faces(g))
  {
    os << degree(f, g);
    for(vertex_descriptor v : vertices_around_face(halfedge(f, g), g))
      os << " " << reindex[v];
    os << '\n';
  }

  return os.good();
}

/*!
\ingroup PkgBGLIOFct

 writes the graph `g` in the OFF format into a file named `fname`.

 \sa Overloads of this function for specific models of the concept `FaceGraph`.

 \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const char* fname,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream out(fname);
  return write_OFF(out, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_OFF(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool write_OFF(std::ostream& os, const FaceGraph& g) { return write_OFF(os, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_OFF(const char* fname, const FaceGraph& g) { return write_OFF(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_OFF(const std::string& fname, const FaceGraph& g) { return write_OFF(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_OFF_H
