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

#include <CGAL/IO/VRML.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>

namespace CGAL {

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the wrl format (VRML 2.0).

  \tparam FaceGraph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
  \param g the graph to be output
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`
    \cgalParamEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful.

  \see \ref IOStreamWRL
*/
template <typename FaceGraph, typename NamedParameters>
bool write_WRL(std::ostream& os,
               const FaceGraph& g,
               const NamedParameters& np)
{
  IO::internal::Generic_facegraph_printer<std::ostream, FaceGraph, CGAL::File_writer_VRML_2> printer(os);
  return printer(g, np);
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
