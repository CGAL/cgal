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

#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <string>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
  \ingroup PkgBGLIOFct

  \brief Writes the graph `g` in the \ref IOStreamWRL (VRML 2.0).

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
  \param g the graph to be output
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
  \cgalNamedParamsEnd

  \returns `true` if writing was successful.
*/
template <typename Graph, typename NamedParameters>
bool write_WRL(std::ostream& os,
               const Graph& g,
               const NamedParameters& np)
{
  IO::internal::Generic_facegraph_printer<std::ostream, Graph, CGAL::File_writer_VRML_2> printer(os);
  return printer(g, np);
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_WRL(const std::string& fname, const Graph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_WRL(fname.c_str(), g, np);
}

template <typename Graph>
bool write_WRL(std::ostream& os, const Graph& g) { return write_WRL(os, g, parameters::all_default()); }
template <typename Graph>
bool write_WRL(const char* fname, const Graph& g) { return write_WRL(fname, g, parameters::all_default()); }
template <typename Graph>
bool write_WRL(const std::string& fname, const Graph& g) { return write_WRL(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_WRL_H
