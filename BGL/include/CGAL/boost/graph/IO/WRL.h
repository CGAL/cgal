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

#ifndef CGAL_BGL_IO_WRL_H
#define CGAL_BGL_IO_WRL_H

#include <CGAL/IO/VRML.h>

#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <string>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
  \ingroup PkgBGLIoFuncsWRL

  \brief writes the graph `g` into the output stream, using the \ref IOStreamWRL (VRML 2.0).

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
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

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`the precision of the stream `os``}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful, `false` otherwise.
*/
template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_WRL(std::ostream& os,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  CGAL::VRML_2_ostream vos(os);
  internal::Generic_facegraph_printer<CGAL::VRML_2_ostream, Graph, CGAL::File_writer_VRML_2> printer(vos);
  return printer(g, np);
}

/*!
  \ingroup PkgBGLIoFuncsWRL

  \brief writes the graph `g` into the output file, using the \ref IOStreamWRL (VRML 2.0).

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the output file
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

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful, `false` otherwise.
*/
template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_WRL(const std::string& fname, const Graph& g, const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream os(fname);
  return write_WRL(os, g, np);
}

template <typename Graph>
bool write_WRL(std::ostream& os, const Graph& g) { return write_WRL(os, g, parameters::all_default()); }
template <typename Graph>
bool write_WRL(const std::string& fname, const Graph& g) { return write_WRL(fname, g, parameters::all_default()); }

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
 \ingroup PkgBGLIOFctDeprecated

 \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_WRL()` should be used instead.
*/
template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool write_wrl(std::ostream& os, const Graph& g, const CGAL_BGL_NP_CLASS& np)
{
  return IO::write_WRL(os, g, np);
}

template <typename Graph>
CGAL_DEPRECATED bool write_wrl(std::ostream& os, const Graph& g)
{
  return write_wrl(os, g, parameters::all_default());
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_BGL_IO_WRL_H
