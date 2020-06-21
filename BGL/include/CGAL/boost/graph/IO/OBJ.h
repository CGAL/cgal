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

#ifndef CGAL_BGL_IO_OBJ_H
#define CGAL_BGL_IO_OBJ_H

#include <CGAL/IO/OBJ.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename Graph, typename Point>
class OBJ_builder
  : public Generic_facegraph_builder<Graph, Point, OBJ_builder<Graph, Point> >
{
  typedef OBJ_builder<Graph, Point>                                         Self;
  typedef Generic_facegraph_builder<Graph, Point, Self>                     Base;

  typedef typename Base::Point_container                                    Point_container;
  typedef typename Base::Face                                               Face;
  typedef typename Base::Face_container                                     Face_container;

public:
  OBJ_builder(std::istream& is, bool verbose) : Base(is, verbose) { }

  template <typename NamedParameters>
  bool read(std::istream& is,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np,
            bool verbose)
  {
    return read_OBJ(is, points, faces, np, verbose);
  }
};

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from the stream `in` in the OBJ format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param is the input stream
  \param g the graph to be built from the input data
  \param verbose whether extra information is printed when an incident occurs during reading
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream is added.

  \returns `true` if the resulting mesh is valid.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOBJ
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(std::istream& is,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                  Point;

  IO::internal::OBJ_builder<Graph, Point> builder(is, verbose);
  return builder(g, np);
}

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OBJ format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the input file
  \param g the graph to be built from the input data
  \param verbose whether extra information is printed when an incident occurs during reading
  \param np optional \ref bgl_namedparameters "Named Parameters" described below

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `Graph`.}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream is added.

  \returns `true` if the resulting mesh is valid.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOBJ
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const char* fname,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  std::ifstream in(fname);
  return read_OBJ(in, g, np, verbose);
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const std::string& fname,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  return read_OBJ(fname.c_str(), g, np, verbose);
}

template <typename Graph>
bool read_OBJ(std::istream& is, Graph& g) { return read_OBJ(is, g, parameters::all_default()); }
template <typename Graph>
bool read_OBJ(const char* fname, Graph& g) { return read_OBJ(fname, g, parameters::all_default()); }
template <typename Graph>
bool read_OBJ(const std::string& fname, Graph& g) { return read_OBJ(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 \ingroup PkgBGLIOFct

  writes the graph `g` in the OBJ format.

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

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOBJ
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(std::ostream& os,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  IO::internal::Generic_facegraph_printer<std::ostream, Graph, CGAL::File_writer_wavefront> printer(os);
  return printer(g, np);
}

/*!
\ingroup PkgBGLIOFct

  writes the graph `g` in the OBJ format into a file named `fname`.

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the output file
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

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOBJ
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const char* fname,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream os(fname);
  return write_OBJ(os, g, np);
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& fname, const Graph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_OBJ(fname.c_str(), g, np);
}

template <typename Graph>
bool write_OBJ(std::ostream& os, const Graph& g) { return write_OBJ(os, g, parameters::all_default()); }
template <typename Graph>
bool write_OBJ(const char* fname, const Graph& g) { return write_OBJ(fname, g, parameters::all_default()); }
template <typename Graph>
bool write_OBJ(const std::string& fname, const Graph& g) { return write_OBJ(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_OBJ_H
