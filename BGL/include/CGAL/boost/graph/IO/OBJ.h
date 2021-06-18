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
  OBJ_builder(std::istream& is) : Base(is) { }

  template <typename NamedParameters>
  bool read(std::istream& is,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np)
  {
    return read_OBJ(is, points, faces, np);
  }
};

} // namespace internal

/*!
  \ingroup PkgBGLIoFuncsOBJ

  \brief reads the graph `g` from the stream `in`, using the \ref IOStreamOBJ.

  The data is expected to represent a 2-manifold (possibly with borders).

  Ignores comment lines which start with a hash, and lines with whitespace.

  \attention The graph `g` is not cleared, and the data from the stream are appended.

  \tparam Graph a model of `MutableFaceGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param is the input stream
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

    \cgalParamNBegin{verbose}
      \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
      \cgalParamType{Boolean}
      \cgalParamDefault{`false`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if reading was successful and the resulting mesh is valid, `false` otherwise.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(std::istream& is,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
              )
{
  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                  Point;

  internal::OBJ_builder<Graph, Point> builder(is);
  return builder(g, np);
}

/// \cond SKIP_IN_MANUAL

template <typename Graph>
bool read_OBJ(std::istream& is, Graph& g,
              typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return read_OBJ(is, g, parameters::all_default());
}

/// \endcond

/*!
  \ingroup PkgBGLIoFuncsOBJ

  \brief reads the graph `g` from the file `fname`, using the \ref IOStreamOBJ.

  The data is expected to represent a 2-manifold (possibly with borders).

  Ignores comment lines which start with a hash, and lines with whitespace.

  \attention The graph `g` is not cleared, and the data from the file are appended.

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

    \cgalParamNBegin{verbose}
      \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
      \cgalParamType{Boolean}
      \cgalParamDefault{`false`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if reading was successful and the resulting mesh is valid, `false` otherwise.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const std::string& fname,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
              )
{
  std::ifstream is(fname);
  CGAL::IO::set_mode(is, CGAL::IO::ASCII);
  return read_OBJ(is, g, np);
}

/// \cond SKIP_IN_MANUAL

template <typename Graph>
bool read_OBJ(const std::string& fname, Graph& g,
              typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return read_OBJ(fname, g, parameters::all_default());
}

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 \ingroup PkgBGLIoFuncsOBJ

  \brief writes the graph `g` into the output stream, using the \ref IOStreamOBJ.

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

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(std::ostream& os,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
               )
{
  internal::Generic_facegraph_printer<std::ostream, Graph, CGAL::File_writer_wavefront> printer(os);
  return printer(g, np);
}

/// \cond SKIP_IN_MANUAL

template <typename Graph>
bool write_OBJ(std::ostream& os, const Graph& g,
               typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_OBJ(os, g, parameters::all_default());
}

/// \endcond

/*!
\ingroup PkgBGLIoFuncsOBJ

  \brief writes the graph `g` into a file named `fname`, using the \ref IOStreamOBJ.

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
    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful, `false` otherwise.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& fname,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
               )
{
  std::ofstream os(fname);
  CGAL::IO::set_mode(os, CGAL::IO::ASCII);
  return write_OBJ(os, g, np);
}

/// \cond SKIP_IN_MANUAL

template <typename Graph>
bool write_OBJ(const std::string& fname, const Graph& g,
               typename boost::disable_if<internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_OBJ(fname, g, parameters::all_default());
}

/// \endcond

}} // namespace CGAL::IO

#endif // CGAL_BGL_IO_OBJ_H
