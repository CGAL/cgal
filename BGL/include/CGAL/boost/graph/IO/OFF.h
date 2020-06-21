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

#include <CGAL/IO/OFF.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/utility/enable_if.hpp>

#include <fstream>
#include <iostream>
#include <string>

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
class OFF_builder
  : public Generic_facegraph_builder<Graph, Point, OFF_builder<Graph, Point> >
{
  typedef OFF_builder<Graph, Point>                                         Self;
  typedef Generic_facegraph_builder<Graph, Point, Self>                     Base;

  typedef typename Base::Point_container                                    Point_container;
  typedef typename Base::Face                                               Face;
  typedef typename Base::Face_container                                     Face_container;

public:
  OFF_builder(std::istream& is, bool verbose) : Base(is, verbose) { }

  template <typename NamedParameters>
  bool read(std::istream& is,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np,
            bool verbose)
  {
    return read_OFF(is, points, faces, np, verbose);
  }
};

// Because some packages can provide overloads with the same signature to automatically initialize
// property maps (see Surface_mesh/IO/ for example)
template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF_BGL(std::istream& is,
                  Graph& g,
                  const CGAL_BGL_NP_CLASS& np,
                  bool verbose = true)
{
  typedef typename CGAL::GetVertexPointMap<Graph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                  Point;

  IO::internal::OFF_builder<Graph, Point> builder(is, verbose);
  return builder(g, np);
}

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash,
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

    \cgalParamNBegin{vertex_normal_map}
      \cgalParamDescription{a property map associating normals to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Vector_3` as value type}
      \cgalParamDefault{vertex normals that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{vertex colors that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_texture_map}
      \cgalParamDescription{a property map associating textures to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_2` as value type}
      \cgalParamDefault{vertex textures that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{face colors that may exist in the input will be ignored}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream is added.

  \returns `true` if the resulting mesh is valid.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOFF
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
              )
{
  return IO::internal::read_OFF_BGL(is, g, np, verbose);
}

template <typename Graph>
bool read_OFF(std::istream& is, Graph& g,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return read_OFF(is, g, parameters::all_default());
}

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from `fname`, a file in the OFF format. Ignores comment lines which start with a hash,
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

    \cgalParamNBegin{vertex_normal_map}
      \cgalParamDescription{a property map associating normals to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Vector_3` as value type}
      \cgalParamDefault{vertex normals that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{vertex colors that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_texture_map}
      \cgalParamDescription{a property map associating textures to the vertices of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_2` as value type}
      \cgalParamDefault{vertex textures that may exist in the input will be ignored}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{face colors that may exist in the input will be ignored}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \attention The graph `g` is not cleared, and the data from the stream is added.

  \returns `true` if the resulting mesh is valid.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOFF
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const char* fname,
              Graph& g,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true
#ifndef DOXYGEN_RUNNING
              , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
              )
{
  std::ifstream is(fname);
  return read_OFF(is, g, np, verbose);
}

template <typename Graph>
bool read_OFF(const char* fname, Graph& g,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return read_OFF(fname, g, parameters::all_default());
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname, Graph& g, const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  return read_OFF(fname.c_str(), g, np, verbose);
}

template <typename Graph>
bool read_OFF(const std::string& fname, Graph& g,
              typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return read_OFF(fname, g, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {
namespace internal {

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_BGL(std::ostream& os,
                   const Graph& g,
                   const CGAL_BGL_NP_CLASS& np)
{
  IO::internal::Generic_facegraph_printer<std::ostream, Graph, CGAL::File_writer_OFF> printer(os);
  return printer(g, np);
}

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the OFF format.

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

    \cgalParamNBegin{vertex_normal_map}
      \cgalParamDescription{a property map associating normals to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Vector_3` as value type}
      \cgalParamDefault{no vertex normals in the output}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{no vertex colors in the output}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_texture_map}
      \cgalParamDescription{a property map associating textures to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_2` as value type}
      \cgalParamDefault{no vertex textures in the output}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{no face colors in the output}
    \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOFF
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
               )
{
  return IO::internal::write_OFF_BGL(os, g, np);
}

template <typename Graph>
bool write_OFF(std::ostream& os, const Graph& g,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_OFF(os, g, parameters::all_default());
}

/*!
  \ingroup PkgBGLIOFct

  writes the graph `g` in the file `fname`, in the OFF format.

  \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the name of the output file
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

    \cgalParamNBegin{vertex_normal_map}
      \cgalParamDescription{a property map associating normals to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Vector_3` as value type}
      \cgalParamDefault{no vertex normals in the output}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_color_map}
      \cgalParamDescription{a property map associating colors to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{no vertex colors in the output}
    \cgalParamNEnd

    \cgalParamNBegin{vertex_texture_map}
      \cgalParamDescription{a property map associating textures to the vertices of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
                     as key type and `%Point_2` as value type}
      \cgalParamDefault{no vertex textures in the output}
    \cgalParamNEnd

    \cgalParamNBegin{face_color_map}
      \cgalParamDescription{a property map associating colors to the faces of `g`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%face_descriptor`
                     as key type and `CGAL::Color` as value type}
      \cgalParamDefault{no face colors in the output}
    \cgalParamNEnd

    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \returns `true` if writing was successful.

  \sa Overloads of this function for specific models of the concept `FaceGraph`.
  \see \ref IOStreamOFF
*/
template <typename Graph,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const char* fname,
               const Graph& g,
               const CGAL_BGL_NP_CLASS& np
#ifndef DOXYGEN_RUNNING
               , typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr
#endif
               )
{
  std::ofstream os(fname);
  if(!os)
  {
    std::cerr<<"Could not create file.";
    return false;
  }
  return write_OFF(os, g, np);
}

template <typename Graph>
bool write_OFF(const char* fname, const Graph& g,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_OFF(fname, g, parameters::all_default());
}

template <typename Graph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& fname, const Graph& g, const CGAL_BGL_NP_CLASS& np,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_OFF(fname.c_str(), g, np);
}

template <typename Graph>
bool write_OFF(const std::string& fname, const Graph& g,
               typename boost::disable_if<IO::internal::is_Point_set_or_Range_or_Iterator<Graph> >::type* = nullptr)
{
  return write_OFF(fname.c_str(), g, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_BGL_IO_OFF_H
