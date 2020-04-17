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
#include <CGAL/boost/graph/IO/Generic_facegraph_builder.h>
#include <CGAL/boost/graph/IO/Generic_facegraph_printer.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

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

// Use CRTP to gain access to the protected members without getters/setters.
template <typename FaceGraph, typename Point>
class OFF_builder
  : public Generic_facegraph_builder<FaceGraph, Point, OFF_builder<FaceGraph, Point> >
{
  typedef OFF_builder<FaceGraph, Point>                                         Self;
  typedef Generic_facegraph_builder<FaceGraph, Point, Self>                     Base;

  typedef typename Base::Point_container                                        Point_container;
  typedef typename Base::Face                                                   Face;
  typedef typename Base::Face_container                                         Face_container;

public:
  OFF_builder(std::istream& is_) : Base(is_) { }

  template <typename NamedParameters>
  bool read(std::istream& input,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np)
  {
    return read_OFF(input, points, faces, np);
  }
};

// Because some packages can provide overloads with the same signature to automatically initialize
// property maps (see Surface_mesh/IO/ for example)
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF_BGL(std::istream& in,
                  FaceGraph& g,
                  const CGAL_BGL_NP_CLASS& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  IO::internal::OFF_builder<FaceGraph, Point> builder(in);
  return builder(g, np);
}

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from data in the OFF format. Ignores comment lines which start with a hash,
  and lines with whitespace.

  \cgalNamedParamsBegin
    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
      If this parameter is omitted, an internal property map for
      `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
    \cgalNamedParamsEnd

  \pre The data must represent a 2-manifold

  \sa Overloads of this function for specific models of the concept `FaceGraph`.

  \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& in, FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return IO::internal::read_OFF_BGL(in, g, np);
}

// document that too
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const char* fname, FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
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
bool read_OFF(std::istream& is, FaceGraph& g,
              typename boost::disable_if<
              typename boost::has_range_const_iterator<FaceGraph>::type
              >::type* =0)
{
  return read_OFF(is, g, parameters::all_default());
}
template <typename FaceGraph>
bool read_OFF(const char* fname, FaceGraph& g) { return read_OFF(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_OFF(const std::string& fname, FaceGraph& g) { return read_OFF(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

namespace IO {
namespace internal {

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_BGL(std::ostream& os,
                   const FaceGraph& g,
                   const CGAL_BGL_NP_CLASS& np)
{
  IO::internal::Generic_facegraph_printer<std::ostream, FaceGraph, CGAL::File_writer_OFF> printer(os);
  return printer(g, np);
}

} // namespace internal
} // namespace IO

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
bool write_OFF(std::ostream& os, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return IO::internal::write_OFF_BGL(os, g, np);
}

// document that too
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const char* fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
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
bool write_OFF(std::ostream& os, const FaceGraph& g
               ,typename boost::disable_if<
               typename boost::has_range_const_iterator<FaceGraph>::type
               >::type* =0)
{
  return write_OFF(os, g, parameters::all_default());
}
template <typename FaceGraph>
bool write_OFF(const char* fname, const FaceGraph& g) { return write_OFF(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_OFF(const std::string& fname, const FaceGraph& g) { return write_OFF(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_OFF_H
