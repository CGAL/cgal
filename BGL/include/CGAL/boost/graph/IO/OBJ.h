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

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

namespace IO {
namespace internal {

// Use CRTP to gain access to the protected members without getters/setters.
template <typename FaceGraph, typename Point>
class OBJ_builder
  : public Generic_facegraph_builder<FaceGraph, Point, OBJ_builder<FaceGraph, Point> >
{
  typedef OBJ_builder<FaceGraph, Point>                                         Self;
  typedef Generic_facegraph_builder<FaceGraph, Point, Self>                     Base;

  typedef typename Base::Point_container                                        Point_container;
  typedef typename Base::Face                                                   Face;
  typedef typename Base::Face_container                                         Face_container;

public:
  OBJ_builder(std::istream& is_) : Base(is_) { }

  template <typename NamedParameters>
  bool read(std::istream& input,
            Point_container& points,
            Face_container& faces,
            const NamedParameters& np)
  {
    return read_OBJ(input, points, faces, np);
  }
};

} // namespace internal
} // namespace IO

/*!
  \ingroup PkgBGLIOFct

  reads the graph `g` from the stream `in` in the OBJ format.

  \returns `true` if the resulting mesh is valid.

  \pre The data must represent a 2-manifold

  \see \ref IOStreamOBJ
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(std::istream& in,
              FaceGraph& g,
              const CGAL_BGL_NP_CLASS& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  IO::internal::OBJ_builder<FaceGraph, Point> builder(in);
  return builder(g, np);
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
bool read_OBJ(const char* fname,
              FaceGraph& g,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  return read_OBJ(in, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OBJ(const std::string& fname,
              FaceGraph& g,
              const CGAL_BGL_NP_CLASS& np)
{
  return read_OBJ(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool read_OBJ(std::istream& is, FaceGraph& g) { return read_OBJ(is, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_OBJ(const char* fname, FaceGraph& g) { return read_OBJ(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool read_OBJ(const std::string& fname, FaceGraph& g) { return read_OBJ(fname, g, parameters::all_default()); }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 \ingroup PkgBGLIOFct

  writes the graph `g` in the OBJ format.

  \returns `true` if writing was successful.

  \see \ref IOStreamOBJ
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(std::ostream& os,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  IO::internal::Generic_facegraph_printer<std::ostream, FaceGraph, CGAL::File_writer_wavefront> printer(os);
  return printer(g, np);
}

/*!
\ingroup PkgBGLIOFct

 writes the graph `g` in the OFF format into a file named `fname`.

 \sa Overloads of this function for specific models of the concept `FaceGraph`.

 \see \ref IOStreamOFF
*/
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const char* fname,
               const FaceGraph& g,
               const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream out(fname);
  return write_OBJ(out, g, np);
}

template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OBJ(const std::string& fname, const FaceGraph& g, const CGAL_BGL_NP_CLASS& np)
{
  return write_OBJ(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool write_OBJ(std::ostream& os, const FaceGraph& g) { return write_OBJ(os, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_OBJ(const char* fname, const FaceGraph& g) { return write_OBJ(fname, g, parameters::all_default()); }
template <typename FaceGraph>
bool write_OBJ(const std::string& fname, const FaceGraph& g) { return write_OBJ(fname, g, parameters::all_default()); }

} // namespace CGAL

#endif // CGAL_BGL_IO_OBJ_H
