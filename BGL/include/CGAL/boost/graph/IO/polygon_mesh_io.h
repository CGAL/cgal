// Copyright (c) 2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno
//                 Mael Rouxel-Labbé

#ifndef CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H
#define CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H

#include <CGAL/boost/graph/IO/3MF.h>
#include <CGAL/boost/graph/IO/GOCAD.h>
#include <CGAL/boost/graph/IO/INP.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/IO/PLY.h>
#include <CGAL/boost/graph/IO/STL.h>
#include <CGAL/boost/graph/IO/VTK.h>
#include <CGAL/boost/graph/IO/WRL.h>
#include <CGAL/IO/helpers.h>

#include <fstream>
#include <string>

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

//not for now : some readers will return "ok" despite not managing to read anything
/*
template <class Graph, typename NamedParameters>
bool read_polygon_mesh(std::istream& is,
                       Graph& g,
                       const NamedParameters& np)
{
  bool ok = false;
  ok = read_OFF(is, g, np, false);
  if(ok)
    return true;
  g.clear();
  is.clear();//reset the error state
  is.seekg (0, is.beg);
  ok = read_OBJ(is, g, np, false);
  if(ok)
    return true;
  g.clear();
  is.clear();
  is.seekg (0, is.beg);
  ok = read_PLY(is, g, np, false);
  if(ok)
    return true;
  g.clear();
  is.clear();
  is.seekg (0, is.beg);
  ok = read_STL(is, g, np, false);
  if(ok)
    return true;
  g.clear();
  is.clear();
  is.seekg (0, is.beg);
  ok = read_GOCAD(is, g, np, false);
  return ok;
}

template <class Graph>
bool read_polygon_mesh(std::istream& is,
                       Graph& g)
{
  return read_polygon_mesh(is, g, parameters::all_default());
}
*/

/*!
 * \ingroup PkgBGLIOFct
 *
 * \brief reads a polygon mesh from a file.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * The format is detected from the filename extension (letter case is not important).
 *
 * The data is expected to represent a 2-manifold (possibly with borders).
 *
 * \tparam Graph a model of `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the file
 * \param g the mesh
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * Other named parameters may be used according to the file extension, see \ref PkgBGLIOFct for an exhaustive list.
 *
 * \return `true` if reading was successful, `false` otherwise.
 *
 * \sa \link PMP_IO_grp `CGAL::Polygon_mesh_processing::IO::read_polygon_mesh()`\endlink if the data is not 2-manifold
*/
template <class Graph, typename NamedParameters>
bool read_polygon_mesh(const std::string& fname,
                       Graph& g,
                       const NamedParameters& np)
{
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  const std::string ext = internal::get_file_extension(fname);
  if(ext == std::string())
  {
    if(verbose)
      std::cerr << "Error: cannot read from file without extension" << std::endl;
    return false;
  }

  if(ext == "obj")
    return read_OBJ(fname, g, np);
  else if(ext == "off")
    return read_OFF(fname, g, np);
  else if(ext == "ply")
    return read_PLY(fname, g, np);
  else if(ext == "stl")
    return read_STL(fname, g, np);
  else if(ext == "ts")
    return read_GOCAD(fname, g, np);
#ifdef CGAL_USE_VTK
  else if(ext == "vtp")
    return read_VTP(fname, g, np);
#endif

  if(verbose)
  {
    std::cerr << "Error: unknown input file extension: " << ext << "\n"
              << "Please refer to the documentation for the list of supported file formats" << std::endl;
  }

  return false;
}

/// \cond SKIP_IN_MANUAL

template <class Graph>
bool read_polygon_mesh(const std::string& fname, Graph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}

/// \endcond

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgBGLIOFct
 *
 * \brief writes a polygon mesh in a file.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * The format is detected from the filename extension (letter case is not important).
 *
 * \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the file
 * \param g the mesh to be output
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{stream_precision}
 *     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
 *     \cgalParamType{int}
 *     \cgalParamDefault{`6`}
 *     \cgalParamExtra{This parameter is only meaningful while using ASCII encoding.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * Other named parameters may be used according to the file extension, see \ref PkgBGLIOFct for an exhaustive list.
 *
 * \return `true` if writing was successful, `false` otherwise.
 */
template <class Graph, typename NamedParameters>
bool write_polygon_mesh(const std::string& fname,
                        Graph& g,
                        const NamedParameters& np)
{
  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  const std::string ext = internal::get_file_extension(fname);
  if(ext == std::string())
  {
    if(verbose)
      std::cerr << "Error: trying to output to file without extension" << std::endl;
    return false;
  }

  if(ext == "obj")
    return write_OBJ(fname, g, np);
  else if(ext == "off")
    return write_OFF(fname, g, np);
  else if(ext == "ply")
    return write_PLY(fname, g, np);
  else if(ext == "stl")
    return write_STL(fname, g, np);
  else if(ext == "ts")
    return write_GOCAD(fname, g, np);
#ifdef CGAL_USE_VTK
  else if(ext == "vtp")
    return write_VTP(fname, g, np);
#endif

  if(verbose)
  {
    std::cerr << "Error: unknown output file extension: " << ext << "\n"
              << "Please refer to the documentation for the list of supported file formats" << std::endl;
  }

  return false;
}

/// \cond SKIP_IN_MANUAL

template <class Graph>
bool write_polygon_mesh(const std::string& fname, Graph& g)
{
  return write_polygon_mesh(fname, g, parameters::all_default());
}

/// \endcond

}} // namespace CGAL::IO

#endif // CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H
