// Copyright (c) 2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno
#ifndef CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H
#define CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H

#include <CGAL/boost/graph/io.h>
#include <CGAL/IO/polygon_soup_io.h>

namespace CGAL {

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
 * \tparam Graph a model of `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format"),
 * `.vtp` (\ref IOStreamVTK "VTP file format")  or `.ts` (\ref IOStreamGocad "GOCAD file format").
 * \param g the mesh
 * \param verbose whether extra information is printed when an incident occurs during reading
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * Other named parameters may be used according to the file extension, see \ref PkgBGLIOFct for an exhaustive list.
 *
 * \pre The data must represent a 2-manifold
 *
 * \return `true` if reading was successful, `false` otherwise.
*/
template <class Graph, typename NamedParameters>
bool read_polygon_mesh(const std::string& fname,
                       Graph& g,
                       const NamedParameters& np,
                       const bool verbose = false)
{
  const std::string ext = IO::internal::get_file_extension(fname);
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

template <class Graph>
bool read_polygon_mesh(const std::string& fname, Graph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}

template <class Graph, typename NamedParameters>
bool read_polygon_mesh(const char* fname, Graph& g, const NamedParameters& np)
{
  return read_polygon_mesh(std::string(fname), g, np);
}

template <class Graph>
bool read_polygon_mesh(const char* fname, Graph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgBGLIOFct
 *
 * \brief writes a polygon mesh in a file.
 *
 * \tparam Graph a model of `FaceListGraph` and `HalfedgeListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format"),
 * `.vtp` (\ref IOStreamVTK "VTP file format")  or `.ts` (\ref IOStreamGocad "GOCAD file format").
 * \param g the mesh to be output
 * \param verbose whether extra information is printed when an incident occurs during writing
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
 * \cgalNamedParamsEnd
 *
 * Other named parameters may be used according to the file extension, see \ref PkgBGLIOFct for an exhaustive list.
 *
 * \return `true` if writing was successful, `false` otherwise.
 */
template <class Graph, typename NamedParameters>
bool write_polygon_mesh(const std::string& fname,
                        Graph& g,
                        const NamedParameters& np,
                        const bool verbose = false)
{
  const std::string ext = IO::internal::get_file_extension(fname);
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

template <class Graph>
bool write_polygon_mesh(const std::string& fname, Graph& g)
{
  return write_polygon_mesh(fname, g, parameters::all_default());
}

template <class Graph, typename NamedParameters>
bool write_polygon_mesh(const char* fname, Graph& g, const NamedParameters& np)
{
  return write_polygon_mesh(std::string(fname), g, np);
}

template <class Graph>
bool write_polygon_mesh(const char* fname, Graph& g)
{
  return write_polygon_mesh(fname, g, parameters::all_default());
}

} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H
