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

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

//not for now.
/*
template <class FaceGraph, typename NamedParameters>
bool read_polygon_mesh(std::istream& is,
                       FaceGraph& g,
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

template <class FaceGraph>
bool read_polygon_mesh(std::istream& is,
                       FaceGraph& g)
{
  return read_polygon_mesh(is, g, parameters::all_default());
}
*/

/*!
 * \ingroup PkgBGLIOFct
 * \brief reads a polygon mesh from a file.
 * \tparam FaceGraph a model of `FaceGraph`
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format"),
 * `.vtp`(\ref IOStreamVTK "VTP file format")  or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param g the mesh
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 * \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
 * \cgalNamedParamsEnd
 * Other named parameters may be used according to the file extension.
 * See `PkgBGLIOFct` for an exhaustive list.

 * \return `true` if the reading worked, `false` otherwise.
 *
 * \pre The data must represent a 2-manifold
 *
 * \see \ref IOStreamOFF
 *
 */
template <class FaceGraph, typename NamedParameters>
bool read_polygon_mesh(const std::string& fname,
                       FaceGraph& g,
                       const NamedParameters& np)
{

  if (fname.find(".obj") != std::string::npos) {
    return read_OBJ(fname, g, np);
  }

  if (fname.find(".off") != std::string::npos) {
    return read_OFF(fname, g, np);
  }

  if (fname.find(".ply") != std::string::npos) {
    return read_PLY(fname, g, np);
  }

  if (fname.find(".stl") != std::string::npos) {
    return read_STL(fname, g, np);
  }

#ifdef CGAL_USE_VTK
  if (fname.find(".vtp") != std::string::npos) {
    return read_VTP(fname, g, np);
  }
#endif

  if (fname.find(".ts") != std::string::npos) {
    return read_GOCAD(fname, g, np);
  }
  return false;
}

template <class FaceGraph>
bool read_polygon_mesh(const std::string& fname,
                       FaceGraph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}



template <class FaceGraph, typename NamedParameters>
bool read_polygon_mesh(const char* fname,
                       FaceGraph& g,
                       const NamedParameters& np)
{
  return read_polygon_mesh(std::string(fname), g, np);
}

template <class FaceGraph>
bool read_polygon_mesh(const char* fname,
                       FaceGraph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

/*!
 * \ingroup PkgBGLIOFct
 * \brief writes a polygon mesh in a file.
 * \tparam FaceGraph a model of `FaceGraph`
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param fname the name of the file. Its extension must be one of the following :
 * `.off` (\ref IOStreamOFF "OFF file format") , `.obj` (\ref IOStreamOBJ "OBJ file format"),
 * `.stl` (\ref IOStreamSTL "STL file format"), `.ply` (\ref IOStreamPLY "PLY file format"),
 * `.vtp`(\ref IOStreamVTK "VTP file format")  or `.ts`(\ref IOStreamGocad "GOCAD file format").
 * \param g the mesh
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 * \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
 * \cgalNamedParamsEnd
 * Other named parameters may be used according to the file extension.
 * See `PkgBGLIOFct`  for an exhaustive list.
 * \return `true` if the writing worked, `false` otherwise.
 *
 * \see \ref IOStreamOFF
 */
template <class FaceGraph, typename NamedParameters>
bool write_polygon_mesh(const std::string& fname,
                       FaceGraph& g,
                       const NamedParameters& np)
{
  if (fname.find(".ts") != std::string::npos) {
    return write_GOCAD(fname, g, np);
  }

  if (fname.find(".obj") != std::string::npos) {
    return write_OBJ(fname, g, np);
  }

  if (fname.find(".off") != std::string::npos) {
    return write_OFF(fname, g, np);
  }

  if (fname.find(".ply") != std::string::npos) {
    return write_PLY(fname, g, np);
  }

  if (fname.find(".stl") != std::string::npos) {
    return write_STL(fname, g, np);
  }
#ifdef CGAL_USE_VTK
  if (fname.find(".vtp") != std::string::npos) {
    return write_VTP(fname, g, np);
  }
#endif

  return false;
}

template <class FaceGraph>
bool write_polygon_mesh(const std::string& fname,
                       FaceGraph& g)
{
  return write_polygon_mesh(fname, g, parameters::all_default());
}

template <class FaceGraph, typename NamedParameters>
bool write_polygon_mesh(const char* fname,
                       FaceGraph& g,
                       const NamedParameters& np)
{
  return write_polygon_mesh(std::string(fname), g, np);
}

template <class FaceGraph>
bool write_polygon_mesh(const char* fname,
                       FaceGraph& g)
{
  return write_polygon_mesh(fname, g, parameters::all_default());
}

}//end CGAL

#endif // CGAL_BOOST_GRAPH_POLYGON_MESH_IO_H
