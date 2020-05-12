// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_POLYGON_MESH_PROCESSING_READ_POLYGON_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_READ_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <fstream>
#include <vector>

namespace CGAL{
namespace Polygon_mesh_processing{

/*!
 * \brief read a file as a polygon soup, and then repair and orient it before trying to convert it
 * into a `FaceGraph`.
 * \param fname the name of the input file.
 * \param g the FaceGraph.
 * \param np sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamBegin{vertex_point_map}
*     a model of `WritablePropertyMap`, the property map with the points associated to the vertices of `out`.
*     If this parameter is omitted, an internal property map for
*     `CGAL::vertex_point_t` must be available in `PolygonMesh`.
*   \cgalParamEnd
* \cgalNamedParamsEnd
 * \return `true` if the reading and conversion worked, `false` otherwise.
 */
template <typename FaceGraph, typename NamedParameter>
bool read_polygon_mesh(const char* fname,
                       FaceGraph& g,
                       const NamedParameter& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameter>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                   Point;
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;
  if(!CGAL::read_polygon_soup(fname, points, faces))
  {
    std::cerr << "Error: cannot read file\n";
    return false;
  }

  std::cout << "Cleaning polygon soup..." << std::endl;
  repair_polygon_soup(points, faces, np);

  if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
  {
    std::cerr << "W: File does not describe a polygon mesh" << std::endl;
  }

  CGAL::Polygon_mesh_processing::
      polygon_soup_to_polygon_mesh(points, faces, g, parameters::all_default(), np);
  return true;
}

template <typename FaceGraph>
bool read_polygon_mesh(const char* fname, FaceGraph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}

template <typename FaceGraph, typename NamedParameter>
bool read_polygon_mesh(const std::string& fname,
                       FaceGraph& g,
                       const NamedParameter& np)
{
  return read_polygon_mesh(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool read_polygon_mesh(const std::string& fname, FaceGraph& g)
{
  return read_polygon_mesh(fname, g, parameters::all_default());
}

}//end PMP
}//end CGAL
#endif // CGAL_POLYGON_MESH_PROCESSING_READ_POLYGON_MESH_H
