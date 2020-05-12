// Copyright (c) 2020  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PMP_IO_POLYGON_MESH_IO_H
#define CGAL_PMP_IO_POLYGON_MESH_IO_H

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace IO {


/*!
  \ingroup pmp_io_grp

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
 *  `repair_polygon_soup` a boolean that decides if the soup should be repaired or not. Default is `true`; \n
 *  named parameters used for `CGAL::Polygon_mesh_processing::repair_polygon_soup()` can also be used with this function.
 *
 * \return `true` if the reading and conversion worked, `false` otherwise.
 */
template <typename FaceGraph, typename NamedParameter>
bool read_polygon_mesh(const char* fname,
                       FaceGraph& g,
                       const NamedParameter& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, NamedParameter>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                   Point;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool ok = ::CGAL::read_polygon_mesh(fname, g, np);

  if(!ok)
  {
    std::vector<Point> points;
    std::vector<std::vector<std::size_t> > faces;
    if(!CGAL::read_polygon_soup(fname, points, faces))
    {
      std::cerr << "Error: cannot read file\n";
      return false;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    const bool do_repair = choose_parameter(get_parameter(np, internal_np::repair_polygon_soup), true);
    if(do_repair)
      CGAL::Polygon_mesh_processing::repair_polygon_soup(points, faces, np);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
    {
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;
    }

    CGAL::Polygon_mesh_processing::
        polygon_soup_to_polygon_mesh(points, faces, g, parameters::all_default(), np);
  }
  return true;
}

template <typename FaceGraph>
bool read_polygon_mesh(const char* fname, FaceGraph& g)
{
  return CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(fname, g, parameters::all_default());
}

template <typename FaceGraph, typename NamedParameter>
bool read_polygon_mesh(const std::string& fname,
                       FaceGraph& g,
                       const NamedParameter& np)
{
  return CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(fname.c_str(), g, np);
}

template <typename FaceGraph>
bool read_polygon_mesh(const std::string& fname, FaceGraph& g)
{
  return CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(fname, g, parameters::all_default());
}
} // namespace IO
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_PMP_IO_POLYGON_MESH_IO_H
