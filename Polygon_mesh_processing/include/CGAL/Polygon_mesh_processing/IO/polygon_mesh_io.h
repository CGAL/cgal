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

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace IO {

// @todo also need named parameters overload
template <typename FaceGraph, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_polygon_mesh(const std::string& filename,
                       FaceGraph& g,
                       const CGAL_BGL_NP_CLASS& np)
{
  typedef typename CGAL::GetVertexPointMap<FaceGraph, CGAL_BGL_NP_CLASS>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  bool ok = ::CGAL::read_polygon_mesh(filename, g, np);

  if(!ok)
  {
    std::vector<Point> points;
    std::vector<std::vector<std::size_t> > faces;
    ok = CGAL::read_polygon_soup(in, points, faces, np);
    if(!ok)
      return false;

    ok = CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
    if(!ok)
      return false;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, g);
  }

  return true;
}

} // namespace IO
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_PMP_IO_POLYGON_MESH_IO_H
