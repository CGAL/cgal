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
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PMP_IO_POLYGON_MESH_IO_H
#define CGAL_PMP_IO_POLYGON_MESH_IO_H

#include <CGAL/license/Polygon_mesh_processing.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace IO {
/*!
  \ingroup PMP_IO_grp

 * \brief reads the file as a polygon soup, repairs (using `repair_polygon_soup()`),
 * and orients it (using `orient_polygon_soup()`) as to obtain a polygon mesh.
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
 * If repairing and orientation are known to not be required, one can use
 * \link PkgBGLIOFct `CGAL::IO::read_polygon_mesh()` \endlink directly.
 *
 * \tparam PolygonMesh a model of `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the input file
 * \param g the polygon mesh
 * \param np sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{repair_polygon_soup}
 *     \cgalParamDescription{a parameter used indicate whether `CGAL::Polygon_mesh_processing::repair_polygon_soup()`
 *                           should be called on the intermediate polygon soup.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the reading, repairing, and orientation operations were successful, `false` otherwise.
 *
 * \sa \link PkgBGLIOFct `CGAL::IO::read_polygon_mesh()` \endlink
 */
template <typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
bool read_polygon_mesh(const std::string& fname,
                       PolygonMesh& g,
                       const NamedParameters& np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type                     Point;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const bool verbose = parameters::choose_parameter(parameters::get_parameter(np, internal_np::verbose), false);

  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;
  if(!CGAL::IO::read_polygon_soup(fname, points, faces, CGAL::parameters::verbose(verbose)))
  {
    if(verbose)
      std::cerr << "Warning: cannot read polygon soup" << std::endl;
    return false;
  }

  const bool do_repair = choose_parameter(get_parameter(np, internal_np::repair_polygon_soup), true);
  if(do_repair)
    PMP::repair_polygon_soup(points, faces, np);

  if(!PMP::orient_polygon_soup(points, faces))
  {
    if(verbose)
      std::cerr << "Some duplication happened during polygon soup orientation" << std::endl;
  }

  if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
  {
    if(verbose)
      std::cerr << "Warning: polygon soup does not describe a polygon mesh" << std::endl;
    return false;
  }

  PMP::polygon_soup_to_polygon_mesh(points, faces, g, parameters::default_values(), np);

  return true;
}

} // namespace IO
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_PMP_IO_POLYGON_MESH_IO_H
