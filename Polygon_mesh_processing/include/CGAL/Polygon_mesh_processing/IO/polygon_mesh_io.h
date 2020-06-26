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

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
  \ingroup PMP_IO_grp

 * \brief Attempts to read a file as a polygon mesh; in case of failure, reads the file as a polygon soup,
 * repairs and orients it to obtain a polygon mesh.
 *
 * \tparam Graph a model of `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the input file
 * \param g the graph
 * \param np sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `Graph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{The traits class must provide the nested functors `Less_xyz_3` and `Equal_3`
 *                    to respectivelycompare lexicographically two points and to check if two points
 *                    are identical. For each functor `Foo`, a function `Foo foo_object()` must be provided.}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{repair_polygon_soup}
 *     \cgalParamDescription{a parameter used indicate whether `CGAL::Polygon_mesh_processing::repair_polygon_soup()`
 *                           should be called on the soup in case of issues in the input.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{erase_all_duplicates}
 *     \cgalParamDescription{Parameter to indicate, when multiple polygons are duplicates,
 *                           whether all the duplicate polygons should be removed
 *                           or if one (arbitrarily chosen) face should be kept.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{require_same_orientation}
 *     \cgalParamDescription{Parameter to indicate if polygon orientation should be taken
 *                           into account when determining whether two polygons are duplicates,
 *                           that is, whether e.g. the triangles `0,1,2` and `0,2,1` are duplicates.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the reading and conversion worked, `false` otherwise.
 *
 * \sa \link PkgBGLIOFct `CGAL::write_polygon_mesh()` \endlink
 */
template <typename Graph, typename NamedParameter>
bool read_polygon_mesh(const char* fname,
                       Graph& g,
                       const NamedParameter& np)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename CGAL::GetVertexPointMap<Graph, NamedParameter>::type      VPM;
  typedef typename boost::property_traits<VPM>::value_type                   Point;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool ok = CGAL::read_polygon_mesh(fname, g, np);

  if(ok)
    return true;

  clear(g);

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
    PMP::repair_polygon_soup(points, faces, np);

  if(!PMP::orient_polygon_soup(points, faces))
    std::cerr << "W: File does not describe a polygon mesh" << std::endl;

  if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    return false;

  PMP::polygon_soup_to_polygon_mesh(points, faces, g, parameters::all_default(), np);

  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename Graph>
bool read_polygon_mesh(const char* fname, Graph& g)
{
  return CGAL::Polygon_mesh_processing::read_polygon_mesh(fname, g, parameters::all_default());
}

template <typename Graph, typename NamedParameter>
bool read_polygon_mesh(const std::string& fname, Graph& g, const NamedParameter& np)
{
  return CGAL::Polygon_mesh_processing::read_polygon_mesh(fname.c_str(), g, np);
}

template <typename Graph>
bool read_polygon_mesh(const std::string& fname, Graph& g)
{
  return CGAL::Polygon_mesh_processing::read_polygon_mesh(fname, g, parameters::all_default());
}

/// \endcond

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_PMP_IO_POLYGON_MESH_IO_H
