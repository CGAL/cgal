// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Maxime Gimeno
#ifndef CGAL_POLYGON_MESH_PROCESSING_TRANSFORM_H
#define CGAL_POLYGON_MESH_PROCESSING_TRANSFORM_H
#include <CGAL/license/Polygon_mesh_processing/core.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL{
namespace Polygon_mesh_processing{
/**
 * \ingroup PkgPolygonMeshProcessing
 * applies a transformation to every vertex of a `PolygonMesh`.
 * 
 * @tparam Transformation a functor that has an `operator()(Point_3)`, with `Point_3`
 * the `value_type` of `vertex_point_map` (see below). Such a functor can be
 * `CGAL::Aff_transformation_3` for example.
 * @tparam PolygonMesh a model of `VertexListGraph`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 * 
 * @param transformation the transformation functor to apply to  the points of `mesh`.
 * @param mesh the `PolygonMesh` to transform.
 * @param np optional sequence of \ref pmp_namedparameters for `mesh`, among the ones listed below
 * 
 * * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `mesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` must be available in `PolygonMesh`\cgalParamEnd
 * \cgalNamedParamsEnd
 * 
 */
template<class Transformation, class PolygonMesh,class NamedParameters>
void transform(const Transformation& transformation, 
               PolygonMesh& mesh,
               const NamedParameters& np)
{
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VPMap;
  VPMap vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                           get_property_map(vertex_point, mesh));
  
  BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd, vertices(mesh))
  {
    put(vpm, vd, transformation(get(vpm, vd)));
  }
}

/// \cond SKIP_IN_MANUAL
template<class Transformation, class PolygonMesh>
void transform(const Transformation& transformation,
               PolygonMesh& mesh)
{
  transform(transformation, mesh, parameters::all_default());
}
/// \endcond
}
}

#endif
