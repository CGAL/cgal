// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
 * \ingroup PkgPolygonMeshProcessingRef
 * applies a transformation to every vertex of a `PolygonMesh`.
 *
 * @tparam Transformation a functor that has an `operator()(Point_3)`, with `Point_3`
 * the `value_type` of `vertex_point_map` (see below). Such a functor can be
 * `CGAL::Aff_transformation_3` for example.
 * @tparam PolygonMesh a model of `VertexListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param transformation the transformation functor to apply to  the points of `mesh`.
 * @param mesh the `PolygonMesh` to transform.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
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

  for(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd : vertices(mesh))
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
