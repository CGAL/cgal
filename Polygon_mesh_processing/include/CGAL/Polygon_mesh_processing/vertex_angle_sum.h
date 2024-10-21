// Copyright (c) 2024 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_VERTEX_ANGLE_SUM_H
#define CGAL_POLYGON_MESH_PROCESSING_VERTEX_ANGLE_SUM_H

#include <CGAL/license/Polygon_mesh_processing/measure.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/utils_classes.h>

#include <CGAL/Lazy.h> // needed for CGAL::exact(FT)/CGAL::exact(Lazy_exact_nt<T>)

#include <boost/graph/graph_traits.hpp>

#include <utility>

namespace CGAL {



namespace Polygon_mesh_processing {

/**
  * \ingroup PMP_vertex_angle_grp
  *
  * computes the sum of the angles around a vertex
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param v the vertex whose sum of angles is computed
  * @param pmesh the polygon mesh to which `v` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the sum of anglesof `v`. The return type `FT` is a number type either deduced
  * from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `pmesh`.
  *
  * \warning This function involves trigonometry.
  *
  */
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
#endif
vertex_angle_sum(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
                 const PolygonMesh& pmesh,
                 const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Geom_traits;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(is_valid_vertex_descriptor(v, pmesh));

  typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
  using Point_3 = typename Geom_traits::Point_3;  // AF: Should it be  the value type of vpm instead"
  using FT = typename Geom_traits::FT;

  FT res(0);
  for(auto h : halfedges_around_source(v,pmesh)){
    res += gt.compute_approximate_angle_3_object()(get(vpm, target(h, pmesh)),
                                                   get(vpm, source(h, pmesh)),
                                                   get(vpm, source(prev(h,pmesh), pmesh)));
  }
  return  res;
}


/**
* \ingroup PMP_vertex_normal_grp
*
* computes the outward unit vector normal for all vertices of the polygon mesh.
*
* @tparam PolygonMesh a model of `FaceListGraph`
* @tparam VertexAngleSumMap a model of `WritablePropertyMap` with
*                         `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
*                         the field type of the geometric traits class as value type.
*
* @param pmesh the polygon mesh
* @param vertex_normals the property map in which the sums of angles are written
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \warning This function involves trigonometry.
*
* @see `compute_vertex_normal()`
*/
template <typename PolygonMesh, typename VertexAngleSumMap, typename NamedParameters = parameters::Default_named_parameters>
void vertex_angle_sums(const PolygonMesh& pmesh,
                       VertexAngleSumMap vertex_angle_sums,
                       const NamedParameters& np = parameters::default_values())
{
   for(auto v : vertices(pmesh)){
      put(vertex_angle_sums, v, vertex_angle_sum(v,pmesh, np));
   }

}



} // namespace Polygon_mesh_processing
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_VERTEX_ANGLE_SUM_H
