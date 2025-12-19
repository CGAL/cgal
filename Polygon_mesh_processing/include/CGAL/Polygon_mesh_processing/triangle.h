// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_H

#include <CGAL/license/Polygon_mesh_processing/core.h>


#include <boost/graph/graph_traits.hpp>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {

  namespace Polygon_mesh_processing {

  /*!
   * \ingroup PkgPolygonMeshProcessingRef
   *
   * computes a triangle for a face descriptor of a triangle mesh.
   *
   * @tparam TriangleMesh a model of `FaceGraph`
   * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * @param tmesh the triangulated surface mesh
   * @param fd the descriptor of the face to construct the triangle from
   * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{vertex_point_map}
   *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
   *                    as key type and `%Point_3` as value type}
   *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
   *   \cgalParamNEnd
   *
   *   \cgalParamNBegin{geom_traits}
   *     \cgalParamDescription{an instance of a geometric traits class}
   *     \cgalParamType{a class model of `Kernel`}
   *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
   *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   */
template<typename TriangleMesh, typename  CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
Triangle_3
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::Triangle_3
#endif
triangle(typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
         const TriangleMesh& tmesh,
         const CGAL_NP_CLASS& np = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  CGAL_precondition(is_valid_face_descriptor(fd, tmesh));

  typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
    vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                           get_const_property_map(CGAL::vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  typename GT::Construct_triangle_3 construct_triangle = gt.construct_triangle_3_object();

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = halfedge(fd, tmesh);
  vertex_descriptor v1  = target(h, tmesh);
  vertex_descriptor v2  = target(next(h, tmesh), tmesh);
  vertex_descriptor v3  = target(next(next(h, tmesh), tmesh), tmesh);
  return construct_triangle(get(vpm, v1), get(vpm, v2), get(vpm, v3));
}

}
}
#endif //CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_H
