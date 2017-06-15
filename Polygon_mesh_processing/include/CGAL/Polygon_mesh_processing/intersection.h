// Copyright (c) 2016 GeometryFactory (France).
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
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace Polygon_mesh_processing{

/**
 * \ingroup PMP_corefinement_grp
 * computes the intersection of triangles of `tm1` and `tm2`. The output is a
 * set of polylines with all vertices but endpoints being of degree 2.
 *
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref namedparameters
 * @tparam NamedParameters2 a sequence of \ref namedparameters
 * @tparam OutputIterator an output iterator in which `std::vector` of points
 *                        can be put. The point type is the one from the
 *                        vertex property map
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param polyline_output output iterator of polylines. Each polyline will be
 *        given as a vector of points
 * @param throw_on_self_intersection if `true`, for each input triangle mesh,
 *        the set of triangles closed to the intersection of `tm1` and `tm2` will be
 *        checked for self-intersection and `CGAL::Corefinement::Self_intersection_exception`
 *        will be thrown if at least one is found.
 * @param np1 optional sequence of \ref namedparameters among the ones listed below
 * @param np2 optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *    a property map with the points associated to the vertices of `tm1`
 *    (`tm2`). The two property map types must be the same.
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
template <class OutputIterator,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2 >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const NamedParameters1& np1,
                     const NamedParameters2& np2,
                     const bool throw_on_self_intersection=false)
{
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::const_type Vpm;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters2>::const_type Vpm2;
  CGAL_USE_TYPE(Vpm2);
  CGAL_assertion_code(
    static const bool same_vpm = (boost::is_same<Vpm,Vpm2>::value);)
  CGAL_static_assertion(same_vpm);

  Vpm vpm1 = choose_const_pmap(get_param(np1, internal_np::vertex_point),
                               tm1,
                               vertex_point);
  Vpm vpm2 = choose_const_pmap(get_param(np2, internal_np::vertex_point),
                               tm2,
                               vertex_point);

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh,Vpm>
    functor(tm1, tm2, vpm1, vpm2);
  functor(polyline_output, throw_on_self_intersection, true);
  return polyline_output;
}

template <class OutputIterator,
          class TriangleMesh >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const bool throw_on_self_intersection=false)
{
  return surface_intersection(tm1, tm2, polyline_output,
    CGAL::Polygon_mesh_processing::parameters::all_default(),
    CGAL::Polygon_mesh_processing::parameters::all_default(),
    throw_on_self_intersection
  );
}

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
