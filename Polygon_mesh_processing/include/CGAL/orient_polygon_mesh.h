// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Ilker O. Yaz


#ifndef CGAL_ORIENT_POLYGON_MESH_H
#define CGAL_ORIENT_POLYGON_MESH_H

#include <algorithm>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template<unsigned int axis>
struct Axis_compare {
  template<class Vertex>
  bool operator()(const Vertex& v0, const Vertex& v1) const
  { return v0.point()[axis] < v1.point()[axis]; }
};

} // namespace internal

/**
 * \ingroup PkgPolygonMeshProcessing
 * Tests whether a closed surface polygon mesh has a positive orientation.
 * A polygon mesh is considered to have positive orientation if the normal vectors
 * of the facets point outside of the polygon mesh. For each facet, its normal vector
 * is considered to point on the side of the facet where the sequence of vertices of
 * the facet is seen counterclockwise.
 * @pre @a `pmesh`.is_closed()
 * @pre @a `pmesh` is consistently oriented
 *
 * @tparam PolygonMesh a %CGAL polyhedron
 *
 * @param pmesh a closed polygon mesh to be tested
 *
 * \todo The following only handle polyhedron with one connected component
 *       the code, the sample example and the plugin must be updated.
 * \todo PolygonMesh should be a model of `FaceListGraph`
 * \todo find a better name for the function
 * @code
 * // use inside out to reverse orientation
 * if(!is_oriented(pmesh)) {
 *   pmesh.inside_out();
 * }
 * @endcode
 */
template<class PolygonMesh>
bool is_oriented(const PolygonMesh& pmesh)
{
  CGAL_precondition(pmesh.is_closed());
  const unsigned int axis = 0;

  typename Polyhedron::Vertex_const_iterator v_min
    = std::min_element(pmesh.vertices_begin(),
                       pmesh.vertices_end(),
                       internal::Axis_compare<axis>());

  typedef typename Polyhedron::Traits K;
  const typename K::Vector_3& normal_v_min
    = compute_vertex_normal<K>(*v_min);

  CGAL_warning(normal_v_min[axis] != 0);
  return normal_v_min[axis] < 0;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H

