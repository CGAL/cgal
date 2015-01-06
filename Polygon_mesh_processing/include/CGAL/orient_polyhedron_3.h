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


#ifndef CGAL_ORIENT_POLYHEDRON_3
#define CGAL_ORIENT_POLYHEDRON_3

#include <algorithm>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

namespace CGAL {
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
 * Tests whether a closed polyhedron has a positive orientation.
 * A polyhedron is considered to have positive orientation if the normal vectors
 * of the facets point outside of the polyhedron. For each facet, its normal vector
 * is considered to point on the side of the facet where the sequence of vertices of
 * the facet is seen counterclockwise.
 * @pre @a `polyhedron`.is_closed()
 * @pre @a `polyhedron` is consistently oriented
 *
 * @tparam Polyhedron a %CGAL polyhedron
 *
 * @param polyhedron a closed polyhedron to be tested
 *
 * \todo The following only handle polyhedron with one connected component
 *       the code, the sample example and the plugin must be updated.
 * @code
 * // use inside out to fix orientation
 * if(!is_oriented(polyhedron)) {
 *   polyhedron.inside_out();
 * }
 * @endcode
 */
template<class Polyhedron>
bool is_oriented(const Polyhedron& polyhedron) {
  CGAL_precondition(polyhedron.is_closed());
  const unsigned int axis = 0;

  typename Polyhedron::Vertex_const_iterator v_min
    = std::min_element(polyhedron.vertices_begin(), polyhedron.vertices_end(), internal::Axis_compare<axis>());

  typedef typename Polyhedron::Traits K;
  const typename K::Vector_3& normal_v_min = Polygon_mesh_processing::compute_vertex_normal<K>(*v_min);

  CGAL_warning(normal_v_min[axis] != 0);
  return normal_v_min[axis] < 0;
}
} // namespace CGAL
#endif // CGAL_ORIENT_POLYHEDRON_3
