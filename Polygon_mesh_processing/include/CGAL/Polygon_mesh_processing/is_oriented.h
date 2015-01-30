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
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Kernel_traits.h>

#include <boost/foreach.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template<unsigned int axis>
struct Axis_compare {
  template<class Point>
  bool operator()(const Point& p0, const Point& p1) const
  { return p0[axis] < p1[axis]; }
};

} // namespace internal

/**
 * \ingroup PkgPolygonMeshProcessing
 * Tests whether a closed surface polygon mesh has a positive orientation.
 * A polygon mesh is considered to have positive orientation if the normal vectors
 * of the facets point outside of the polygon mesh. For each facet, its normal vector
 * is considered to point on the side of the facet where the sequence of vertices of
 * the facet is seen counterclockwise.
 * @pre @a `pmesh` is closed
 * @pre @a `pmesh` is consistently oriented
 *
 * @tparam PolygonMesh a model of `FaceListGraph`, possibly a %CGAL polyhedron
 *
 * @param pmesh a closed polygon mesh to be tested
 *
 * \todo The following only handles polyhedron with one connected component
 *       the code, the sample example and the plugin must be updated.
 * @code
 * // use inside_out operator to reverse orientation
 * if(!is_outward_oriented(pmesh)) {
 *   pmesh.inside_out();
 * }
 * @endcode
 */
template<
  class PolygonMesh,
  typename Kernel
    = typename CGAL::Kernel_traits<typename PolygonMesh::Point>::Kernel >
bool is_outward_oriented(const PolygonMesh& pmesh,
                         const Kernel& = Kernel())
{
  CGAL_precondition(CGAL::is_closed(pmesh));
  CGAL_precondition(CGAL::is_valid(pmesh));

  const unsigned int axis = 0;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type
    ppmap = get(boost::vertex_point, pmesh);

  vertex_descriptor v_min = *vertices(pmesh).first;
  BOOST_FOREACH(vertex_descriptor vd, vertices(pmesh)) {
    if(internal::Axis_compare<axis>()(ppmap[vd], ppmap[v_min]))
      v_min = vd;
  }

  const typename Kernel::Vector_3&
    normal_v_min = compute_vertex_normal<Kernel>(v_min, pmesh);

  CGAL_warning(normal_v_min[axis] != 0);
  return normal_v_min[axis] < 0;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H

