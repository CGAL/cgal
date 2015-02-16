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

namespace internal{
  template <class Less_xyz, class VPmap>
  struct Compare_vertex_points_xyz_3{
    Less_xyz less;
    VPmap vpmap;

   Compare_vertex_points_xyz_3(VPmap vpmap)
	: vpmap(vpmap){}

    typedef bool result_type;
    template <class vertex_descriptor>
    bool operator()(vertex_descriptor v1, vertex_descriptor v2) const
    {
      return less(get(vpmap, v1), get(vpmap, v2));
    }

  };
} // end of namespace internal
#ifndef DOXYGEN_RUNNING
template <class PolygonMesh>
bool
is_outward_oriented(
  const PolygonMesh& pmesh)
{
  typedef typename Kernel_traits<typename boost::property_traits<typename boost::property_map<PolygonMesh,CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  return is_outward_oriented(pmesh, Kernel());
}
#endif
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
  template<class Kernel, class PolygonMesh >
  bool is_outward_oriented(const PolygonMesh& pmesh, const Kernel&k)
{
  CGAL_warning(CGAL::is_closed(pmesh));
  CGAL_precondition(CGAL::is_valid(pmesh));

  typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type VPmap;
  VPmap ppmap = get(boost::vertex_point, pmesh);

  internal::Compare_vertex_points_xyz_3< typename Kernel::Less_xyz_3, VPmap > less_xyz(ppmap);

  typename boost::graph_traits<PolygonMesh>::vertex_iterator vbegin, vend;
  cpp11::tie(vbegin, vend) = vertices(pmesh);
  typename boost::graph_traits<PolygonMesh>::vertex_iterator v_min
    = std::min_element(vbegin, vend, less_xyz);

  const typename Kernel::Vector_3&
    normal_v_min = compute_vertex_normal(*v_min, pmesh, k);

  return normal_v_min[0] < 0 || (
            normal_v_min[0] == 0 && (
              normal_v_min[1] < 0  ||
              ( normal_v_min[1]==0  && normal_v_min[2] < 0 )
            )
         );
}

} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H

