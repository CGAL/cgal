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

/**
 * \ingroup PkgPolygonMeshProcessing
 * Tests whether a closed surface polygon mesh has a positive orientation.
 * A polygon mesh is considered to have positive orientation if the normal vectors
 * of the facets point outside the domain described by the polygon mesh. For each facet, its normal vector
 * is considered to point on the side of the facet where the sequence of vertices of
 * the facet is seen counterclockwise.
 * @pre @a `pmesh` is closed
 * @pre @a `pmesh` is consistently oriented
 *
 * @tparam PolygonMesh a model of `FaceListGraph`
 * @tparam VertexPointMap a model of `ReadablePropertyMap` with
    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
    `Kernel::Vector_3` as value type
 * @tparam Kernel Geometric traits class. It can be omitted and deduced automatically from the point type of `PolygonMesh`.
 *
 * @param pmesh a closed polygon mesh to be tested
 * @param vpmap the property map with the points associated to the vertices of `pmesh`
 * @param k a traits class instance, can be omitted
 *
 * \todo The following only handles polyhedron with one connected component
 *       the code, the sample example and the plugin must be updated.
 * @code
 * if(!is_outward_oriented(pmesh)) {
 *   reverse_face_orientations(pmesh);
 * }
 * @endcode
 */
template<typename PolygonMesh
       , typename VertexPointMap
#ifdef DOXYGEN_RUNNING
       = typename boost::property_map<PolygonMesh, vertex_point_t>::type
#endif
       , typename Kernel
#ifdef DOXYGEN_RUNNING
       = typename Kernel_traits<
           typename boost::property_traits<VertexPointMap>::value_type>::Kernel
#endif
         >
bool is_outward_oriented(const PolygonMesh& pmesh
                       , VertexPointMap vpmap
#ifdef DOXYGEN_RUNNING
                       = get(vertex_point, pmesh)
#endif
                       , const Kernel& k
#ifdef DOXYGEN_RUNNING
                       = Kernel()
#endif
                       )
{
  CGAL_warning(CGAL::is_closed(pmesh));
  CGAL_precondition(CGAL::is_valid(pmesh));

  internal::Compare_vertex_points_xyz_3<typename Kernel::Less_xyz_3,
                                        VertexPointMap > less_xyz(vpmap);

  typename boost::graph_traits<PolygonMesh>::vertex_iterator vbegin, vend;
  cpp11::tie(vbegin, vend) = vertices(pmesh);
  typename boost::graph_traits<PolygonMesh>::vertex_iterator v_min
    = std::min_element(vbegin, vend, less_xyz);

  const typename Kernel::Vector_3&
    normal_v_min = compute_vertex_normal(*v_min, pmesh, vpmap, k);

  return normal_v_min[0] < 0 || (
            normal_v_min[0] == 0 && (
              normal_v_min[1] < 0  ||
              ( normal_v_min[1]==0  && normal_v_min[2] < 0 )
            )
         );
}

///\cond SKIP_IN_MANUAL

template<typename PolygonMesh
       , typename VertexPointMap>
bool is_outward_oriented(const PolygonMesh& pmesh
                       , VertexPointMap vpmap)
{
  typedef typename Kernel_traits<
      typename boost::property_traits<
        typename boost::property_map<PolygonMesh, 
          CGAL::vertex_point_t>::type>::value_type>::Kernel Kernel;
  return is_outward_oriented(pmesh, vpmap, Kernel());
}

template<typename PolygonMesh>
bool is_outward_oriented(const PolygonMesh& pmesh)
{
  return is_outward_oriented(pmesh, get(boost::vertex_point, pmesh));
}

/// \endcond

/**
* \ingroup PkgPolygonMeshProcessing
* reverses faces orientations.

* @pre @a `pmesh` is closed
* @pre @a `pmesh` is consistently oriented
*
* @tparam PolygonMesh a model of `FaceListGraph`
*
* @param pmesh a closed polygon mesh to be tested
* @param vpmap the property map with the points associated to the vertices of `pmesh`
* @param k a traits class instance, can be omitted
*
* \todo write the code. Can be copy-pasted and BGLized from function
* inside_out in Polyhedron_3
*/
template<typename PolygonMesh>
void reverse_face_orientations(PolygonMesh& pmesh)
{
  ;
}

} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H

