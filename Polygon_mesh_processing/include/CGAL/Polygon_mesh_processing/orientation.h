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
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
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

   Compare_vertex_points_xyz_3(VPmap const& vpmap)
	: vpmap(vpmap){}

    typedef bool result_type;
    template <class vertex_descriptor>
    bool operator()(vertex_descriptor v1, vertex_descriptor v2) const
    {
      return less(get(vpmap, v1), get(vpmap, v2));
    }

  };
} // end of namespace internal

void toto(int i)
{
  int j = i;
}
/**
 * \ingroup PkgPolygonMeshProcessing
 * tests whether a closed polygon mesh has a positive orientation.
 * A closed polygon mesh is considered to have a positive orientation if the normal vectors
 * of its faces point outside the domain bounded by the polygon mesh. For each face, its normal vector
 * is considered to point on the side of the face where the sequence of vertices of
 * the facet is seen counterclockwise.
 * @pre `CGAL::is_closed(pmesh)`
 * @pre If `pmesh` contains several connected components they are oriented consistentl,
 *      that is the answer to this predicate would be the same if called on each
 *       connected component isolated.
 *
 * @tparam PolygonMesh a model of `FaceListGraph` that has a property map
*         for `CGAL::vertex_point_t`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh the closed polygon mesh to be tested
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 \b Named \b parameters
 <ul>
 <li>\b vertex_point_map the property map with the points associated to the vertices of `pmesh`
 <li>\b kernel a geometric traits class instance
 </ul>
 *
 * \todo code : The following only handles polyhedron with one connected component
 *       the code, the sample example and the plugin must be updated.
 *
 * \sa `CGAL::Polygon_mesh_processing::reverse_face_orientations()`
 */
template<typename PolygonMesh, typename NamedParameters>
bool is_outward_oriented(const PolygonMesh& pmesh,
                         const NamedParameters& np)
{
  CGAL_warning(CGAL::is_closed(pmesh));
  CGAL_precondition(CGAL::is_valid(pmesh));

  using boost::choose_const_pmap;
  using boost::get_param;

  //VertexPointMap
  typedef typename boost::lookup_named_param_def <
    boost::vertex_point_t,
    NamedParameters,
    boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type//default
  > ::type  VPMap;
  VPMap vpmap = choose_const_pmap(get_param(np, boost::vertex_point),
                                  pmesh,
                                  boost::vertex_point);
  //Kernel
  typedef typename CGAL::Kernel_traits <
    typename property_map_value<PolygonMesh, CGAL::vertex_point_t>::type
  > ::Kernel DefaultKernel;
  typedef typename boost::lookup_named_param_def <
    CGAL::geom_traits_t,
    NamedParameters,
    DefaultKernel //default
  > ::type  Kernel;

  internal::Compare_vertex_points_xyz_3<typename Kernel::Less_xyz_3, VPMap >
    less_xyz(vpmap);

  typename boost::graph_traits<PolygonMesh>::vertex_iterator vbegin, vend;
  cpp11::tie(vbegin, vend) = vertices(pmesh);
  typename boost::graph_traits<PolygonMesh>::vertex_iterator v_min
    = std::min_element(vbegin, vend, less_xyz);

  const typename Kernel::Vector_3&
    normal_v_min = compute_vertex_normal(*v_min, pmesh, np);

  return normal_v_min[0] < 0 || (
            normal_v_min[0] == 0 && (
              normal_v_min[1] < 0  ||
              ( normal_v_min[1]==0  && normal_v_min[2] < 0 )
            )
         );
}

///\cond SKIP_IN_MANUAL

template<typename PolygonMesh>
bool is_outward_oriented(const PolygonMesh& pmesh)
{
  return is_outward_oriented(pmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

template<typename PolygonMesh>
void reverse_orientation(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor first, PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    if ( first == halfedge_descriptor())
        return;
    halfedge_descriptor last  = first;
    halfedge_descriptor prev  = first;
    halfedge_descriptor start = first;
    first = next(first, pmesh);
    vertex_descriptor  new_v = target( start, pmesh);
    while (first != last) {
      vertex_descriptor  tmp_v = target( first, pmesh);
      set_target( first, new_v, pmesh);
      set_halfedge(new_v, first, pmesh);
        new_v = tmp_v;
        halfedge_descriptor n = next(first, pmesh);
        set_next(first, prev, pmesh);
        prev  = first;
        first = n;
    }
    set_target( start, new_v, pmesh);
    set_halfedge( new_v, start, pmesh);
    set_next(start, prev,pmesh);
}

/**
* \ingroup PkgPolygonMeshProcessing
* reverses for each face the order of the vertices along the face boundary.
*
* @tparam PolygonMesh a model of `FaceListGraph`
*/
template<typename PolygonMesh>
void reverse_face_orientations(PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(face_descriptor fd, faces(pmesh)){
    reverse_orientation(halfedge(fd,pmesh),pmesh);
  } 
  // Note: A border edge is now parallel to its opposite edge.
  // We scan all border edges for this property. If it holds, we
  // reorient the associated hole and search again until no border
  // edge with that property exists any longer. Then, all holes are
  // reoriented.
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh)){
    if ( is_border(h,pmesh) &&
         target(h,pmesh) == target(opposite(h,pmesh),pmesh)){
      reverse_orientation(h, pmesh);
    }
  }
}

} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H

