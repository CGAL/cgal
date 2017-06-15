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

#include <CGAL/license/Polygon_mesh_processing/orientation.h>


#include <algorithm>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal{

  template <class GT, class VPmap>
  struct Compare_vertex_points_z_3
  {
    VPmap vpmap;
    typename GT::Compare_z_3 compare_z;

    Compare_vertex_points_z_3(VPmap const& vpmap, const GT& gt)
      : vpmap(vpmap)
      , compare_z(gt.compare_z_3_object())
    {}

    typedef bool result_type;
    template <class vertex_descriptor>
    bool operator()(vertex_descriptor v1, vertex_descriptor v2) const
    {
      return CGAL::SMALLER == compare_z(get(vpmap, v1), get(vpmap, v2));
    }
  };


  template<typename PolygonMesh, typename NamedParameters>
  bool is_outward_oriented(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_max,
                           const PolygonMesh& pmesh,
                           const NamedParameters& np)
  {
    using boost::choose_param;
    using boost::get_param;

    //VertexPointMap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPMap;
    VPMap vpmap = choose_param(get_param(np, vertex_point),
                               get_const_property_map(vertex_point, pmesh));
    //Kernel
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    GT gt = choose_param(get_param(np, internal_np::geom_traits), GT());

    //among the incoming edges of `v_max`, find one edge `e` with the minimal slope
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
    halfedge_descriptor min_slope_he = halfedge(v_max, pmesh);
    CGAL_assertion(v_max == target(min_slope_he, pmesh));

    typename GT::Compare_slope_3 compare_slope = gt.compare_slope_3_object();
    BOOST_FOREACH(halfedge_descriptor he, halfedges_around_target(v_max, pmesh))
    {
      CGAL_assertion(v_max == target(min_slope_he, pmesh));
      CGAL_assertion(v_max == target(he, pmesh));

      if(CGAL::SMALLER == compare_slope(get(vpmap, source(he, pmesh)),
                                        get(vpmap, v_max),
                                        get(vpmap, source(min_slope_he, pmesh)),
                                        get(vpmap, v_max)))
      {
        min_slope_he = he;
      }
    }

    // We compute the orientations of the two triangles incident to the edge
    // of `min_slope_he` projected in the xy-plane. We can conclude using
    // the 2D orientation of the 3D triangle that is the top one along the z-axis
    // in the neighborhood of `min_slope_he`.
    Projection_traits_xy_3<GT> p_gt;
    typename Projection_traits_xy_3<GT>::Orientation_2 orientation_2 = p_gt.orientation_2_object();

    typename boost::property_traits<VPMap>::reference p1 = get(vpmap, source(min_slope_he, pmesh));
    typename boost::property_traits<VPMap>::reference p2 = get(vpmap, target(min_slope_he, pmesh));
    typename boost::property_traits<VPMap>::reference p3 = get(vpmap, target(next(min_slope_he, pmesh), pmesh));
    typename boost::property_traits<VPMap>::reference p4 = get(vpmap, target(next(opposite(min_slope_he, pmesh), pmesh), pmesh));

    Orientation p1p2p3_2d = orientation_2(p1, p2, p3);
    Orientation p2p1p4_2d = orientation_2(p2, p1, p4);

    CGAL_assertion( p1p2p3_2d!=COLLINEAR || p2p1p4_2d!=COLLINEAR ); // no self-intersection

    if ( p1p2p3_2d == COLLINEAR)
      return p2p1p4_2d == LEFT_TURN;
    if (p2p1p4_2d ==COLLINEAR)
      return p1p2p3_2d == LEFT_TURN;

    // if the local dihedral angle is strictly larger that PI/2, we can conclude with any of two triangles
    if (p1p2p3_2d==p2p1p4_2d)
      return p1p2p3_2d == LEFT_TURN;

    typename GT::Orientation_3 orientation_3 = gt.orientation_3_object();

    CGAL_assertion( orientation_3(p1, p2, p3, p4) != COPLANAR ); // same side of min_slope_he and no self-intersection

    // if p1p2p3_2d is left turn, then it must be the top face so that the orientation is outward oriented
    if (p1p2p3_2d == LEFT_TURN)
      return orientation_3(p1, p2, p3, p4) == NEGATIVE;

    // same test with the other face
    CGAL_assertion(p2p1p4_2d == LEFT_TURN);
    return orientation_3(p2, p1, p4, p3) == NEGATIVE;
  }
} // end of namespace internal

/**
 * \ingroup PMP_orientation_grp
 * tests whether a closed polygon mesh has a positive orientation.
 * A closed polygon mesh is considered to have a positive orientation if the normal vectors
 * to all its faces point outside the domain bounded by the polygon mesh.
 * The normal vector to each face is chosen pointing on the side of the face
 * where its sequence of vertices is seen counterclockwise.
 * @pre `CGAL::is_closed(pmesh)`
 * @pre If `pmesh` contains several connected components, they are oriented consistently.
 *      In other words, the answer to this predicate would be the same for each
 *      isolated connected component.
 *
 * @tparam PolygonMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param pmesh the closed polygon mesh to be tested
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
 * \cgalNamedParamsEnd
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

  //check for empty pmesh
  CGAL_warning(faces(pmesh).first != faces(pmesh).second);
  if (faces(pmesh).first == faces(pmesh).second)
    return true;

  using boost::choose_param;
  using boost::get_param;

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));
  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  GT gt = choose_param(get_param(np, internal_np::geom_traits), GT());

  //find the vertex with maximal z coordinate
  typename boost::graph_traits<PolygonMesh>::vertex_iterator vbegin, vend;
  cpp11::tie(vbegin, vend) = vertices(pmesh);

  internal::Compare_vertex_points_z_3<GT, VPMap> less_z(vpmap, gt);
  typename boost::graph_traits<PolygonMesh>::vertex_iterator v_max_it
    = std::max_element(vbegin, vend, less_z);
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_max = *v_max_it;

  return internal::is_outward_oriented(v_max, pmesh, np);
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
* \ingroup PMP_orientation_grp
* reverses for each face the order of the vertices along the face boundary.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
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

// Do the same thing as `reverse_face_orientations()` except that for
// the reversal of the border cycles (last step in the aforementioned function),
// this function guarantees that each cycle is reversed only once. This is
// particularly useful if you mesh contains polylines (i.e. edge which halfedges
// are both border halfedges).
template<typename PolygonMesh>
void reverse_face_orientations_of_mesh_with_polylines(PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  // reverse the orientation of each face
  BOOST_FOREACH(face_descriptor fd, faces(pmesh))
    reverse_orientation(halfedge(fd,pmesh),pmesh);

  //extract all border cycles
  boost::unordered_set<halfedge_descriptor> already_seen;
  std::vector<halfedge_descriptor> border_cycles;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
    if ( is_border(h,pmesh) && already_seen.insert(h).second )
    {
      border_cycles.push_back(h);
      BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_face(h,pmesh))
        already_seen.insert(h2);
    }

  // now reverse the border cycles
  BOOST_FOREACH(halfedge_descriptor h, border_cycles)
    reverse_orientation(h, pmesh);
}

/**
* \ingroup PMP_orientation_grp
* reverses for each face in `face_range` the order of the vertices along the face boundary.
* The function does not perform any control and if the orientation change of the faces
* makes the polygon mesh invalid, the behavior is undefined.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam FaceRange range of face descriptors, model of `Range`.
*         Its iterator type is `InputIterator`.
*/
template<typename PolygonMesh, typename FaceRange>
void reverse_face_orientations(const FaceRange& face_range, PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  BOOST_FOREACH(face_descriptor fd, face_range){
    reverse_orientation(halfedge(fd,pmesh),pmesh);
  }

  // Note: A border edge is now parallel to its opposite edge.
  // We scan all border edges for this property. If it holds, we
  // reorient the associated hole and search again until no border
  // edge with that property exists any longer. Then, all holes are
  // reoriented.
  BOOST_FOREACH(face_descriptor fd, face_range)
    BOOST_FOREACH(halfedge_descriptor hd,
                  halfedges_around_face(halfedge(fd, pmesh), pmesh))
    {
      halfedge_descriptor ohd = opposite(hd, pmesh);
      if ( is_border(ohd, pmesh) &&
         target(hd,pmesh) == target(ohd,pmesh))
      {
        reverse_orientation(ohd, pmesh);
      }
    }
}

} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H
