// Copyright (c) 2013 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_PMP_ORIENTATION_BASIC_H
#define CGAL_PMP_ORIENTATION_BASIC_H

#include <CGAL/license/Polygon_mesh_processing/orientation.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/utility.h>
#include <unordered_set>

#include <algorithm>

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
  template <class vertex_descriptor1, class vertex_descriptor2>
  bool operator()(vertex_descriptor1 v1, vertex_descriptor2 v2) const
  {
    return CGAL::SMALLER == compare_z(get(vpmap, v1), get(vpmap, v2));
  }
};


template<typename PolygonMesh, typename NamedParameters>
bool is_outward_oriented(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_max,
                         const PolygonMesh& pmesh,
                         const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(halfedge(v_max, pmesh)!=boost::graph_traits<PolygonMesh>::null_halfedge());

  //VertexPointMap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_const_property_map(vertex_point, pmesh));
  //Kernel
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  //among the incoming edges of `v_max`, find one edge `e` with the minimal slope
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor min_slope_he = halfedge(v_max, pmesh);
  CGAL_assertion(v_max == target(min_slope_he, pmesh));

  typename GT::Compare_slope_3 compare_slope = gt.compare_slope_3_object();
  for(halfedge_descriptor he : halfedges_around_target(v_max, pmesh))
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
 *
 * \brief tests whether a closed triangle mesh has a positive orientation.
 *
 * A closed triangle mesh is considered to have a positive orientation if the normal vectors
 * to all its faces point outside the domain bounded by the triangle mesh.
 * The normal vector to each face is chosen pointing on the side of the face
 * where its sequence of vertices is seen counterclockwise.
 *
 * @pre \link CGAL::is_closed `CGAL::is_closed(tm)` \endlink
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
 * @pre If `tm` contains several connected components, they are oriented consistently.
 *      In other words, the answer to this predicate would be the same for each
 *      isolated connected component.
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm the closed triangle mesh free from self-intersections to be tested
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
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
 * \note This function is only doing an orientation test for one connected component of `tm`.
 *       For performance reasons, it is left to the user to call the function `does_bound_a_volume()`
 *       on a triangulated version of `tm` to ensure the result returned is relevant.
 *       For advanced usages, the function `volume_connected_components()` should be used instead.
 *
 * \sa `CGAL::Polygon_mesh_processing::reverse_face_orientations()`
 */
template<typename TriangleMesh,
         typename NamedParameters = parameters::Default_named_parameters>
bool is_outward_oriented(const TriangleMesh& tm,
                         const NamedParameters& np = parameters::default_values())
{
  CGAL_warning(CGAL::is_closed(tm));
  CGAL_warning(CGAL::is_triangle_mesh(tm));
  CGAL_precondition(CGAL::is_valid_polygon_mesh(tm));

#ifdef CGAL_PMP_DEBUG_CODE
  //check for empty tm
  CGAL_warning(faces(tm).first != faces(tm).second);
#endif

  if (faces(tm).first == faces(tm).second)
    return true;


  using parameters::choose_parameter;
  using parameters::get_parameter;

  //VertexPointMap
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_const_property_map(vertex_point, tm));
  //Kernel
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  //find the vertex with maximal z coordinate
  internal::Compare_vertex_points_z_3<GT, VPMap> less_z(vpmap, gt);
  typename boost::graph_traits<TriangleMesh>::vertex_descriptor v_max = *(vertices(tm).first);
  for (typename boost::graph_traits<TriangleMesh>::vertex_iterator
          vit=std::next(vertices(tm).first), vit_end = vertices(tm).second;
          vit!=vit_end; ++vit)
  {
    // skip isolated vertices
    if (halfedge(*vit, tm)==boost::graph_traits<TriangleMesh>::null_halfedge())
      continue;
    if( less_z(v_max, *vit) )
      v_max=*vit;
  }

  // only isolated vertices
  if (halfedge(v_max, tm)==boost::graph_traits<TriangleMesh>::null_halfedge())
    return true;

  return internal::is_outward_oriented(v_max, tm, np);
}

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
*
* reverses for each face the order of the vertices along the face boundary.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
*
* @sa `is_outward_oriented()`
*/
template<typename PolygonMesh>
void reverse_face_orientations(PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  for(face_descriptor fd : faces(pmesh)){
    reverse_orientation(halfedge(fd,pmesh),pmesh);
  }
  // Note: A border edge is now parallel to its opposite edge.
  // We scan all border edges for this property. If it holds, we
  // reorient the associated hole and search again until no border
  // edge with that property exists any longer. Then, all holes are
  // reoriented.
  for(halfedge_descriptor h : halfedges(pmesh)){
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
  for(face_descriptor fd : faces(pmesh))
    reverse_orientation(halfedge(fd,pmesh),pmesh);

  //extract all border cycles
  std::unordered_set<halfedge_descriptor> already_seen;
  std::vector<halfedge_descriptor> border_cycles;
  for(halfedge_descriptor h : halfedges(pmesh))
    if ( is_border(h,pmesh) && already_seen.insert(h).second )
    {
      border_cycles.push_back(h);
      for(halfedge_descriptor h2 : halfedges_around_face(h,pmesh))
        already_seen.insert(h2);
    }

  // now reverse the border cycles
  for(halfedge_descriptor h : border_cycles)
    reverse_orientation(h, pmesh);
}

/**
* \ingroup PMP_orientation_grp
*
* reverses for each face in `face_range` the order of the vertices along the face boundary.
* The function does not perform any control and if the change of orientation of the faces
* makes the polygon mesh invalid, the behavior is undefined.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam FaceRange range of face descriptors, model of `Range`.
*         Its iterator type is `InputIterator`.
*
* @sa `is_outward_oriented()`
*/
template<typename PolygonMesh, typename FaceRange>
void reverse_face_orientations(const FaceRange& face_range, PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  for(face_descriptor fd : face_range){
    reverse_orientation(halfedge(fd,pmesh),pmesh);
  }

  // Note: A border edge is now parallel to its opposite edge.
  // We scan all border edges for this property. If it holds, we
  // reorient the associated hole and search again until no border
  // edge with that property exists any longer. Then, all holes are
  // reoriented.
  for(face_descriptor fd : face_range)
    for(halfedge_descriptor hd :
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

/**
* \ingroup PMP_orientation_grp
* makes each closed connected component of a triangulated surface mesh
* inward or outward oriented. If a connected component is not closed,
* the orientation may or may not be changed or not is not guaranteed.
*
* @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters
*
* @param tm a closed triangulated surface mesh
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* @pre \link CGAL::is_closed `CGAL::is_closed(tm)` \endlink
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `tm`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{face_index_map}
*     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
*                    as key type and `std::size_t` as value type}
*     \cgalParamDefault{an automatically indexed internal map}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{outward_orientation}
*     \cgalParamDescription{If `true`, each connected component will be outward oriented (and inward oriented if `false`).}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template<class TriangleMesh,
         class NamedParameters = parameters::Default_named_parameters>
void orient(TriangleMesh& tm,
            const NamedParameters& np = parameters::default_values())
{
  typedef boost::graph_traits<TriangleMesh>                                        Graph_traits;
  typedef typename Graph_traits::vertex_descriptor                                 vertex_descriptor;
  typedef typename Graph_traits::face_descriptor                                   face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor                               halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type    Vpm;
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::type FaceIndexMap;

  CGAL_precondition(is_triangle_mesh(tm));
  CGAL_precondition(is_valid_polygon_mesh(tm));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool orient_outward = choose_parameter(get_parameter(np, internal_np::outward_orientation),true);

  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  FaceIndexMap fid_map = CGAL::get_initialized_face_index_map(tm, np);

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

  // set the connected component id of each face
  std::size_t nb_cc = connected_components(tm,
                                           make_compose_property_map(fid_map,make_property_map(face_cc)),
                                           parameters::face_index_map(fid_map));

  // extract a vertex with max z coordinate for each connected component
  std::vector<vertex_descriptor> xtrm_vertices(nb_cc, Graph_traits::null_vertex());
  for(vertex_descriptor vd : vertices(tm))
  {
    halfedge_descriptor test_hd = halfedge(vd, tm);
    if(test_hd == Graph_traits::null_halfedge())
      continue;
    face_descriptor test_face = face(halfedge(vd, tm), tm);
    if(test_face == Graph_traits::null_face())
      test_face = face(opposite(halfedge(vd, tm), tm), tm);
    CGAL_assertion(test_face != Graph_traits::null_face());
    std::size_t cc_id = face_cc[get(fid_map,test_face )];
    if (xtrm_vertices[cc_id]==Graph_traits::null_vertex())
      xtrm_vertices[cc_id]=vd;
    else
      if (get(vpm, vd).z()>get(vpm,xtrm_vertices[cc_id]).z())
        xtrm_vertices[cc_id]=vd;
  }
  std::vector<std::vector<face_descriptor> > ccs(nb_cc);
  for(face_descriptor fd : faces(tm))
  {
    ccs[face_cc[get(fid_map,fd)]].push_back(fd);
  }

  //orient ccs outward
  for(std::size_t id=0; id<nb_cc; ++id)
  {
    // skip it if the vertex is on the boundary
    bool v_is_border = false;
    for(halfedge_descriptor h : halfedges_around_target(xtrm_vertices[id], tm))
      if (is_border(h, tm))
      {
        v_is_border = true;
        break;
      }

    if(!v_is_border && (internal::is_outward_oriented(xtrm_vertices[id], tm, np)
                        != orient_outward))
    {
      reverse_face_orientations(ccs[id], tm);
    }
  }
}

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_PMP_ORIENTATION_BASIC_H
