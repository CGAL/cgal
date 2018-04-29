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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Ilker O. Yaz


#ifndef CGAL_ORIENT_POLYGON_MESH_H
#define CGAL_ORIENT_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/orientation.h>


#include <algorithm>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>
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

    CGAL_assertion(halfedge(v_max, pmesh)!=boost::graph_traits<PolygonMesh>::null_halfedge());

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
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param pmesh the closed polygon mesh to be tested
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
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
  internal::Compare_vertex_points_z_3<GT, VPMap> less_z(vpmap, gt);
  typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_max = *(vertices(pmesh).first);
  for (typename boost::graph_traits<PolygonMesh>::vertex_iterator
          vit=cpp11::next(vertices(pmesh).first), vit_end = vertices(pmesh).second;
          vit!=vit_end; ++vit)
  {
    // skip isolated vertices
    if (halfedge(*vit, pmesh)==boost::graph_traits<PolygonMesh>::null_halfedge())
      continue;
    if( less_z(v_max, *vit) )
      v_max=*vit;
  }

  // only isolated vertices
  if (halfedge(v_max, pmesh)==boost::graph_traits<PolygonMesh>::null_halfedge())
    return true;

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
namespace internal {

template <class Kernel, class TriangleMesh, class VD, class Fid_map, class Vpm>
void recursive_orient_volume_ccs( TriangleMesh& tm,
                                  Vpm& vpm,
                                  Fid_map& fid_map,
                                  const std::vector<VD>& xtrm_vertices,
                                  boost::dynamic_bitset<>& cc_handled,
                                  const std::vector<std::size_t>& face_cc,
                                  std::size_t xtrm_cc_id,
                                  bool is_parent_outward_oriented)
{
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef Side_of_triangle_mesh<TriangleMesh, Kernel, Vpm> Side_of_tm;
  std::vector<face_descriptor> cc_faces;
  BOOST_FOREACH(face_descriptor fd, faces(tm))
  {
    if(face_cc[get(fid_map, fd)]==xtrm_cc_id)
      cc_faces.push_back(fd);
  }
// first check that the orientation of the current cc is consistant with its
// parent cc containing it
  bool new_is_parent_outward_oriented = internal::is_outward_oriented(
         xtrm_vertices[xtrm_cc_id], tm, parameters::vertex_point_map(vpm));
  if (new_is_parent_outward_oriented==is_parent_outward_oriented)
  {
    Polygon_mesh_processing::reverse_face_orientations(cc_faces, tm);
    new_is_parent_outward_oriented = !new_is_parent_outward_oriented;
  }
  cc_handled.set(xtrm_cc_id);

  std::size_t nb_cc = cc_handled.size();

// get all cc that are inside xtrm_cc_id

  typename Side_of_tm::AABB_tree aabb_tree(cc_faces.begin(), cc_faces.end(),
                                           tm, vpm);
  Side_of_tm side_of_cc(aabb_tree);

  std::vector<std::size_t> cc_inside;
  for(std::size_t id=0; id<nb_cc; ++id)
  {
    if (cc_handled.test(id)) continue;
    if (side_of_cc(get(vpm,xtrm_vertices[id]))==ON_BOUNDED_SIDE)
      cc_inside.push_back(id);
  }

// check whether we need another recursion for cc inside xtrm_cc_id
  if (!cc_inside.empty())
  {
    std::size_t new_xtrm_cc_id = cc_inside.front();
    boost::dynamic_bitset<> new_cc_handled(nb_cc,0);
    new_cc_handled.set();
    new_cc_handled.reset(new_xtrm_cc_id);
    cc_handled.set(new_xtrm_cc_id);

    std::size_t nb_candidates = cc_inside.size();
    for (std::size_t i=1;i<nb_candidates;++i)
    {
      std::size_t candidate = cc_inside[i];
      if(get(vpm,xtrm_vertices[candidate]).z() >
         get(vpm,xtrm_vertices[new_xtrm_cc_id]).z()) new_xtrm_cc_id=candidate;
      new_cc_handled.reset(candidate);
      cc_handled.set(candidate);
    }

    internal::recursive_orient_volume_ccs<Kernel>(
           tm, vpm, fid_map, xtrm_vertices, new_cc_handled, face_cc,
           new_xtrm_cc_id, new_is_parent_outward_oriented);
  }

// now explore remaining cc included in the same cc as xtrm_cc_id
  boost::dynamic_bitset<> cc_not_handled = ~cc_handled;
  std::size_t new_xtrm_cc_id = cc_not_handled.find_first();
  if (new_xtrm_cc_id == cc_not_handled.npos) return ;

  for (std::size_t candidate = cc_not_handled.find_next(new_xtrm_cc_id);
                   candidate < cc_not_handled.npos;
                   candidate = cc_not_handled.find_next(candidate))
  {
     if(get(vpm,xtrm_vertices[candidate]).z() > get(vpm,xtrm_vertices[new_xtrm_cc_id]).z())
        new_xtrm_cc_id = candidate;
  }

  internal::recursive_orient_volume_ccs<Kernel>(
            tm, vpm, fid_map, xtrm_vertices, cc_handled, face_cc,
            new_xtrm_cc_id, is_parent_outward_oriented);
}

}//end internal

/**
* \ingroup PMP_orientation_grp
* makes each connected component of a closed triangulated surface mesh
* inward or outward oriented.
*
* @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph` .
*                      If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
*                      as a named parameter, then it must be initialized.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters
*
* @param tm a closed triangulated surface mesh
* @param np optional sequence of \ref pmp_namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamBegin{vertex_point_map}
*     the property map with the points associated to the vertices of `tm`.
*     If this parameter is omitted, an internal property map for
*     `CGAL::vertex_point_t` should be available in `TriangleMesh`
*   \cgalParamEnd
*   \cgalParamBegin{face_index_map}
*     a property map containing the index of each face of `tm`.
*   \cgalParamEnd
*   \cgalParamBegin{outward_orientation}
*     if set to `true` (default) indicates that each connected component will be outward oriented,
*     (inward oriented if `false`).
*   \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<class TriangleMesh, class NamedParameters>
void orient(TriangleMesh& tm, const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh,
      NamedParameters>::const_type Vpm;
  typedef typename GetFaceIndexMap<TriangleMesh,
      NamedParameters>::const_type Fid_map;

  CGAL_assertion(is_triangle_mesh(tm));
  CGAL_assertion(is_valid(tm));
  CGAL_assertion(is_closed(tm));

  using boost::choose_param;
  using boost::get_param;

  bool orient_outward = choose_param(
                          get_param(np, internal_np::outward_orientation),true);

  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(boost::vertex_point, tm));

  Fid_map fid_map = choose_param(get_param(np, internal_np::face_index),
                                 get_const_property_map(boost::face_index, tm));

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

  // set the connected component id of each face
  std::size_t nb_cc = connected_components(tm,
                                           bind_property_maps(fid_map,make_property_map(face_cc)),
                                           parameters::face_index_map(fid_map));

  // extract a vertex with max z coordinate for each connected component
  std::vector<vertex_descriptor> xtrm_vertices(nb_cc, Graph_traits::null_vertex());
  BOOST_FOREACH(vertex_descriptor vd, vertices(tm))
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
  BOOST_FOREACH(face_descriptor fd, faces(tm))
  {
    ccs[face_cc[get(fid_map,fd)]].push_back(fd);
  }

  //orient ccs outward
  for(std::size_t id=0; id<nb_cc; ++id)
  {
    if(internal::is_outward_oriented(xtrm_vertices[id], tm, np)
        != orient_outward)
    {
      reverse_face_orientations(ccs[id], tm);
    }
  }
}

template<class TriangleMesh>
void orient(TriangleMesh& tm)
{
  orient(tm, parameters::all_default());
}


/** \ingroup PMP_orientation_grp
 *
 * orients the connected components of `tm` to make it bound a volume.
 * See \ref coref_def_subsec for a precise definition.
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
 *                      If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
 *                      as a named parameter, then it must be initialized.
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters
 *
 * @param tm a closed triangulated surface mesh
 * @param np optional sequence of \ref pmp_namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *     the property map with the points associated to the vertices of `tm`.
 *     If this parameter is omitted, an internal property map for
 *     `CGAL::vertex_point_t` should be available in `TriangleMesh`
 *   \cgalParamEnd
 *   \cgalParamBegin{face_index_map}
 *     a property map containing the index of each face of `tm`.
 *   \cgalParamEnd
 *   \cgalParamBegin{outward_orientation}
 *     if set to `true` (default) the outer connected components will be outward oriented (inward oriented if set to `false`).
 *     If the outer connected components are inward oriented, it means that the infinity will be considered
 *     as part of the volume bounded by `tm`.
 *   \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \see `CGAL::Polygon_mesh_processing::does_bound_a_volume()`
 */
template <class TriangleMesh, class NamedParameters>
void orient_to_bound_a_volume(TriangleMesh& tm,
                                        const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh,
      NamedParameters>::const_type Vpm;
  typedef typename GetFaceIndexMap<TriangleMesh,
      NamedParameters>::const_type Fid_map;
  typedef typename Kernel_traits<
      typename boost::property_traits<Vpm>::value_type >::Kernel Kernel;
  if (!is_closed(tm)) return;
  if (!is_triangle_mesh(tm)) return;

  using boost::choose_param;
  using boost::get_param;

  bool orient_outward = choose_param(
                          get_param(np, internal_np::outward_orientation),true);

  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(boost::vertex_point, tm));

  Fid_map fid_map = choose_param(get_param(np, internal_np::face_index),
                                 get_const_property_map(boost::face_index, tm));

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));


  // set the connected component id of each face
  std::size_t nb_cc = connected_components(tm,
                                           bind_property_maps(fid_map,make_property_map(face_cc)),
                                           parameters::face_index_map(fid_map));

  if (nb_cc == 1)
  {
    if( orient_outward != is_outward_oriented(tm))
      reverse_face_orientations(faces(tm), tm);
    return ;
  }


  boost::dynamic_bitset<> cc_handled(nb_cc, 0);

  // extract a vertex with max z coordinate for each connected component
  std::vector<vertex_descriptor> xtrm_vertices(nb_cc, Graph_traits::null_vertex());
  BOOST_FOREACH(vertex_descriptor vd, vertices(tm))
  {
    std::size_t cc_id = face_cc[get(fid_map, face(halfedge(vd, tm), tm))];
    if (xtrm_vertices[cc_id]==Graph_traits::null_vertex())
      xtrm_vertices[cc_id]=vd;
    else
      if (get(vpm, vd).z()>get(vpm,xtrm_vertices[cc_id]).z())
        xtrm_vertices[cc_id]=vd;
  }

  //extract a vertex with max z amongst all components
  std::size_t xtrm_cc_id=0;
  for(std::size_t id=1; id<nb_cc; ++id)
    if (get(vpm, xtrm_vertices[id]).z()>get(vpm,xtrm_vertices[xtrm_cc_id]).z())
      xtrm_cc_id=id;

  bool is_parent_outward_oriented =
      ! orient_outward;

  internal::recursive_orient_volume_ccs<Kernel>(tm, vpm, fid_map,
                                                xtrm_vertices,
                                                cc_handled,
                                                face_cc,
                                                xtrm_cc_id,
                                                is_parent_outward_oriented);
}

template <class TriangleMesh>
void orient_to_bound_a_volume(TriangleMesh& tm)
{
  orient_to_bound_a_volume(tm, parameters::all_default());
}
} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H
