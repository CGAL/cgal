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
// Author(s)     : Sebastien Loriot and Ilker O. Yaz


#ifndef CGAL_ORIENT_POLYGON_MESH_H
#define CGAL_ORIENT_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/orientation.h>


#include <algorithm>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/utility.h>

#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/ref.hpp>

#include <functional>
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
    using parameters::choose_parameter;
    using parameters::get_parameter;

    CGAL_assertion(halfedge(v_max, pmesh)!=boost::graph_traits<PolygonMesh>::null_halfedge());

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
 * tests whether a closed triangle mesh has a positive orientation.
 * A closed triangle mesh is considered to have a positive orientation if the normal vectors
 * to all its faces point outside the domain bounded by the triangle mesh.
 * The normal vector to each face is chosen pointing on the side of the face
 * where its sequence of vertices is seen counterclockwise.
 * @pre `CGAL::is_closed(tm)`
 * @pre `CGAL::is_triangle_mesh(tm)`
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
template<typename TriangleMesh, typename NamedParameters>
bool is_outward_oriented(const TriangleMesh& tm,
                         const NamedParameters& np)
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

///\cond SKIP_IN_MANUAL

template<typename TriangleMesh>
bool is_outward_oriented(const TriangleMesh& tm)
{
  return is_outward_oriented(tm,
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
  boost::unordered_set<halfedge_descriptor> already_seen;
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
* makes each connected component of a closed triangulated surface mesh
* inward or outward oriented.
*
* @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph` .
* @tparam NamedParameters a sequence of \ref bgl_namedparameters
*
* @param tm a closed triangulated surface mesh
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \pre `CGAL::is_closed(tm)`
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
template<class TriangleMesh, class NamedParameters>
void orient(TriangleMesh& tm,
            const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh>                                        Graph_traits;
  typedef typename Graph_traits::vertex_descriptor                                 vertex_descriptor;
  typedef typename Graph_traits::face_descriptor                                   face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor                               halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type    Vpm;
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::type FaceIndexMap;

  CGAL_assertion(is_triangle_mesh(tm));
  CGAL_assertion(is_valid_polygon_mesh(tm));
  CGAL_assertion(is_closed(tm));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool orient_outward = choose_parameter(get_parameter(np, internal_np::outward_orientation),true);

  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  FaceIndexMap fid_map = CGAL::get_initialized_face_index_map(tm, np);

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

  // set the connected component id of each face
  std::size_t nb_cc = connected_components(tm,
                                           bind_property_maps(fid_map,make_property_map(face_cc)),
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

/*!
 * \ingroup PMP_orientation_grp
 * Enumeration type used to indicate the status of a set of faces
 * classified by the function `volume_connected_components()`.
 * The set of faces defines either a volume connected connected component
 * in the case of `VALID_VOLUME` or a surface connected component otherwise.
 */
enum Volume_error_code { VALID_VOLUME, ///< The set of faces bounds a volume
                         SURFACE_WITH_SELF_INTERSECTIONS, ///< The set of faces is self-intersecting
                         VOLUME_INTERSECTION, ///< The set of faces intersects another surface connected component
                         INCOMPATIBLE_ORIENTATION ///< The set of faces is included in a volume but has an incompatible orientation
                       };
namespace internal {

// helper function to copy data
template<class T, class RefToContainer>
void copy_container_content(
  const std::vector<T>& vec,
  RefToContainer ref_wrapper)
{
  ref_wrapper.get().reserve(vec.size());
  for(const T& t : vec)
  {
   ref_wrapper.get().push_back(t);
  }
}

template<class T>
void copy_container_content(
  std::vector<T>& vec,
  std::reference_wrapper<std::vector<T> > ref_wrapper)
{
  vec.swap(ref_wrapper.get());
}

template<class T>
void copy_container_content(
  std::vector<T>& vec,
  boost::reference_wrapper<std::vector<T> > ref_wrapper)
{
  vec.swap(ref_wrapper.get());
}

template<class T>
inline
void copy_container_content(
  std::vector<T>&,
  internal_np::Param_not_found)
{}

template <class RefToContainer>
void copy_nested_parents(
  const std::vector< std::vector<std::size_t> >& nested_parents,
  RefToContainer ref_to_vector)
{
  typedef typename RefToContainer::type Container;
  typedef typename Container::value_type Container_value;

  ref_to_vector.get().reserve(nested_parents.size());
  for(const auto& t : nested_parents)
  {
    Container_value c;
    c.reserve(t.size());
    for(const std::size_t& val : t)
      c.push_back(val);
    ref_to_vector.get().push_back(c);
  }
}

inline
void copy_nested_parents(
  std::vector< std::vector<std::size_t> >& nested_parents,
  std::reference_wrapper<std::vector< std::vector<std::size_t> > > ref_to_vector)
{
  nested_parents.swap(ref_to_vector.get());
}

inline
void copy_nested_parents(
  std::vector< std::vector<std::size_t> >& nested_parents,
  boost::reference_wrapper<std::vector< std::vector<std::size_t> > > ref_to_vector)
{
  nested_parents.swap(ref_to_vector.get());
}

inline
void copy_nested_parents(
  std::vector< std::vector<std::size_t> >&,
  internal_np::Param_not_found)
{}


// helper function for setting id maps
template <class TriangleMesh, class FaceIndexMap, class FaceCCIdMap>
void set_f_cc_id(
  const std::vector<std::size_t>& f_cc,
  FaceIndexMap face_index_map,
  FaceCCIdMap face_cc_map,
  const TriangleMesh& tm)
{
  for(typename boost::graph_traits<TriangleMesh>::face_descriptor fd : faces(tm))
    put(face_cc_map, fd, f_cc[ get(face_index_map, fd) ]);
}

template <class TriangleMesh, class FaceIndexMap>
void set_f_cc_id(
  const std::vector<std::size_t>&,
  FaceIndexMap,
  internal_np::Param_not_found,
  const TriangleMesh&)
{}

template <class RefToVector>
void copy_cc_to_volume_id(
  std::vector<std::size_t>& cc_volume_ids,
  RefToVector ref_to_vector)
{
  ref_to_vector.get().swap( cc_volume_ids );
}

inline
void copy_cc_to_volume_id(
  std::vector<std::size_t>&,
  internal_np::Param_not_found)
{}

template <class RefToVector>
void copy_nesting_levels(
  std::vector<std::size_t>& nesting_levels,
  RefToVector ref_to_vector)
{
  ref_to_vector.get().swap( nesting_levels );
}

inline
void copy_nesting_levels(
  std::vector<std::size_t>&,
  internal_np::Param_not_found)
{}

template <class RefToBitset>
void copy_orientation_bitset(
  const std::vector<bool>& is_cc_outward_oriented,
  RefToBitset ref_to_bs)
{
  ref_to_bs.get() = is_cc_outward_oriented;
}

inline
void copy_orientation_bitset(
  const std::vector<bool>&,
  internal_np::Param_not_found)
{}

template <class OutputIterator>
void set_cc_intersecting_pairs(
  const std::set< std::pair<std::size_t, std::size_t> >& self_intersecting_cc,
  OutputIterator out)
{
  for (const std::pair<std::size_t, std::size_t>& p : self_intersecting_cc)
    *out++=p;
}

inline
void set_cc_intersecting_pairs(
  const std::set< std::pair<std::size_t, std::size_t> >&,
  internal_np::Param_not_found)
{}

} // internal


/*!
 * \ingroup PMP_orientation_grp
 * assigns to each face of `tm` an id corresponding to the volume connected component
 * it contributes to.
 *
 * Using the adjacency relation of two faces along an edge, a triangle mesh can be split
 * into connected components (*surface components* in the following).
 * A surface component without boundary separates the 3D space into an infinite and
 * a finite volume. We say that the finite volume is <i>enclosed</i> by this surface
 * component.
 *
 * The volume connected components (*volume components* in the following) are defined as follows:
 * Each surface component `S` that is outside any volume enclosed by
 * another surface component defines the *outer boundary* of a volume component.
 * Each surface component that is inside the volume enclosed by `S`
 * defines a *hole* if it is included in no other volume enclosed by a surface component
 * but `S`. Ignoring the identified volume component, the same procedure is recursively
 * repeated for all surface components in each hole.
 *
 * There are some special cases:
 * - a non-closed surface component is reported as a volume component ignoring any inclusion test
 * - a self-intersecting surface component is reported as a volume component ignoring any inclusion test
 * - a surface component intersecting another surface component
 *   is reported as a volume component, and so are all the surface components inside its
 *   enclosed volume
 * - if `do_orientation_tests` is set to `true`, if the holes are not all equally oriented
 *   (all inward or all outward) or if the holes and the outer boundary are equally
 *   oriented, each surface component is reported as a volume component,
 *   and so are all the surface components inside the corresponding enclosed volumes
 * - If `do_orientation_tests` is set to `true` and the surface components that are
 *   outside all enclosed volumes are inward oriented, they are then considered as holes
 *   of the unbounded volume (that has no outer boundary)
 *
 * A property map for `CGAL::vertex_point_t` must be either available as an internal property map
 * of `tm` or provided as one of the \ref bgl_namedparameters "Named Parameters".
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam VolumeFaceIndexMap a model of `WritablePropertyMap` with
 *                      `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type and
 *                      `boost::graph_traits<TriangleMesh>::%faces_size_type` as value type.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm the input triangle mesh
 * @param volume_id_map the property map filled by this function with indices of volume components associated to the faces of `tm`
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * @pre `CGAL::is_closed(tm)`
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
 *   \cgalParamNBegin{face_connected_component_map}
 *     \cgalParamDescription{a property map filled by this function and that will contain for each face the id
 *                           of its surface component in the range `[0, number of surface components - 1]`}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{volume_inclusions}
 *     \cgalParamDescription{a container, which contains at position `k` the ids of all the
 *                           surface components that are the first intersected by any ray with source on
 *                           the surface component `k` and directed outside the volume enclosed by the
 *                           surface component `k`. There is only one such id but when some surface components intersect.
 *                           The size of the container is exactly the number of surface components of `tm`.}
 *     \cgalParamType{a `reference_wrapper` (either from `boost` or the standard library) containing
 *                    a reference to an object that must be a model of the `BackInsertionSequence` concept,
 *                    with a value type being a model of `BackInsertionSequence` of `std::size_t`,
 *                    both types having the functions `reserve()` and `push_back()`}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_orientation_tests}
 *     \cgalParamDescription{If `true`, the orientation of the faces of each surface components
 *                           will be taken into account for the definition of the volume.
 *                           If `false`, the face orientation is ignored and the volumes are defined
 *                           only by the nesting of surface components.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{error_codes}
 *     \cgalParamDescription{a container which indicates the status of a volume assigned to a set of faces.
 *                           The description of the value type is given in the documentation of the enumeration type.
 *                           The size of the container is exactly the number of volume components.}
 *     \cgalParamType{a `reference_wrapper` (either from `boost` or the standard library) containing
 *                    a reference to a container that must be a model of the `BackInsertionSequence` concept,
 *                    with value type \link PMP_orientation_grp `Volume_error_code` \endlink}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{do_self_intersection_tests}
 *     \cgalParamDescription{If `false`, it is assumed that `tm` does not contains any self-intersections.
 *                           If `true`, the input might contain some self-intersections and a test is done prior to the volume decomposition.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{connected_component_id_to_volume_id}
 *     \cgalParamDescription{For each connected component identified using its id `ccid`, the id
 *                           of the volume it contributes to describe is the value at the position
 *                           `ccid` in the container. The size of the container is exactly the number
 *                           of connected components.}
 *     \cgalParamType{a `reference_wrapper` (either from `boost` or the standard library) containing
 *                    a reference to a container that must be a model of the `BackInsertionSequence` concept,
 *                    with value type `std::size_t`}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *
 *
 *   \cgalParamNBegin{nesting_levels}
 *     \cgalParamDescription{For each connected component identified using its id `ccid`, the container contains the number of
 *                           connected components containing on its bounded side this component.
 *                           The size of the container is exactly the number of connected components.}
 *     \cgalParamType{a `reference_wrapper` (either from `boost` or the standard library) containing
 *                    a reference to a container that must be a model of the `BackInsertionSequence` concept,
 *                    with value type `std::size_t`}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{is_cc_outward_oriented}
 *     \cgalParamDescription{For each connected component identified using its id `ccid`, the output of `is_outward_oriented()`
 *                           called on the triangle mesh corresponding to this connected component
 *                           is the value at the position `ccid` in the container.
 *                           The size of the container is exactly the number of connected components.}
 *     \cgalParamType{a `reference_wrapper` (either from `boost` or the standard library) containing
 *                    a reference to a container that must be a model of the `BackInsertionSequence` concept,
 *                    with value type `bool`}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{intersecting_volume_pairs_output_iterator}
 *     \cgalParamDescription{Output iterator into which pairs of ids (id must be convertible to `std::size_t`) can be put.
 *                           Each pair of connected components intersecting will be reported using their ids.
 *                           If `do_self_intersection_tests` named parameter is set to `false`, nothing will be reported.}
 *     \cgalParamType{a model of `OutputIterator`}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return the number of volume components defined by `tm`
 */
template <class TriangleMesh, class VolumeFaceIndexMap, class NamedParameters>
std::size_t
volume_connected_components(const TriangleMesh& tm,
                                  VolumeFaceIndexMap volume_id_map,
                            const NamedParameters& np)
{
  CGAL_precondition(is_triangle_mesh(tm));
  CGAL_precondition(is_closed(tm));

  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::const_type Vpm;
  typedef typename GetInitializedFaceIndexMap<TriangleMesh,
                                              NamedParameters>::type FaceIndexMap;

  typedef typename Kernel_traits<
    typename boost::property_traits<Vpm>::value_type >::Kernel Kernel;

  Vpm vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                get_const_property_map(boost::vertex_point, tm));

  FaceIndexMap fid_map = CGAL::get_initialized_face_index_map(tm, np);

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

// set the connected component id of each face
  const std::size_t nb_cc = connected_components(tm,
                              bind_property_maps(fid_map,make_property_map(face_cc)),
                              parameters::face_index_map(fid_map));

  // contains for each CC the CC that are in its bounded side
  std::vector<std::vector<std::size_t> > nested_cc_per_cc(nb_cc);

  // copy cc-id info
  internal::set_f_cc_id(face_cc, fid_map, parameters::get_parameter(np, internal_np::face_connected_component_map), tm);

  const bool do_self_intersection_tests =
        parameters::choose_parameter(parameters::get_parameter(np, internal_np::do_self_intersection_tests),
                            false);
  const bool ignore_orientation_of_cc =
    !parameters::choose_parameter(parameters::get_parameter(np, internal_np::do_orientation_tests),
                         true);

  const bool used_as_a_predicate =
    parameters::choose_parameter(parameters::get_parameter(np, internal_np::i_used_as_a_predicate),
                        false); // indicate if the function is called by does_bound_a_volume

  CGAL_assertion(!used_as_a_predicate || !ignore_orientation_of_cc);

  std::vector<Volume_error_code> error_codes;
  std::vector<bool> is_cc_outward_oriented;
  std::vector<std::size_t> cc_volume_ids(nb_cc, -1);
  std::vector < std::size_t > nesting_levels(nb_cc, 0); // indicates for each CC its nesting level

  boost::dynamic_bitset<> cc_handled(nb_cc, 0);
  std::size_t next_volume_id = 0;

// Handle self-intersecting connected components
  typedef std::pair<face_descriptor, face_descriptor> Face_pair;
  std::vector<Face_pair> si_faces;
  std::set< std::pair<std::size_t, std::size_t> > self_intersecting_cc; // due to self-intersections
  if (do_self_intersection_tests)
    self_intersections(tm, std::back_inserter(si_faces));
  std::vector<bool> is_involved_in_self_intersection(nb_cc, false);

  if (!si_faces.empty() && used_as_a_predicate)
    return 0;

  for(const Face_pair& fp : si_faces)
  {
    std::size_t first_cc_id = face_cc[ get(fid_map, fp.first) ];
    std::size_t second_cc_id = face_cc[ get(fid_map, fp.second) ];

    if (first_cc_id==second_cc_id)
    {
      if ( !cc_handled.test(first_cc_id) )
      {
        cc_handled.set(first_cc_id);
        cc_volume_ids[first_cc_id]=next_volume_id++;
        error_codes.push_back(SURFACE_WITH_SELF_INTERSECTIONS);
      }
    }
    else
    {
      is_involved_in_self_intersection[first_cc_id] = true;
      is_involved_in_self_intersection[second_cc_id] = true;
      self_intersecting_cc.insert( make_sorted_pair(first_cc_id, second_cc_id) );
    }
  }

  std::vector< std::vector<std::size_t> > nesting_parents(nb_cc);
  if (!cc_handled.all())
  {
  // extract a vertex with max z coordinate for each connected component
    std::vector<vertex_descriptor> xtrm_vertices(nb_cc, GT::null_vertex());
    for(vertex_descriptor vd : vertices(tm))
    {
      halfedge_descriptor h = halfedge(vd, tm);
      if (is_border(h, tm)) h = opposite(h, tm);
      std::size_t cc_id = face_cc[get(fid_map, face(h, tm))];
      if (cc_handled.test(cc_id)) continue;
      if (xtrm_vertices[cc_id]==GT::null_vertex())
        xtrm_vertices[cc_id]=vd;
      else
        if (get(vpm, vd).z()>get(vpm,xtrm_vertices[cc_id]).z())
          xtrm_vertices[cc_id]=vd;
    }

  // fill orientation vector for each surface CC
    if (!ignore_orientation_of_cc)
    {
      is_cc_outward_oriented.resize(nb_cc);
      for(std::size_t cc_id=0; cc_id<nb_cc; ++cc_id)
      {
        if (cc_handled.test(cc_id)) continue;
        is_cc_outward_oriented[cc_id] = internal::is_outward_oriented(xtrm_vertices[cc_id], tm, np);
      }
    }

  //collect faces per CC
    std::vector< std::vector<face_descriptor> > faces_per_cc(nb_cc);
    std::vector< std::size_t > nb_faces_per_cc(nb_cc, 0);
    for(face_descriptor fd : faces(tm))
    {
      std::size_t cc_id = face_cc[ get(fid_map, fd) ];
      ++nb_faces_per_cc[ cc_id ];
    }
    for (std::size_t i=0; i<nb_cc; ++i)
      if (!cc_handled.test(i))
        faces_per_cc[i].reserve( nb_faces_per_cc[i] );
    for(face_descriptor fd : faces(tm))
    {
      std::size_t cc_id = face_cc[ get(fid_map, fd) ];
      faces_per_cc[ cc_id ].push_back(fd);
    }

  // init the main loop
    // similar as above but exclusively contains cc ids included by more that one CC.
    // The result will be then merged with nested_cc_per_cc but temporarilly we need
    // another container to not more than once the inclusion testing (in case a CC is
    // included by more than 2 CC) + associate such CC to only one volume
    std::vector<std::vector<std::size_t> > nested_cc_per_cc_shared(nb_cc);
    std::vector < boost::dynamic_bitset<> > level_k_nestings; // container containing CCs in the same volume (one bitset per volume) at level k
    level_k_nestings.push_back( ~cc_handled );

  // the following loop is exploring the nesting level by level (0 -> max_level)
    std::size_t k = 0;
    while (!level_k_nestings.empty())
    {
      std::vector < boost::dynamic_bitset<> > level_k_plus_1_nestings;
      for(boost::dynamic_bitset<> cc_to_handle : level_k_nestings)
      {
        CGAL_assertion( cc_to_handle.any() );
        while(cc_to_handle.any())
        {
        //extract a vertex with max z amongst all components
          std::size_t xtrm_cc_id=cc_to_handle.find_first();
          for(std::size_t id  = cc_to_handle.find_next(xtrm_cc_id);
                          id != cc_to_handle.npos;
                          id  = cc_to_handle.find_next(id))
          {
            if (get(vpm, xtrm_vertices[id]).z()>get(vpm,xtrm_vertices[xtrm_cc_id]).z())
              xtrm_cc_id=id;
          }
          cc_to_handle.reset(xtrm_cc_id);
          nesting_levels[xtrm_cc_id] = k;

        // collect id inside xtrm_cc_id CC
          typedef Side_of_triangle_mesh<TriangleMesh, Kernel, Vpm> Side_of_tm;
          typename Side_of_tm::AABB_tree aabb_tree(faces_per_cc[xtrm_cc_id].begin(),
                                                   faces_per_cc[xtrm_cc_id].end(),
                                                   tm, vpm);
          Side_of_tm side_of_cc(aabb_tree);

          std::vector<std::size_t> cc_intersecting; // contains id of CC intersecting xtrm_cc_id

          boost::dynamic_bitset<> nested_cc_to_handle(nb_cc, 0);
          for(std::size_t id  = cc_to_handle.find_first();
                          id != cc_to_handle.npos;
                          id  = cc_to_handle.find_next(id))
          {
            if (self_intersecting_cc.count( make_sorted_pair(xtrm_cc_id, id) )!= 0)
            {
              cc_intersecting.push_back(id);
              continue; // to not dot inclusion test for intersecting CCs
            }

            if (side_of_cc(get(vpm,xtrm_vertices[id]))==ON_BOUNDED_SIDE)
            {
              nested_cc_per_cc[xtrm_cc_id].push_back(id);
              // mark nested CC as handle and collect them for the handling of the next level
              nested_cc_to_handle.set(id);
              cc_to_handle.reset(id);
            }
          }

        //for each CC intersecting xtrm_cc_id, find the CCs included in both
          for(std::size_t id : cc_intersecting)
          {
            typename Side_of_tm::AABB_tree aabb_tree(faces_per_cc[id].begin(),
                                                     faces_per_cc[id].end(),
                                                     tm, vpm);
            Side_of_tm side_of_cc(aabb_tree);
            for(std::size_t ncc_id : nested_cc_per_cc[xtrm_cc_id])
            {
              if (self_intersecting_cc.count( make_sorted_pair(ncc_id, id) )!= 0)
                continue;
              if (side_of_cc(get(vpm,xtrm_vertices[ncc_id]))==ON_BOUNDED_SIDE)
                nested_cc_per_cc_shared[id].push_back(ncc_id);
            }
          }

          if ( nested_cc_per_cc[xtrm_cc_id].empty() ) continue;
          level_k_plus_1_nestings.push_back(nested_cc_to_handle);
        }
      }
      ++k;
      level_k_nestings.swap(level_k_plus_1_nestings);
    }

    // early return for orient_to_bound_a_volume
    if (parameters::choose_parameter(parameters::get_parameter(np, internal_np::i_used_for_volume_orientation),false))
    {
      internal::copy_container_content(nesting_levels, parameters::get_parameter(np, internal_np::nesting_levels));
      internal::copy_container_content(is_cc_outward_oriented, parameters::get_parameter(np, internal_np::is_cc_outward_oriented));
      return 0;
    }

  // detect inconsistencies of the orientation at the level 0
  // and check if all CC at level 0 are in the same volume
    std::size_t ref_cc_id = nb_cc;
    std::size_t FIRST_LEVEL = 0; // used to know if even or odd nesting is the top level
    if(!ignore_orientation_of_cc)
    {
      for(std::size_t cc_id=0; cc_id<nb_cc; ++cc_id)
      {
        if (cc_handled.test(cc_id)) continue;
        if( nesting_levels[cc_id]==0 )
        {
          if(ref_cc_id==nb_cc)
            ref_cc_id=cc_id;
          else
            if( is_cc_outward_oriented[cc_id] != is_cc_outward_oriented[ref_cc_id] )
            {
              if (used_as_a_predicate) return 0;
              // all is indefinite
              for(std::size_t id=0; id<nb_cc; ++id)
              {
                if (cc_handled.test(cc_id)) continue;
                cc_volume_ids[id] = next_volume_id++;
                error_codes.push_back(INCOMPATIBLE_ORIENTATION);
              }
              cc_handled.set();
              break;
            }
        }
      }

      if (!cc_handled.all() && !is_cc_outward_oriented[ref_cc_id])
      {
        // all level 0 CC are in the same volume
        for(std::size_t cc_id=0; cc_id<nb_cc; ++cc_id)
        {
          if (cc_handled.test(cc_id)) continue;
          if( nesting_levels[cc_id]==0 )
          {
            cc_handled.set(cc_id);
            cc_volume_ids[cc_id]=next_volume_id;
          }
        }
        ++next_volume_id;
        error_codes.push_back(VALID_VOLUME);
        FIRST_LEVEL = 1;
      }
    }

  // apply volume classification using level 0 nesting
    for(std::size_t cc_id=0; (!cc_handled.all()) && cc_id<nb_cc; ++cc_id)
    {
      if (cc_handled.test(cc_id)) continue;
      CGAL_assertion( nesting_levels[cc_id]!=0 || ignore_orientation_of_cc || is_cc_outward_oriented[cc_id] );
      if( nesting_levels[cc_id]%2 != FIRST_LEVEL ) continue; // we look for outer boundaries of volume only
      cc_handled.set(cc_id);
      cc_volume_ids[cc_id] = next_volume_id++;

      //if the CC is involved in a self-intersection all nested CC are put in a separate volumes
      if (is_involved_in_self_intersection[cc_id])
      {
        error_codes.push_back(VOLUME_INTERSECTION);
        for(std::size_t ncc_id : nested_cc_per_cc[cc_id])
        {
          cc_handled.set(ncc_id);
          cc_volume_ids[ncc_id] = next_volume_id++;
          error_codes.push_back(VOLUME_INTERSECTION);
        }
        continue;
      }
      else
      {
        if (!ignore_orientation_of_cc && !is_cc_outward_oriented[cc_id])
        {
          // invalid orientation, all children are marked as incorrectly oriented
          if (used_as_a_predicate) return 0;
          cc_handled.set(cc_id);
          cc_volume_ids[cc_id] = next_volume_id++;
          error_codes.push_back(INCOMPATIBLE_ORIENTATION);
          for(std::size_t ncc_id : nested_cc_per_cc[cc_id])
          {
            cc_handled.set(ncc_id);
            cc_volume_ids[ncc_id] = next_volume_id++;
            error_codes.push_back(INCOMPATIBLE_ORIENTATION);
          }
          continue;
        }
        else
          error_codes.push_back(VALID_VOLUME);
      }

      for(std::size_t ncc_id : nested_cc_per_cc[cc_id])
      {
        if (nesting_levels[ncc_id]==nesting_levels[cc_id]+1)
        {
          cc_handled.set(ncc_id);
          if (!ignore_orientation_of_cc)
          {
            if (is_cc_outward_oriented[cc_id]==is_cc_outward_oriented[ncc_id])
            {
              // the surface component has an incorrect orientation wrt to its parent:
              // we dump it and all included surface components as independant volumes.
              cc_volume_ids[ncc_id] = next_volume_id++;
              error_codes.push_back(INCOMPATIBLE_ORIENTATION);
              if (used_as_a_predicate) return 0;
              for(std::size_t nncc_id : nested_cc_per_cc[ncc_id])
              {
                cc_handled.set(nncc_id);
                cc_volume_ids[nncc_id] = next_volume_id++;
                error_codes.push_back(INCOMPATIBLE_ORIENTATION);
              }
              continue;
            }
          }
          if (is_involved_in_self_intersection[ncc_id])
          {
            cc_volume_ids[ncc_id] = next_volume_id++;
            error_codes.push_back(VOLUME_INTERSECTION);
            for(std::size_t nncc_id : nested_cc_per_cc[ncc_id])
            {
              if (cc_handled.test(nncc_id))
              {
                error_codes[ cc_volume_ids[nncc_id] ] = VOLUME_INTERSECTION;
                continue;
              }
              cc_handled.set(nncc_id);
              cc_volume_ids[nncc_id] = next_volume_id++;
              error_codes.push_back(VOLUME_INTERSECTION);
            }
            continue;
          }
          cc_volume_ids[ncc_id] = cc_volume_ids[cc_id];
        }
      }
    }
    if (used_as_a_predicate)
    {
      internal::copy_container_content(is_cc_outward_oriented, parameters::get_parameter(np, internal_np::is_cc_outward_oriented));
      return 1;
    }

  // merge nested_cc_per_cc and nested_cc_per_cc_shared
  // (done after the volume creation to assign a CC to a unique volume)
    for(std::size_t id=0; id<nb_cc; ++id)
    {
      if (!nested_cc_per_cc_shared[id].empty())
        nested_cc_per_cc[id].insert(nested_cc_per_cc[id].end(),
                                    nested_cc_per_cc_shared[id].begin(),
                                    nested_cc_per_cc_shared[id].end());
    }


  // extract direct nested parent (more than one in case of self-intersection)
    for(std::size_t cc_id=0; cc_id<nb_cc; ++cc_id)
    {
      for(std::size_t ncc_id : nested_cc_per_cc[cc_id])
      {
        if (nesting_levels[cc_id]+1 == nesting_levels[ncc_id])
          nesting_parents[ncc_id].push_back(cc_id);
      }
    }

  // update volume id map
    for(std::size_t cc_id=0; cc_id<nb_cc; ++cc_id)
    {
      for(face_descriptor fd : faces_per_cc[cc_id])
        put(volume_id_map, fd, cc_volume_ids[cc_id]);
    }
  }
  else
  {
    for(face_descriptor fd : faces(tm))
    {
      std::size_t cc_id = face_cc[ get(fid_map, fd) ];
      put(volume_id_map, fd, cc_volume_ids[cc_id]);
    }
  }

  CGAL_assertion(next_volume_id == error_codes.size());
  internal::copy_container_content(error_codes, parameters::get_parameter(np, internal_np::error_codes));
  internal::copy_nested_parents(nesting_parents, parameters::get_parameter(np, internal_np::volume_inclusions));
  internal::copy_container_content(nesting_levels, parameters::get_parameter(np, internal_np::nesting_levels));
  internal::copy_container_content(cc_volume_ids, parameters::get_parameter(np, internal_np::connected_component_id_to_volume_id));
  internal::copy_container_content(is_cc_outward_oriented, parameters::get_parameter(np, internal_np::is_cc_outward_oriented));
  internal::set_cc_intersecting_pairs(self_intersecting_cc, parameters::get_parameter(np, internal_np::intersecting_volume_pairs_output_iterator));

  return next_volume_id;
}

/** \ingroup PMP_orientation_grp
 *
 * indicates if `tm` bounds a volume.
 * See \ref coref_def_subsec for details.
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm a closed triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * @pre `CGAL::is_closed(tm)`
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{is_cc_outward_oriented}
 *     \cgalParamDescription{For each connected component identified using its id `ccid`, the output of `is_outward_oriented()`
 *                           called on the triangle mesh corresponding to this connected component
 *                           is the value at the position `ccid` in the container.
 *                           The size of the container is exactly the number of connected components.}
 *     \cgalParamType{a `reference_wrapper` (either from `boost` or the standard library) containing
 *                    a reference to a container that must be a model of the `BackInsertionSequence` concept,
 *                    with value type `bool`}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see `CGAL::Polygon_mesh_processing::orient_to_bound_a_volume()`
 */
template <class TriangleMesh, class NamedParameters>
bool does_bound_a_volume(const TriangleMesh& tm, const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;

  CGAL_precondition(is_closed(tm));
  CGAL_precondition(is_triangle_mesh(tm));

  Constant_property_map<face_descriptor, std::size_t> vidmap(0); // dummy map not used
  std::size_t res =
    volume_connected_components(tm, vidmap, np.do_orientation_tests(true)
                                              .i_used_as_a_predicate(true));
  CGAL_assertion(res==0 || res==1);

  return res!=0;
}

/// \cond SKIP_IN_MANUAL
template <class TriangleMesh>
bool does_bound_a_volume(const TriangleMesh& tm)
{
  return does_bound_a_volume(tm, parameters::all_default());
}

template <class TriangleMesh, class VolumeFaceIndexMap>
std::size_t volume_connected_components(const TriangleMesh& tm, VolumeFaceIndexMap volume_id_map)
{
  return volume_connected_components(tm, volume_id_map, parameters::all_default());
}
/// \endcond


/** \ingroup PMP_orientation_grp
 *
 * orients the connected components of `tm` to make it bound a volume.
 * See \ref coref_def_subsec for a precise definition.
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters
 *
 * @param tm a closed triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * @pre `CGAL::is_closed(tm)`
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
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
 *     \cgalParamExtra{If the outer connected components are inward oriented,
 *                     it means that the infinity will be considered as part of the volume bounded by `tm`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see `CGAL::Polygon_mesh_processing::does_bound_a_volume()`
 */
template <class TriangleMesh, class NamedParameters>
void orient_to_bound_a_volume(TriangleMesh& tm,
                              const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh>                                        Graph_traits;
  typedef typename Graph_traits::face_descriptor                                   face_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type    Vpm;
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::type FaceIndexMap;

  CGAL_assertion(is_closed(tm));
  CGAL_assertion(is_triangle_mesh(tm));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool orient_outward = choose_parameter(get_parameter(np, internal_np::outward_orientation),true);

  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  FaceIndexMap fid_map = CGAL::get_initialized_face_index_map(tm, np);

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));
  std::vector<std::size_t> nesting_levels;
  std::vector<bool> is_cc_outward_oriented;
  Constant_property_map<face_descriptor, std::size_t> vidmap(0); // dummy map not used

  volume_connected_components(tm, vidmap,
                              parameters::vertex_point_map(vpm)
                                          .nesting_levels(boost::ref(nesting_levels))
                                          .face_connected_component_map(bind_property_maps(fid_map,make_property_map(face_cc)))
                                          .i_used_for_volume_orientation(true)
                                          .do_orientation_tests(true)
                                          .is_cc_outward_oriented(boost::ref(is_cc_outward_oriented))
                                          );

  // set the connected component id of each face


  if (nesting_levels.empty()) //case 1 cc
  {
    if( orient_outward != is_cc_outward_oriented[0])
      reverse_face_orientations(faces(tm), tm);
    return ;
  }
  std::size_t nb_cc = nesting_levels.size();
  boost::dynamic_bitset<> cc_to_reverse(nb_cc, 0);
  for(std::size_t i=0; i<nb_cc; ++i)
  {
    if ( ((nesting_levels[i]%2==0) == orient_outward) != is_cc_outward_oriented[i] )
    {
      cc_to_reverse.set(i);
    }
  }

  std::vector<face_descriptor> faces_to_reverse;
  for (face_descriptor f : faces(tm))
    if ( cc_to_reverse.test( face_cc[get(fid_map, f)] ) )
      faces_to_reverse.push_back(f);

  reverse_face_orientations(faces_to_reverse, tm);
}

template <class TriangleMesh>
void orient_to_bound_a_volume(TriangleMesh& tm)
{
  orient_to_bound_a_volume(tm, parameters::all_default());
}

/*!
 * \ingroup PMP_orientation_grp
 * reverses the connected components of `tm` having compatible boundary cycles
 * that could be merged if their orientation were made compatible, and stitches them.
 * Connected components are examined by increasing number of faces.
 *
 * @tparam PolygonMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters
 *
 * @param pm a surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `pm`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `pm` a unique index
 *                           between `0` and `num_faces(pm) - 1`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_number_of_faces}
 *     \cgalParamDescription{If not `0`, a connected component is considered reversible only
 *                           if it has no more faces than the value given.
 *                           Otherwise, it is always considered reversible.}
 *     \cgalParamType{`std::size_t`}
 *     \cgalParamDefault{`0`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template <class PolygonMesh, class NamedParameters>
void merge_reversible_connected_components(PolygonMesh& pm,
                                           const NamedParameters& np)
{
  typedef boost::graph_traits<PolygonMesh> GrT;
  typedef typename GrT::face_descriptor face_descriptor;
  typedef typename GrT::halfedge_descriptor halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type Vpm;

  typedef typename boost::property_traits<Vpm>::value_type Point_3;
  Vpm vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                         get_const_property_map(vertex_point, pm));

  typedef std::size_t F_cc_id;
  typedef std::size_t B_cc_id;


  typedef typename CGAL::GetInitializedFaceIndexMap<PolygonMesh, NamedParameters>::const_type Fidmap;
  Fidmap fim = CGAL::get_initialized_face_index_map(pm, np);

  typedef dynamic_face_property_t<F_cc_id>                   Face_property_tag;
  typedef typename boost::property_map<PolygonMesh, Face_property_tag>::type   Face_cc_map;
  Face_cc_map f_cc_ids  = get(Face_property_tag(), pm);
  F_cc_id nb_cc = connected_components(pm, f_cc_ids, parameters::face_index_map(fim));

  std::vector<std::size_t> nb_faces_per_cc(nb_cc, 0);
  for (face_descriptor f : faces(pm))
    nb_faces_per_cc[ get(f_cc_ids, f) ]+=1;

  std::map< std::pair<Point_3, Point_3>, std::vector<halfedge_descriptor> > border_hedges_map;
  std::vector<halfedge_descriptor> border_hedges;
  typedef typename boost::property_map<PolygonMesh, dynamic_halfedge_property_t<B_cc_id> >::type H_to_bcc_id;
  H_to_bcc_id h_bcc_ids = get(dynamic_halfedge_property_t<B_cc_id>(), pm);
  const B_cc_id base_value(-1);
  const B_cc_id FILTERED_OUT(-2);

  // collect border halfedges
  for (halfedge_descriptor h : halfedges(pm))
    if ( is_border(h, pm) )
    {
      put(h_bcc_ids, h, base_value);
      border_hedges.push_back(h);
    }

  // compute the border cc id of all halfedges and mark those duplicated in their own cycle
  B_cc_id bcc_id=0;
  for (halfedge_descriptor h : border_hedges)
  {
    if (get(h_bcc_ids,h) == base_value)
    {
      typedef std::map< std::pair<Point_3, Point_3>, halfedge_descriptor> Hmap;
      Hmap hmap;
      for (halfedge_descriptor hh : halfedges_around_face(h, pm))
      {
        std::pair< typename Hmap::iterator, bool > insert_res =
            hmap.insert(
              std::make_pair(
                make_sorted_pair(get(vpm, source(hh, pm)),
                                 get(vpm, target(hh,pm))), hh) );
        if (insert_res.second)
          put(h_bcc_ids, hh, bcc_id);
        else
        {
          put(h_bcc_ids, hh, FILTERED_OUT);
          put(h_bcc_ids, insert_res.first->second, FILTERED_OUT);
        }
      }
      ++bcc_id;
    }
  }

  // fill endpoints -> hedges
  for (halfedge_descriptor h : border_hedges)
  {
    if ( get(h_bcc_ids, h) != FILTERED_OUT)
      border_hedges_map[std::make_pair(get(vpm, source(h, pm)), get(vpm, target(h, pm)))].push_back(h);
  }

  // max nb of faces for a CC to be reversed
  const std::size_t threshold =
    parameters::choose_parameter( parameters::get_parameter(np, internal_np::maximum_number_of_faces), 0);

  std::vector<bool> border_cycle_to_ignore(bcc_id, false);
  std::vector<F_cc_id> cycle_f_cc_id(bcc_id);
  std::vector< std::vector<F_cc_id> > patch_neighbors(nb_cc);

  for (const auto& p : border_hedges_map)
  {
    const std::vector<halfedge_descriptor>& hedges = p.second;
    switch(hedges.size())
    {
      case 1:
        // isolated border hedge nothing to do
      break;
      case 2:
      {
        F_cc_id cc_id_0 = get(f_cc_ids, face(opposite(hedges[0], pm), pm)),
                cc_id_1 = get(f_cc_ids, face(opposite(hedges[1], pm), pm));

        if (cc_id_0!=cc_id_1)
        {
          cycle_f_cc_id[ get(h_bcc_ids, hedges[0]) ] = cc_id_0;
          cycle_f_cc_id[ get(h_bcc_ids, hedges[1]) ] = cc_id_1;
          // WARNING: we might have duplicates here but it is not important for our usage
          patch_neighbors[cc_id_0].push_back(cc_id_1);
          patch_neighbors[cc_id_1].push_back(cc_id_0);
          break;
        }
        CGAL_FALLTHROUGH;
      }
      default:
        for (halfedge_descriptor h : hedges)
          border_cycle_to_ignore[get(h_bcc_ids, h)]=true;
    }
  }

  // sort the connected components with potential matches using their number
  // of faces (sorted by decreasing number of faces)
  std::set<F_cc_id> ccs_to_reverse;
  std::vector<bool> reversible(nb_cc, false);
  std::set< F_cc_id, std::function<bool(F_cc_id,F_cc_id)> > queue(
    [&nb_faces_per_cc](F_cc_id i, F_cc_id j)
    {return nb_faces_per_cc[i]==nb_faces_per_cc[j] ? i<j : nb_faces_per_cc[i]>nb_faces_per_cc[j];}
  );

  for (B_cc_id i=0; i<bcc_id; ++i)
  {
    if ( !border_cycle_to_ignore[i] )
    {
      reversible[ cycle_f_cc_id[i] ] = true;
      queue.insert(cycle_f_cc_id[i]);
    }
  }

  // consider largest CC selected and reverse the neighbor patches if
  // not already reversed or not marked as reversible
  while( !queue.empty() )
  {
    F_cc_id f_cc_id = *queue.begin();
    queue.erase( queue.begin() );
    CGAL_assertion( reversible[f_cc_id] );
    for (F_cc_id id : patch_neighbors[f_cc_id])
    {
      if (reversible[id] && (threshold==0 || threshold >= nb_faces_per_cc[id]))
      {
        CGAL_assertion( nb_faces_per_cc[f_cc_id] >= nb_faces_per_cc[id] );
        ccs_to_reverse.insert(id);
        reversible[id]=false;
        queue.erase(id);
      }
    }
  }

  // reverse ccs and stitches boundaries
  std::vector<face_descriptor> faces_to_reverse;
  for (face_descriptor f : faces(pm))
    if ( ccs_to_reverse.count( get(f_cc_ids, f) ) != 0 )
      faces_to_reverse.push_back(f);

  if ( !faces_to_reverse.empty() )
  {
    reverse_face_orientations(faces_to_reverse, pm);
    stitch_borders(pm, np);
  }
}

template <class PolygonMesh>
void merge_reversible_connected_components(PolygonMesh& pm)
{
  merge_reversible_connected_components(pm, parameters::all_default());
}
} // namespace Polygon_mesh_processing
} // namespace CGAL
#endif // CGAL_ORIENT_POLYGON_MESH_H
