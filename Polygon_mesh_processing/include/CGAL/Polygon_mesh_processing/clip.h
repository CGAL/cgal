// Copyright (c) 2016 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/convex_hull_3.h>

namespace CGAL{
namespace Polygon_mesh_processing {

namespace internal
{
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip_open_impl(      TriangleMesh& tm,
                     TriangleMesh& clipper,
               const NamedParameters1& np_tm,
               const NamedParameters2& np_c)
{
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::type Vpm;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters2>::type GeomTraits;
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

// First build an AABB-tree of the clipper triangles as it will be modified
  typedef std::vector<typename GeomTraits::Triangle_3> Clipper_triangles;
  typedef typename Clipper_triangles::iterator Tr_iterator;
  typedef CGAL::AABB_triangle_primitive<GeomTraits, Tr_iterator> Primitive;
  typedef CGAL::AABB_traits<GeomTraits, Primitive> AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits> Clipper_tree;

  // vector of clipper triangles
  Clipper_triangles clipper_triangles;
  clipper_triangles.reserve( num_faces(clipper) );
  Vpm vpm_c = boost::choose_param(boost::get_param(np_c, internal_np::vertex_point),
                                  get_property_map(vertex_point, clipper));
  BOOST_FOREACH(face_descriptor f, faces(clipper))
  {
    halfedge_descriptor h = halfedge(f, clipper);
    clipper_triangles.push_back( typename GeomTraits::Triangle_3(
      get(vpm_c, source(h, clipper)),
      get(vpm_c, target(h, clipper)),
      get(vpm_c, target(next(h, clipper), clipper)) ) );
  }
  // tree
  Clipper_tree clipper_tree(clipper_triangles.begin(), clipper_triangles.end());
  // predicate functor
  Side_of_triangle_mesh<TriangleMesh, GeomTraits, Vpm, Clipper_tree> side_of(clipper_tree);

// Second corefine the meshes
  typedef CGAL::dynamic_edge_property_t<bool> Ecm_tag;
  typedef typename boost::property_map<TriangleMesh, Ecm_tag>::type Ecm;
  Ecm ecm = get(Ecm_tag(), tm);

  corefine(tm, clipper, np_tm.edge_is_constrained_map(ecm), np_c);

// Extract connected components
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters1>::type Fid_map;

  Fid_map fid_map = boost::choose_param(boost::get_param(np_tm, internal_np::face_index),
                                        get_property_map(boost::face_index, tm));
  Vpm vpm1 = boost::choose_param(boost::get_param(np_tm, internal_np::vertex_point),
                                 get_property_map(vertex_point, tm));

  typedef CGAL::dynamic_vertex_property_t<std::size_t> Vid_tag;
  typedef typename boost::property_map<TriangleMesh, Vid_tag>::type Vid_map;
  Vid_map vid_map = get(Vid_tag(), tm);

  // init indices if needed
  helpers::init_face_indices(tm, fid_map);
  helpers::init_vertex_indices(tm, vid_map);

  // set the connected component id of each face
  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));
  std::size_t nb_cc =
    connected_components(tm,
                         bind_property_maps(fid_map, make_property_map(face_cc)),
                         parameters::face_index_map(fid_map).
                         edge_is_constrained_map(ecm));


  boost::dynamic_bitset<> cc_not_handled(nb_cc);
  cc_not_handled.set();
  std::vector <std::size_t> ccs_to_remove;

  BOOST_FOREACH(face_descriptor f, faces(tm))
  {
    std::size_t cc_id = face_cc[ get(fid_map, f) ];
    if ( !cc_not_handled.test(cc_id) ) continue;

    halfedge_descriptor h=halfedge(f, tm);
    for(int i=0;i<3;++i)
    {
      // look for a vertex not on a constrained edge
      bool no_marked_edge=true;
      BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_target(h, tm))
        if ( get(ecm, edge(h2, tm)) ){
          no_marked_edge=false;
          break;
        }
      if (no_marked_edge){
        if ( side_of( get(vpm1, target(h, tm) ) ) == ON_UNBOUNDED_SIDE )
          ccs_to_remove.push_back(cc_id);
        cc_not_handled.reset(cc_id);
        break;
      }
      h=next(h, tm);
    }
    if (!cc_not_handled.any()) break;
  }

  if (cc_not_handled.any())
  {
    // A patch without no vertex incident to a non-constrained edges
    //  is a coplanar patch: drop it or keep it!
    if (!boost::choose_param(boost::get_param(np_tm, internal_np::use_compact_clipper), true))
    {
      for (std::size_t cc_id = cc_not_handled.find_first();
                       cc_id < cc_not_handled.npos;
                       cc_id = cc_not_handled.find_next(cc_id))
      {
        ccs_to_remove.push_back(cc_id);
      }
    }
  }
// Filter out the cc
  remove_connected_components(tm,
    ccs_to_remove,
    bind_property_maps(fid_map, make_property_map(face_cc)),
    parameters::vertex_index_map(vid_map));

  return true;
}

/// \todo remove convex_hull_3
template <class Plane_3,
          class TriangleMesh,
          class NamedParameters>
Oriented_side
clip_to_bbox(const Plane_3& plane,
             const Bbox_3& bbox,
                   TriangleMesh& tm_out,
             const NamedParameters& np )
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;

  Vpm vpm_out = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                    get_property_map(boost::vertex_point, tm_out));


  cpp11::array<Point_3,8> corners= {{
    Point_3(bbox.xmin(),bbox.ymin(),bbox.zmin()),
    Point_3(bbox.xmin(),bbox.ymax(),bbox.zmin()),
    Point_3(bbox.xmax(),bbox.ymax(),bbox.zmin()),
    Point_3(bbox.xmax(),bbox.ymin(),bbox.zmin()),
    Point_3(bbox.xmin(),bbox.ymin(),bbox.zmax()),
    Point_3(bbox.xmin(),bbox.ymax(),bbox.zmax()),
    Point_3(bbox.xmax(),bbox.ymax(),bbox.zmax()),
    Point_3(bbox.xmax(),bbox.ymin(),bbox.zmax())
  }};

  cpp11::array<CGAL::Oriented_side,8> orientations = {{
    plane.oriented_side(corners[0]),
    plane.oriented_side(corners[1]),
    plane.oriented_side(corners[2]),
    plane.oriented_side(corners[3]),
    plane.oriented_side(corners[4]),
    plane.oriented_side(corners[5]),
    plane.oriented_side(corners[6]),
    plane.oriented_side(corners[7])
  }};

  std::vector<Point_3> points;

  // look for intersections on edges
  cpp11::array<int,24> edge_indices = {{ // 2 *12 edges
    0,1, 1,2, 2,3, 3,0, // bottom face edges
    4,5, 5,6, 6,7, 7,4, // top face edges
    0,4, 1,5, 2,6, 3,7
  }};

  for (int i=0; i<12; ++i)
  {
    int i1=edge_indices[2*i], i2=edge_indices[2*i+1];
    if (orientations[i1]==ON_ORIENTED_BOUNDARY) continue;
    if (orientations[i2]==ON_ORIENTED_BOUNDARY) continue;
    if (orientations[i1]!=orientations[i2])
      points.push_back(
        boost::get<Point_3>(
          *intersection(plane, Segment_3(corners[i1], corners[i2]) )
        )
      );
  }


  Oriented_side last_os = ON_ORIENTED_BOUNDARY;
  for (int i=0; i<8; ++i)
    if (orientations[i]!=ON_ORIENTED_BOUNDARY)
    {
      if (last_os==ON_ORIENTED_BOUNDARY)
        last_os=orientations[i];
      else
      {
        if(last_os!=orientations[i])
        {
          last_os=ON_ORIENTED_BOUNDARY;
          break;
        }
      }
    }

  // the intersection is the full bbox
  if (last_os!=ON_ORIENTED_BOUNDARY)
    return last_os;

  //add points on negative side and on the plane
  for (int i=0; i<8; ++i)
    if (orientations[i]!=ON_POSITIVE_SIDE)
      points.push_back(corners[i]);

  // take the convex hull of the points on the negative side+intersection points
  // overkill...
  TriangleMesh ch_tm;
  CGAL::convex_hull_3(points.begin(), points.end(), ch_tm);
  copy_face_graph(ch_tm, tm_out,
                  Emptyset_iterator(), Emptyset_iterator(), Emptyset_iterator(),
                  get(vertex_point, ch_tm), vpm_out);
  return ON_ORIENTED_BOUNDARY;
}

} // end of internal namespace

/**
  * \ingroup PMP_corefinement_grp
  * clips `tm` by keeping the part that is inside the volume \link coref_def_subsec bounded \endlink
  * by `clipper`.
  * If `tm` is closed, the clipped part can be closed too if the named parameter `clip_volumes` is set to `true`.
  * \attention With the current implementation, `clipper` will be modified (refined with the intersection with `tm`).
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(clipper)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(clipper)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
  *                      as a named parameter, then it must be initialized.
  *
  * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param clipper triangulated surface mesh used to clip `tm`
  * @param np_tm optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  * @param np_c optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map}
  *     the property map with the points associated to the vertices of `tm` (`clipper`).
  *     If this parameter is omitted, an internal property map for
  *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
  *   \cgalParamEnd
  *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm` (`clipper`).
  *     Note that if the property map is writable, the indices of the faces
  *     of `tm` and `clipper` will be set after the refining `tm` with the intersection with `plane`.
  *   \cgalParamEnd
  *   \cgalParamBegin{new_face_visitor} a class model of `PMPCorefinementNewFaceVisitor`
  *                                     that is used to track the creation of new faces.
  *   \cgalParamEnd
  *   \cgalParamBegin{throw_on_self_intersection} if `true`,
  *      the set of triangles closed to the intersection of `tm` and `clipper` will be
  *      checked for self-intersections and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *      will be thrown if at least one is found.
  *   \cgalParamEnd
  *   \cgalParamBegin{clip_volumes} if `true` and `tm` is closed, the clipping will be done on
  *      the volume \link coref_def_subsec bounded \endlink by `tm` rather than on its surface
  *      (i.e. `tm` will be kept closed).
  *   \cgalParamEnd
  *   \cgalParamBegin{use_compact_clipper} if `false` and `clip_volumes` is `false` and `tm` is open, the parts of `tm` coplanar with `clipper`
  *                                        will not be part of the output.
  *   \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold.
  *         If `false` is returned `tm` and `clipper` are only corefined.
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip(      TriangleMesh& tm,
           TriangleMesh& clipper,
     const NamedParameters1& np_tm,
     const NamedParameters2& np_c)
{
  const bool close =
    boost::choose_param(boost::get_param(np_tm, internal_np::clip_volumes), false);

  if (close && is_closed(tm))
    return corefine_and_compute_intersection(tm, clipper, tm, np_tm, np_c);

  return internal::clip_open_impl(tm, clipper, np_tm, np_c);
}

/**
  * \ingroup PMP_corefinement_grp
  * clips `tm` by keeping the part that is on the negative side of `plane` (side opposite to its normal vector).
  * If `tm` is closed, the clipped part can be closed too if the named parameter `clip_volumes` is set to `true`.
  *
  * \note In the current implementation it is not possible to set the vertex point map and the default will be used.
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      If `TriangleMesh` has an internal property map for `CGAL::face_index_t`,
  *                      as a named parameter, then it must be initialized.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param plane plane whose negative side defines the half-space to intersect `tm` with.
  *              `Plane_3` is the plane type for the same CGAL kernel as the point of the vertex point map of `tm`.
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm`.
  *     Note that if the property map is writable, the indices of the faces
  *     of `tm` will be set after the refining of `tm` with the intersection with `plane`.
  *   \cgalParamEnd
  *   \cgalParamBegin{new_face_visitor} a class model of `PMPCorefinementNewFaceVisitor`
  *                                     that is used to track the creation of new faces.
  *   \cgalParamEnd
  *   \cgalParamBegin{throw_on_self_intersection} if `true`,
  *      the set of triangles closed to the intersection of `tm` and `plane` will be
  *      checked for self-intersections and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *      will be thrown if at least one is found.
  *   \cgalParamEnd
  *   \cgalParamBegin{clip_volumes} if `true` and `tm` is closed, the clipping will be done on
  *      the volume \link coref_def_subsec bounded \endlink by `tm` rather than on its surface
  *      (i.e. `tm` will be kept closed).
  *   \cgalParamEnd
  *   \cgalParamBegin{use_compact_clipper} if `false` and `clip_volumes` is `false` and `tm` is open, the parts of `tm` coplanar with `plane`
  *                                        will not be part of the output.
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold.
  *         If `false` is returned `tm` is only refined by the intersection with `plane`.
  */
template <class TriangleMesh,
          class NamedParameters>
bool clip(      TriangleMesh& tm,
          #ifdef DOXYGEN_RUNNING
          const Plane_3& plane,
          #else
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::Plane_3& plane,
          #endif
          const NamedParameters& np)
{
  if( boost::begin(faces(tm))==boost::end(faces(tm)) ) return true;

  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm);

  //extend the bbox a bit to avoid border cases
  double xd=(bbox.xmax()-bbox.xmin())/100;
  double yd=(bbox.ymax()-bbox.ymin())/100;
  double zd=(bbox.zmax()-bbox.zmin())/100;
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
  TriangleMesh clipper;
  Oriented_side os = internal::clip_to_bbox(plane, bbox, clipper, parameters::all_default());

  switch(os)
  {
    case ON_NEGATIVE_SIDE:
      return true; // nothing to clip, the full mesh is on the negative side
    case ON_POSITIVE_SIDE:
      clear(tm); // clear the mesh that is fully on the positive side
      return true;
    default:
      break;
  }
  return clip(tm, clipper, np, parameters::all_default());
}

/// \cond SKIP_IN_MANUAL
// convenience overloads
template <class TriangleMesh,
          class NamedParameters>
bool clip(      TriangleMesh& tm,
          const typename GetGeomTraits<TriangleMesh,
                                       NamedParameters>::type::Plane_3& plane)
{
  return clip(tm, plane, parameters::all_default());
}

// convenience overload
template <class TriangleMesh,
          class NamedParameters1>
bool
clip(      TriangleMesh& tm,
           TriangleMesh& clipper,
     const NamedParameters1& np_tm)
{
  return clip(tm, clipper, np_tm, parameters::all_default());
}

// convenience overload
template <class TriangleMesh>
bool
clip(      TriangleMesh& tm,
           TriangleMesh& clipper)
{
  return clip(tm, clipper, parameters::all_default());
}
/// \endcond

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_H
