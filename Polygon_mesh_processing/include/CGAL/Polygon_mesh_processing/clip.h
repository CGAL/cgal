// Copyright (c) 2016 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/iterator.h>

#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <boost/property_map/property_map.hpp>

#include <unordered_map>

namespace CGAL{
namespace Polygon_mesh_processing {

namespace internal
{

template <class Geom_traits, class Plane_3, class Point_3>
int
inter_pt_index(int i, int j,
               const Plane_3& plane,
               std::vector<Point_3>& points,
               std::map<std::pair<int,int>, int>& id_map)
{
  std::pair<std::map<std::pair<int,int>, int>::iterator, bool> res =
    id_map.insert(std::make_pair(make_sorted_pair(i,j),
                  static_cast<int> (points.size())));
  if(res.second)
    points.push_back(
      typename Geom_traits::Construct_plane_line_intersection_point_3()
        (plane, points[i], points[j]));

  return res.first->second;
}

template <class Plane_3,
          class TriangleMesh,
          class NamedParameters>
Oriented_side
clip_to_bbox(const Plane_3& plane,
             const Bbox_3& bbox,
                   TriangleMesh& tm_out,
             const NamedParameters& np)
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;

  Vpm vpm_out = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                             get_property_map(boost::vertex_point, tm_out));


  std::vector<Point_3> corners(8);
  corners[0] = Point_3(bbox.xmin(),bbox.ymin(),bbox.zmin());
  corners[1] = Point_3(bbox.xmin(),bbox.ymax(),bbox.zmin());
  corners[2] = Point_3(bbox.xmax(),bbox.ymax(),bbox.zmin());
  corners[3] = Point_3(bbox.xmax(),bbox.ymin(),bbox.zmin());
  corners[4] = Point_3(bbox.xmin(),bbox.ymin(),bbox.zmax());
  corners[5] = Point_3(bbox.xmin(),bbox.ymax(),bbox.zmax());
  corners[6] = Point_3(bbox.xmax(),bbox.ymax(),bbox.zmax());
  corners[7] = Point_3(bbox.xmax(),bbox.ymin(),bbox.zmax());

  std::array<CGAL::Oriented_side,8> orientations = {{
    plane.oriented_side(corners[0]),
    plane.oriented_side(corners[1]),
    plane.oriented_side(corners[2]),
    plane.oriented_side(corners[3]),
    plane.oriented_side(corners[4]),
    plane.oriented_side(corners[5]),
    plane.oriented_side(corners[6]),
    plane.oriented_side(corners[7])
  }};

  // description of faces of the bbox
  std::array<int, 24> face_indices =
    {{ 0, 1, 2, 3,
       2, 1, 5, 6,
       3, 2, 6, 7,
       1, 0, 4, 5,
       4, 0, 3, 7,
       6, 5, 4, 7 }};

  std::map<std::pair<int,int>, int> id_map;
  std::vector< std::vector<int> > output_faces(6);
  bool all_in = true;
  bool all_out = true;
  std::set<int> in_point_ids; // to collect the set of points in the clipped bbox

  // for each face of the bbox, we look for intersection of the plane with its edges
  for(int i=0; i<6; ++i)
  {
    for(int k=0; k< 4; ++k)
    {
      int current_id = face_indices[4*i + k];
      int next_id = face_indices[4*i + (k+1)%4];

      switch(orientations[ current_id ])
      {
        case ON_NEGATIVE_SIDE:
        {
          all_out=false;
          // point on or on the negative side
          output_faces[i].push_back(current_id);
          in_point_ids.insert(output_faces[i].back());
          // check for intersection of the edge
          if(orientations[ next_id ] == ON_POSITIVE_SIDE)
          {
            output_faces[i].push_back(
              inter_pt_index<Geom_traits>(current_id, next_id, plane, corners, id_map));
            in_point_ids.insert(output_faces[i].back());
          }
          break;
        }
        case ON_POSITIVE_SIDE:
        {
          all_in = false;
          // check for intersection of the edge
          if(orientations[ next_id ] == ON_NEGATIVE_SIDE)
          {
            output_faces[i].push_back(
              inter_pt_index<Geom_traits>(current_id, next_id, plane, corners, id_map));
            in_point_ids.insert(output_faces[i].back());
          }
          break;
        }
        case ON_ORIENTED_BOUNDARY:
        {
          output_faces[i].push_back(current_id);
          in_point_ids.insert(output_faces[i].back());
        }
      }
    }
    if(output_faces[i].size() < 3){
      CGAL_assertion(output_faces[i].empty() ||
                     (output_faces[i].front()<8 && output_faces[i].back()<8));
      output_faces[i].clear(); // edge of the bbox included in the plane
    }
  }

  // the intersection is the full bbox
  if(all_in) return ON_NEGATIVE_SIDE;
  if(all_out) return ON_POSITIVE_SIDE;

  // build the clipped bbox
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;

  std::map<int, vertex_descriptor> out_vertices;
  for(int i : in_point_ids)
  {
    vertex_descriptor v = add_vertex(tm_out);
    out_vertices.insert(std::make_pair(i, v));
    put(vpm_out, v, corners[i]);
  }

  std::map< std::pair<int,int>, halfedge_descriptor> hedge_map;
  const halfedge_descriptor null_hedge = graph_traits::null_halfedge();
  const face_descriptor null_fd = graph_traits::null_face();
  for(const std::vector<int>& findices : output_faces)
  {
    if(findices.empty()) continue;
    const face_descriptor fd=add_face(tm_out);
    int prev_id = findices.back();

    // create of recover face boundary halfedges
    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(findices.size());
    for(int current_id : findices)
    {
      vertex_descriptor src = out_vertices[prev_id], tgt = out_vertices[current_id];

      std::pair<typename std::map< std::pair<int,int>,
                halfedge_descriptor>::iterator, bool> res =
        hedge_map.insert(std::make_pair(std::make_pair(prev_id, current_id), null_hedge));
      if(res.second)
      {
        res.first->second = halfedge(add_edge(tm_out), tm_out);
        hedge_map.insert(std::make_pair(std::make_pair(current_id, prev_id),
                            opposite(res.first->second, tm_out)));
        set_face(opposite(res.first->second, tm_out), null_fd, tm_out);

      }
      hedges.push_back(res.first->second);
      // set edge source and target
      set_target(hedges.back(), tgt, tm_out);
      set_target(opposite(hedges.back(), tm_out), src, tm_out);
      // set face pointer of halfedges
      set_face(hedges.back(), fd, tm_out);
      // set vertex halfedge
      set_halfedge(src, opposite(hedges.back(), tm_out), tm_out);
      set_halfedge(tgt, hedges.back(), tm_out);

      if(current_id==findices.front())
        set_halfedge(fd, hedges.back(), tm_out);

      prev_id = current_id;
    }
    CGAL_assertion(hedges.size() == findices.size());

    // set next/prev relationship
    halfedge_descriptor prev_h=hedges.back();
    for(halfedge_descriptor h : hedges)
    {
      set_next(prev_h, h, tm_out);
      prev_h = h;
    }
  }

  // handle the face of the plane:
  // look for a border halfedge and reconstruct the face of the plane
  // by turning around vertices inside the mesh constructed above
  // until we reach another border halfedge
  for(halfedge_descriptor h : halfedges(tm_out))
  {
    if(face(h, tm_out) == null_fd)
    {
      face_descriptor fd = add_face(tm_out);
      set_halfedge(fd, h, tm_out);

      halfedge_descriptor h_prev=h;
      halfedge_descriptor h_curr=h;
      do{
        h_curr=opposite(h_curr, tm_out);
        do{
          h_curr=opposite(prev(h_curr, tm_out), tm_out);
        } while(face(h_curr, tm_out) != null_fd && h_curr!=h);
        set_face(h_prev, fd, tm_out);
        set_next(h_prev, h_curr, tm_out);
        if(h_curr==h)
          break;
        h_prev=h_curr;
      } while(true);
      break;
    }
  }
  CGAL_assertion(is_valid_polygon_mesh(tm_out));

  // triangulate the faces
  CGAL::Polygon_mesh_processing::triangulate_faces(tm_out, np);

  return ON_ORIENTED_BOUNDARY;
}

template <class TriangleMesh, class Ecm, class VPM>
void split_along_edges(TriangleMesh& tm,
                       Ecm ecm,
                       VPM vpm)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  std::vector<edge_descriptor> shared_edges;
  for(edge_descriptor e : edges(tm))
    if(get(ecm, e))
      shared_edges.push_back(e);

  std::size_t nb_shared_edges = shared_edges.size();
  std::vector<halfedge_descriptor> hedges_to_update;

  typedef CGAL::dynamic_halfedge_property_t<bool> H_tag;
  typename boost::property_map<TriangleMesh, H_tag>::type
    no_target_update = get(H_tag(), tm);

  std::vector< std::pair<halfedge_descriptor, vertex_descriptor> > vertices_to_duplicate;

  //collect border halfedges having as target one of the edge endpoints
  std::set<halfedge_descriptor> extra_border_hedges;
  for(std::size_t k=0; k<nb_shared_edges; ++k)
  {
    for(halfedge_descriptor h : halfedges_around_target(target(shared_edges[k], tm), tm))
      if(is_border(h, tm))
        extra_border_hedges.insert(h);
    for(halfedge_descriptor h : halfedges_around_target(source(shared_edges[k], tm), tm))
      if(is_border(h, tm))
        extra_border_hedges.insert(h);
  }

  for(halfedge_descriptor h : extra_border_hedges)
  {
    put(no_target_update, h, true);
    set_halfedge(target(h, tm), h, tm);
    hedges_to_update.push_back(h);
  }

  // now duplicate the edge and set its pointers
  for(std::size_t k=0; k<nb_shared_edges; ++k)
  {
    halfedge_descriptor h    = halfedge(shared_edges[k], tm);
    face_descriptor fh = face(h, tm);
    //add edge
    halfedge_descriptor new_hedge = halfedge(add_edge(tm), tm),
                        new_opp   = opposite(new_hedge,tm);

    vertex_descriptor vt = target(h, tm);
    vertex_descriptor vs = source(h, tm);

    //replace h with new_hedge
    set_next(new_hedge, next(h, tm), tm);
    set_next(prev(h, tm), new_hedge, tm);
    set_face(new_hedge, fh, tm);
    set_halfedge(fh, new_hedge, tm);

    set_target(new_hedge, vt, tm);
    set_target(new_opp, vs, tm);

    set_face(new_opp, GT::null_face(), tm);
    set_face(h, GT::null_face(), tm);

    // handle vertices to duplicate
    halfedge_descriptor h_vt = halfedge(vt, tm);
    if(get(no_target_update, h_vt))
      vertices_to_duplicate.push_back(std::make_pair(h, vt));
    else
      set_halfedge(vt, h, tm);
    halfedge_descriptor h_vs = halfedge(vs, tm);
    if(get(no_target_update, h_vs))
      vertices_to_duplicate.push_back(std::make_pair(new_opp, vs));
    else
      set_halfedge(vs, new_opp, tm);

    hedges_to_update.push_back(h);
    put(no_target_update, h, true);
    hedges_to_update.push_back(new_opp);
    put(no_target_update, new_opp, true);

    CGAL_assertion(next(prev(new_hedge, tm), tm) == new_hedge);
    CGAL_assertion(prev(next(new_hedge, tm), tm) == new_hedge);
  }

  // update next/prev relationship
  for(halfedge_descriptor h : hedges_to_update)
  {
    CGAL_assertion(is_border(h, tm));
    halfedge_descriptor h_opp = opposite(h, tm);

    // set next pointer of h, visiting faces inside the patch we consider
    halfedge_descriptor candidate = opposite(prev(h_opp, tm), tm);
    while (!is_border(candidate, tm))
      candidate = opposite(prev(candidate, tm), tm);
    set_next(h, candidate, tm);
    CGAL_assertion(prev(next(h_opp, tm), tm)==h_opp);

    CGAL_assertion(prev(next(h, tm), tm) == h);
    CGAL_assertion(is_border(next(h, tm), tm));
  }

  for(const std::pair<halfedge_descriptor, vertex_descriptor>& p : vertices_to_duplicate)
  {
    vertex_descriptor nv = add_vertex(tm);
    put(vpm, nv, get(vpm, p.second));
    for(halfedge_descriptor h : halfedges_around_target(p.first, tm))
      set_target(h, nv, tm);
    set_halfedge(nv, p.first, tm);
   }

  CGAL_assertion_code(for(halfedge_descriptor h : hedges_to_update))
  {
    CGAL_assertion(next(prev(h, tm), tm) == h);
    CGAL_assertion(prev(next(h, tm), tm) == h);
  }

  for(halfedge_descriptor h : hedges_to_update)
  {
    for(halfedge_descriptor hh : halfedges_around_target(h, tm))
        if(h!=hh)
          set_target(hh, target(h, tm), tm);
  }

  CGAL_assertion(is_valid_polygon_mesh(tm));
}

} // end of internal namespace

/**
  * \ingroup PMP_corefinement_grp
  *
  * clips `tm` by keeping the part that is inside the volume \link coref_def_subsec bounded \endlink
  * by `clipper`.
  * If `tm` is closed, the clipped part can be closed too if the named parameter `clip_volume` is set to `true`.
  * See Subsection \ref coref_clip for more details.
  * \attention With the current implementation, `clipper` will be modified (refined with the intersection with `tm`).
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(clipper)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(clipper)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param clipper triangulated surface mesh used to clip `tm`
  * @param np_tm optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  * @param np_c optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm` (resp. `clipper`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm (resp. clipper))`}
  *     \cgalParamExtra{if this parameter is omitted an internal property map for `CGAL::vertex_point_t` must be available in `TriangleMesh`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `tm` (`clipper`) a unique index between `0` and `num_faces(tm (resp. clipper)) - 1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *     \cgalParamExtra{if the property map is writable, the indices of the faces of `tm` and `clipper` will be set after refining `tm` with the intersection with `clipper`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{if `true` the set of triangles closed to the intersection of `tm` and `clipper` will be
  *                           checked for self-intersections and `Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{clip_volume}
  *     \cgalParamDescription{if `true` and `tm` is closed, the clipping will be done on the volume
  *                           \link coref_def_subsec bounded \endlink by `tm` rather than on its surface
  *                           (i.e., `tm` will be kept closed)}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{use_compact_clipper}
  *     \cgalParamDescription{if `false`, the parts of `tm` coplanar with `clipper` will not be part of the output.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *     \cgalParamExtra{This option has an effect only if a surface and not a volume is clipped,
  *                     (i.e., if `clip_volume` is `false` or if `tm` is open).}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold.
  *         If `false` is returned `tm` and `clipper` are only corefined.
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
clip(TriangleMesh& tm,
     TriangleMesh& clipper,
     const NamedParameters1& np_tm,
     const NamedParameters2& np_c)
{
  const bool clip_volume =
    parameters::choose_parameter(parameters::get_parameter(np_tm, internal_np::clip_volume), false);

  if(clip_volume && is_closed(tm))
    return corefine_and_compute_intersection(tm, clipper, tm, np_tm, np_c);
  return corefine_and_compute_intersection(tm, clipper, tm,
                                           np_tm.use_bool_op_to_clip_surface(true),
                                           np_c);
}

/**
  * \ingroup PMP_corefinement_grp
  * clips `tm` by keeping the part that is on the negative side of `plane` (side opposite to its normal vector).
  * If `tm` is closed, the clipped part can be closed too if the named parameter `clip_volume` is set to `true`.
  * See Subsection \ref coref_clip for more details.
  *
  * \note In the current implementation it is not possible to set the vertex point map and the default will be used. `Plane_3` must be
  * from the same %Kernel as the point of the vertex point map.
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
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
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{if `true`, the set of triangles closed to the intersection of `tm`
  *                           and `plane` will be checked for self-intersections
  *                           and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{clip_volume}
  *     \cgalParamDescription{if `true` and `tm` is closed, the clipping will be done on
  *                           the volume \link coref_def_subsec bounded \endlink by `tm`
  *                           rather than on its surface (i.e., `tm` will be kept closed)}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{use_compact_clipper}
  *     \cgalParamDescription{if `false` the parts of `tm` coplanar with `plane` will not be part of the output}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold.
  *         If `false` is returned `tm` is only refined by the intersection with `plane`.
  */
template <class TriangleMesh,
          class NamedParameters>
bool clip(TriangleMesh& tm,
#ifdef DOXYGEN_RUNNING
          const Plane_3& plane,
#else
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::Plane_3& plane,
#endif
          const NamedParameters& np)
{
  if(boost::begin(faces(tm))==boost::end(faces(tm))) return true;

  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm);

  //extend the bbox a bit to avoid border cases
  double xd=(std::max)(1.,(bbox.xmax()-bbox.xmin())/100);
  double yd=(std::max)(1.,(bbox.ymax()-bbox.ymin())/100);
  double zd=(std::max)(1.,(bbox.zmax()-bbox.zmin())/100);
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


/*!
  * \ingroup PMP_corefinement_grp
  * corefines `tm` and `splitter` and duplicates edges in `tm` that are on the intersection with `splitter`.
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(splitter)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *
  * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param splitter triangulated surface mesh used to split `tm`
  * @param np_tm optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  * @param np_s optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm` (`splitter`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *     \cgalParamDefault{If this parameter is omitted, an internal property map for
  *                       `CGAL::vertex_point_t` must be available in `TriangleMesh`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{if `true`, the set of triangles closed to the intersection of `tm`
  *                           and `splitter` will be checked for self-intersections
  *                           and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
*/
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
void split(TriangleMesh& tm,
           TriangleMesh& splitter,
           const NamedParameters1& np_tm,
           const NamedParameters2& np_s)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type VPM1;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::type VPM2;

  typedef typename boost::template property_map<TriangleMesh, CGAL::dynamic_edge_property_t<bool> >::type Ecm;

  VPM1 vpm_tm = choose_parameter(get_parameter(np_tm, internal_np::vertex_point),
                                 get_property_map(vertex_point, tm));
  VPM2 vpm_s = choose_parameter(get_parameter(np_s, internal_np::vertex_point),
                                get_property_map(vertex_point, splitter));

  Ecm ecm  = get(CGAL::dynamic_edge_property_t<bool>(), tm);

  // create a constrained edge map and corefine input mesh with the splitter,
  // and mark edges

  PMP::corefine(tm, splitter,
                CGAL::parameters::vertex_point_map(vpm_tm).edge_is_constrained_map(ecm),
                CGAL::parameters::vertex_point_map(vpm_s));

  //split mesh along marked edges
  internal::split_along_edges(tm, ecm, vpm_tm);
}

/**
  * \ingroup PMP_corefinement_grp
  * adds intersection edges of `plane` and `tm` in `tm` and duplicates those edges.
  *
  * \note In the current implementation it is not possible to set the vertex point map and the default will be used.
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param plane the plane that will be used to split `tm`.
  *              `Plane_3` is the plane type for the same CGAL kernel as the point of the vertex point map of `tm`.
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{if `true`, the set of triangles closed to the intersection of `tm`
  *                           and `plane` will be checked for self-intersections
  *                           and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  */
template <class TriangleMesh,
          class NamedParameters>
void split(TriangleMesh& tm,
#ifdef DOXYGEN_RUNNING
           const Plane_3& plane,
#else
           const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::Plane_3& plane,
#endif
           const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace params = PMP::parameters;

  // create a splitter mesh for the splitting plane using an internal CGAL function
  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm, np);
  double xd = (std::max)(1., 0.01 * (bbox.xmax() - bbox.xmin()));
  double yd = (std::max)(1., 0.01 * (bbox.ymax() - bbox.ymin()));
  double zd = (std::max)(1., 0.01 * (bbox.zmax() - bbox.zmin()));
  bbox = CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                      bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);

  TriangleMesh splitter;
  CGAL::Oriented_side os = PMP::internal::clip_to_bbox(plane, bbox, splitter, PMP::parameters::all_default());

  if(os == CGAL::ON_ORIENTED_BOUNDARY)
    return split(tm, splitter, np, params::all_default());

  //else nothing to do, no intersection.
}

/**
  * \ingroup PMP_corefinement_grp
  * clips `tm` by keeping the part that is inside `iso_cuboid`.
  * If `tm` is closed, the clipped part can be closed too if the named parameter `clip_volume` is set to `true`.
  * See Subsection \ref coref_clip for more details.
  *
  * \note In the current implementation it is not possible to set the vertex point map and the default will be used. `Iso_cuboid_3` must be
  * from the same %Kernel as the point of the vertex point map.
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh
  * @param iso_cuboid iso-cuboid used to clip `tm`.
  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{if `true`, the set of triangles closed to the intersection of `tm`
  *                           and `iso_cuboid` will be checked for self-intersections
  *                           and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{clip_volume}
  *     \cgalParamDescription{if `true` and `tm` is closed, the clipping will be done on
  *                           the volume \link coref_def_subsec bounded \endlink by `tm`
  *                           rather than on its surface (i.e., `tm` will be kept closed)}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{use_compact_clipper}
  *     \cgalParamDescription{if `false` the parts of `tm` coplanar with `plane` will not be part of the output}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold.
  *         If `false` is returned `tm` is only refined by the intersection with `iso_cuboid`.
  */
template <class TriangleMesh,
          class NamedParameters>
bool clip(TriangleMesh& tm,
#ifdef DOXYGEN_RUNNING
          const Iso_cuboid_3& iso_cuboid,
#else
          const typename GetGeomTraits<TriangleMesh, NamedParameters>::type::Iso_cuboid_3& iso_cuboid,
#endif
          const NamedParameters& np)
{
  if(boost::begin(faces(tm))==boost::end(faces(tm))) return true;
  TriangleMesh clipper;

  make_hexahedron(iso_cuboid[0], iso_cuboid[1], iso_cuboid[2], iso_cuboid[3],
                  iso_cuboid[4], iso_cuboid[5], iso_cuboid[6], iso_cuboid[7],
                  clipper);
  triangulate_faces(clipper);

  return clip(tm, clipper, np, parameters::all_default());
}

/// \cond SKIP_IN_MANUAL

// convenience overloads
template <class TriangleMesh>
bool clip(TriangleMesh& tm,
          const typename GetGeomTraits<TriangleMesh>::type::Plane_3& plane)
{
  return clip(tm, plane, parameters::all_default());
}

// convenience overloads
template <class TriangleMesh>
bool clip(TriangleMesh& tm,
          const typename GetGeomTraits<TriangleMesh>::type::Iso_cuboid_3& iso_cuboid)
{
  return clip(tm, iso_cuboid, parameters::all_default());
}

// convenience overload
template <class TriangleMesh,
          class NamedParameters1>
bool
clip(TriangleMesh& tm,
     TriangleMesh& clipper,
     const NamedParameters1& np_tm)
{
  return clip(tm, clipper, np_tm, parameters::all_default());
}

// convenience overload
template <class TriangleMesh>
bool
clip(TriangleMesh& tm,
     TriangleMesh& clipper)
{
  return clip(tm, clipper, parameters::all_default());
}


// convenience overload
template <class TriangleMesh,
          class NamedParameters1>
void
split(TriangleMesh& tm,
      TriangleMesh& splitter,
      const NamedParameters1& np_tm)
{
  split(tm, splitter, np_tm, parameters::all_default());
}

// convenience overload
template <class TriangleMesh>
void
split(TriangleMesh& tm,
      TriangleMesh& splitter)
{
  split(tm, splitter, parameters::all_default());
}

template <class TriangleMesh>
void split(TriangleMesh& tm,
           const typename GetGeomTraits<TriangleMesh>::type::Plane_3& plane)
{
   split(tm, plane, parameters::all_default());
}

/// \endcond

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_H
