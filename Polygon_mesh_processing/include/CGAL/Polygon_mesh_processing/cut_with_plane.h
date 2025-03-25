// Copyright (c) 2024 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_CUT_WITH_PLANE_H
#define CGAL_POLYGON_MESH_PROCESSING_CUT_WITH_PLANE_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#ifndef CGAL_PLANE_CLIP_DO_NOT_USE_TRIANGULATION
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#endif
#ifndef CGAL_PLANE_CLIP_DO_NOT_USE_BOX_INTERSECTION_D
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

template <class PolygonMesh>
struct Default_cut_visitor
{
  void before_edge_split(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor, PolygonMesh&) {}
  void edge_split(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor, PolygonMesh&) {}
  void vertices_on_cut(const std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>&, PolygonMesh&){}
};

// TODO: doc me or hide me in the np
template <class Kernel>
struct Orthogonal_cut_plane_traits
{
  using FT = typename Kernel::FT;
  using Plane_3 = std::pair<int, FT>;
  using Point_3 = typename Kernel::Point_3;

  struct Oriented_side_3
  {
    Oriented_side operator()(const Plane_3& plane, const Point_3& p)  const
    {
      if (p[plane.first]==plane.second) return ON_ORIENTED_BOUNDARY;
      return p[plane.first]<plane.second ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
    }
  };

  // TODO here we should integrate the epsilon to reuse existing points (moving them on the grid as an option?)
  struct Construct_plane_line_intersection_point_3
  {
    Point_3 operator()(const Plane_3& plane, const Point_3& p, const Point_3& q)
    {
      //TODO: divide by largest value?
      FT alpha = (plane.second - q[plane.first]) / (p[plane.first] - q[plane.first]);
      std::array<FT,3> coords;
      for (int i=0;i<3;++i)
      {
        if (i==plane.first)
          coords[i] = plane.second;
        else
          coords[i] = (p[i]==q[i])?p[i]:p[i]*alpha +(1-alpha)*q[i];
      }
      return Point_3(coords[0], coords[1], coords[2]);
    }
  };

  Oriented_side_3 oriented_side_3_object() const
  {
    return Oriented_side_3();
  }

  Construct_plane_line_intersection_point_3 construct_plane_line_intersection_point_3_object() const
  {
    return Construct_plane_line_intersection_point_3();
  }

#ifndef CGAL_PLANE_CLIP_DO_NOT_USE_BOX_INTERSECTION_D
// for does self-intersect
  using Segment_3 = typename Kernel::Segment_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Construct_segment_3 = typename Kernel::Construct_segment_3;
  using Construct_triangle_3 =typename  Kernel::Construct_triangle_3;
  using Do_intersect_3 = typename Kernel::Do_intersect_3;
  Construct_segment_3 construct_segment_3_object() const { return Construct_segment_3(); }
  Construct_triangle_3 construct_triangle_3_object() const { return Construct_triangle_3(); }
  Do_intersect_3 do_intersect_3_object() const { return Do_intersect_3(); }
#endif
};

/*!
 *  \ingroup PMP_corefinement_grp
 *
 *  refines `pm` by inserting the intersection points of `plane` with its edges and faces.
 *
 *  \tparam PolygonMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 *  \tparam Plane_3 plane type, equal to `GeomTraits::Plane_3`, `GeomTraits` being the type of the parameter `geom_traits`.
 *  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 *  \param pm input mesh to be refined
 *  \param plane the plane used to refine the mesh
 *  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 *  \cgalNamedParamsBegin
 *
 *    \cgalParamNBegin{concurrency_tag}
 *      \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *      \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *      \cgalParamDefault{`CGAL::Sequential_tag`}
 *    \cgalParamNEnd
 *
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `pm`.
 *                           If an edge marked as constrained is split, the two resulting edges will be marked as constrained.}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{edge_is_marked_map}
 *     \cgalParamDescription{a property map filled by this function that puts `true` for all intersection edge of faces
 *                           of `pm` and `plane`, and `false` for all other edges.}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_oriented_side_map}
 *     \cgalParamDescription{a property map filled by this function containing the position
 *                           of each vertex relative to the oriented plane `plane`.}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `Oriented_side` as value type}
 *     \cgalParamDefault{Dynamic vertex property map}
 *     \cgalParamExtra{If the concurrenty tag is set to `Parallel_tag`, the property map might be filled by several thread at the same time.}
 *   \cgalParamNEnd
 *
 *    \cgalParamNBegin{visitor}
 *      \cgalParamDescription{TODO add concept}
 *      \cgalParamType{reference wrapper recommeded if it has state TODO}
 *      \cgalParamDefault{None}
 *    \cgalParamNEnd
 *
 *    \cgalParamNBegin{do_not_triangulate_faces}
 *      \cgalParamDescription{If the input mesh is triangulated and this parameter is set to `false`, the mesh will be kept triangulated.}
 *      \cgalParamType{`bool`}
 *      \cgalParamDefault{`true`}
 *    \cgalParamNEnd
 *
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `pm`}
 *      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in `PolygonMesh`.}
 *    \cgalParamNEnd
 *
 *    \cgalParamNBegin{geom_traits}
 *      \cgalParamDescription{an instance of a geometric traits class}
 *      \cgalParamType{a class model of `Kernel`}
 *      \cgalParamDefault{a \cgal kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *    \cgalParamNEnd
 *
 *  \cgalNamedParamsEnd
 */
template <class PolygonMesh, class Plane_3, class NamedParameters =  parameters::Default_named_parameters>
void cut_with_plane(PolygonMesh& pm,
                    const Plane_3& plane,
                    const NamedParameters& np = parameters::default_values())
{
  // TODO: concurrency tag
  // TODO: if you want to clip with many planes (**Kernel**),
  //       it might be interesting to first classify all vertices with all planes
  //       to limit the number of intersection points computed: several classify done lazily
  //       on points not already eliminated (verices all out with adjacent vertices out too)
  // actually might be a global classifier to filter out edges, testing all planes at once per vertex and stop as soon as one is out
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::is_default_parameter;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  using face_descriptor = typename BGT::face_descriptor;
  using edge_descriptor = typename BGT::edge_descriptor;
  using halfedge_descriptor = typename BGT::halfedge_descriptor;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  // np typedefs
  using Default_ecm = Static_boolean_property_map<edge_descriptor, false>;
  using Default_visitor = Default_cut_visitor<PolygonMesh>;
  using Visitor_ref = typename internal_np::Lookup_named_param_def<internal_np::visitor_t, NamedParameters, Default_visitor>::reference;
  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  static_assert(std::is_same_v<typename GT::Plane_3,Plane_3>);

  auto ecm = choose_parameter<Default_ecm>(get_parameter(np, internal_np::edge_is_constrained));
  auto edge_is_marked = choose_parameter<Default_ecm>(get_parameter(np, internal_np::edge_is_marked_map));

  Default_visitor default_visitor;
  Visitor_ref visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_property_map(vertex_point, pm));

  bool triangulate = !choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), true);
  if (triangulate && !is_triangle_mesh(pm))
    triangulate = false;

  bool throw_on_self_intersection = choose_parameter(get_parameter(np, internal_np::use_compact_clipper), false);
  if (throw_on_self_intersection && !is_triangle_mesh(pm))
    throw_on_self_intersection = false;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  // constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

  auto oriented_side = traits.oriented_side_3_object();
  auto intersection_point = traits.construct_plane_line_intersection_point_3_object();

  // TODO: the default is not thread-safe for example for Polyhedron
  using V_os_tag = dynamic_vertex_property_t<Oriented_side>;
  static constexpr bool use_default_vosm =
    is_default_parameter<NamedParameters, internal_np::vertex_oriented_side_map_t>::value;

  using Vertex_oriented_side_map =
    std::conditional_t<use_default_vosm,
                       typename boost::property_map<PolygonMesh, V_os_tag>::type,
                       typename internal_np::Get_param<typename NamedParameters::base,
                                                       internal_np::vertex_oriented_side_map_t>::type>;


  Vertex_oriented_side_map vertex_os;
  if constexpr (use_default_vosm)
    vertex_os = get(V_os_tag(), pm);
  else
    vertex_os = get_parameter(np, internal_np::vertex_oriented_side_map);

  std::vector<edge_descriptor> inters;

  bool all_in = true;
  bool all_out = true;
  bool at_least_one_on = false;
  std::vector<vertex_descriptor> on_obnd;
  //TODO: parallel for
  for (vertex_descriptor v : vertices(pm))
  {
    Oriented_side os = oriented_side(plane,  get(vpm, v));
    put(vertex_os,v,os);
    switch(os)
    {
      case ON_POSITIVE_SIDE:
        all_in = false;
      break;
      case ON_NEGATIVE_SIDE:
        all_out = false;
      break;
      case ON_ORIENTED_BOUNDARY:
        at_least_one_on=true;
        on_obnd.push_back(v);
    }
  }

  if (at_least_one_on || (!all_in && !all_out))
  {
    //TODO: parallel for
    for(edge_descriptor e : edges(pm))
    {
      vertex_descriptor src = source(e,pm), tgt = target(e,pm);
      if (get(vertex_os, src)==CGAL::ON_ORIENTED_BOUNDARY)
      {
        if (get(vertex_os, tgt)==CGAL::ON_ORIENTED_BOUNDARY)
        {
          bool pure_coplanar=true;
          if (!is_border(e, pm))
          {
            halfedge_descriptor he=halfedge(e, pm);
            for (halfedge_descriptor h : halfedges_around_face(he,pm))
              if (get(vertex_os,target(h, pm))!=CGAL::ON_ORIENTED_BOUNDARY)
              {
                pure_coplanar=false;
                break;
              }
            if (pure_coplanar)
            {
              he=opposite(he, pm);
              for (halfedge_descriptor h : halfedges_around_face(he,pm))
                if (get(vertex_os, target(h, pm))!=CGAL::ON_ORIENTED_BOUNDARY)
                {
                  pure_coplanar=false;
                  break;
                }
            }
          }
          if (!pure_coplanar)
            put(edge_is_marked, e, true);
        }
      }
      else
        if (get(vertex_os, tgt)!=CGAL::ON_ORIENTED_BOUNDARY &&
            get(vertex_os, src)!=get(vertex_os, tgt))
          inters.push_back(e);
    }
  }

  if (all_in || all_out)
  {
    visitor.vertices_on_cut(on_obnd, pm);
    return;
  }

  std::unordered_map<face_descriptor, std::vector<halfedge_descriptor> > splitted_faces;

  if (throw_on_self_intersection)
  {
    std::vector<face_descriptor> test_faces;
    for (edge_descriptor e : inters)
    {
      halfedge_descriptor h = halfedge(e, pm);
      if (!is_border(h,pm))
        test_faces.push_back(face(h,pm));
      h=opposite(h, pm);
      if (!is_border(h,pm))
        test_faces.push_back(face(h,pm));
    }
    std::sort(test_faces.begin(), test_faces.end());
    auto last = std::unique(test_faces.begin(), test_faces.end());
    test_faces.erase(last, test_faces.end());
    if (does_self_intersect<Concurrency_tag>(test_faces, pm, np))
      throw std::runtime_error("TODO Corefinement::Self_intersection_exception");
  }

  //TODO: parallel for
  for (edge_descriptor e : inters)
  {
    halfedge_descriptor h = halfedge(e, pm);
    auto pts = CGAL::make_sorted_pair(get(vpm, source(h,pm)), get(vpm, target(h, pm)));
    typename GT::Point_3 ip = intersection_point(plane, pts.first, pts.second);

    bool was_marked = get(ecm, edge(h, pm));
    visitor.before_edge_split(h, pm);
    h = CGAL::Euler::split_edge(h, pm);
    put(vpm, target(h, pm), ip);
    put(vertex_os, target(h, pm), ON_ORIENTED_BOUNDARY);
    visitor.edge_split(h, pm);
    if (was_marked)
      put(ecm, edge(h, pm), true);

    if (!is_border(h, pm))
      splitted_faces[face(h, pm)].push_back(h);
    h=prev(opposite(h,pm),pm);
    if (!is_border(h, pm))
      splitted_faces[face(h, pm)].push_back(h);
  }

  visitor.vertices_on_cut(on_obnd, pm);

  // collect faces to be cut that have one vertex on the cut plane
  for (vertex_descriptor v : on_obnd)
  {
    halfedge_descriptor hv = halfedge(v, pm);
    for (halfedge_descriptor h : halfedges_around_target(hv, pm))
    {
      if (is_border(h, pm)) continue;
      Oriented_side prev_ori = get(vertex_os, source(h, pm)),
                    next_ori = get(vertex_os, target(next(h, pm), pm));
      if ( prev_ori == ON_ORIENTED_BOUNDARY || next_ori == ON_ORIENTED_BOUNDARY) continue; // skip full edge
      if (prev_ori!=next_ori) splitted_faces[face(h, pm)].push_back(h); // skip tangency point
    }
  }

  //TODO: parallel for
  for (std::pair<const face_descriptor, std::vector<halfedge_descriptor>>& f_and_hs : splitted_faces)
  {
    std::size_t nb_hedges = f_and_hs.second.size();

    CGAL_assertion( nb_hedges%2 ==0 );

    if (nb_hedges==2)
    {
      halfedge_descriptor h1=f_and_hs.second[0], h2=f_and_hs.second[1];
      CGAL_assertion(next(h1,pm)!=h2 && next(h2,pm)!=h1); // the edge does not already exist
      halfedge_descriptor res = CGAL::Euler::split_face(h1, h2, pm);
      put(edge_is_marked, edge(res, pm), true);

      if (triangulate)
      {
        if (!is_triangle(res, pm))
        {
          // TODO: take the criteria in triangulate_faces ?
          halfedge_descriptor newh =
            CGAL::Euler::split_face(res, next(next(res, pm), pm), pm);
          put(edge_is_marked, edge(newh, pm), false);
        }
        else
        {
          res = opposite(res, pm);
          if (!is_triangle(res, pm))
          {
            // TODO: take the criteria in triangulate_faces ?
            halfedge_descriptor newh =
              CGAL::Euler::split_face(res, next(next(res, pm), pm), pm);
            put(edge_is_marked, edge(newh, pm), false);
          }
        }
      }
    }
    else
    {
      // sort hedges to make them match
      CGAL_assertion(!triangulate);
      // TODO: need mechanism to make it robust even with EPICK
      auto less_hedge = [&pm, vpm](halfedge_descriptor h1, halfedge_descriptor h2)
      {
        return lexicographically_xyz_smaller(get(vpm,target(h1,pm)), get(vpm,target(h2,pm)));
      };
      std::sort(f_and_hs.second.begin(), f_and_hs.second.end(), less_hedge);

      for (std::size_t i=0; i<nb_hedges; i+=2)
      {
        halfedge_descriptor h1=f_and_hs.second[i], h2=f_and_hs.second[i+1];
        CGAL_assertion(next(h1,pm)!=h2 && next(h2,pm)!=h1); // the edge does not already exist
        halfedge_descriptor res = CGAL::Euler::split_face(h1, h2, pm);
        put(edge_is_marked, edge(res, pm), true);
      }
    }
  }
}

} } // CGAL::Polygon_mesh_processing


#endif // CGAL_POLYGON_MESH_PROCESSING_CUT_WITH_PLANE_H
