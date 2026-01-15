// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque


#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_CONVEX_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_CONVEX_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <boost/property_map/function_property_map.hpp>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal{

/**
  * Given a convex mesh and a plane, return an halfedge crossing the plane from a vertex on positive side to a vertex on the plane or the negative side.
  * Return a null halfedge if the mesh and the plane does not intersect.
  */
template <class PolygonMesh, class NamedParameters =  parameters::Default_named_parameters>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
find_crossing_edge(PolygonMesh& pm,
                   const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Plane_3& plane,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  // np typedefs
  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_property_map(vertex_point, pm));

  auto oriented_side = traits.oriented_side_3_object();
  auto sq = traits.compute_squared_distance_3_object();

  // ____________________ Find a crossing edge _____________________

  vertex_descriptor start = choose_parameter(get_parameter(np, internal_np::starting_vertex_descriptor), *vertices(pm).begin());
  vertex_descriptor src = start;
  FT sp_src = sq(plane, get(vpm, src));
  Sign direction_to_zero = sign(sp_src);

  vertex_descriptor trg;
  FT sp_trg;

  bool is_crossing_edge=false;
  if(direction_to_zero!=EQUAL){
    do{
      bool is_local_max=true;
      for(auto v: vertices_around_target(src ,pm)){
        sp_trg = sq(plane, get(vpm, v));
        CGAL_assertion(sq(plane, get(vpm, v)) == sp_trg);
        // Check if v in the direction to the plane
        if(compare(sp_src, sp_trg)==direction_to_zero){
          if(sign(sp_trg)!=direction_to_zero){
            // Fund a crossing edge
            trg = v;
            is_crossing_edge=true;
          } else {
            // Continue from v
            sp_src = sp_trg;
            src = v;
          }
          is_local_max=false; // repeat with the new vertex
          break;
        }
      }
      // No intersection with the plane, kernel is either empty or full
      if(is_local_max){
        if(direction_to_zero==POSITIVE)
          clear(pm); // The result is empty
        return BGT::null_halfedge();
      }
    } while(!is_crossing_edge);
  } else {
    // src on the plane
    trg=src;
    sp_trg = sp_src;
  }

  if(sign(sp_trg)==EQUAL && direction_to_zero!=POSITIVE){
    // Search a vertex around trg coming from positive side
    bool no_positive_side = true;
    for(auto v: vertices_around_target(trg ,pm)){
      Oriented_side side_v = oriented_side(plane, get(vpm, v));
      if(side_v==ON_POSITIVE_SIDE){
        src = v;
        sp_src = sq(plane, get(vpm, src));
        no_positive_side = false;
        break;
      }
    }
    // Nothing on negative side
    if(no_positive_side)
      return BGT::null_halfedge();
  } else if(direction_to_zero==NEGATIVE){
    // Orient the edge from positive to negative
    std::swap(src, trg);
  }

  CGAL_assertion(oriented_side(plane, get(vpm, src)) == ON_POSITIVE_SIDE);
  CGAL_assertion(oriented_side(plane, get(vpm, trg)) != ON_POSITIVE_SIDE);

  return halfedge(src, trg, pm).first;
}

/**
  * Given a convex mesh, a plane and an halfedge crossing the plane from positive side, refine the mesh with the plane
  */
template <class PolygonMesh, class NamedParameters =  parameters::Default_named_parameters>
std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>
refine_convex_with_plane(PolygonMesh& pm,
            const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Plane_3& plane,
            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
            const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  using halfedge_descriptor = typename BGT::halfedge_descriptor;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  // np typedefs
  // using Default_ecm = Static_boolean_property_map<edge_descriptor, false>;
  using Default_visitor = Corefinement::Default_visitor<PolygonMesh>;
  using Visitor_ref = typename internal_np::Lookup_named_param_def<internal_np::visitor_t, NamedParameters, Default_visitor>::reference;
  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  struct Default_Bbox{
    // Needed to compile with MSVC
    vertex_descriptor operator[](std::size_t /*i*/){ return vertex_descriptor(); }
    Default_Bbox operator*(){ return *this; }
  };
  using Bbox = typename internal_np::Lookup_named_param_def<internal_np::bounding_box_t, NamedParameters, Default_Bbox*>::type;
  constexpr bool update_bbox = !std::is_same_v< std::remove_reference_t<Bbox>, Default_Bbox*>;

  Default_Bbox* default_bbox;
  Bbox bbox_pointer = choose_parameter(get_parameter_reference(np, internal_np::bounding_box), default_bbox);

  Default_visitor default_visitor;
  Visitor_ref visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);
  // constexpr bool has_visitor = !std::is_same_v<Default_visitor, std::remove_cv_t<std::remove_reference_t<Visitor_ref>>>;

  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_property_map(vertex_point, pm));

  // config flags
  bool triangulate = !choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), true);

  auto oriented_side = traits.oriented_side_3_object();
  auto intersection_point = traits.construct_plane_line_intersection_point_3_object();

  // ____ Cut the convex along the plane by marching along crossing edges starting from the previous edge _____
  std::vector<halfedge_descriptor> boundaries;

  if(oriented_side(plane, get(vpm, target(h,pm)))!=ON_ORIENTED_BOUNDARY)
  {
    //split the first edge
    auto pts = make_sorted_pair(get(vpm, source(h, pm)), get(vpm, target(h, pm)));
    typename GT::Point_3 ip = intersection_point(plane, pts.first, pts.second);
    visitor.before_edge_split(h, pm);
    h = CGAL::Euler::split_edge(h, pm);
    put(vpm, target(h, pm), ip);
    // visitor.new_vertex_added(vpm.size()-1, target(h,pm), pm);
    visitor.edge_split(h, pm);
    visitor.after_edge_split();
  }

  vertex_descriptor v_start = target(h, pm);
  halfedge_descriptor h_start=h;
  do{
    halfedge_descriptor h_previous = h;
    CGAL_assertion(oriented_side(plane, get(vpm, source(h,pm)))==ON_POSITIVE_SIDE);
    CGAL_assertion(oriented_side(plane, get(vpm, target(h,pm)))==ON_ORIENTED_BOUNDARY);

    h = next(h, pm);
    Oriented_side side_trg = oriented_side(plane, get(vpm, target(h,pm)));

    if(side_trg != ON_NEGATIVE_SIDE){
      // The face does not cross the plane, we look for a face incident to the same vertex than h that crosses the plane.
      while(side_trg == ON_ORIENTED_BOUNDARY){
        // The edge is along the plane, add it to boundaries
        boundaries.emplace_back(h);
        set_halfedge(target(h, pm), h, pm);
        h = next(h, pm);
        side_trg=oriented_side(plane, get(vpm, target(h,pm)));
      }
      // continue on next face
      h = opposite(h, pm);
      continue;
    }

    // Search a crossing edge
    h = next(h, pm);
    side_trg=oriented_side(plane, get(vpm, target(h,pm)));
    while(side_trg == ON_NEGATIVE_SIDE){
      h = next(h,pm);
      side_trg=oriented_side(plane, get(vpm, target(h,pm)));
    }

    if(side_trg != ON_ORIENTED_BOUNDARY){
      // Split the edge
      auto pts = make_sorted_pair(get(vpm, source(h,pm)), get(vpm, target(h,pm)));
      typename GT::Point_3 ip = intersection_point(plane, pts.first, pts.second);
      visitor.before_edge_split(h, pm);
      h = CGAL::Euler::split_edge(h, pm);
      put(vpm, target(h, pm), ip);
      // visitor.new_vertex_added(vpm.size()-1, target(h,pm), pm);
      visitor.edge_split(h, pm);
      visitor.after_edge_split();
    }

    // Split the face
    visitor.before_subface_created(pm);
    halfedge_descriptor sh = CGAL::Euler::split_face(h_previous, h, pm);
    visitor.after_subface_created(face(h, pm), pm);
    boundaries.emplace_back(sh);
    set_halfedge(target(sh, pm), sh, pm);
    visitor.add_retriangulation_edge(sh, pm);

    if (triangulate){
      halfedge_descriptor sh_opp = opposite(sh, pm);
      if(!is_triangle(sh_opp, pm)){
        visitor.before_subface_created(pm);
        halfedge_descriptor newh = CGAL::Euler::split_face(sh_opp, next(next(sh_opp, pm), pm), pm);
        visitor.after_subface_created(face(opposite(newh, pm), pm), pm);
        visitor.add_retriangulation_edge(newh, pm);
      }
    }

    CGAL_assertion(target(sh, pm) == target(h, pm));
    h = opposite(next(sh,pm), pm);

    // During the loop, if a mesh vertex lies on the plane, we look for a face incident to that vertex that crosses the plane.
    // The second part of the while-condition ensures we don't exit prematurely
  } while(target(h, pm)!=v_start || (boundaries.empty() && h!=h_start));

  CGAL_assertion(is_valid_polygon_mesh(pm));
  CGAL_assertion(!boundaries.empty());

  return boundaries;
}

/**
  * Given a mesh, a cycle of halfedges forming a boundary and a vertex, remove all connected elements to the vertex without crossing the boundary
  */
template <class PolygonMesh, class NamedParameters =  parameters::Default_named_parameters>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
remove_bounded_region_and_fill(PolygonMesh& pm,
                               const std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor> &boundaries,
                               typename boost::graph_traits<PolygonMesh>::vertex_descriptor inside_vertex,
                               const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  using face_descriptor = typename BGT::face_descriptor;
  using edge_descriptor = typename BGT::edge_descriptor;
  using halfedge_descriptor = typename BGT::halfedge_descriptor;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  // np typedefs
  // using Default_ecm = Static_boolean_property_map<edge_descriptor, false>;
  using Default_visitor = Corefinement::Default_visitor<PolygonMesh>;
  using Visitor_ref = typename internal_np::Lookup_named_param_def<internal_np::visitor_t, NamedParameters, Default_visitor>::reference;
  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;

  struct Default_Bbox{
    // Needed to compile with MSVC
    vertex_descriptor operator[](std::size_t /*i*/){ return vertex_descriptor(); }
    Default_Bbox operator*(){ return *this; }
  };
  using Bbox = typename internal_np::Lookup_named_param_def<internal_np::bounding_box_t, NamedParameters, Default_Bbox*>::type;
  constexpr bool update_bbox = !std::is_same_v< std::remove_reference_t<Bbox>, Default_Bbox*>;

  Default_Bbox* default_bbox;
  Bbox bbox_pointer = choose_parameter(get_parameter_reference(np, internal_np::bounding_box), default_bbox);

  Default_visitor default_visitor;
  Visitor_ref visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);
  // constexpr bool has_visitor = !std::is_same_v<Default_visitor, std::remove_cv_t<std::remove_reference_t<Visitor_ref>>>;

  // Used only if do_triangulate_faces or bounding_box check
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_property_map(vertex_point, pm));

  // config flags
  bool clip_volume =  choose_parameter(get_parameter(np, internal_np::clip_volume), true);
  bool triangulate = !choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), true);

  // ________________________ Remove the negative side _________________________
  std::vector<vertex_descriptor> vertices_to_remove;
  std::vector<edge_descriptor> edges_to_remove;
  std::vector<face_descriptor> faces_to_remove;

  std::vector<halfedge_descriptor> boundary_edges;
  std::vector<vertex_descriptor> boundary_vertices;
  boundary_edges.reserve(boundaries.size());
  boundary_vertices.reserve(boundaries.size());
  for(halfedge_descriptor h: boundaries){
    boundary_edges.push_back(h);
    boundary_vertices.push_back(target(h, pm));
  }
  std::sort(boundary_vertices.begin(), boundary_vertices.end());
  std::sort(boundary_edges.begin(), boundary_edges.end());
  auto is_boundary=[&](edge_descriptor e){
    halfedge_descriptor h = halfedge(e, pm);
    return std::binary_search(boundary_edges.begin(), boundary_edges.end(), h) ||
           std::binary_search(boundary_edges.begin(), boundary_edges.end(), opposite(h, pm));
  };
  connected_component(face(halfedge(inside_vertex, pm), pm), pm, std::back_inserter(faces_to_remove),
                      parameters::edge_is_constrained_map(boost::make_function_property_map<edge_descriptor>(is_boundary)));
  for(face_descriptor f: faces_to_remove){
    halfedge_descriptor h_start = halfedge(f, pm);
    halfedge_descriptor h = h_start;
    do {
      if(!is_boundary(edge(h, pm))){
        if(h < opposite(h, pm)) // To avoid multiple assertions of a same edge
          edges_to_remove.push_back(edge(h, pm));
        if(halfedge(target(h, pm), pm) == h && // To avoid multiple assertions of a same vertex
           !std::binary_search(boundary_vertices.begin(), boundary_vertices.end(), target(h, pm)))
          vertices_to_remove.push_back(target(h, pm));
      }
      h = next(h, pm);
    } while (h != h_start);
  }

  if constexpr(update_bbox){
    auto &bbox = *bbox_pointer;
    struct BBoxEntry {
      std::size_t index;
      std::function<double(const Point_3&)> bound;
      bool take_min;
    };

    std::array<BBoxEntry,6> entries {{
        {0, [](const Point_3& p){ return to_interval(p.x()).first;  }, true},  // xmin
        {1, [](const Point_3& p){ return to_interval(p.x()).second; }, false}, // xmax
        {2, [](const Point_3& p){ return to_interval(p.y()).first;  }, true},  // ymin
        {3, [](const Point_3& p){ return to_interval(p.y()).second; }, false}, // ymax
        {4, [](const Point_3& p){ return to_interval(p.z()).first;  }, true},  // zmin
        {5, [](const Point_3& p){ return to_interval(p.z()).second; }, false}  // zmax
    }};

    for (vertex_descriptor v : vertices_to_remove){
      for (const auto& e : entries){
        // Only update this bbox entry if the removed vertex *was* this extremum
        std::size_t i = e.index;
        if (v == bbox[i]){
          auto it = boundary_vertices.begin();
          bbox[i] = *it;
          double bound = e.bound(get(vpm, bbox[i]));
          for (++it; it != boundary_vertices.end(); ++it){
            vertex_descriptor v = *it;
            double val = e.bound(get(vpm, v));
            if ((e.take_min && val < bound) ||
              (!e.take_min && val > bound)){
              bound = val;
              bbox[i] = v;
            }
          }
        }
      }
    }
  }

  for (vertex_descriptor v : vertices_to_remove)
    remove_vertex(v, pm);
  for (edge_descriptor e : edges_to_remove)
    remove_edge(e, pm);
  for (face_descriptor f : faces_to_remove)
    remove_face(f, pm);

  // Reorder halfedges of the hole
  for(size_t i=1; i<boundaries.size(); ++i)
    set_next(boundaries[i-1], boundaries[i], pm);
  set_next(boundaries.back(), boundaries[0], pm);

  // Fill the hole
  if (clip_volume){
    face_descriptor f=pm.add_face();
    for(auto h: boundaries){
      set_face(h, f, pm);
    }
    set_halfedge(f, boundaries[0], pm);
    if(triangulate)
      triangulate_face(f, pm, parameters::vertex_point_map(vpm));
  } else {
    for(auto h: boundaries){
      set_face(h, BGT::null_face(), pm);
    }
  }

  if( pm.number_of_vertices() < 3) //Degenerate to a segment
    clear(pm);
  CGAL_assertion(is_valid_polygon_mesh(pm));
  // std::ofstream("clip.off") << pm;

  return *boundary_vertices.begin();
}

  /**
  * \ingroup PMP_corefinement_grp
  *
  * \brief clips `pm` by keeping the part that is on the negative side of `plane` (the side opposite to its normal vector).
  * The input mesh must be convex to guarantee a correct execution and results.
  *
  * By default, the clipped part is kept closed. Set the named parameter `clip_volume` to `false` to disable this.
  *
  * @tparam PolygonMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pm input surface mesh
  * @param plane plane whose negative side defines the halfspace to intersect `pm` with.
  *              `Plane_3` is the plane type for the same \cgal kernel as the point of the vertex point map of `pm`.
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces, edges, and faces.
  *                           Note that as there are no mesh associated with `plane`,
  *                           `boost::graph_traits<PolygonMesh>::null_halfedge()` and `boost::graph_traits<PolygonMesh>::null_face()` will be used when calling
  *                           functions of the visitor expecting a halfedge or a face from `plane`. Similarly, `pm` will be used as the mesh of `plane`.}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<PolygonMesh>`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{clip_volume}
  *     \cgalParamDescription{If `true`, the clipping will be done on
  *                           the volume \link coref_def_subsec bounded \endlink by `pm`
  *                           rather than on its surface (i.e., `pm` will remain closed).}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  *
  *    \cgalParamNBegin{do_not_triangulate_faces}
  *      \cgalParamDescription{If this parameter is set to `false`, the face added or modified by the algorithm will be triangulated.}
  *      \cgalParamType{Boolean}
  *      \cgalParamDefault{`true`}
  *    \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @return `true`
  *
  * @see `clip()`
  */
template <class PolygonMesh, class NamedParameters =  parameters::Default_named_parameters>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
clip_convex(PolygonMesh& pm,
#ifdef DOXYGEN_RUNNING
            const Plane_3& plane,
#else
            const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Plane_3& plane,
#endif
            const NamedParameters& np = parameters::default_values()){
  auto he = find_crossing_edge(pm, plane, np);
  // Early exit
  if(he == boost::graph_traits<PolygonMesh>::null_halfedge())
    return parameters::choose_parameter(parameters::get_parameter(np, internal_np::starting_vertex_descriptor), *vertices(pm).begin());
  const auto &boundaries =refine_convex_with_plane(pm, plane, he, np);
  return remove_bounded_region_and_fill(pm, boundaries, source(he, pm), np);
}

} // end of namespace internal


} } // CGAL::Polygon_mesh_processing


#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_CONVEX_H
