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

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/property_map/property_map.hpp>

namespace CGAL {
namespace Polygon_mesh_processing {

template <class PolygonMesh, class NamedParameters =  parameters::Default_named_parameters>
void clip_convex(PolygonMesh& pm,
                 const typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Plane_3& plane,
                 const NamedParameters& np = parameters::default_values())
{
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
  // using Default_ecm = Static_boolean_property_map<edge_descriptor, false>;
  // using Default_visitor = Default_cut_visitor<PolygonMesh>;
  // using Visitor_ref = typename internal_np::Lookup_named_param_def<internal_np::visitor_t, NamedParameters, Default_visitor>::reference;
  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using FT = typename GT::FT;
  // using Point_3 = typename GT::Point_3;

  struct Default_Bbox{
    // Needed to compile with MSVC
    vertex_descriptor operator[](std::size_t /*i*/){ return vertex_descriptor(); }
    Default_Bbox operator*(){ return *this; }
  };
  using Bbox = typename internal_np::Lookup_named_param_def<internal_np::bounding_box_t, NamedParameters, Default_Bbox*>::type;
  constexpr bool update_bbox = !std::is_same_v< std::remove_reference_t<Bbox>, Default_Bbox*>;

  Default_Bbox* default_bbox;
  Bbox bbox_pointer = choose_parameter(get_parameter_reference(np, internal_np::bounding_box), default_bbox);

  // Default_visitor default_visitor;
  // Visitor_ref visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);
  // constexpr bool has_visitor = !std::is_same_v<Default_visitor, std::remove_cv_t<std::remove_reference_t<Visitor_ref>>>;

  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_property_map(vertex_point, pm));

  bool triangulate = !choose_parameter(get_parameter(np, internal_np::do_not_triangulate_faces), true);
  if (triangulate && !is_triangle_mesh(pm))
    triangulate = false;

  bool throw_on_self_intersection = choose_parameter(get_parameter(np, internal_np::use_compact_clipper), false);
  if (throw_on_self_intersection && !is_triangle_mesh(pm))
    throw_on_self_intersection = false;

  auto oriented_side = traits.oriented_side_3_object();
  auto intersection_point = traits.construct_plane_line_intersection_point_3_object();
  auto sq = traits.compute_squared_distance_3_object();
  // auto csq = traits.compare_squared_distance_3_object();
  // auto vector_3 = traits.construct_vector_3_object();

  // ____________________ Find a crossing edge _____________________

  vertex_descriptor src=*vertices(pm).begin();
  FT sp_src = sq(plane, get(vpm, src)); // Not normalized distance
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
        // TODO with EPICK, use compare_distance(plane, src, plane, trg) (But no possibility to memorize some computations for the next)
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
        return;
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
    // Nothing to clip
    if(no_positive_side)
      return;
  } else if(direction_to_zero==NEGATIVE){
    // Orient the edge from negative to positive
    std::swap(src, trg);
  }

  CGAL_assertion(oriented_side(plane, get(vpm, src)) == ON_POSITIVE_SIDE);
  CGAL_assertion(oriented_side(plane, get(vpm, trg)) != ON_POSITIVE_SIDE);

  // Cut the convex along the plane by marching along crossing edges starting from the previous edge
  std::vector<halfedge_descriptor> boundaries;
  std::vector<vertex_descriptor> boundary_vertices;

  halfedge_descriptor h = halfedge(src, trg, pm).first;
  if(sign(sp_trg)!=EQUAL)
  {
    //split the first edge
    auto pts = make_sorted_pair(get(vpm, src), get(vpm, trg));
    typename GT::Point_3 ip = intersection_point(plane, pts.first, pts.second);
    //visitor.before_edge_split(h, pm);
    h = CGAL::Euler::split_edge(h, pm);
    put(vpm, target(h, pm), ip);
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
        boundary_vertices.emplace_back(target(h, pm));
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
      // visitor.before_edge_split(h, pm);
      h = CGAL::Euler::split_edge(h, pm);
      put(vpm, target(h, pm), ip);
      // visitor.new_vertex_added(vid, target(h,pm), pm);
      // visitor.edge_split(h, pm);
      // visitor.after_edge_split();

    }

    // Split the face
    halfedge_descriptor sh = CGAL::Euler::split_face(h_previous, h, pm);
    boundaries.emplace_back(sh);
    boundary_vertices.emplace_back(target(sh, pm));
    set_halfedge(target(sh, pm), sh, pm);

    CGAL_assertion(target(sh, pm) == target(h, pm));
    h = opposite(next(sh,pm), pm);

    // During the loop, if a mesh vertex lies on the plane, we look for a face incident to that vertex that crosses the plane.
    // The second part of the while-condition ensures we don't exit prematurely in this looking.
  } while(target(h, pm)!=v_start || (boundaries.empty() && h!=h_start));

  CGAL_assertion(is_valid_polygon_mesh(pm));
  CGAL_assertion(!boundaries.empty());

  // ________________________ Remove the negative side _________________________
  std::vector<vertex_descriptor> vertices_to_remove;
  std::vector<edge_descriptor> edges_to_remove;
  std::vector<face_descriptor> faces_to_remove;

  std::vector<halfedge_descriptor> boundary_edges;
  for(halfedge_descriptor h: boundaries)
    boundary_edges.push_back(h);
  std::sort(boundary_vertices.begin(), boundary_vertices.end());
  std::sort(boundary_edges.begin(), boundary_edges.end());
  auto is_boundary=[&](edge_descriptor e){
    halfedge_descriptor h = halfedge(e, pm);
    return std::binary_search(boundary_edges.begin(), boundary_edges.end(), h) ||
           std::binary_search(boundary_edges.begin(), boundary_edges.end(), opposite(h, pm));
  };
  connected_component(face(halfedge(src, pm), pm), pm, std::back_inserter(faces_to_remove),
                      parameters::edge_is_constrained_map(boost::make_function_property_map<edge_descriptor>(is_boundary)));
  for(face_descriptor f: faces_to_remove){
    halfedge_descriptor h_start = halfedge(f, pm);
    halfedge_descriptor h = h_start;
    do {
      if(!is_boundary(edge(h, pm))){ // TODO avoid this linear operation
        if(h < opposite(h, pm)) // To avoid multiple assertions of a same edge
          edges_to_remove.push_back(edge(h, pm));
        if(halfedge(target(h, pm), pm) == h && // To avoid multiple assertions of a same vertex
           !std::binary_search(boundary_vertices.begin(), boundary_vertices.end(), target(h, pm)))
          vertices_to_remove.push_back(target(h, pm));
        CGAL_assertion(oriented_side(plane, get(vpm, target(h, pm)))!=ON_NEGATIVE_SIDE);
      }
      h = next(h, pm);
    } while (h != h_start);
  }

  if constexpr(update_bbox){
    auto &bbox = *bbox_pointer;
    for (vertex_descriptor v : vertices_to_remove){
      // TODO find a way to factorize

      // Update xmin
      if(v==bbox[0]){
        bbox[0] = *(boundary_vertices.begin());
        double xmin = to_interval(get(vpm, bbox[0]).x()).first;
        for(auto it=++boundary_vertices.begin(); it!=boundary_vertices.end(); ++it){
          vertex_descriptor v = *it;
          double x = to_interval(get(vpm, v).x()).first;
          if(x < xmin){
            xmin = x;
            bbox[0] = v;
          }
        }
      }

      // Update xmax
      if(v==bbox[1]){
        bbox[1] = *(boundary_vertices.begin());
        double xmax = to_interval(get(vpm, bbox[1]).x()).second;
        for(auto it=++boundary_vertices.begin(); it!=boundary_vertices.end(); ++it){
          vertex_descriptor v = *it;
          double x = to_interval(get(vpm, v).x()).second;
          if(x > xmax){
            xmax = x;
            bbox[1] = v;
          }
        }
      }

      // Update ymin
      if(v==bbox[2]){
        bbox[2] = *(boundary_vertices.begin());
        double ymin = to_interval(get(vpm, bbox[2]).y()).first;
        for(auto it=++boundary_vertices.begin(); it!=boundary_vertices.end(); ++it){
          vertex_descriptor v = *it;
          double y = to_interval(get(vpm, v).y()).first;
          if(y < ymin){
            ymin = y;
            bbox[2] = v;
          }
        }
      }

      // Update ymax
      if(v==bbox[3]){
        bbox[3] = *(boundary_vertices.begin());
        double ymax = to_interval(get(vpm, bbox[3]).y()).second;
        for(auto it=++boundary_vertices.begin(); it!=boundary_vertices.end(); ++it){
          vertex_descriptor v = *it;
          double y = to_interval(get(vpm, v).y()).second;
          if(y > ymax){
            ymax = y;
            bbox[3] = v;
          }
        }
      }

      // Update zmin
      if(v==bbox[4]){
        bbox[4] = *(boundary_vertices.begin());
        double zmin = to_interval(get(vpm, bbox[4]).z()).first;
        for(auto it=++boundary_vertices.begin(); it!=boundary_vertices.end(); ++it){
          vertex_descriptor v = *it;
          double z = to_interval(get(vpm, v).z()).first;
          if(z < zmin){
            zmin = z;
            bbox[4] = v;
          }
        }
      }

      // Update zmax
      if(v==bbox[5]){
        bbox[5] = *(boundary_vertices.begin());
        double zmax = to_interval(get(vpm, bbox[5]).z()).second;
        for(auto it=++boundary_vertices.begin(); it!=boundary_vertices.end(); ++it){
          vertex_descriptor v = *it;
          double z = to_interval(get(vpm, v).z()).second;
          if(z > zmax){
            zmax = z;
            bbox[5] = v;
          }
        }
      }
    }
  }

  for (vertex_descriptor v : vertices_to_remove){
    remove_vertex(v, pm);
    CGAL_assertion(oriented_side(plane, get(vpm, v))==ON_POSITIVE_SIDE);
  }
  for (edge_descriptor e : edges_to_remove)
    remove_edge(e, pm);
  for (face_descriptor f : faces_to_remove)
    remove_face(f, pm);

  // Reorder halfedges of the hole
  for(size_t i=1; i<boundaries.size(); ++i)
    set_next(boundaries[i-1], boundaries[i], pm);
  set_next(boundaries.back(), boundaries[0], pm);

  // Fill the hole
  face_descriptor f=pm.add_face();
  for(auto h: boundaries){
    set_face(h, f, pm);
  }
  set_halfedge(f, boundaries[0], pm);

  if( pm.number_of_vertices() < 3) //Degenerate to a segment
    clear(pm);
  CGAL_assertion(is_valid_polygon_mesh(pm));
  // std::ofstream("clip.off") << pm;
}

} } // CGAL::Polygon_mesh_processing


#endif // CGAL_POLYGON_MESH_PROCESSING_CLIP_CONVEX_H
