// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Leo Valque

#ifndef CGAL_SURFACE_MESH_BOOLEAN_OPERATIONS_UTILS_H
#define CGAL_SURFACE_MESH_BOOLEAN_OPERATIONS_UTILS_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace Corefinement {

template<class Point>
struct is_surface_mesh<Surface_mesh<Point>> : std::true_type {};

// Specialization of the function for Surface_mesh to exploit vectors structure of Surface_mesh.
template < bool reverse_patch_orientation,
           class Point,
           class PatchDescription,
           class VertexPointMap,
           class VertexPointMapOut,
           class EdgeMarkMapOut,
           class EdgeMarkMapIn,
           class EdgeToEdgeMap,
           class VertexToVertexMap,
           class UserVisitor>
void append_patch(
  const Surface_mesh<Point>& tm,
  Surface_mesh<Point>& output,
  PatchDescription& patch,
  const VertexPointMapOut& vpm_out,
  const VertexPointMap& vpm_tm,
  EdgeMarkMapOut& edge_mark_map_out,
  const EdgeMarkMapIn& edge_mark_map_in,
  EdgeToEdgeMap& tm_to_output_edges,
  VertexToVertexMap& tm_to_output_vertices,
  UserVisitor& user_visitor,
  std::size_t vertices_idx_begin,
  std::size_t edges_idx_begin,
  std::size_t faces_idx_begin)
{
  using SM = Surface_mesh<Point>;
  using vertex_descriptor = typename SM::Vertex_index;
  using edge_descriptor = typename SM::Edge_index;
  using halfedge_descriptor = typename SM::Halfedge_index;
  using face_descriptor = typename SM::Face_index;

  // Fill Vertex to vertex map and vpm for interior vertices
  auto fill_vertex = [&](std::size_t i){
    vertex_descriptor v = patch.interior_vertices[i];
    vertex_descriptor new_v(vertices_idx_begin + i);

    put(tm_to_output_vertices, v, new_v);
    put(vpm_out, new_v, get(vpm_tm, v));
    output.set_halfedge(new_v, SM::null_halfedge());
  };

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(true){
    tbb::parallel_for(std::size_t(0), patch.interior_vertices.size(), fill_vertex);
  }
  else
#endif
  {
    for(std::size_t i=0; i<patch.interior_vertices.size(); ++i) fill_vertex(i);
  }

  // Fill Edge to edge map for interior edges
  // Also fill target of halfedges and halfedge of vertices
  auto fill_edge = [&](std::size_t i){
    halfedge_descriptor h = patch.interior_edges[i];
    edge_descriptor new_edge(i+edges_idx_begin);
    halfedge_descriptor new_h = output.halfedge(new_edge);

    vertex_descriptor src = get(tm_to_output_vertices, tm.source(h));
    vertex_descriptor tgt = get(tm_to_output_vertices, tm.target(h));
    user_visitor.before_edge_copy(h, tm, output);
    put(tm_to_output_edges, tm.edge(h), new_edge);
    output.set_target(output.opposite(new_h), src);
    output.set_target(new_h, tgt);
    user_visitor.after_edge_copy(h, tm, new_h, output);

    // copy the mark on input edge to the output edge
    copy_edge_mark<SM>(tm.edge(h), new_edge,
                       edge_mark_map_in, edge_mark_map_out);

    CGAL_assertion(is_border(new_h, output));
    CGAL_assertion(is_border(opposite(new_h, output), output));

    // Fill the halfedge of the vertices if it is not set yet
    if (tm.halfedge(tm.target(h)) == h &&
        output.halfedge(output.target(new_h)) == SM::null_halfedge()){
      user_visitor.before_vertex_copy(tm.target(h), tm, output); // Call visitor when the vertex would be complete
      output.set_halfedge(output.target(new_h), new_h);
      user_visitor.after_vertex_copy(tm.target(h), tm, output.target(new_h), output);
    }

    if (tm.halfedge(tm.source(h)) == tm.opposite(h) &&
        output.halfedge(output.source(new_h)) == SM::null_halfedge()){
      user_visitor.before_vertex_copy(tm.source(h), tm, output); // Call visitor when the vertex would be complete
      output.set_halfedge(output.source(new_h), output.opposite(new_h));
      user_visitor.after_vertex_copy(tm.source(h), tm, output.target(new_h), output);
    }
  };
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(true){
    tbb::parallel_for(std::size_t(0), patch.interior_edges.size(), fill_edge);
  }
  else
#endif
  {
    for(std::size_t i=0; i<patch.interior_edges.size(); ++i) fill_edge(i);
  }

  auto get_halfedge = [&](halfedge_descriptor h){
    edge_descriptor e = get(tm_to_output_edges, tm.edge(h));
    halfedge_descriptor h_out = output.halfedge(e);
    if( output.target(h_out) == get(tm_to_output_vertices, tm.target(h))){
      CGAL_assertion( output.source(h_out) == get(tm_to_output_vertices, tm.source(h)) );
      if constexpr(reverse_patch_orientation)
        return output.opposite(h_out);
      else
        return h_out;
    }
    CGAL_assertion( output.target(h_out) == get(tm_to_output_vertices, output.source(h, tm)) );
    CGAL_assertion( output.source(h_out) == get(tm_to_output_vertices, output.target(h, tm)) );
    if constexpr(reverse_patch_orientation)
      return h_out;
    else
      return output.opposite(h_out);
  };
  //create faces and connect halfedges
  auto fill_face = [&](std::size_t i){
    face_descriptor f = patch.faces[i];
    halfedge_descriptor h_in_1 = tm.halfedge(f);
    halfedge_descriptor h_in_2 = tm.next(h_in_1);
    halfedge_descriptor h_in_3 = tm.next(h_in_2);
    CGAL_assertion(tm.next(h_in_3) == h_in_1);

    std::array<halfedge_descriptor, 3> hedges = { get_halfedge(h_in_1), get_halfedge(h_in_2), get_halfedge(h_in_3) };

    user_visitor.before_face_copy(f, tm, output);
    SM_Face_index new_f(faces_idx_begin + i);
    user_visitor.after_face_copy(f, tm, new_f, output);
    output.set_halfedge(new_f, hedges[0]);

    for (int i=0;i<3;++i)
    {
      CGAL_assertion(hedges[i] != null_halfedge());
      if(reverse_patch_orientation)
        output.set_next(hedges[i], hedges[(i+2)%3]);
      else
        output.set_next(hedges[i], hedges[(i+1)%3]);
      output.set_face(hedges[i], new_f);
    }
  };
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(true){
    tbb::parallel_for(std::size_t(0), patch.faces.size(), fill_face);
  }
  else
#endif
  {
    for(std::size_t i=0; i<patch.faces.size(); ++i) fill_face(i);
  }
}

// Specialization of the function for Surface_mesh to exploit vectors structure of Surface_mesh.
template < bool reverse_patch_orientation,
           class Point,
           class PatchContainer,
           class VertexPointMap,
           class VertexPointMapOut,
           class EdgeMarkMapOut,
           class EdgeMarkMapIn,
           class EdgeToEdgeMap,
           class VertexToVertexMap,
           class UserVisitor>
void append_patches_to_triangle_mesh(
  Surface_mesh<Point>& output,
  const boost::dynamic_bitset<>& patches_to_append,
  PatchContainer& patches,
  const VertexPointMapOut& vpm_out,
  const VertexPointMap& vpm_tm,
  EdgeMarkMapOut& edge_mark_map_out,
  const EdgeMarkMapIn& edge_mark_map_in,
  EdgeToEdgeMap& tm_to_output_edges,
  VertexToVertexMap& tm_to_output_vertices,
  UserVisitor& user_visitor)
{
  using SM = Surface_mesh<Point>;
  const SM& tm = patches.pm;

  std::size_t vertices_idx_begin = output.number_of_vertices();
  std::size_t edges_idx_begin = output.number_of_edges();
  std::size_t faces_idx_begin = output.number_of_faces();
  std::size_t tnv = vertices_idx_begin, tne = edges_idx_begin, tnf = faces_idx_begin;
  for (std::size_t i= patches_to_append.find_first();
                   i < patches_to_append.npos;
                   i = patches_to_append.find_next(i))
  {
    tnv += patches[i].interior_vertices.size();
    tne += patches[i].interior_edges.size();
    tnf += patches[i].faces.size();
  }
  output.resize(vertices_idx_begin, edges_idx_begin, faces_idx_begin);
  for (std::size_t i= patches_to_append.find_first();
                   i < patches_to_append.npos;
                   i = patches_to_append.find_next(i))
  {
    Patch_description<SM>& patch=patches[i];
    append_patch<reverse_patch_orientation>(tm, output, patch,
                                            vpm_out, vpm_tm,
                                            edge_mark_map_out, edge_mark_map_in,
                                            tm_to_output_edges, tm_to_output_vertices,
                                            user_visitor,
                                            vertices_idx_begin, edges_idx_begin, faces_idx_begin);
    vertices_idx_begin += patch.interior_vertices.size();
    edges_idx_begin += patch.interior_edges.size();
    faces_idx_begin += patch.faces.size();
  }

  // TODO post process borders
}

// Specialization of fill_new_triangle_mesh for Surface_mesh to exploit vectors structure of Surface_mesh.
template < class Point,
           class IntersectionEdgeMap,
           class VertexPointMap1,
           class VertexPointMap2,
           class VertexPointMapOut,
           class EdgeMarkMap1,
           class EdgeMarkMap2,
           class EdgeMarkMapOut,
           class IntersectionPolylines,
           class PatchContainer1,
           class PatchContainer2,
           class UserVisitor>
void fill_new_triangle_mesh(
  Surface_mesh<Point>& output,
  const boost::dynamic_bitset<>& patches_of_tm1_to_import,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer1& patches_of_tm1,
  PatchContainer2& patches_of_tm2,
  bool reverse_orientation_of_patches_from_tm1,
  bool reverse_orientation_of_patches_from_tm2,
  const IntersectionPolylines& polylines,
  const IntersectionEdgeMap& intersection_edges1,
  const IntersectionEdgeMap& intersection_edges2,
  const VertexPointMap1& vpm1,
  const VertexPointMap2& vpm2,
  const VertexPointMapOut& vpm_out,
  const EdgeMarkMap1& edge_mark_map1,
  const EdgeMarkMap2& edge_mark_map2,
        EdgeMarkMapOut& edge_mark_map_out,
  std::vector< typename Surface_mesh<Point>::Edge_index >& output_shared_edges,
  UserVisitor& user_visitor)
{
  using SM = Surface_mesh<Point>;
  using vertex_descriptor = typename SM::Vertex_index;
  using edge_descriptor = typename SM::Edge_index;
  using halfedge_descriptor = typename SM::Halfedge_index;
  using face_descriptor = typename SM::Face_index;

  using V2V_tag = typename CGAL::dynamic_vertex_property_t<vertex_descriptor>;
  using Vertex_to_vertex_map = typename boost::property_map<SM, V2V_tag>::const_type;

  using E2E_tag = typename CGAL::dynamic_edge_property_t<edge_descriptor>;
  using Edge_to_edge_map = typename boost::property_map<SM, E2E_tag>::const_type;

  const SM &tm1 = patches_of_tm1.pm;
  const SM &tm2 = patches_of_tm2.pm;

  Edge_to_edge_map tm1_to_output_edges = get(E2E_tag(), tm1, SM::null_edge()),
                   tm2_to_output_edges = get(E2E_tag(), tm2, SM::null_edge());
  Vertex_to_vertex_map tm1_to_output_vertices = get(V2V_tag(), tm1, SM::null_vertex()),
                       tm2_to_output_vertices = get(V2V_tag(), tm2, SM::null_vertex());

  std::size_t vertices_idx_begin = output.number_of_vertices();
  std::size_t edges_idx_begin = output.number_of_edges();
  std::size_t faces_idx_begin = output.number_of_faces();

  output_shared_edges.reserve( std::accumulate(polylines.lengths.begin(), polylines.lengths.end(), std::size_t(0)) );
  std::size_t nb_polylines = polylines.lengths.size();
  for (std::size_t i=0; i < nb_polylines; ++i)
    if (!polylines.to_skip.test(i))
      import_polyline(output,
                      polylines.tm1[i], polylines.tm2[i],
                      tm1, tm2,
                      polylines.lengths[i],
                      tm1_to_output_edges, tm2_to_output_edges,
                      tm1_to_output_vertices, tm2_to_output_vertices,
                      intersection_edges1, intersection_edges2,
                      vpm1, vpm2, vpm_out,
                      output_shared_edges,
                      user_visitor);

  // Get ids of patch to append
  std::vector<std::size_t> ids_of_patches_to_append_from_tm1;
  std::vector<std::size_t> ids_of_patches_to_append_from_tm2;
  ids_of_patches_to_append_from_tm1.reserve(patches_of_tm1_to_import.count());
  for (std::size_t i= patches_of_tm1_to_import.find_first();
                  i < patches_of_tm1_to_import.npos;
                  i = patches_of_tm1_to_import.find_next(i)){
    ids_of_patches_to_append_from_tm1.push_back(i);
  }
  ids_of_patches_to_append_from_tm2.reserve(patches_of_tm2_to_import.count());
  for (std::size_t i= patches_of_tm2_to_import.find_first();
                  i < patches_of_tm2_to_import.npos;
                  i = patches_of_tm2_to_import.find_next(i)){
    ids_of_patches_to_append_from_tm2.push_back(i);
  }

  // Compute final sizes
  std::size_t nv = output.number_of_vertices(), ne = output.number_of_edges(), nf = output.number_of_faces();
  std::size_t tnv = output.number_of_vertices(), tne = output.number_of_edges(), tnf = output.number_of_faces();
  for (std::size_t i : ids_of_patches_to_append_from_tm1){
    tnv += patches_of_tm1[i].interior_vertices.size();
    tne += patches_of_tm1[i].interior_edges.size();
    tnf += patches_of_tm1[i].faces.size();
  }
  for (std::size_t i : ids_of_patches_to_append_from_tm2){
    tnv += patches_of_tm2[i].interior_vertices.size();
    tne += patches_of_tm2[i].interior_edges.size();
    tnf += patches_of_tm2[i].faces.size();
  }
  output.resize(tnv, tne, tnf);

  // Append patches
  for (std::size_t i : ids_of_patches_to_append_from_tm1){
    if(reverse_orientation_of_patches_from_tm1)
      append_patch<true>(tm1, output,
                         patches_of_tm1[i],
                         vpm_out, vpm1,
                         edge_mark_map_out,
                         edge_mark_map1,
                         tm1_to_output_edges,
                         tm1_to_output_vertices,
                         user_visitor,
                         nv, ne, nf);
    else
      append_patch<false>(tm1, output,
                          patches_of_tm1[i],
                          vpm_out, vpm1,
                          edge_mark_map_out,
                          edge_mark_map1,
                          tm1_to_output_edges,
                          tm1_to_output_vertices,
                          user_visitor,
                          nv, ne, nf);
    nv += patches_of_tm1[i].interior_vertices.size();
    ne += patches_of_tm1[i].interior_edges.size();
    nf += patches_of_tm1[i].faces.size();
  }

  // TODO post process borders of patches from tm1

  for(std::size_t i : ids_of_patches_to_append_from_tm2){
    if(reverse_orientation_of_patches_from_tm2)
      append_patch<true>(tm2, output,
                         patches_of_tm2[i],
                         vpm_out, vpm2,
                         edge_mark_map_out,
                         edge_mark_map2,
                         tm2_to_output_edges,
                         tm2_to_output_vertices,
                         user_visitor,
                         nv, ne, nf);
    else
      append_patch<false>(tm2, output,
                          patches_of_tm2[i],
                          vpm_out, vpm2,
                          edge_mark_map_out,
                          edge_mark_map2,
                          tm2_to_output_edges,
                          tm2_to_output_vertices,
                          user_visitor,
                          nv, ne, nf);
    nv += patches_of_tm2[i].interior_vertices.size();
    ne += patches_of_tm2[i].interior_edges.size();
    nf += patches_of_tm2[i].faces.size();
  }

  // TODO post process borders of patches from tm2
}

} // namespace Corefinement
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_BOOLEAN_OPERATIONS_UTILS_H