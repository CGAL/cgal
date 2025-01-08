// Copyright (c) 2015-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTION_NO_POLYHEDRAL_ENVELOPE_CHECK
#include <CGAL/Polyhedral_envelope.h>
#endif

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/box_intersection_d.h>
#ifdef CGAL_PMP_REPAIR_SI_USE_OBB_IN_COMPACTIFICATION
#include <CGAL/Optimal_bounding_box/oriented_bounding_box.h>
#endif
#include <CGAL/utility.h>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// #define CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING
// #define CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING
// #define CGAL_PMP_REMOVE_SELF_INTERSECTION_NO_POLYHEDRAL_ENVELOPE_CHECK

// Self-intersection removal is done by making a big-enough hole and filling it
//
// Local self-intersection removal is more subtle and only considers self-intersections
// within a connected component. It then tries to fix those by trying successively:
// - smoothing with the sharp edges in the area being constrained
// - smoothing without the sharp edges in the area being constrained
// - hole-filling with the sharp edges in the area being constrained
// - hole-filling without the sharp edges in the area being constrained
//
// The working area grows as long as we haven't been able to fix the self-intersection,
// up to a user-defined number of times.

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
static int unsolved_self_intersections = 0;
static int self_intersections_solved_by_constrained_smoothing = 0;
static int self_intersections_solved_by_unconstrained_smoothing = 0;
static int self_intersections_solved_by_constrained_hole_filling = 0;
static int self_intersections_solved_by_unconstrained_hole_filling = 0;
#endif

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename Point, typename PolygonMesh, typename VertexPointMap, typename FaceOutputIterator>
FaceOutputIterator replace_faces_with_patch_without_reuse(const std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& border_vertices,
                                                          const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                                                          const std::vector<std::vector<Point> >& patch,
                                                          PolygonMesh& pmesh,
                                                          VertexPointMap vpm,
                                                          FaceOutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor        face_descriptor;

  typedef std::vector<Point>                                                Point_face;
  typedef std::vector<vertex_descriptor>                                    Vertex_face;

  std::map<Point, vertex_descriptor> point_to_vd;

  // First, add those for which the vertex will not change
  for(const vertex_descriptor v : border_vertices)
  {
    // In this version, remove_face() will get rid of isolated vertices so only vertices incident
    // to at least one face that is not going to be removed will be kept
    bool kept_vertex = false;
    for(face_descriptor f : faces_around_target(halfedge(v, pmesh), pmesh))
    {
      if(f != boost::graph_traits<PolygonMesh>::null_face() && faces.count(f) == 0)
      {
        kept_vertex = true;
        break;
      }
    }

    if(kept_vertex)
      point_to_vd[get(vpm, v)] = v;
  }

  for(face_descriptor f : faces)
    Euler::remove_face(halfedge(f, pmesh), pmesh);

  CGAL_assertion(is_valid_polygon_mesh(pmesh));

  // now build a correspondence map and the faces with vertices
  const vertex_descriptor null_v = boost::graph_traits<PolygonMesh>::null_vertex();
  for(const Point_face& face : patch)
  {
    Vertex_face vface;
    vface.reserve(face.size());

    for(const Point& p : face)
    {
      bool success;
      typename std::map<Point, vertex_descriptor>::iterator it;
      std::tie(it, success) = point_to_vd.emplace(p, null_v);
      vertex_descriptor& v = it->second;

      if(success)
      {
        // first time we meet that point, means it's an interior point and we need to make a new vertex
        v = add_vertex(pmesh);
        put(vpm, v, p);
      }

      vface.push_back(v);
    }

    face_descriptor new_f = boost::graph_traits<PolygonMesh>::null_face();
    if(Euler::can_add_face(vface, pmesh))
       new_f = Euler::add_face(vface, pmesh);

    if(new_f == boost::graph_traits<PolygonMesh>::null_face())
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cerr << "Error: failed to insert patch face??" << std::endl;
#endif
      return out;
    }

    out++ = new_f;
  }

  return out;
}


// @todo these could be extracted to somewhere else, it's useful in itself
template <typename Point, typename PolygonMesh, typename VertexPointMap, typename FaceOutputIterator>
FaceOutputIterator replace_faces_with_patch(const std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& border_vertices,
                                            const std::set<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& interior_vertices,
                                            const std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>& border_hedges,
                                            const std::set<typename boost::graph_traits<PolygonMesh>::edge_descriptor>& interior_edges,
                                            const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                                            const std::vector<std::vector<Point> >& patch,
                                            PolygonMesh& pmesh,
                                            VertexPointMap vpm,
                                            FaceOutputIterator out)
{
  static_assert(std::is_same<typename boost::property_traits<VertexPointMap>::value_type, Point>::value);

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor        edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor        face_descriptor;

  typedef std::vector<Point>                                                Point_face;
  typedef std::vector<vertex_descriptor>                                    Vertex_face;

  CGAL_precondition(is_valid_polygon_mesh(pmesh));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Replacing range with patch: ";
  std::cout << faces.size() << " triangles removed, " << patch.size() << " created\n";
#endif

  // To be used to create new elements
  std::vector<vertex_descriptor> vertex_stack(interior_vertices.begin(), interior_vertices.end());
  std::vector<edge_descriptor> edge_stack(interior_edges.begin(), interior_edges.end());
  std::vector<face_descriptor> face_stack(faces.begin(), faces.end());

  // Introduce new vertices, convert the patch in vertex patches
  std::vector<Vertex_face> patch_with_vertices;
  patch_with_vertices.reserve(patch.size());

  std::map<Point, vertex_descriptor> point_to_vd;

  // first, add those for which the vertex will not change
  for(const vertex_descriptor v : border_vertices)
    point_to_vd[get(vpm, v)] = v;

  // now build a correspondence map and the faces with vertices
  const vertex_descriptor null_v = boost::graph_traits<PolygonMesh>::null_vertex();
  for(const Point_face& face : patch)
  {
    Vertex_face vface;
    vface.reserve(face.size());

    for(const Point& p : face)
    {
      bool success;
      typename std::map<Point, vertex_descriptor>::iterator it;
      std::tie(it, success) = point_to_vd.emplace(p, null_v);
      vertex_descriptor& v = it->second;

      if(success) // first time we meet that point, means it's an interior point and we need to make a new vertex
      {
        if(vertex_stack.empty())
        {
          v = add_vertex(pmesh);
        }
        else
        {
          v = vertex_stack.back();
          vertex_stack.pop_back();
        }

        put(vpm, v, p);
      }

      vface.push_back(v);
    }

    patch_with_vertices.push_back(vface);
  }

  typedef std::pair<vertex_descriptor, vertex_descriptor>                        Vertex_pair;
  typedef std::map<Vertex_pair, halfedge_descriptor>                             Vertex_pair_halfedge_map;

  Vertex_pair_halfedge_map halfedge_map;

  // register border halfedges
  for(halfedge_descriptor h : border_hedges)
  {
    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    halfedge_map.emplace(std::make_pair(vs, vt), h);

    set_halfedge(target(h, pmesh), h, pmesh); // update vertex halfedge pointer
  }

  face_descriptor f = boost::graph_traits<PolygonMesh>::null_face();

  for(const Vertex_face& vface : patch_with_vertices)
  {
    if(face_stack.empty())
    {
      f = add_face(pmesh);
    }
    else
    {
      f = face_stack.back();
      face_stack.pop_back();
    }

    CGAL_assertion(f != boost::graph_traits<PolygonMesh>::null_face());
    *out++ = f;

    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(vface.size());

    for(std::size_t i=0, n=vface.size(); i<n; ++i)
    {
      const vertex_descriptor vi = vface[i];
      const vertex_descriptor vj = vface[(i+1)%n];

      // get the corresponding halfedge (either a new one or an already created)
      bool success;
      typename Vertex_pair_halfedge_map::iterator it;
      std::tie(it, success) = halfedge_map.emplace(std::make_pair(vi, vj),
                                                   boost::graph_traits<PolygonMesh>::null_halfedge());
      halfedge_descriptor& h = it->second;

      if(success) // this halfedge is an interior halfedge
      {
        if(edge_stack.empty())
        {
          h = halfedge(add_edge(pmesh), pmesh);
        }
        else
        {
          h = halfedge(edge_stack.back(), pmesh);
          edge_stack.pop_back();
        }

        halfedge_map[std::make_pair(vj, vi)] = opposite(h, pmesh);
      }

      hedges.push_back(h);
    }

    CGAL_assertion(vface.size() == hedges.size());

    // update halfedge connections + face pointers
    for(std::size_t i=0, n=vface.size(); i<n; ++i)
    {
      set_next(hedges[i], hedges[(i+1)%n], pmesh);
      set_face(hedges[i], f, pmesh);

      set_target(hedges[i], vface[(i+1)%n], pmesh);
      set_halfedge(vface[(i+1)%n], hedges[i], pmesh);
    }

    set_halfedge(f, hedges[0], pmesh);
  }

  // now remove the remaining superfluous vertices, edges, faces
  for(vertex_descriptor v : vertex_stack)
    remove_vertex(v, pmesh);
  for(edge_descriptor e : edge_stack)
    remove_edge(e, pmesh);
  for(face_descriptor f : face_stack)
    remove_face(f, pmesh);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT_INTERMEDIATE_FULL_MESH
  CGAL::IO::write_polygon_mesh("results/last_patch_replacement.off", pmesh, CGAL::parameters::stream_precision(17));
#endif

  CGAL_postcondition(is_valid_polygon_mesh(pmesh));

  return out;
}

template <typename Point, typename PolygonMesh, typename VertexPointMap, typename FaceOutputIterator>
FaceOutputIterator replace_faces_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& face_range,
                                            const std::vector<std::vector<Point> >& patch,
                                            PolygonMesh& pmesh,
                                            VertexPointMap vpm,
                                            FaceOutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor       face_descriptor;

  std::vector<vertex_descriptor> border_vertices;
  std::set<vertex_descriptor> interior_vertices;
  std::vector<halfedge_descriptor> border_hedges;
  std::set<edge_descriptor> interior_edges;

  for(face_descriptor fh : face_range)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(fh, pmesh), pmesh))
    {
      if(halfedge(target(h, pmesh), pmesh) == h) // limit the number of insertions
        interior_vertices.insert(target(h, pmesh));
    }
  }

  for(face_descriptor fh : face_range)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(fh, pmesh), pmesh))
    {
      CGAL_assertion(!is_border(h, pmesh));

      const edge_descriptor e = edge(h, pmesh);
      const halfedge_descriptor opp_h = opposite(h, pmesh);
      const face_descriptor opp_f = face(opp_h, pmesh);

      if(is_border(opp_h, pmesh) || face_range.count(opp_f) == 0)
      {
        vertex_descriptor v = target(h, pmesh);
        interior_vertices.erase(v);
        border_hedges.push_back(h);
        border_vertices.push_back(v);
      }
      else
      {
        interior_edges.insert(e);
      }
    }
  }

  return replace_faces_with_patch(border_vertices, interior_vertices,
                                  border_hedges, interior_edges, face_range, patch,
                                  pmesh, vpm, out);
}

template <typename Point, typename PolygonMesh, typename VertexPointMap>
void replace_faces_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                              const std::vector<std::vector<Point> >& patch,
                              PolygonMesh& pmesh,
                              VertexPointMap vpm)
{
  CGAL::Emptyset_iterator out;
  replace_faces_with_patch(faces, patch, pmesh, vpm, out);
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename FaceRange, typename EdgeConstrainMap,
          typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
void constrain_edges(const FaceRange& faces,
                     TriangleMesh& tmesh,
                     const bool constrain_border_edges,
                     const bool constrain_sharp_edges,
                     const double dihedral_angle,
                     const double /*weak_DA*/,
                     EdgeConstrainMap& eif,
                     VertexPointMap vpm,
                     const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename GeomTraits::FT                                           FT;
  typedef typename GeomTraits::Vector_3                                     Vector;

  std::unordered_map<edge_descriptor, bool> is_border_of_selection;
  for(face_descriptor f : faces)
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      // Default initialization is guaranteed to be `false`. Thus, meet it once will switch
      // the value to `true` and meeting it twice will switch back to `false`.
      const edge_descriptor e = edge(h, tmesh);
      if(constrain_sharp_edges)
        is_border_of_selection[e] = !(is_border_of_selection[e]);
      else
        is_border_of_selection[e] = false;
    }
  }

#if 0 // Until detect_features++ is integrated
  CGAL::Polygon_mesh_processing::experimental::detect_sharp_edges_pp(faces, tmesh, dihedral_angle, eif,
                                                                     parameters::weak_dihedral_angle(weak_DA));

  // ...
#else
  // this is basically the code that is in detect_features (at the very bottom)
  // but we do not want a folding to be marked as a sharp feature so the dihedral angle is also
  // bounded from above
  const double bound = dihedral_angle;
  const double cos_angle = std::cos(bound * CGAL_PI / 180.);

  for(const auto& ep : is_border_of_selection)
  {
    bool flag = ep.second;
    if(!constrain_border_edges)
      flag = false;

    if(constrain_sharp_edges && !flag)
    {
      const halfedge_descriptor h = halfedge(ep.first, tmesh);
      CGAL_assertion(!is_border(edge(h, tmesh), tmesh));

      const face_descriptor f1 = face(h, tmesh);
      const face_descriptor f2 = face(opposite(h, tmesh), tmesh);

      // @speed cache normals
      const Vector n1 = compute_face_normal(f1, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
      const Vector n2 = compute_face_normal(f2, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
      if(n1 != CGAL::NULL_VECTOR && n2 != CGAL::NULL_VECTOR)
      {
        const FT c = gt.compute_scalar_product_3_object()(n1, n2);

        // Do not mark as sharp edges with a dihedral angle that is almost `pi` because this is likely
        // due to a fold on the mesh rather than a sharp edge that we would like to preserve
        // (Ideally this would be pre-treated as part of the flatness treatment)
        flag = (c <= cos_angle && c >= -cos_angle);
      }
    }

    is_border_of_selection[ep.first] = flag; // Only needed for output, really
    put(eif, ep.first, flag);
  }
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/constrained_edges.polylines.txt");
  out << std::setprecision(17);
  for(edge_descriptor e : edges(tmesh))
    if(get(eif, e))
       out << "2 " << tmesh.point(source(e, tmesh)) << " " << tmesh.point(target(e, tmesh)) << std::endl;
  out.close();
#endif
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename PolyhedralEnvelope>
bool remove_self_intersections_with_smoothing(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& face_range,
                                              TriangleMesh& tmesh,
                                              const bool constrain_sharp_edges,
                                              const double dihedral_angle,
                                              const double weak_DA,
                                              const PolyhedralEnvelope& cc_envelope,
                                              VertexPointMap vpm,
                                              const GeomTraits& gt)
{
  namespace CP = CGAL::parameters;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: repair with smoothing... (constraining sharp edges: ";
  std::cout << std::boolalpha << constrain_sharp_edges << ")" << std::endl;
#endif

  CGAL_precondition(does_self_intersect(face_range, tmesh));

  // Rather than working directly on the mesh, copy a range and work on this instead
  const CGAL::Face_filtered_graph<TriangleMesh> ffg(tmesh, face_range);
  TriangleMesh local_mesh;
  CGAL::copy_face_graph(ffg, local_mesh, CP::vertex_point_map(vpm));

  // smoothing cannot be applied if the input has degenerate faces
  for(face_descriptor fd : faces(local_mesh))
    if(is_degenerate_triangle_face(fd, local_mesh))
      return false;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  CGAL::IO::write_polygon_mesh("results/local_mesh.off", local_mesh, CGAL::parameters::stream_precision(17));
#endif

  // Constrain sharp and border edges
  typedef CGAL::dynamic_edge_property_t<bool>                                 Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type EIFMap;
  EIFMap eif = get(Edge_property_tag(), local_mesh);

  VertexPointMap local_vpm = get_property_map(vertex_point, local_mesh);

  constrain_edges(faces(local_mesh), local_mesh, true /*constrain_borders*/,
                  constrain_sharp_edges, dihedral_angle, weak_DA, eif, local_vpm, gt);

  // @todo choice of number of iterations? Till convergence && max of 100?
  Polygon_mesh_processing::angle_and_area_smoothing(faces(local_mesh),
                                                    local_mesh,
                                                    CP::edge_is_constrained_map(eif)
                                                    .number_of_iterations(100)
                                                    .use_safety_constraints(false)
#ifndef CGAL_PMP_USE_CERES_SOLVER
                                                    .use_area_smoothing(false)
#endif
                                                    );

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  CGAL::IO::write_polygon_mesh("results/post_smoothing_local_mesh.off", local_mesh, CGAL::parameters::stream_precision(17));
#endif

  if(does_self_intersect(local_mesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: patch still self-intersecting after smoothing\n";
#endif
    return false;
  }
  if (!cc_envelope.is_empty() && !cc_envelope(local_mesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: patch is not in the input polyhedral envelope\n";
#endif
    return false;
  }

  // Patch is acceptable, swap it in
  std::vector<std::vector<Point> > patch;
  for(const face_descriptor f : faces(local_mesh))
  {
    halfedge_descriptor h = halfedge(f, local_mesh);
    patch.emplace_back(std::initializer_list<Point>{get(local_vpm, target(h, local_mesh)),
                                                    get(local_vpm, target(next(h, local_mesh), local_mesh)),
                                                    get(local_vpm, target(prev(h, local_mesh), local_mesh))});
  }

  std::set<face_descriptor> new_faces;
  replace_faces_with_patch(face_range, patch, tmesh, vpm, std::inserter(new_faces, new_faces.end()));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(constrain_sharp_edges)
    ++self_intersections_solved_by_constrained_smoothing;
  else
    ++self_intersections_solved_by_unconstrained_smoothing;
#endif

  return true;
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename TriangleMesh>
bool order_border_halfedge_range(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& hrange,
                                 const TriangleMesh& tmesh)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;

  CGAL_precondition(hrange.size() > 2);

  for(std::size_t i=0; i<hrange.size()-2; ++i)
  {
    const vertex_descriptor tgt = target(hrange[i], tmesh);
    for(std::size_t j=i+1; j<hrange.size(); ++j)
    {
      if(tgt == source(hrange[j], tmesh))
      {
        std::swap(hrange[i+1], hrange[j]);
        break;
      }

      // something went wrong while ordering halfedge (e.g. hole has more than one boundary cycle)
      if(j == hrange.size() - 1)
        return false;
    }
  }

  CGAL_postcondition(source(hrange.front(), tmesh) == target(hrange.back(), tmesh));
  return true;
}

template <class Box, class TM, class VPM, class GT, class OutputIterator>
struct Strict_intersect_edges // "strict" as in "not sharing a vertex"
{
  typedef typename boost::graph_traits<TM>::halfedge_descriptor               halfedge_descriptor;
  typedef typename GT::Segment_3                                              Segment;

  mutable OutputIterator m_iterator;
  const TM& m_tmesh;
  const VPM m_vpmap;

  typename GT::Construct_segment_3 m_construct_segment;
  typename GT::Do_intersect_3 m_do_intersect;

  Strict_intersect_edges(const TM& tmesh, VPM vpmap, const GT& gt, OutputIterator it)
    :
      m_iterator(it),
      m_tmesh(tmesh),
      m_vpmap(vpmap),
      m_construct_segment(gt.construct_segment_3_object()),
      m_do_intersect(gt.do_intersect_3_object())
  {}

  void operator()(const Box* b, const Box* c) const
  {
    const halfedge_descriptor h = b->info();
    const halfedge_descriptor g = c->info();

    if(source(h, m_tmesh) == target(g, m_tmesh) || target(h, m_tmesh) == source(g, m_tmesh))
      return;

    const Segment s1 = m_construct_segment(get(m_vpmap, source(h, m_tmesh)), get(m_vpmap, target(h, m_tmesh)));
    const Segment s2 = m_construct_segment(get(m_vpmap, source(g, m_tmesh)), get(m_vpmap, target(g, m_tmesh)));

    if(m_do_intersect(s1, s2))
      *m_iterator++ = std::make_pair(b->info(), c->info());
  }
};

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool is_simple_3(const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                 const TriangleMesh& tmesh,
                 VertexPointMap vpm,
                 const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor                       halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::reference                            Point_ref;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS                                         Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor, Box_policy> Box;

  std::vector<Box> boxes;
  boxes.reserve(cc_border_hedges.size());

  for(halfedge_descriptor h : cc_border_hedges)
  {
    const Point_ref p = get(vpm, source(h, tmesh));
    const Point_ref q = get(vpm, target(h, tmesh));
    CGAL_assertion(!gt.equal_3_object()(p, q));

    boxes.emplace_back(p.bbox() + q.bbox(), h);
  }

  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(boxes.size());

  for(Box& b : boxes)
    box_ptr.push_back(&b);

  typedef boost::function_output_iterator<CGAL::internal::Throw_at_output>          Throwing_output_iterator;
  typedef internal::Strict_intersect_edges<Box, TriangleMesh, VertexPointMap,
                                           GeomTraits, Throwing_output_iterator>    Throwing_filter;
  Throwing_filter throwing_filter(tmesh, vpm, gt, Throwing_output_iterator());

  try
  {
    const std::ptrdiff_t cutoff = 2000;
    CGAL::box_self_intersection_d<Parallel_if_available_tag>(box_ptr.begin(), box_ptr.end(), throwing_filter, cutoff);
  }
  catch(CGAL::internal::Throw_at_output_exception&)
  {
    return false;
  }

  return true;
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
template <typename Point>
void dump_patch(const std::string filename,
                std::vector<std::vector<Point> >& point_patch)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Writing " << point_patch.size() << " face(s) into " << filename << std::endl;
#endif

  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;

  std::map<Point, int> unique_points_with_id;
  for(const std::vector<Point>& face : point_patch)
    for(const Point& p : face)
      unique_points_with_id.emplace(p, 0);

  out << "OFF\n";
  out << unique_points_with_id.size() << " " << point_patch.size() << " 0\n";

  int unique_id = 0;
  for(auto& e : unique_points_with_id)
  {
    e.second = unique_id++;
    out << e.first << "\n";
  }

  for(const std::vector<Point>& face : point_patch)
  {
    out << face.size();
    for(const Point& p : face)
      out << " " << unique_points_with_id.at(p);
    out << "\n";
  }

  out << std::endl;
  out.close();
}

template <typename FaceContainer, typename PolygonMesh, typename VertexPointMap>
void dump_cc(const std::string filename,
             const FaceContainer& cc_faces,
             const PolygonMesh& pmesh,
             const VertexPointMap vpm)
{
  typedef typename boost::property_traits<VertexPointMap>::value_type      Point;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor       face_descriptor;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Writing " << cc_faces.size() << " face(s) into " << filename << std::endl;
#endif

  std::unordered_map<vertex_descriptor, std::size_t> vertex_ids;
  std::stringstream vss, fss;
  std::size_t id = 0;
  for(face_descriptor f : cc_faces)
  {
    fss << degree(f, pmesh);
    for(vertex_descriptor v : vertices_around_face(halfedge(f, pmesh), pmesh))
    {
      auto res = vertex_ids.emplace(v, id);
      if(res.second) // insert was successful (first time seeing this vertex)
      {
        ++id;
        vss << get(vpm, v) << "\n";
      }

      fss << " " << res.first->second /*id*/;
    }
    fss << "\n";
  }

  std::ofstream out(filename);
  out << std::setprecision(17);

  out << "OFF\n";
  out << id << " " << cc_faces.size() << " 0\n";
  out << vss.str() << "\n" << fss.str() << std::endl;
}
#endif // CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

// Hole filling can be influenced by setting a third point associated to an edge on the border of the hole.
// This third point is supposed to represent how the mesh continues on the other side of the hole.
// If that edge is a border edge, there is no third point (since the opposite face is the null face).
// Similarly if the edge is an internal sharp edge, we don't really want to use the opposite face because
// there is by definition a strong discontinuity and it might thus mislead the hole filling algorithm.
//
// Rather, we construct an artificial third point that is in the same plane as the face incident to `h`,
// defined as the third point of the imaginary equilateral triangle incident to opp(h, tmesh)
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
typename boost::property_traits<VertexPointMap>::value_type
construct_artificial_third_point(const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                                 const TriangleMesh& tmesh,
                                 const VertexPointMap vpm,
                                 const GeomTraits& gt)
{
  typedef typename GeomTraits::FT                                           FT;
  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;
  typedef typename boost::property_traits<VertexPointMap>::reference        Point_ref;
  typedef typename GeomTraits::Vector_3                                     Vector;

  const Point_ref p1 = get(vpm, source(h, tmesh));
  const Point_ref p2 = get(vpm, target(h, tmesh));
  const Point_ref opp_p = get(vpm, target(next(h, tmesh), tmesh));

  // sqrt(3)/2 to have an equilateral triangle with p1, p2, and third_point
  const FT dist = 0.5 * CGAL::sqrt(3.) * CGAL::approximate_sqrt(gt.compute_squared_distance_3_object()(p1, p2));

  const Vector ve1 = gt.construct_vector_3_object()(p1, p2);
  const Vector ve2 = gt.construct_vector_3_object()(p1, opp_p);

  // gram schmidt
  const FT e1e2_sp = gt.compute_scalar_product_3_object()(ve1, ve2);
  Vector orthogonalized_ve2 = gt.construct_sum_of_vectors_3_object()(
                                ve2, gt.construct_scaled_vector_3_object()(ve1, - e1e2_sp));
  Polygon_mesh_processing::internal::normalize(orthogonalized_ve2, gt);

  const Point mid_p1p2 = gt.construct_midpoint_3_object()(p1, p2);
  const Point third_p = gt.construct_translated_point_3_object()(
                          mid_p1p2, gt.construct_scaled_vector_3_object()(orthogonalized_ve2, -dist));

  return third_p;
}

template <typename Point, typename TriangleMesh, typename VertexPointMap>
bool check_patch_compatibility(const std::vector<std::vector<Point> >& patch,
                               const std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& border_vertices,
                               const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& border_hedges,
                               const std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor>& interior_edges,
                               const TriangleMesh& tmesh,
                               const VertexPointMap vpm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;

  if(patch.empty())
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Empty patch" << std::endl;
#endif
    return false;
  }

  std::map<Point, vertex_descriptor> point_to_vd;
  for(vertex_descriptor v : border_vertices)
    point_to_vd[get(vpm, v)] = v;

  // make sure that the hole filling is valid: check that no edge
  // already in the mesh is present in hole_faces.
  bool non_manifold_edge_found = false;
  for(const std::vector<Point>& f : patch)
  {
    for(int i=0; i<3; ++i)
    {
      const Point& p0 = f[i];
      const Point& p1 = f[(i+1)%3];

      auto p0_it = point_to_vd.find(p0);
      auto p1_it = point_to_vd.find(p1);

      // @fixme
      // If any of the vertices is an inner point created through refine(), we don't have an easy way
      //to know the possible correspondency with an existing vertex of the mesh: it might be a vertex
      // part of a completely different CC. Unfortunately, a nm edge could be created with this vertex,
      // but the complexity to check all vertices of the mesh is horrible (even spatially filtered,
      // this needs to be updated, ...)
      if(p0_it == point_to_vd.end() || p1_it == point_to_vd.end())
        continue;

      const vertex_descriptor v0 = p0_it->second;
      const vertex_descriptor v1 = p1_it->second;

      halfedge_descriptor h = halfedge(v0, v1, tmesh).first; // null halfedge if not found
      if(h != boost::graph_traits<TriangleMesh>::null_halfedge())
      {
        if(std::find(border_hedges.begin(), border_hedges.end(), h) == border_hedges.end() &&
           interior_edges.count(edge(h, tmesh)) == 0)
        {
          non_manifold_edge_found = true;
          break;
        }
      }
    }

    if(non_manifold_edge_found)
      break;
  }

  if(non_manifold_edge_found)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Triangulation produced is non-manifold when plugged into the mesh.\n";
#endif

    return false;
  }

  return true;
}

// Patch is not valid if:
// - we insert the same face more than once
// - insert (geometric) non-manifold edges
template <typename TriangleMesh, typename Point>
bool check_patch_sanity(const std::vector<std::vector<Point> >& patch)
{
  std::set<std::set<Point> > unique_faces;
  std::map<std::set<Point>, int> unique_edges;

  for(const std::vector<Point>& face : patch)
  {
    if(!unique_faces.emplace(face.begin(), face.end()).second) // this face had already been found
      return false;

    int i = (unique_edges.insert(std::make_pair(std::set<Point> { face[0], face[1] }, 0)).first->second)++;
    if(i == 2) // non-manifold edge
      return false;

    i = (unique_edges.insert(std::make_pair(std::set<Point> { face[1], face[2] }, 0)).first->second)++;
    if(i == 2) // non-manifold edge
      return false;

    i = (unique_edges.insert(std::make_pair(std::set<Point> { face[2], face[0] }, 0)).first->second)++;
    if(i == 2) // non-manifold edge
      return false;
  }

  // Check for self-intersections within the patch
  // @todo something better than just making a mesh out of the soup?
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;
  std::map<Point, std::size_t> ids;

  std::size_t c = 0;
  for(const std::vector<Point>& face : patch)
  {
    std::vector<std::size_t> ps_f;
    for(const Point& pt : face)
    {
      std::size_t id = c;
      auto is_insert_successful = ids.emplace(pt, c);
      if(is_insert_successful.second) // first time we've seen that point
      {
        ++c;
        points.push_back(pt);
      }
      else // already seen that point
      {
        id = is_insert_successful.first->second;
      }

      CGAL_assertion(id < points.size());
      ps_f.push_back(id);
    }

    faces.push_back(ps_f);
  }

  TriangleMesh patch_mesh;
  if(is_polygon_soup_a_polygon_mesh(faces))
    polygon_soup_to_polygon_mesh(points, faces, patch_mesh);
  else
    return false;

  if(does_self_intersect(patch_mesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Tentative patch has self-intersections." << std::endl;
#endif

    return false;
  }

  return true;
}

template <typename Point, typename GeomTraits>
bool construct_hole_patch(std::vector<CGAL::Triple<int, int, int> >& hole_faces,
                          const std::vector<Point>& hole_points,
                          const std::vector<Point>& third_points,
                          const GeomTraits& gt)
{
  if(hole_points.size() > 3)
  {
    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(hole_faces),
                              parameters::geom_traits(gt));
  }
  else
  {
    hole_faces.emplace_back(0, 1, 2); // trivial hole filling
  }

  if(hole_faces.empty())
  {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to fill a hole using Delaunay search space.\n";
 #endif

    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(hole_faces),
                              parameters::use_delaunay_triangulation(false).geom_traits(gt));
#endif
    if(hole_faces.empty())
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Failed to fill a hole using the whole search space.\n";
#endif
      return false;
    }
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::cout << "  DEBUG: " << hole_faces.size() << " faces in the patch" << std::endl;
  std::vector<std::vector<Point> > to_dump;
  for(const auto& face : hole_faces)
  {
    to_dump.emplace_back(std::initializer_list<Point>{ hole_points[face.first],
                                                       hole_points[face.second],
                                                       hole_points[face.third] });
  }

  CGAL_assertion(to_dump.size() == hole_faces.size());

  static int patch_id = 0;
  std::stringstream oss;
  oss << "results/raw_patch_" << patch_id++ << ".off" << std::ends;
  const std::string filename = oss.str().c_str();

  dump_patch(filename, to_dump);
#endif

  return true;
}

template <typename GeomTraits>
struct Mesh_projection_functor
{
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Triangle_3 Triangle_3;

  typedef std::vector<Triangle_3> Triangle_container;
  typedef CGAL::AABB_triangle_primitive_3<GeomTraits, typename Triangle_container::const_iterator> Primitive;
  typedef CGAL::AABB_traits_3<GeomTraits, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  template <typename TriangleMesh, typename VPM>
  Mesh_projection_functor(const TriangleMesh& mesh,
                          const VPM vpm)
  {
    triangles.reserve(num_faces(mesh));
    for(auto f : faces(mesh))
      triangles.emplace_back(get(vpm, target(halfedge(f, mesh), mesh)),
                             get(vpm, target(next(halfedge(f, mesh), mesh), mesh)),
                             get(vpm, source(halfedge(f, mesh), mesh)));

    tree.insert(std::cbegin(triangles), std::cend(triangles));
  }

  Point_3 operator()(const Point_3& p) const { return tree.closest_point(p); }

private:
  Triangle_container triangles;
  Tree tree;
};

template <typename Point, typename Projector, typename TriangleMesh, typename GeomTraits>
bool adapt_patch(std::vector<std::vector<Point> >& point_patch,
                 const Projector& projector,
                 const TriangleMesh&,
                 const GeomTraits&)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  dump_patch("results/pre-adapt.off", point_patch);
#endif

  CGAL_precondition(!point_patch.empty());

  std::vector<Point> soup_points;
  std::vector<std::array<std::size_t, 3> > soup_faces;

  std::size_t pid = 0;
  std::map<Point, std::size_t> point_ids;
  for(const auto& fp : point_patch)
  {
    CGAL_assertion(fp.size() == 3);
    std::array<std::size_t, 3> f;
    for(std::size_t i=0; i<3; ++i)
    {
      auto res = point_ids.emplace(fp[i], pid);
      if(res.second)
      {
        soup_points.push_back(fp[i]);
        ++pid;
      }
      f[i] = res.first->second;
    }
    soup_faces.push_back(f);
  }

  CGAL_assertion(is_polygon_soup_a_polygon_mesh(soup_faces));

  TriangleMesh local_mesh;
  auto local_vpm = get(vertex_point, local_mesh);

  polygon_soup_to_polygon_mesh(soup_points, soup_faces, local_mesh);
  bool has_SI = does_self_intersect(local_mesh);

  std::vector<halfedge_descriptor> border_hedges;
  border_halfedges(faces(local_mesh), local_mesh, std::back_inserter(border_hedges));

  std::vector<vertex_descriptor> new_vertices;
  refine(local_mesh, faces(local_mesh), CGAL::Emptyset_iterator(), std::back_inserter(new_vertices));

  for(vertex_descriptor v : new_vertices)
    put(local_vpm, v, projector(get(local_vpm, v)));

  // The projector can create degenerate faces
  for (halfedge_descriptor h : border_hedges)
    if (is_degenerate_triangle_face(face(opposite(h, local_mesh), local_mesh), local_mesh))
      return !has_SI;
  if(!remove_degenerate_faces(local_mesh))
    return !has_SI;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  static int adapted_patch_id = 0;
  std::stringstream oss;
  oss << "results/adapted_patch_" << adapted_patch_id++ << ".off" << std::ends;
  const std::string filename = oss.str().c_str();
  std::cout << "  DEBUG: Writing " << point_patch.size() << " faces into " << filename << std::endl;
  IO::write_polygon_mesh(filename, local_mesh);
#endif

  // If the adapted tentative patch has SI, revert back to the base patch
  if(does_self_intersect(local_mesh))
    return !has_SI; // if the base patch also has self-intersections, we are done

  // Replace the tentative patch with the new, self-intersection-less, adapted patch
  point_patch.clear();
  point_patch.reserve(num_faces(local_mesh));

  for(face_descriptor f : faces(local_mesh))
  {
    std::vector<Point> fp { get(local_vpm, target(halfedge(f, local_mesh), local_mesh)),
                            get(local_vpm, target(next(halfedge(f, local_mesh), local_mesh), local_mesh)),
                            get(local_vpm, source(halfedge(f, local_mesh), local_mesh)) };
    point_patch.push_back(fp);
  }

  return true;
}

// This overload uses hole filling to construct a patch and tests the manifoldness of the patch
template <typename Point, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool construct_manifold_hole_patch(std::vector<std::vector<Point> >& point_patch,
                                   const std::vector<Point>& hole_points,
                                   const std::vector<Point>& third_points,
                                   const std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_border_vertices,
                                   const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                   const std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor>& cc_interior_edges,
                                   const Projector& projector,
                                   const TriangleMesh& tmesh,
                                   const VertexPointMap vpm,
                                   const GeomTraits& gt)
{
  typedef CGAL::Triple<int, int, int>                                          Face_indices;

  // Try to triangulate the hole using default parameters
  // (using Delaunay search space if CGAL_HOLE_FILLING_DO_NOT_USE_DT3 is not defined)
  std::vector<Face_indices> hole_faces;
  construct_hole_patch(hole_faces, hole_points, third_points, gt);

  std::vector<std::vector<Point> > local_point_patch;
  local_point_patch.reserve(hole_faces.size());
  for(const Face_indices& face : hole_faces)
  {
    local_point_patch.emplace_back(std::initializer_list<Point>{hole_points[face.first],
                                                                hole_points[face.second],
                                                                hole_points[face.third]});
  }

  if(!adapt_patch(local_point_patch, projector, tmesh, gt))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to adapt the patch..." << std::endl;
#endif
    return false;
  }

  // Check manifoldness compatibility with the rest of the mesh
  if(!check_patch_compatibility(local_point_patch, cc_border_vertices, cc_border_hedges, cc_interior_edges, tmesh, vpm))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Incompatible patch" << std::endl;
#endif
    return false;
  }

  point_patch.reserve(point_patch.size() + local_point_patch.size());
  std::move(std::begin(local_point_patch), std::end(local_point_patch), std::back_inserter(point_patch));

  bool is_sane = check_patch_sanity<TriangleMesh>(point_patch);
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(is_sane)
    std::cout << "  DEBUG: Found sane hole-filling patch (" << point_patch.size() << " faces)\n";
  else
    std::cout << "  DEBUG: Insane hole-filling patch\n";
#endif

  return is_sane;
}

// This overloads fill the containers `cc_interior_vertices` and `cc_interior_edges`
template <typename Point, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool construct_tentative_hole_patch_with_border(std::vector<std::vector<Point> >& point_patch,
                                                const std::vector<Point>& hole_points,
                                                const std::vector<Point>& third_points,
                                                const std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_border_vertices,
                                                const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                                std::set<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_interior_vertices,
                                                std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor>& cc_interior_edges,
                                                const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                                const Projector& projector,
                                                const TriangleMesh& tmesh,
                                                const VertexPointMap vpm,
                                                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  CGAL_assertion(hole_points.size() == third_points.size());

  // Collect vertices and edges inside the current selection cc: first collect all vertices and
  // edges incident to the faces to remove...
  for(const face_descriptor f : cc_faces)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      if(halfedge(target(h, tmesh), tmesh) == h) // to limit the number of insertions
        cc_interior_vertices.insert(target(h, tmesh));

      cc_interior_edges.insert(edge(h, tmesh));
    }
  }

  // ... and then remove those on the boundary
  for(halfedge_descriptor h : cc_border_hedges)
  {
    cc_interior_vertices.erase(target(h, tmesh));
    cc_interior_edges.erase(edge(h, tmesh));
  }

  return construct_manifold_hole_patch(point_patch, hole_points, third_points,
                                       cc_border_vertices, cc_border_hedges, cc_interior_edges,
                                       projector, tmesh, vpm, gt);
}

// This function constructs the ranges `hole_points` and `third_points`. Note that for a sub-hole,
// these two ranges are constructed in another function because we don't want to set 'third_points'
// for edges that are on the border of the sub-hole but not on the border of the (full) hole.
template <typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool construct_tentative_hole_patch(std::vector<std::vector<typename boost::property_traits<VertexPointMap>::value_type> >& patch,
                                    std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_border_vertices,
                                    std::set<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_interior_vertices,
                                    std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor>& cc_interior_edges,
                                    const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                    const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                    const Projector& projector,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap vpm,
                                    const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

  cc_border_vertices.reserve(cc_border_hedges.size());

  std::vector<Point> hole_points, third_points;
  hole_points.reserve(cc_border_hedges.size());
  third_points.reserve(cc_border_hedges.size());

  for(const halfedge_descriptor h : cc_border_hedges)
  {
    const vertex_descriptor v = source(h, tmesh);
    hole_points.push_back(get(vpm, v));
    cc_border_vertices.push_back(v);

    CGAL_assertion(!is_border(h, tmesh));

    if(is_border_edge(h, tmesh))
      third_points.push_back(construct_artificial_third_point(h, tmesh, vpm, gt));
    else
      third_points.push_back(get(vpm, target(next(opposite(h, tmesh), tmesh), tmesh)));
  }

  CGAL_postcondition(hole_points.size() >= 3);

  return construct_tentative_hole_patch_with_border(patch, hole_points, third_points,
                                                    cc_border_vertices, cc_border_hedges,
                                                    cc_interior_vertices, cc_interior_edges,
                                                    cc_faces, projector, tmesh, vpm, gt);
}

// In this overload, we don't know the border of the patch because the face range is a sub-region
// of the hole. We also construct `hole_points` and `third_points`, but with no third point for internal
// sharp edges because a local self-intersection is usually caused by folding and thus we do not want
// a third point resulting from folding to wrongly influence the hole filling process.
template <typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool construct_tentative_sub_hole_patch(std::vector<std::vector<typename boost::property_traits<VertexPointMap>::value_type> >& patch,
                                        const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& sub_cc_faces,
                                        const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                        const Projector& projector,
                                        TriangleMesh& tmesh,
                                        VertexPointMap vpm,
                                        const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

  // Collect halfedges on the boundary of the region to be selected
  // (pointing inside the domain to be remeshed)
  std::set<halfedge_descriptor> internal_hedges;
  std::vector<halfedge_descriptor> cc_border_hedges;
  for(const face_descriptor fd : sub_cc_faces)
  {
    halfedge_descriptor h = halfedge(fd, tmesh);
    for(int i=0; i<3;++i)
    {
      if(is_border(opposite(h, tmesh), tmesh))
      {
         cc_border_hedges.push_back(h);
      }
      else
      {
        const face_descriptor opp_f = face(opposite(h, tmesh), tmesh);
        if(sub_cc_faces.count(opp_f) == 0)
        {
          cc_border_hedges.push_back(h);
          if(cc_faces.count(opp_f) != 0)
            internal_hedges.insert(h);
        }
      }

      h = next(h, tmesh);
    }
  }

  // Sort halfedges so that they describe the sequence of halfedges of the hole to be made
  if(!order_border_halfedge_range(cc_border_hedges, tmesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: More than one border in sub-hole. Not currently handled." << std::endl;
#endif

    return false;
  }

  if(!is_simple_3(cc_border_hedges, tmesh, vpm, gt))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Hole filling cannot handle non-simple sub border" << std::endl;
#endif
    return false;
  }

  // @todo we don't care about these sets, so instead there could be a system of output iterators
  // in construct_tentative_hole_patch instead (and here would be emptyset iterators).
  std::set<vertex_descriptor> cc_interior_vertices;
  std::set<edge_descriptor> cc_interior_edges;

  std::vector<vertex_descriptor> cc_border_vertices;
  cc_border_vertices.reserve(cc_border_hedges.size());

  std::vector<Point> hole_points, third_points;
  hole_points.reserve(cc_border_hedges.size());
  third_points.reserve(cc_border_hedges.size());

  for(const halfedge_descriptor h : cc_border_hedges)
  {
    const vertex_descriptor v = source(h, tmesh);
    hole_points.push_back(get(vpm, v));
    cc_border_vertices.push_back(v);

    CGAL_assertion(!is_border(h, tmesh));

    if(internal_hedges.count(h) == 0 && // `h` is on the border of the full CC
       !is_border_edge(h, tmesh))
    {
      third_points.push_back(get(vpm, target(next(opposite(h, tmesh), tmesh), tmesh)));
    }
    else // `h` is on the border of the sub CC but not on the border of the full CC
    {
      const Point tp = construct_artificial_third_point(h, tmesh, vpm, gt);
      third_points.push_back(tp);
    }
  }

  return construct_tentative_hole_patch_with_border(patch, hole_points, third_points,
                                                    cc_border_vertices, cc_border_hedges,
                                                    cc_interior_vertices, cc_interior_edges,
                                                    sub_cc_faces, projector, tmesh, vpm, gt);
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

// This function is only called when the hole is NOT subdivided into smaller holes
template <typename PolyhedralEnvelope, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
               const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
               const PolyhedralEnvelope& cc_envelope,
               const Projector& projector,
               TriangleMesh& tmesh,
               VertexPointMap vpm,
               const GeomTraits& gt,
               bool reuse_faces = true)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Attempting hole-filling (no constraints), " << cc_faces.size() << " faces\n";
#endif

  std::set<vertex_descriptor> cc_interior_vertices;
  std::set<edge_descriptor> cc_interior_edges;

  std::vector<vertex_descriptor> cc_border_vertices;
  cc_border_vertices.reserve(cc_border_hedges.size());

  std::vector<std::vector<Point> > patch;
  if(!construct_tentative_hole_patch(patch, cc_border_vertices, cc_interior_vertices, cc_interior_edges,
                                     cc_border_hedges, cc_faces, projector, tmesh, vpm, gt))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to find acceptable hole patch\n";
#endif

    return false;
  }

  if(!cc_envelope.is_empty() && !cc_envelope(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Patch is not inside the input polyhedral envelope\n";
#endif
    return false;
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Found acceptable hole-filling patch (" << patch.size() << " faces)\n";
#endif

  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  // Plug the new triangles in the mesh, reusing previous edges and faces
  if(reuse_faces)
  {
    replace_faces_with_patch(cc_border_vertices, cc_interior_vertices,
                             cc_border_hedges, cc_interior_edges,
                             cc_faces, patch, tmesh, vpm,
                             std::inserter(working_face_range, working_face_range.end()));
  }
  else
  {
    replace_faces_with_patch_without_reuse(cc_border_vertices, cc_faces, patch, tmesh, vpm,
                                           std::inserter(working_face_range, working_face_range.end()));
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT_INTERMEDIATE_FULL_MESH
  static int filed_hole_id = 0;
  std::stringstream oss;
  oss << "results/filled_basic_" << filed_hole_id++ << ".off" << std::ends;
  CGAL::IO::write_polygon_mesh(oss.str().c_str(), tmesh, CGAL::parameters::stream_precision(17));
#endif

  CGAL_postcondition(is_valid_polygon_mesh(tmesh));

  return true;
}

// Same function as above but border of the hole is not known
template <typename PolyhedralEnvelope, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole(const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
               const PolyhedralEnvelope& cc_envelope,
               const Projector& projector,
               TriangleMesh& tmesh,
               VertexPointMap vpm,
               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  std::vector<halfedge_descriptor> cc_border_hedges;
  for(face_descriptor fd : cc_faces)
  {
    halfedge_descriptor h = halfedge(fd, tmesh);
    for(int i=0; i<3; ++i)
    {
      if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
        cc_border_hedges.push_back(h);

      h = next(h, tmesh);
    }
  }

  if(order_border_halfedge_range(cc_border_hedges, tmesh))
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, cc_envelope, projector, tmesh, vpm, gt);
  else
    return false;
}

template <typename PolyhedralEnvelope, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole_with_constraints(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                TriangleMesh& tmesh,
                                const double dihedral_angle,
                                const double weak_DA,
                                const PolyhedralEnvelope& cc_envelope,
                                const Projector& projector,
                                VertexPointMap vpm,
                                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Attempting local hole-filling with constrained sharp edges..." << std::endl;
#endif

  // If we are treating self intersections locally, first try to constrain sharp edges in the hole
  typedef CGAL::dynamic_edge_property_t<bool>                                 Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type EIFMap;
  EIFMap eif = get(Edge_property_tag(), tmesh);

  constrain_edges(cc_faces, tmesh, true /*constrain_border_edges*/, true /*constrain_sharp_edges*/,
                  dihedral_angle, weak_DA, eif, vpm, gt);

  // Partition the hole using these constrained edges
  std::set<face_descriptor> visited_faces;
  std::vector<std::vector<Point> > patch;

  for(face_descriptor f : cc_faces)
  {
    if(!visited_faces.insert(f).second) // already visited that face
      continue;

    // gather the faces making a sub-hole
    std::set<face_descriptor> sub_cc;
    Polygon_mesh_processing::connected_component(f, tmesh, std::inserter(sub_cc, sub_cc.end()),
                                                 CGAL::parameters::edge_is_constrained_map(eif));

    visited_faces.insert(sub_cc.begin(), sub_cc.end());

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    dump_cc("results/current_cc.off", sub_cc, tmesh, vpm);
#endif

    // The mesh is not modified, but 'patch' gets filled
    if(!construct_tentative_sub_hole_patch(patch, sub_cc, cc_faces, projector, tmesh, vpm, gt))
    {
      // Something went wrong while finding a potential cover for the a sub-hole --> use basic hole-filling
      return fill_hole(cc_border_hedges, cc_faces, working_face_range, cc_envelope, projector, tmesh, vpm, gt);
    }
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  dump_patch("results/hole_fillers.off", patch);
#endif

  // We're assembling multiple patches so we could have the same face appearing multiple times...
  if(!check_patch_sanity<TriangleMesh>(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Unhealthy patch, defaulting to basic fill_hole" << std::endl;
#endif
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, cc_envelope, projector, tmesh, vpm, gt);
  }

  // check if the patch is inside the input polyhedral envelope
  if(!cc_envelope.is_empty() && !cc_envelope(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Patch is not entirely inside the input polyhedral envelope, defaulting to basic fill_hole" << std::endl;
#endif
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, cc_envelope, projector, tmesh, vpm, gt);
  }

  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  // Plug the hole-filling patch in the mesh
  replace_faces_with_patch(cc_faces, patch, tmesh, vpm,
                           std::inserter(working_face_range, working_face_range.end()));

  return true;
}

template <typename PolyhedralEnvelope, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool remove_self_intersections_with_hole_filling(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                                 const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                                 std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                                 TriangleMesh& tmesh,
                                                 const double strong_dihedral_angle,
                                                 const double weak_dihedral_angle,
                                                 const PolyhedralEnvelope& cc_envelope,
                                                 const Projector& projector,
                                                 VertexPointMap vpm,
                                                 const GeomTraits& gt)
{
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/zone_border.polylines.txt");
  out << std::setprecision(17);
  for(const auto& h : cc_border_hedges)
    out << "2 " << tmesh.point(source(h, tmesh)) << " " << tmesh.point(target(h, tmesh)) << std::endl;
  out.close();
#endif

  if(!order_border_halfedge_range(cc_border_hedges, tmesh))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to orient the boundary??\n";
#endif

    CGAL_assertion(false); // we shouldn't fail to orient the boundary cycle of the complete hole
    return false;
  }

  if(!is_simple_3(cc_border_hedges, tmesh, vpm, gt))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Hole filling cannot handle non-simple border" << std::endl;
#endif
    return false;
  }

  bool fixed_by_hole_filling = false;

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING
  fixed_by_hole_filling = fill_hole_with_constraints(cc_border_hedges, cc_faces, working_face_range,
                                                     tmesh, strong_dihedral_angle, weak_dihedral_angle,
                                                     cc_envelope, projector, vpm, gt);
  if(fixed_by_hole_filling)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    ++self_intersections_solved_by_constrained_hole_filling;
#endif

    return true;
  }
#endif // CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING

  fixed_by_hole_filling = fill_hole(cc_border_hedges, cc_faces, working_face_range, cc_envelope,
                                    projector, tmesh, vpm, gt);
  if(fixed_by_hole_filling)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    ++self_intersections_solved_by_constrained_hole_filling;
#endif

    return true;
  }

  return false;
}

template <typename PolyhedralEnvelope, typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool handle_CC_with_complex_topology(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                     const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                     std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                     TriangleMesh& tmesh,
                                     const double strong_dihedral_angle,
                                     const double weak_dihedral_angle,
                                     const bool preserve_genus,
                                     const PolyhedralEnvelope& cc_envelope,
                                     const Projector& projector,
                                     VertexPointMap vpm,
                                     const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor          edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename boost::property_traits<VertexPointMap>::value_type          Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: CC with Euler_chi != 1" << std::endl;
#endif

  if(preserve_genus)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: CC not handled, selection is not a topological disk (preserve_genus=true)\n";
#endif
    return false;
  }

  const CGAL::Face_filtered_graph<TriangleMesh> ccmesh(tmesh, cc_faces);
  if(!ccmesh.is_selection_valid())
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: invalid FFG selection\n";
#endif
    return false;
  }

  std::vector<halfedge_descriptor> boundary_reps;
  extract_boundary_cycles(ccmesh, std::back_inserter(boundary_reps));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: " << boundary_reps.size() << " borders in the CC\n";
#endif

  if(boundary_reps.size() == 1)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Complex topology but single border --> standard hole filling\n";
#endif

    // If there is a single border, fill the hole as if it were a topological disk.
    // This will lose some information since chi != -1, but preserve_genus = false here
    return remove_self_intersections_with_hole_filling(cc_border_hedges, cc_faces, working_face_range,
                                                       tmesh, strong_dihedral_angle, weak_dihedral_angle,
                                                       cc_envelope, projector, vpm, gt);
  }

  // From there on, there is more than one border

  std::vector<bool> is_hole_incident_to_patch(boundary_reps.size());
  std::vector<FT> hole_lengths(boundary_reps.size());

  int holes_incident_to_patches_n = 0;
  for(std::size_t hole_id = 0; hole_id<boundary_reps.size(); ++hole_id)
  {
    FT border_length = 0;
    bool is_incident_to_patch = false;
    bool count_once = true;
    halfedge_descriptor bh = boundary_reps[hole_id], end = bh;

    // check whether the patch is incident to a face of the input mesh that is not part of the CC
    do
    {
      border_length += edge_length(edge(bh, tmesh), tmesh, CGAL::parameters::vertex_point_map(vpm)
                                                                            .geom_traits(gt));
      if(!is_border(bh, tmesh)) // note the 'tmesh'
      {
        is_incident_to_patch = true;

        if(count_once)
        {
          count_once = false;
          ++holes_incident_to_patches_n;
        }
      }

      bh = next(bh, ccmesh);
    }
    while(bh != end);

    is_hole_incident_to_patch[hole_id] = is_incident_to_patch;
    hole_lengths[hole_id] = border_length;
  }

  // If all border halfedges are "real" border halfedges (i.e., they are border halfedges
  // when looked at in tmesh), then fill only the longest hole
  // @todo when islands can be handled, something better could be attempted
  if(holes_incident_to_patches_n == 0)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Complex topology, multiple borders, hole filling the longest border\n";
#endif
    const auto longest_border_id = std::distance(hole_lengths.begin(),
                                                 std::max_element(hole_lengths.begin(), hole_lengths.end()));
    std::vector<halfedge_descriptor> longest_border_hedges;
    halfedge_descriptor bh = boundary_reps[longest_border_id], end = bh;
    do
    {
      longest_border_hedges.push_back(opposite(bh, tmesh));
      bh = prev(bh, ccmesh); // prev because we insert the opposite
    }
    while(bh != end);

    // 'false' because we can't do on-the-fly patching due to multiple boundary cycles
    // @todo this currently doesn't attempt to constrain sharp edges
    return fill_hole(longest_border_hedges, cc_faces, working_face_range, cc_envelope, projector,
                     tmesh, vpm, gt, false /*reuse*/);
  }

  // If there exists some boundary cycles with "fake" border halfedges, hole-fill those
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Complex topology, some fake borders (" << holes_incident_to_patches_n << ")\n";
#endif

  // This is needed for the patch insertion process at the end
  std::vector<vertex_descriptor> all_border_vertices; // border vertices for all the borders to be filled

  // The patch is built iteratively and made of as many CCs as there are holes being filled
  std::vector<std::vector<Point> > patch;

  for(std::size_t hole_id=0; hole_id<boundary_reps.size(); ++hole_id)
  {
    if(!is_hole_incident_to_patch[hole_id])
      continue;

    std::vector<halfedge_descriptor> border_hedges;
    halfedge_descriptor bh = boundary_reps[hole_id], end = bh;
    do
    {
      border_hedges.push_back(opposite(bh, tmesh));
      bh = prev(bh, ccmesh); // prev because we insert the opposite
    }
    while(bh != end);

    std::vector<vertex_descriptor> border_vertices;
    border_vertices.reserve(border_hedges.size());
    all_border_vertices.reserve(all_border_vertices.size() + border_hedges.size());

    std::vector<Point> hole_points, third_points;
    hole_points.reserve(border_hedges.size());
    third_points.reserve(border_hedges.size());

    for(const halfedge_descriptor h : border_hedges)
    {
      CGAL_assertion(!is_border(h, tmesh));

      const vertex_descriptor v = source(h, tmesh);
      hole_points.push_back(get(vpm, v));

      border_vertices.push_back(v);
      all_border_vertices.push_back(v);

      if(is_border_edge(h, tmesh)) // h is incident to a real face
        third_points.push_back(construct_artificial_third_point(h, tmesh, vpm, gt));
      else
        third_points.push_back(get(vpm, target(next(opposite(h, tmesh), tmesh), tmesh)));
    }

    std::set<vertex_descriptor> interior_vertices;
    std::set<edge_descriptor> interior_edges;

    if(!construct_tentative_hole_patch_with_border(patch, hole_points, third_points,
                                                   border_vertices, border_hedges,
                                                   interior_vertices, interior_edges,
                                                   cc_faces, projector, tmesh, vpm, gt))
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Failed to fill hole #" << hole_id << "\n";
#endif

      return false;
    }
  }

  // Built the patch from all the boundary cycles, put it in
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  dump_patch("results/multiple_real_borders.off", patch);
#endif

  // We're assembling multiple patches so we could have the same face appearing multiple times...
  if(!check_patch_sanity<TriangleMesh>(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Unhealthy patch(s), defaulting to basic fill_hole" << std::endl;
#endif
    return false;
  }

  // check if the patch is inside the input polyhedral envelope
  if(!cc_envelope.is_empty() && !cc_envelope(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Patch is not entirely inside the input polyhedral envelope, defaulting to basic fill_hole" << std::endl;
#endif
    return false;
  }

  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  // Plug the hole-filling patch in the mesh
  replace_faces_with_patch_without_reuse(all_border_vertices, cc_faces, patch, tmesh, vpm,
                                         std::inserter(working_face_range, working_face_range.end()));

  return true;
}

// the parameter `step` controls how many extra layers of faces we take around the range `faces_to_treat`
template <typename Projector, typename TriangleMesh, typename VertexPointMap, typename GeomTraits, typename Visitor>
std::pair<bool, bool>
remove_self_intersections_one_step(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& faces_to_treat,
                                   std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                   TriangleMesh& tmesh,
                                   const int step,
                                   const bool preserve_genus,
                                   const bool treat_all_CCs,
                                   const double strong_dihedral_angle,
                                   const double weak_dihedral_angle,
                                   const bool use_smoothing,
                                   const double containment_epsilon,
                                   const Projector& projector,
                                   VertexPointMap vpm,
                                   const GeomTraits& gt,
                                   Visitor& visitor)
{
  typedef boost::graph_traits<TriangleMesh>                                    graph_traits;
  typedef typename graph_traits::halfedge_descriptor                           halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                               face_descriptor;

#ifdef CGAL_PMP_REPAIR_SI_USE_OBB_IN_COMPACTIFICATION
  typedef typename graph_traits::vertex_descriptor                             vertex_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::value_type          Point;
#endif

  std::set<face_descriptor> faces_to_treat_copy = faces_to_treat;

  bool something_was_done = false; // indicates if a region was successfully remeshed
  bool all_fixed = true; // indicates if all removal went well
  // indicates if a removal was not possible because the region handle has
  // some boundary cycle of halfedges
  bool topology_issue = false;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: is_valid in one_step(tmesh)? " << is_valid_polygon_mesh(tmesh) << std::endl;

  unsolved_self_intersections = 0;
#endif

  CGAL_precondition(is_valid_polygon_mesh(tmesh));

#if defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG) || defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT)
  int cc_id = -1;
#endif

  while(!faces_to_treat.empty())
  {
    if(visitor.stop())
      return std::make_pair(false, false);

    visitor.start_component_handling();
    visitor.status_update(faces_to_treat);

#if defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG) || defined(CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT)
    ++cc_id;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Remaining faces to remove: " << faces_to_treat.size() << "\n";
    std::cout << "  DEBUG: --------------- Considering CC #" << cc_id << " ---------------\n";
    std::cout << "  DEBUG: Initial face " << *faces_to_treat.begin() << "\n";
    std::cout << "  DEBUG: first face: " << get(vpm, target(halfedge(*(faces_to_treat.begin()), tmesh), tmesh)) << " "
                                         << get(vpm, target(next(halfedge(*(faces_to_treat.begin()), tmesh), tmesh), tmesh)) << " "
                                         << get(vpm, source(halfedge(*(faces_to_treat.begin()), tmesh), tmesh)) << "\n";
#endif

    // Collect all the faces from the connected component
    std::set<face_descriptor> cc_faces;
    std::vector<face_descriptor> queue(1, *faces_to_treat.begin()); // temporary queue
    cc_faces.insert(queue.back());
    while(!queue.empty())
    {
      face_descriptor top = queue.back();
      queue.pop_back();
      halfedge_descriptor h = halfedge(top, tmesh);
      for(int i=0; i<3; ++i)
      {
        face_descriptor adjacent_face = face(opposite(h, tmesh), tmesh);
        if(adjacent_face != boost::graph_traits<TriangleMesh>::null_face())
        {
          if(faces_to_treat.count(adjacent_face) != 0 && cc_faces.insert(adjacent_face).second)
            queue.push_back(adjacent_face);
        }

        h = next(h, tmesh);
      }
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in base CC\n";
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    std::string fname = "results/initial_step_"+std::to_string(step)+"_CC_" + std::to_string(cc_id)+".off";
    dump_cc(fname, cc_faces, tmesh, vpm);
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT_INTERMEDIATE_FULL_MESH
    fname = "results/mesh_at_step_"+std::to_string(step)+"_CC_"+std::to_string(cc_id)+".off";
    CGAL::IO::write_polygon_mesh(fname, tmesh, CGAL::parameters::stream_precision(17));
#endif

    // expand the region to be filled
    if(step > 0)
    {
      expand_face_selection(cc_faces, tmesh, step,
                            make_boolean_property_map(cc_faces),
                            Emptyset_iterator());
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in expanded CC\n";

    fname = "results/expanded_step_"+std::to_string(step)+"_CC_"+std::to_string(cc_id)+".off";
    dump_cc(fname, cc_faces, tmesh, vpm);
#endif

    // @todo keep this?
    // try to compactify the selection region by also selecting all the faces included
    // in the bounding box of the initial selection
    std::vector<halfedge_descriptor> stack_for_expension;

#ifdef CGAL_PMP_REPAIR_SI_USE_OBB_IN_COMPACTIFICATION
    std::set<Point> cc_points;
    for(face_descriptor f : cc_faces)
      for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh), tmesh))
          cc_points.insert(get(vpm, v));

    typedef typename GeomTraits::Aff_transformation_3 Aff_transformation;
    Aff_transformation tr{CGAL::Identity_transformation()};

    if(cc_points.size() > 3)
      CGAL::oriented_bounding_box(cc_points, tr, CGAL::parameters::random_seed(0));

    // Construct the rotated OBB
    Bbox_3 bb;
    for(const Point& p : cc_points)
      bb += (tr.transform(p)).bbox();

#else
    Bbox_3 bb;
#endif

    for(face_descriptor fd : cc_faces)
    {
      for(halfedge_descriptor h : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
#ifndef CGAL_PMP_REPAIR_SI_USE_OBB_IN_COMPACTIFICATION
        bb += get(vpm, target(h, tmesh)).bbox();
#endif
        face_descriptor nf = face(opposite(h, tmesh), tmesh);
        if(nf != boost::graph_traits<TriangleMesh>::null_face() && cc_faces.count(nf) == 0)
          stack_for_expension.push_back(opposite(h, tmesh));
      }
    }

    while(!stack_for_expension.empty())
    {
      halfedge_descriptor h = stack_for_expension.back();
      stack_for_expension.pop_back();
      if(cc_faces.count(face(h, tmesh)) == 1)
        continue;

#ifdef CGAL_PMP_REPAIR_SI_USE_OBB_IN_COMPACTIFICATION
      if(do_overlap(bb, tr.transform(get(vpm, target(next(h, tmesh), tmesh))).bbox()))
#else
      if(do_overlap(bb, get(vpm, target(next(h, tmesh), tmesh)).bbox()))
#endif
      {
        cc_faces.insert(face(h, tmesh));
        halfedge_descriptor candidate = opposite(next(h, tmesh), tmesh);
        if(face(candidate, tmesh) != boost::graph_traits<TriangleMesh>::null_face())
          stack_for_expension.push_back(candidate);

        candidate = opposite(prev(h, tmesh), tmesh);
        if(face(candidate, tmesh) != boost::graph_traits<TriangleMesh>::null_face())
          stack_for_expension.push_back(candidate);
      }
    }

    Boolean_property_map<std::set<face_descriptor> > is_selected(cc_faces);
    expand_face_selection_for_removal(cc_faces, tmesh, is_selected);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in expanded and compactified CC\n";
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    fname = "results/expanded_compactified_step_"+std::to_string(step)+"_CC_"+std::to_string(cc_id)+".off";
    dump_cc(fname, cc_faces, tmesh, vpm);
#endif

    // Now, we have a proper selection to work on.

    for(const face_descriptor f : cc_faces)
      faces_to_treat.erase(f);

    if(cc_faces.size() == 1)
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Compactified CC of size 1, moving on\n";
#endif
      visitor.end_component_handling();
      continue;
    }

    bool self_intersects = does_self_intersect(cc_faces, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    if(!self_intersects)
      std::cout << "  DEBUG: No self-intersection within the CC\n";
#endif

    if(!treat_all_CCs && !self_intersects)
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      ++unsolved_self_intersections;
#endif

      all_fixed = false;
      visitor.end_component_handling();
      continue;
    }

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTION_NO_POLYHEDRAL_ENVELOPE_CHECK
    Polyhedral_envelope<GeomTraits> cc_envelope;
    if(containment_epsilon != 0)
      cc_envelope = Polyhedral_envelope<GeomTraits>(cc_faces, tmesh, containment_epsilon);
#else
    struct Return_true
    {
      constexpr bool is_empty() const { return true; }
      bool operator()(const std::vector<std::vector<typename GeomTraits::Point_3> >&) const { return true; }
      bool operator()(const TriangleMesh&) const { return true; }
    };

    Return_true cc_envelope;
    CGAL_USE(containment_epsilon);
#endif

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING
    // First, try to smooth if we only care about local self-intersections
    // Two different approaches:
    // - First, try to constrain edges that are in the zone to smooth and whose dihedral angle is large,
    //   but not too large (we don't want to constrain edges that are foldings);
    // - If that fails, try to smooth without any constraints, but make sure that the deviation from
    //   the first zone is small.
    //
    // If smoothing fails, the face patch is restored to its pre-smoothing state.
    //
    // There is no need to update the working range because smoothing doesn`t change
    // the number of faces (and old faces are reused).
    //
    // Do not smooth if there are no self-intersections within the patch: this means the intersection
    // is with another CC and smoothing is unlikely to move the surface sufficiently
    if(use_smoothing && self_intersects)
    {
      bool fixed_by_smoothing = false;

      fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, tmesh, true /*constrain_sharp_edges*/,
                                                                    strong_dihedral_angle, weak_dihedral_angle,
                                                                    cc_envelope, vpm, gt);

      if(!fixed_by_smoothing)
      {
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: Could not be solved via smoothing with constraints\n";
 #endif

        // try again, but without constraining sharp edges
        fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, tmesh, false /*constrain_sharp_edges*/,
                                                                      strong_dihedral_angle, weak_dihedral_angle,
                                                                      cc_envelope, vpm, gt);
      }

      if(fixed_by_smoothing)
      {
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: Solved with smoothing!\n";
 #endif

        something_was_done = true;
        visitor.end_component_handling();
        continue;
      }
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      else
      {
        std::cout << "  DEBUG: Could not be solved via smoothing\n";
      }
 #endif
    }

#endif // ndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Trying hole-filling based approach...\n";
#endif

    // Collect halfedges on the boundary of the region to be selected
    // (incident to faces that are part of the CC)
    std::vector<halfedge_descriptor> cc_border_hedges;
    for(face_descriptor fd : cc_faces)
    {
      for(halfedge_descriptor h : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
        if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
          cc_border_hedges.push_back(h);
      }
    }

    // Whichever step we are at, no border means no expansion will change this selection
    // This CC was not fixed by smoothing, and there is nothing hole filling can do
    // @todo just remove the CC?
    if(cc_border_hedges.empty())
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: CC is closed!\n"; // @todo wrap?
      ++unsolved_self_intersections;
#endif

      all_fixed = false;
      visitor.end_component_handling();
      continue;
    }

    int selection_chi = euler_characteristic_of_selection(cc_faces, tmesh);
    if(selection_chi != 1) // not a topological disk
    {
      if(!handle_CC_with_complex_topology(cc_border_hedges, cc_faces, working_face_range,
                                          tmesh, strong_dihedral_angle, weak_dihedral_angle,
                                          preserve_genus, cc_envelope, projector, vpm, gt))
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: Failed to handle complex CC\n";
        ++unsolved_self_intersections;
#endif
        topology_issue = true;
        all_fixed = false;
      }
      else
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        ++self_intersections_solved_by_unconstrained_hole_filling;
#endif

        something_was_done = true;
      }

      visitor.end_component_handling();
      continue;
    }

    // From here on, the CC is a topological disk

    if(!remove_self_intersections_with_hole_filling(cc_border_hedges, cc_faces, working_face_range,
                                                    tmesh, strong_dihedral_angle, weak_dihedral_angle,
                                                    cc_envelope, projector, vpm, gt))
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Failed to fill hole\n";
      ++unsolved_self_intersections;
#endif

      all_fixed = false;
    }
    else
    {
      something_was_done = true;
    }

    visitor.end_component_handling();
  }

  if(!something_was_done)
  {
    faces_to_treat.swap(faces_to_treat_copy);
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Nothing was changed during this step, self-intersections won`t be recomputed." << std::endl;
#endif
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: " << unsolved_self_intersections << " unsolved SI" << std::endl;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT_INTERMEDIATE_FULL_MESH
  std::stringstream oss;
  oss << "results/after_step_" << step << ".off" << std::ends;
  CGAL::IO::write_polygon_mesh(oss.str().c_str(), tmesh, CGAL::parameters::stream_precision(17));
#endif

  return std::make_pair(all_fixed, topology_issue);
}

} // namespace internal

namespace experimental {

template <class TriangleMesh>
struct Remove_self_intersection_default_visitor
{
  constexpr bool stop() const { return false; }
  template <class FaceContainer>
  void status_update(const FaceContainer&) {}
  void start_main_loop() {}
  void end_main_loop() {}
  void start_iteration() {}
  void end_iteration() {}
  void start_component_handling() {}
  void end_component_handling() {}
  void parameters_used( bool /* parameters_used(preserve_genus */,
                        bool /* treat_all_CCs */,
                        int /* max_steps */,
                        double /* strong_dihedral_angle */,
                        double /* weak_dihedral_angle */,
                        double /* containment_epsilon */ ) {}
};

template <typename FaceRange, typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
bool remove_self_intersections(const FaceRange& face_range,
                               TriangleMesh& tmesh,
                               const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef boost::graph_traits<TriangleMesh>                                 graph_traits;
  typedef typename graph_traits::face_descriptor                            face_descriptor;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  IO::write_polygon_mesh("results/input.off", tmesh, parameters::stream_precision(17));
#endif

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       GeomTraits;
  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));

  bool preserve_genus = choose_parameter(get_parameter(np, internal_np::preserve_genus), true);

  // @tmp Squatting that named parameter to signify that treatment should be applied within the CC
  // even if there are no self-intersections. For example, two spheres intersecting each other.
  const bool treat_all_CCs = choose_parameter(get_parameter(np, internal_np::apply_per_connected_component), true);

  // When treating intersections locally, we don't want to grow the working range too much as
  // either the solution is found fast, or it's too difficult and neither local smoothing or local
  // hole filling are going to provide nice results.
  const int default_max_step = 7;
  const int max_steps = choose_parameter(get_parameter(np, internal_np::number_of_iterations), default_max_step);

  // @fixme give it its own named parameter rather than abusing 'with_dihedral_angle'?
  const double strong_dihedral_angle = choose_parameter(get_parameter(np, internal_np::with_dihedral_angle), 60.);

  // detect_feature_pp NP (unused for now)
  const double weak_dihedral_angle = 0.; // choose_parameter(get_parameter(np, internal_np::weak_dihedral_angle), 20.);

  const bool use_smoothing = choose_parameter(get_parameter(np, internal_np::use_smoothing), false);

  struct Return_false
  {
    bool operator()(std::pair<face_descriptor, face_descriptor>) const { return false; }
  };

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::filter_t,
    NamedParameters,
    Return_false//default
  > ::type  Output_iterator_predicate;
  Output_iterator_predicate out_it_predicates
    = choose_parameter<Return_false>(get_parameter(np, internal_np::filter));

  // use containment check
  const double containment_epsilon = choose_parameter(get_parameter(np, internal_np::polyhedral_envelope_epsilon), 0.);

  internal::Mesh_projection_functor<GeomTraits> projector(tmesh, vpm);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";
  std::cout << "\tpreserve_genus: " << preserve_genus << std::endl;
  std::cout << "\ttreat_all_CCs: " << treat_all_CCs << std::endl;
  std::cout << "\tmax_steps: " << max_steps << std::endl;
  std::cout << "\tstrong_dihedral_angle: " << strong_dihedral_angle << std::endl;
  std::cout << "\tweak_dihedral_angle: " << weak_dihedral_angle << std::endl;
  std::cout << "\tcontainment_epsilon: " << containment_epsilon << std::endl;
#endif

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Remove_self_intersection_default_visitor<TriangleMesh>//default
  > ::type Visitor;
  Visitor visitor = choose_parameter<Visitor>(get_parameter(np, internal_np::visitor));

  visitor.parameters_used(preserve_genus,
                          treat_all_CCs,
                          max_steps,
                          strong_dihedral_angle,
                          weak_dihedral_angle,
                          containment_epsilon);

  if(!preserve_genus)
    duplicate_non_manifold_vertices(tmesh, np);

  // Look for self-intersections in the mesh and remove them
  int step = -1;
  bool all_fixed = true; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> faces_to_treat;
  std::set<face_descriptor> working_face_range(face_range.begin(), face_range.end());

  visitor.start_main_loop();
  while(++step < max_steps)
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: ========== STEP " << step << " / " << max_steps - 1 << " ==========" << std::endl;
#endif

    if(visitor.stop())
      break;

    visitor.start_iteration();

    if(faces_to_treat.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;

      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(working_face_range, tmesh,
                         filter_output_iterator(std::back_inserter(self_inter), out_it_predicates),
                         parameters::vertex_point_map(vpm).geom_traits(gt));
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: " << self_inter.size() << " intersecting pairs" << std::endl;
#endif
      for(const Face_pair& fp : self_inter)
      {
        faces_to_treat.insert(fp.first);
        faces_to_treat.insert(fp.second);
      }
    }

    if(faces_to_treat.empty())
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "DEBUG: There are no more faces to treat." << std::endl;
#endif
      break;
    }

    visitor.status_update(faces_to_treat);

    std::tie(all_fixed, topology_issue) =
      internal::remove_self_intersections_one_step(
          faces_to_treat, working_face_range, tmesh, step,
          preserve_genus, treat_all_CCs, strong_dihedral_angle, weak_dihedral_angle,
          use_smoothing, containment_epsilon, projector, vpm, gt, visitor);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    if(all_fixed && topology_issue)
        std::cout << "DEBUG: boundary cycles of boundary edges involved in self-intersections.\n";
#endif

    visitor.end_iteration();
  }
  visitor.end_main_loop();

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "solved by constrained smoothing: " << internal::self_intersections_solved_by_constrained_smoothing << std::endl;
  std::cout << "solved by unconstrained smoothing: " << internal::self_intersections_solved_by_unconstrained_smoothing << std::endl;
  std::cout << "solved by constrained hole-filling: " << internal::self_intersections_solved_by_constrained_hole_filling << std::endl;
  std::cout << "solved by unconstrained hole-filling: " << internal::self_intersections_solved_by_unconstrained_hole_filling << std::endl;
  std::cout << "issues during CC treatment: " << internal::unsolved_self_intersections << std::endl;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream("results/final.off") << std::setprecision(17) << tmesh;
#endif

  bool self_intersects = does_self_intersect(working_face_range, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(self_intersects)
    std::cout << "DEBUG: Failed to solve all self-intersections.\n";
#endif

  return !self_intersects;
}

template <typename TriangleMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool remove_self_intersections(TriangleMesh& tmesh, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return remove_self_intersections(faces(tmesh), tmesh, np);
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H
