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

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/utility.h>

#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

// #define CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING
// #define CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING

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

// @todo these could be extracted to somewhere else, it's useful in itself
template <typename PolygonMesh, typename VPM, typename Point, typename FaceOutputIterator>
FaceOutputIterator replace_faces_with_patch(const std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& border_vertices,
                                            const std::set<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& interior_vertices,
                                            const std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>& border_hedges,
                                            const std::set<typename boost::graph_traits<PolygonMesh>::edge_descriptor>& interior_edges,
                                            const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                                            const std::vector<std::vector<Point> >& patch,
                                            PolygonMesh& pmesh,
                                            VPM& vpm,
                                            FaceOutputIterator out)
{
  CGAL_static_assertion((std::is_same<typename boost::property_traits<VPM>::value_type, Point>::value));

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor        edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor        face_descriptor;

  typedef std::vector<Point>                                                Point_face;
  typedef std::vector<vertex_descriptor>                                    Vertex_face;

  CGAL_precondition(is_valid_polygon_mesh(pmesh));

  // To be used to create new elements
  std::vector<vertex_descriptor> vertex_stack(interior_vertices.begin(), interior_vertices.end());
  std::vector<edge_descriptor> edge_stack(interior_edges.begin(), interior_edges.end());
  std::vector<face_descriptor> face_stack(faces.begin(), faces.end());

  // Introduce new vertices, convert the patch in vertex patches
  std::vector<Vertex_face> patch_with_vertices;
  patch_with_vertices.reserve(patch.size());

  std::map<Point, vertex_descriptor> point_to_vs;

  // first, add those for which the vertex will not change
  for(const vertex_descriptor v : border_vertices)
    point_to_vs[get(vpm, v)] = v;

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
      std::tie(it, success) = point_to_vs.insert(std::make_pair(p, null_v));
      vertex_descriptor& v = it->second;

      if(success) // first time we meet that point, means it`s an interior point and we need to make a new vertex
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
  int i = 0;
  for(halfedge_descriptor h : border_hedges)
  {
    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    halfedge_map.insert(std::make_pair(std::make_pair(vs, vt), h));

    set_halfedge(target(h, pmesh), h, pmesh); // update vertex halfedge pointer
    ++i;
  }

  face_descriptor f = boost::graph_traits<PolygonMesh>::null_face();
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::vector<face_descriptor> new_faces;
#endif

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

    *out++ = f;
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    new_faces.push_back(f);
#endif

    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(vface.size());

    for(std::size_t i=0, n=vface.size(); i<n; ++i)
    {
      const vertex_descriptor vi = vface[i];
      const vertex_descriptor vj = vface[(i+1)%n];

      // get the corresponding halfedge (either a new one or an already created)
      bool success;
      typename Vertex_pair_halfedge_map::iterator it;
      std::tie(it, success) = halfedge_map.insert(std::make_pair(std::make_pair(vi, vj),
                                                  boost::graph_traits<PolygonMesh>::null_halfedge()));
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

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream res_out("results/last_patch_replacement.off");
  res_out << std::setprecision(17) << pmesh;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Replacing range with patch: ";
  std::cout << faces.size() << " triangles removed, " << patch.size() << " created\n";
#endif

  CGAL_postcondition(is_valid_polygon_mesh(pmesh));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  CGAL_postcondition(!does_self_intersect(new_faces, pmesh));
#endif

  return out;
}

template <typename PolygonMesh, typename VPM, typename Point, typename FaceOutputIterator>
FaceOutputIterator replace_faces_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& face_range,
                                            const std::vector<std::vector<Point> >& patch,
                                            PolygonMesh& pmesh,
                                            VPM& vpm,
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

template <typename PolygonMesh, typename VPM, typename Point>
void replace_faces_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                              const std::vector<std::vector<Point> >& patch,
                              PolygonMesh& pmesh,
                              VPM& vpm)
{
  CGAL::Emptyset_iterator out;
  replace_faces_with_patch(faces, patch, pmesh, vpm, out);
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename Point, typename FaceRange, typename PolygonMesh, typename VertexPointMap>
void back_up_face_range_as_point_patch(std::vector<std::vector<Point> >& point_patch,
                                       const FaceRange& face_range,
                                       const PolygonMesh& tmesh,
                                       const VertexPointMap vpm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor        face_descriptor;

  point_patch.reserve(face_range.size());

  for(const face_descriptor f : face_range)
  {
    std::vector<Point> face_points;
    for(const halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
      face_points.push_back(get(vpm, target(h, tmesh)));

    point_patch.push_back(face_points);
  }
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

template <typename FaceRange, typename EdgeConstrainMap,
          typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
void constrain_sharp_and_border_edges(const FaceRange& faces,
                                      TriangleMesh& tmesh,
                                      EdgeConstrainMap& eif,
                                      const bool constrain_sharp_edges,
                                      const double dihedral_angle,
                                      const double /*weak_DA*/,
                                      VertexPointMap vpm,
                                      const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename GeomTraits::FT                                           FT;
  typedef typename GeomTraits::Vector_3                                     Vector;

  std::map<edge_descriptor, bool> is_border_of_selection;
  for(face_descriptor f : faces)
  {
    // @fixme what about nm vertices
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      // Default initialization is guaranteed to be `false`. Thus, meet it once will switch
      // the value to `true` and meeting it twice will switch back to `false`.
      const edge_descriptor e = edge(h, tmesh);
      is_border_of_selection[e] = !(is_border_of_selection[e]);
    }
  }

#if 0 // Until detect_features++ is integrated
  CGAL::Polygon_mesh_processing::experimental::detect_sharp_edges_pp(faces, tmesh, dihedral_angle, eif,
                                                                     parameters::weak_dihedral_angle(weak_DA));

  // borders are also constrained
  for(const auto& ep : is_border_of_selection)
    if(ep.second)
      put(eif, ep.first, true);
#else
  // this is basically the code that is in detect_features (at the very bottom)
  // but we do not want a folding to be marked as a sharp feature so the dihedral angle is also
  // bounded from above
  const double bound = dihedral_angle;
  const double cos_angle = std::cos(bound * CGAL_PI / 180.);

  for(const auto& ep : is_border_of_selection)
  {
    bool flag = ep.second;
    if(constrain_sharp_edges && !flag)
    {
      const halfedge_descriptor h = halfedge(ep.first, tmesh);
      CGAL_assertion(!is_border(edge(h, tmesh), tmesh));

      const face_descriptor f1 = face(h, tmesh);
      const face_descriptor f2 = face(opposite(h, tmesh), tmesh);

      // @todo cache normals
      const Vector n1 = compute_face_normal(f1, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
      const Vector n2 = compute_face_normal(f2, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt));
      const FT c = gt.compute_scalar_product_3_object()(n1, n2);

      // Do not mark as sharp edges with a dihedral angle that is almost `pi` because this is likely
      // due to a foldness on the mesh rather than a sharp edge that we wish to preserve
      // (Ideally this would be pre-treated as part of the flatness treatment)
      flag = (c <= cos_angle && c >= -cos_angle);
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

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool remove_self_intersections_with_smoothing(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& face_range,
                                              TriangleMesh& tmesh,
                                              const bool constrain_sharp_edges,
                                              const double dihedral_angle,
                                              const double weak_DA,
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

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out_p("results/local_mesh.off");
  out_p << std::setprecision(17) << local_mesh;
  out_p.close();
#endif

  // Constrain sharp and border edges
  typedef CGAL::dynamic_edge_property_t<bool>                                 Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type EIFMap;
  EIFMap eif = get(Edge_property_tag(), local_mesh);

  VertexPointMap local_vpm = get_property_map(vertex_point, local_mesh);

  constrain_sharp_and_border_edges(faces(local_mesh), local_mesh, eif, constrain_sharp_edges,
                                   dihedral_angle, weak_DA, local_vpm, gt);

  // @todo choice of number of iterations? Till convergence && max of 100?
  Polygon_mesh_processing::smooth_mesh(faces(local_mesh), local_mesh, CP::edge_is_constrained_map(eif)
                                                                         .number_of_iterations(100)
                                                                         .use_safety_constraints(false));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/post_smoothing_local_mesh.off");
  out << std::setprecision(17) << local_mesh;
  out.close();
#endif

  if(does_self_intersect(local_mesh))
    return false;

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
  CGAL_assertion(!does_self_intersect(new_faces, tmesh, parameters::vertex_point_map(vpm)));

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

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT

template <typename FaceContainer, typename TriangleMesh>
void dump_cc(const FaceContainer& cc_faces,
             const TriangleMesh& mesh,
             const std::string filename)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh>::const_type VertexPointMap;
  VertexPointMap vpm =get_const_property_map(vertex_point, mesh);

  std::ofstream out(filename);
  out.precision(17);

  out << "OFF\n";
  out << 3*cc_faces.size() << " " << cc_faces.size() << " 0\n";

  for(const face_descriptor f : cc_faces)
  {
    out << get(vpm, source(halfedge(f, mesh), mesh)) << "\n";
    out << get(vpm, target(halfedge(f, mesh), mesh)) << "\n";
    out << get(vpm, target(next(halfedge(f, mesh), mesh), mesh)) << "\n";
  }

  int id = 0;
  for(const face_descriptor f : cc_faces)
  {
    CGAL_USE(f);
    out << "3 " << id << " " << id+1 << " " << id+2 << "\n";
    id += 3;
  }

  out.close();
}

template <typename Point>
void dump_tentative_hole(std::vector<std::vector<Point> >& point_patch,
                         const std::string filename)
{
  std::ofstream out(filename);
  out << std::setprecision(17);

  std::map<Point, int> unique_points_with_id;
  for(const std::vector<Point>& face : point_patch)
    for(const Point& p : face)
      unique_points_with_id.insert(std::make_pair(p, 0));

  out << "OFF\n";
  out << unique_points_with_id.size() << " " << point_patch.size() << " 0\n";

  int unique_id = 0;
  for(auto& pp : unique_points_with_id)
  {
    out << pp.first << "\n";
    pp.second = unique_id++;
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

#endif // CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

// Hole filling can be influenced by setting a third point associated to an edge on the border of the hole.
// This third point is supposed to represent how the mesh continues on the other side of the hole.
// If that edge is a border edge, there is no third point (since the opposite face is the null face).
// Similarly if the edge is an internal sharp edge, we don't really want to use the opposite face because
// there is by definition a strong discontinuity and it might thus mislead the hole filling algorithm.
//
// Rather, we construct an artifical third point that is in the same plane as the face incident to `h`,
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

template <typename TriangleMesh, typename VertexPointMap, typename Point, typename GeomTraits>
bool construct_tentative_hole_patch(std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_border_vertices,
                                    std::set<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_interior_vertices,
                                    std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor>& cc_interior_edges,
                                    const std::vector<Point>& hole_points,
                                    const std::vector<Point>& third_points,
                                    const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                    const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                    std::vector<std::vector<Point> >& point_patch,
                                    const TriangleMesh& tmesh,
                                    VertexPointMap /*vpm*/,
                                    const GeomTraits& /*gt*/)
{
  CGAL_static_assertion((std::is_same<typename boost::property_traits<VertexPointMap>::value_type, Point>::value));

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef CGAL::Triple<int, int, int>                                       Face_indices;

  CGAL_assertion(cc_border_hedges.size() == cc_border_hedges.size());
  CGAL_assertion(hole_points.size() == third_points.size());

  // Collect vertices and edges inside the current selection cc: first collect all vertices and
  // edges incident to the faces to remove...
  for(const face_descriptor f : cc_faces)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      if(halfedge(target(h, tmesh), tmesh) == h) // limit the number of insertions
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

  // try to triangulate the hole using default parameters
  // (using Delaunay search space if CGAL_HOLE_FILLING_DO_NOT_USE_DT3 is not defined)
  std::vector<Face_indices> hole_faces;
  if(hole_points.size() > 3)
    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(hole_faces));
  else
    hole_faces.emplace_back(0, 1, 2); // trivial hole filling

  if(hole_faces.empty())
  {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to fill a hole using Delaunay search space.\n";
 #endif

    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(hole_faces),
                              parameters::use_delaunay_triangulation(false));
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
  std::cout << hole_faces.size() << " faces in the patch" << std::endl;
  std::vector<std::vector<Point> > to_dump;
  for(const Face_indices& face : hole_faces)
  {
    to_dump.emplace_back(std::initializer_list<Point>{hole_points[face.first],
                                                      hole_points[face.second],
                                                      hole_points[face.third]});
  }

  CGAL_assertion(to_dump.size() == hole_faces.size());

  static int hole_id = 0;
  std::stringstream oss;
  oss << "results/tentative_hole_" << hole_id++ << ".off" << std::ends;
  const std::string filename = oss.str().c_str();

  dump_tentative_hole(to_dump, filename);
#endif

  // make sure that the hole filling is valid, we check that no
  // edge already in the mesh is present in hole_faces.
  bool non_manifold_edge_found = false;
  for(const Face_indices& triangle : hole_faces)
  {
    std::array<int, 6> edges = make_array(triangle.first, triangle.second,
                                          triangle.second, triangle.third,
                                          triangle.third, triangle.first);
    for(int k=0; k<3; ++k)
    {
      const int vi = edges[2*k], vj = edges[2*k+1];

      // ignore boundary edges
      if(vi+1 == vj || (vj == 0 && static_cast<std::size_t>(vi) == cc_border_vertices.size()-1))
        continue;

      halfedge_descriptor h = halfedge(cc_border_vertices[vi], cc_border_vertices[vj], tmesh).first;
      if(h != boost::graph_traits<TriangleMesh>::null_halfedge() &&
         cc_interior_edges.count(edge(h, tmesh)) == 0)
      {
        non_manifold_edge_found = true;
        break;
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

  point_patch.reserve(point_patch.size() + hole_faces.size());
  for(const Face_indices& face : hole_faces)
  {
    point_patch.emplace_back(std::initializer_list<Point>{hole_points[face.first],
                                                          hole_points[face.second],
                                                          hole_points[face.third]});
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Found acceptable hole-filling patch.\n";
#endif

  return true;
}

// This function constructs the ranges `hole_points` and `third_points`. Note that for a sub-hole,
// these two ranges are constructed in another function because we don't want to set 'third_points'
// for edges that are on the border of the sub-hole but not on the border of the (full) hole.
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool construct_tentative_hole_patch(std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_border_vertices,
                                    std::set<typename boost::graph_traits<TriangleMesh>::vertex_descriptor>& cc_interior_vertices,
                                    std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor>& cc_interior_edges,
                                    const std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                    const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                    std::vector<std::vector<typename boost::property_traits<VertexPointMap>::value_type> >& patch,
                                    const TriangleMesh& tmesh,
                                    VertexPointMap vpm,
                                    const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

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

  return construct_tentative_hole_patch(cc_border_vertices, cc_interior_vertices, cc_interior_edges,
                                        hole_points, third_points, cc_border_hedges, cc_faces,
                                        patch, tmesh, vpm, gt);
}

// In that overload, we don't know the border of the patch because the face range is a sub-region
// of the hole. We also construct `hole_points` and `third_points`, but with no third point for internal
// sharp edges because a local self-intersection is usually caused by folding and thus we do not want
// a third point resulting from folding to constrain the way we fill the hole in the wrong way.
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool construct_tentative_sub_hole_patch(std::vector<std::vector<typename boost::property_traits<VertexPointMap>::value_type> >& patch,
                                        const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& sub_cc_faces,
                                        const std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
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

  // @todo we don't care about those sets, so instead there could be a system of output iterators
  // in construct_tentative_hole_patch() instead (and here would be emptyset iterators).
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

  return construct_tentative_hole_patch(cc_border_vertices, cc_interior_vertices, cc_interior_edges,
                                        hole_points, third_points, cc_border_hedges, sub_cc_faces,
                                        patch, tmesh, vpm, gt);
}

// -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

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

  // Check for self-intersections
  // Don't know anything better than just making a mesh out of the soup for now...
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
      std::pair<typename std::map<Point, std::size_t>::iterator, bool> is_insert_successful =
        ids.insert(std::make_pair(pt, c));
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

// This function is only called when the hole is NOT subdivided into smaller holes
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
               TriangleMesh& tmesh,
               VertexPointMap vpm,
               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Attempting hole-filling (no constraints), " << cc_faces.size() << " faces\n";
#endif

  if(!order_border_halfedge_range(cc_border_hedges, tmesh))
  {
    CGAL_assertion(false); // we shouldn't fail to orient the boundary cycle of the complete hole
    return false;
  }

  std::set<vertex_descriptor> cc_interior_vertices;
  std::set<edge_descriptor> cc_interior_edges;

  std::vector<vertex_descriptor> cc_border_vertices;
  cc_border_vertices.reserve(cc_border_hedges.size());

  std::vector<std::vector<Point> > patch;
  if(!construct_tentative_hole_patch(cc_border_vertices, cc_interior_vertices, cc_interior_edges,
                                     cc_border_hedges, cc_faces, patch, tmesh, vpm, gt) ||
     !check_patch_sanity<TriangleMesh>(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Failed to find acceptable hole patch\n";
#endif

    return false;
  }

  // Could renew the range directly within the patch replacement function
  // to avoid erasing and re-adding the same face
  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  // Plug the new triangles in the mesh, reusing previous edges and faces
  replace_faces_with_patch(cc_border_vertices, cc_interior_vertices,
                           cc_border_hedges, cc_interior_edges,
                           cc_faces, patch, tmesh, vpm,
                           std::inserter(working_face_range, working_face_range.end()));

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  static int filed_hole_id = 0;
  std::stringstream oss;
  oss << "results/filled_basic_" << filed_hole_id++ << ".off" << std::ends;
  std::ofstream(oss.str().c_str()) << std::setprecision(17) << tmesh;
#endif

  CGAL_postcondition(is_valid_polygon_mesh(tmesh));

  return true;
}

// Same function as above but border of the hole is not known
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
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
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh, vpm, gt);
  else
    return false;
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole_with_constraints(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                TriangleMesh& tmesh,
                                const double dihedral_angle,
                                const double weak_DA,
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

  constrain_sharp_and_border_edges(cc_faces, tmesh, eif, true /*constrain_sharp_edges*/, dihedral_angle, weak_DA, vpm, gt);

  // Partition the hole using these constrained edges
  std::set<face_descriptor> visited_faces;
  std::vector<std::vector<Point> > patch;

  int cc_counter = 0;
  for(face_descriptor f : cc_faces)
  {
    if(!visited_faces.insert(f).second) // already visited that face
      continue;

    // gather the faces of the sub-hole
    std::set<face_descriptor> sub_cc;
    Polygon_mesh_processing::connected_component(f, tmesh, std::inserter(sub_cc, sub_cc.end()),
                                                 CGAL::parameters::edge_is_constrained_map(eif));

    visited_faces.insert(sub_cc.begin(), sub_cc.end());
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "CC of size " << sub_cc.size() << " (total: " << cc_faces.size() << ")" << std::endl;
#endif
    ++cc_counter;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    dump_cc(sub_cc, tmesh, "results/current_cc.off");
#endif

    // The mesh is not modified, but 'patch' gets filled
    if(!construct_tentative_sub_hole_patch(patch, sub_cc, cc_faces, tmesh, vpm, gt))
    {
      // Something went wrong while finding a potential cover for the a sub-hole --> use basic hole-filling
      return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh, vpm, gt);
    }
  }
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << cc_counter << " independent sub holes" << std::endl;
#endif
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream out("results/hole_fillers.off");
  out.precision(17);

  out << "OFF\n";
  out << 3*patch.size() << " " << patch.size() << " 0\n";

  for(const auto& f : patch)
  {
    for(const auto& pt : f)
      out << pt << "\n";
  }

  int id = 0;
  for(std::size_t i=0; i<patch.size(); ++i)
  {
    out << "3 " << id << " " << id+1 << " " << id+2 << "\n";
    id += 3;
  }
  out.close();
#endif

  // We're assembling multiple patches so we could have the same face appearing multiple times...
  if(!check_patch_sanity<TriangleMesh>(patch))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Unhealthy patch, use base fill_hole" << std::endl;
#endif
    return fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh, vpm, gt);
  }

  // Plug the hole-filling patch in the mesh
  std::set<face_descriptor> new_faces;
  replace_faces_with_patch(cc_faces, patch, tmesh, vpm, std::inserter(new_faces, new_faces.end()));

  // Otherwise it should have failed the sanity check
  CGAL_assertion(!does_self_intersect(new_faces, tmesh, parameters::vertex_point_map(vpm)));

  // Update working range with the new faces
  for(const face_descriptor f : cc_faces)
    working_face_range.erase(f);

  working_face_range.insert(new_faces.begin(), new_faces.end());

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

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool remove_self_intersections_with_hole_filling(std::vector<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>& cc_border_hedges,
                                                 std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
                                                 std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                                 TriangleMesh& tmesh,
                                                 bool local_self_intersection_removal,
                                                 const double strong_dihedral_angle,
                                                 const double weak_dihedral_angle,
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

  if(!is_simple_3(cc_border_hedges, tmesh, vpm, gt))
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "Hole filling cannot handle non-simple border" << std::endl;
#endif
    return false;
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_CONSTRAINTS_IN_HOLE_FILLING
  // Do not try to impose sharp edge constraints if we are not doing local-only self intersections removal
  local_self_intersection_removal = false;
#endif

  bool success = false;
  if(local_self_intersection_removal)
  {
    success = fill_hole_with_constraints(cc_border_hedges, cc_faces, working_face_range, tmesh,
                                         strong_dihedral_angle, weak_dihedral_angle, vpm, gt);
  }
  else
  {
    success = fill_hole(cc_border_hedges, cc_faces, working_face_range, tmesh, vpm, gt);
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(success)
  {
    if(local_self_intersection_removal)
      ++self_intersections_solved_by_constrained_hole_filling;
    else
      ++self_intersections_solved_by_unconstrained_hole_filling;
  }
#endif

  return success;
}

// the parameter `step` controls how many extra layers of faces we take around the range `faces_to_remove`
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
std::pair<bool, bool>
remove_self_intersections_one_step(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& faces_to_remove,
                                   std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                   TriangleMesh& tmesh,
                                   const int step,
                                   const bool preserve_genus,
                                   const bool only_treat_self_intersections_locally,
                                   const double strong_dihedral_angle,
                                   const double weak_dihedral_angle,
                                   VertexPointMap vpm,
                                   const GeomTraits& gt)
{
  typedef boost::graph_traits<TriangleMesh>                               graph_traits;
  typedef typename graph_traits::vertex_descriptor                        vertex_descriptor;
  typedef typename graph_traits::halfedge_descriptor                      halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                          face_descriptor;

  std::set<face_descriptor> faces_to_remove_copy = faces_to_remove;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "##### running remove_self_intersections_one_step, step " << step
            << " with " << faces_to_remove.size() << " intersecting faces\n";
#endif

  bool something_was_done = false; // indicates if a region was successfully remeshed
  bool all_fixed = true; // indicates if all removal went well
  // indicates if a removal was not possible because the region handle has
  // some boundary cycle of halfedges
  bool topology_issue = false;

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: is_valid in one_step(tmesh)? ";
  std::cout.flush();

  unsolved_self_intersections = 0;
#endif

  CGAL_precondition(is_valid_polygon_mesh(tmesh));

  while(!faces_to_remove.empty())
  {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: --------------- " << faces_to_remove.size() << " faces to remove (step: " << step << ")\n";
#endif

    // Process a connected component of faces to remove.
    // collect all the faces from the connected component
    std::set<face_descriptor> cc_faces;
    std::vector<face_descriptor> queue(1, *faces_to_remove.begin()); // temporary queue
    cc_faces.insert(queue.back());
    while(!queue.empty())
    {
      face_descriptor top = queue.back();
      queue.pop_back();
      halfedge_descriptor h = halfedge(top, tmesh);
      for(int i=0; i<3; ++i)
      {
        face_descriptor adjacent_face = face(opposite(h, tmesh), tmesh);
        if(adjacent_face!=boost::graph_traits<TriangleMesh>::null_face())
        {
          if(faces_to_remove.count(adjacent_face) != 0 && cc_faces.insert(adjacent_face).second)
            queue.push_back(adjacent_face);
        }

        h = next(h, tmesh);
      }
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in CC\n";
    std::cout << "  DEBUG: first face: " << get(vpm, source(halfedge(*(cc_faces.begin()), tmesh), tmesh)) << " "
              << get(vpm, target(halfedge(*(cc_faces.begin()), tmesh), tmesh)) << " "
              << get(vpm, target(next(halfedge(*(cc_faces.begin()), tmesh), tmesh), tmesh)) << "\n";
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    static int ini_cc_id = 0;
    std::stringstream ini_oss, mesh_oss;
    std::cout << "Output initial CC #" << ini_cc_id << std::endl;
    ini_oss << "results/initial_cc_" << ini_cc_id << ".off" << std::ends;
    dump_cc(cc_faces, tmesh, ini_oss.str().c_str());

    mesh_oss << "results/mesh_at_cc_ " << ini_cc_id++ << ".off" << std::ends;
    std::ofstream mout(mesh_oss.str().c_str());
    mout << std::setprecision(17) << tmesh;
    mout.close();
#endif

    // expand the region to be filled
    if(step > 0)
    {
      expand_face_selection(cc_faces, tmesh, step,
                            make_boolean_property_map(cc_faces),
                            Emptyset_iterator());
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    static int exp_cc_id = 0;
    std::stringstream oss;
    std::cout << "Output expanded CC #" << exp_cc_id << std::endl;
    oss << "results/expanded_cc_" << exp_cc_id++ << ".off" << std::ends;
    dump_cc(cc_faces, tmesh, oss.str().c_str());
#endif

    // try to compactify the selection region by also selecting all the faces included
    // in the bounding box of the initial selection
    std::vector<halfedge_descriptor> stack_for_expension;
    Bbox_3 bb;
    for(face_descriptor fd : cc_faces)
    {
      for(halfedge_descriptor h : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
        bb += get(vpm, target(h, tmesh)).bbox();
        face_descriptor nf = face(opposite(h, tmesh), tmesh);
        if(nf != boost::graph_traits<TriangleMesh>::null_face() && cc_faces.count(nf) == 0)
        {
          stack_for_expension.push_back(opposite(h, tmesh));
        }
      }
    }

    while(!stack_for_expension.empty())
    {
      halfedge_descriptor h = stack_for_expension.back();
      stack_for_expension.pop_back();
      if(cc_faces.count(face(h, tmesh)) == 1)
        continue;

      if(do_overlap(bb, get(vpm, target(next(h, tmesh), tmesh)).bbox()))
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

    if(only_treat_self_intersections_locally)
    {
      if(!does_self_intersect(cc_faces, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt)))
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: No self-intersection in CC\n";
#endif

        for(const face_descriptor f : cc_faces)
          faces_to_remove.erase(f);

        continue;
      }
    }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: " << cc_faces.size() << " faces in expanded CC\n";
#endif

    // remove faces from the set to process
    for(const face_descriptor f : cc_faces)
      faces_to_remove.erase(f);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
    std::stringstream ex_oss;
    std::cout << "Output FULLY expanded CC #" << exp_cc_id-1 << std::endl;
    ex_oss << "results/fully_expanded_cc_" << exp_cc_id-1 << ".off" << std::ends;
    dump_cc(cc_faces, tmesh, ex_oss.str().c_str());
#endif

    //Check for non-manifold vertices in the selection and remove them by selecting all incident faces:
    //  extract the set of halfedges that is on the boundary of the holes to be
    //  made. In addition, we make sure no hole to be created contains a vertex
    //  visited more than once along a hole border (pinched surface)
    //  We save the size of boundary_hedges to make sur halfedges added
    //  from non_filled_hole are not removed.
    bool non_manifold_vertex_remaining_in_selection = false;
    do
    {
      bool non_manifold_vertex_removed = false; //here non-manifold is for the 1D polyline
      std::vector<halfedge_descriptor> boundary_hedges;
      for(face_descriptor fh : cc_faces)
      {
        halfedge_descriptor h = halfedge(fh, tmesh);
        for(int i=0; i<3; ++i)
        {
          if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
            boundary_hedges.push_back(h);

          h = next(h, tmesh);
        }
      }

      // detect vertices visited more than once along
      // a hole border. We then remove all faces incident
      // to such a vertex to force the removal of the vertex.
      // Actually even if two holes are sharing a vertex, this
      // vertex will be removed. It is not needed but since
      // we do not yet have one halfedge per hole it is simpler
      // and does not harm
      std::set<vertex_descriptor> border_vertices;
      for(halfedge_descriptor h : boundary_hedges)
      {
        if(!border_vertices.insert(target(h, tmesh)).second)
        {
          bool any_face_added = false;
          for(halfedge_descriptor hh : halfedges_around_target(h, tmesh))
          {
            if(!is_border(hh, tmesh))
            {
              // add the face to the current selection
              any_face_added |= cc_faces.insert(face(hh, tmesh)).second;
              faces_to_remove.erase(face(hh, tmesh));
            }
          }

          if(any_face_added)
            non_manifold_vertex_removed = true;
          else
            non_manifold_vertex_remaining_in_selection = true;
        }
      }

      if(!non_manifold_vertex_removed)
        break;
    }
    while(true);

    if(preserve_genus && non_manifold_vertex_remaining_in_selection)
    {
      topology_issue = true;
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: CC not handled due to the presence at least one non-manifold vertex\n";

      ++unsolved_self_intersections;
#endif

      continue; // cannot replace a patch containing a nm vertex by a disk
    }

    // before running this function if preserve_genus=false, we duplicated
    // all of them
    CGAL_assertion(!non_manifold_vertex_remaining_in_selection);

    // Collect halfedges on the boundary of the region to be selected
    // (pointing inside the domain to be remeshed)
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

    if(cc_faces.size() == 1) // it is a triangle nothing better can be done
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
     ++unsolved_self_intersections;
#endif
      continue;
    }

    working_face_range.insert(cc_faces.begin(), cc_faces.end());

    // Now, we have a proper selection that we can work on.

#ifndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING
    // First, try to smooth if we only care about local self-intersections
    // Two different approaches:
    // - First, try to constrain edges that are in the zone to smooth and whose dihedral angle is large,
    //   but not too large (we don't want to constrain edges that are created by foldings)
    // - If that fails, try to smooth without any constraints, but make sure that the deviation from
    //   the first zone is small
    //
    // If smoothing fails, the face patch is restored to its pre-smoothing state.
    //
    // Note that there is no need to update the working range because smoothing doesn`t change
    // the number of faces (and old faces are re-used).
    bool fixed_by_smoothing = false;

    if(only_treat_self_intersections_locally)
    {
      fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, tmesh, true,
                                                                    strong_dihedral_angle,
                                                                    weak_dihedral_angle, vpm, gt);

      if(!fixed_by_smoothing) // try again, but without constraining sharp edges
      {
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: Could not be solved via smoothing with constraints\n";
 #endif

        fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, tmesh, false,
                                                                      strong_dihedral_angle,
                                                                      weak_dihedral_angle, vpm, gt);
      }
    }

    if(fixed_by_smoothing)
    {
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "  DEBUG: Solved with smoothing!\n";
 #endif

      something_was_done = true;
      continue;
    }
 #ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    else
    {
      std::cout << "  DEBUG: Could not be solved via smoothing\n";
    }
 #endif
#endif // ndef CGAL_PMP_REMOVE_SELF_INTERSECTIONS_NO_SMOOTHING

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Trying hole-filling based approach\n";
#endif

    if(!is_selection_a_topological_disk(cc_faces, tmesh))
    {
      // check if the selection contains cycles of border halfedges
      bool only_border_edges = true;
      std::set<halfedge_descriptor> mesh_border_hedge;

      for(halfedge_descriptor h : cc_border_hedges)
      {
        if(!is_border(opposite(h, tmesh), tmesh))
          only_border_edges = false;
        else
          mesh_border_hedge.insert(opposite(h, tmesh));
      }

      int nb_cycles = 0;
      while(!mesh_border_hedge.empty())
      {
        // we must count the number of cycle of boundary edges
        halfedge_descriptor h_b = *mesh_border_hedge.begin(), h=h_b;
        mesh_border_hedge.erase(mesh_border_hedge.begin());
        do
        {
          h = next(h, tmesh);
          if(h == h_b)
          {
            // found a cycle
            ++nb_cycles;
            break;
          }
          else
          {
            typename std::set<halfedge_descriptor>::iterator it = mesh_border_hedge.find(h);
            if(it == mesh_border_hedge.end())
              break; // not a cycle

            mesh_border_hedge.erase(it);
          }
        }
        while(true);
      }

      if(nb_cycles > (only_border_edges ? 1 : 0))
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "  DEBUG: CC not handled due to the presence of "
                  << nb_cycles << " of boundary edges\n";
     ++unsolved_self_intersections;
#endif

        topology_issue = true;
        continue;
      }
      else
      {
        if(preserve_genus)
        {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
          std::cout << "  DEBUG: CC not handled because it is not a topological disk (preserve_genus=true)\n";
          ++unsolved_self_intersections;
#endif

          all_fixed = false;
          continue;
        }

        // count the number of cycles of halfedges of the boundary
        std::map<vertex_descriptor, vertex_descriptor> bhs;
        for(halfedge_descriptor h : cc_border_hedges)
          bhs[source(h, tmesh)] = target(h, tmesh);

        int nbc=0;
        while(!bhs.empty())
        {
          ++nbc;
          std::pair<vertex_descriptor, vertex_descriptor > top=*bhs.begin();
          bhs.erase(bhs.begin());

          do
          {
            typename std::map<vertex_descriptor, vertex_descriptor>::iterator it_find = bhs.find(top.second);
            if(it_find == bhs.end())
              break;

            top = *it_find;
            bhs.erase(it_find);
          }
          while(true);
        }

        if(nbc != 1)
        {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
          std::cout << "  DEBUG: CC not handled because it is not a topological disk("
                    << nbc << " boundary cycles)\n";
          ++unsolved_self_intersections;
#endif

          all_fixed = false;
          continue;
        }
        else
        {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
          std::cout << "  DEBUG: CC that is not a topological disk but has only one boundary cycle(preserve_genus=false)\n";
#endif
        }
      }
    }

    if(!remove_self_intersections_with_hole_filling(cc_border_hedges, cc_faces, working_face_range,
                                                    tmesh, only_treat_self_intersections_locally,
                                                    strong_dihedral_angle, weak_dihedral_angle,
                                                    vpm, gt))
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
  }

  if(!something_was_done)
  {
    faces_to_remove.swap(faces_to_remove_copy);
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    std::cout << "  DEBUG: Nothing was changed during this step, self-intersections won`t be recomputed." << std::endl;
#endif
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::stringstream oss;
  oss << "results/after_step_" << step << ".off" << std::ends;
  std::ofstream(oss.str().c_str()) << std::setprecision(17) << tmesh;
#endif

  return std::make_pair(all_fixed, topology_issue);
}

} // namespace internal

namespace experimental {

template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_self_intersections(const FaceRange& face_range,
                               TriangleMesh& tmesh,
                               const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef boost::graph_traits<TriangleMesh>                                 graph_traits;
  typedef typename graph_traits::face_descriptor                            face_descriptor;

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       GeomTraits;
  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));

  bool preserve_genus = choose_parameter(get_parameter(np, internal_np::preserve_genus), true);
  const bool only_treat_self_intersections_locally = choose_parameter(get_parameter(np, internal_np::apply_per_connected_component), false);

  // When treating intersections locally, we don't want to grow the working range too much as
  // either the solution is found fast, or it's too difficult and neither local smoothing or local
  // hole filling are going to provide nice results.
  const int default_max_step = only_treat_self_intersections_locally ? 2 : 7;
  const int max_steps = choose_parameter(get_parameter(np, internal_np::number_of_iterations), default_max_step);

  // @fixme give it its own named parameter rather than abusing 'with_dihedral_angle'?
  const double strong_dihedral_angle = choose_parameter(get_parameter(np, internal_np::with_dihedral_angle), 60.);

  // detect_feature_pp NP (unused for now)
  const double weak_dihedral_angle = 0.; // choose_parameter(get_parameter(np, internal_np::weak_dihedral_angle), 20.);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";
  std::cout << "\tpreserve_genus: " << preserve_genus << std::endl;
  std::cout << "\tonly_treat_self_intersections_locally: " << only_treat_self_intersections_locally << std::endl;
  std::cout << "\tmax_steps: " << max_steps << std::endl;
  std::cout << "\tstrong_dihedral_angle: " << strong_dihedral_angle << std::endl;
  std::cout << "\tweak_dihedral_angle: " << weak_dihedral_angle << std::endl;
#endif

  if(!preserve_genus)
    duplicate_non_manifold_vertices(tmesh, np);

  // Look for self-intersections in the mesh and remove them
  int step = -1;
  bool all_fixed = true; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> faces_to_remove;
  std::set<face_descriptor> working_face_range(face_range.begin(), face_range.end());

  while(++step < max_steps)
  {
    if(faces_to_remove.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;

      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(working_face_range, tmesh, std::back_inserter(self_inter));
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << self_inter.size() << " intersecting pairs" << std::endl;
#endif
      for(const Face_pair& fp : self_inter)
      {
        faces_to_remove.insert(fp.first);
        faces_to_remove.insert(fp.second);
      }
    }

    if(faces_to_remove.empty() && all_fixed)
    {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
      std::cout << "DEBUG: There are no more faces to remove." << std::endl;
#endif
      break;
    }

    std::tie(all_fixed, topology_issue) =
      internal::remove_self_intersections_one_step(
          faces_to_remove, working_face_range, tmesh,
          step, preserve_genus, only_treat_self_intersections_locally,
          strong_dihedral_angle, weak_dihedral_angle, vpm, gt);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    if(all_fixed && topology_issue)
        std::cout << "DEBUG: boundary cycles of boundary edges involved in self-intersections.\n";
#endif
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "solved by constrained smoothing: " << internal::self_intersections_solved_by_constrained_smoothing << std::endl;
  std::cout << "solved by unconstrained smoothing: " << internal::self_intersections_solved_by_unconstrained_smoothing << std::endl;
  std::cout << "solved by constrained hole-filling: " << internal::self_intersections_solved_by_constrained_hole_filling << std::endl;
  std::cout << "solved by unconstrained hole-filling: " << internal::self_intersections_solved_by_unconstrained_hole_filling << std::endl;
  std::cout << "unsolved: " << internal::unsolved_self_intersections << std::endl;
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream("results/final.off") << std::setprecision(17) << tmesh;
#endif

  return step < max_steps;
}

template <typename FaceRange, typename TriangleMesh>
bool remove_self_intersections(const FaceRange& face_range, TriangleMesh& tmesh)
{
  return remove_self_intersections(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_self_intersections(TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  return remove_self_intersections(faces(tmesh), tmesh, np);
}

template <typename TriangleMesh>
bool remove_self_intersections(TriangleMesh& tmesh)
{
  return remove_self_intersections(tmesh, parameters::all_default());
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_SELF_INTERSECTIONS_H
