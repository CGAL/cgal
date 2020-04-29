// Copyright (c) 2018-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_HELPER_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_HELPER_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/iterator.h>
#include <CGAL/utility.h>

#include <array>
#include <iostream>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

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

// @todo these could be extracted to somewhere else, it's useful in itself
template <typename BorderVerticesContainer,
          typename PolygonMesh, typename VPM, typename Point,
          typename FaceOutputIterator>
bool replace_faces_with_patch(const BorderVerticesContainer& border_vertices,
                              const std::set<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& interior_vertices,
                              const std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>& border_hedges,
                              const std::set<typename boost::graph_traits<PolygonMesh>::edge_descriptor>& interior_edges,
                              const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& face_range,
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

  // To be used to create new elements
  std::vector<vertex_descriptor> vertex_stack(interior_vertices.begin(), interior_vertices.end());
  std::vector<edge_descriptor> edge_stack(interior_edges.begin(), interior_edges.end());
  std::vector<face_descriptor> face_stack(face_range.begin(), face_range.end());

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << vertex_stack.size() << " interior vertices" << std::endl;
  std::cout << border_hedges.size() << " border halfedges" << std::endl;
  std::cout << edge_stack.size() << " interior edges" << std::endl;
  std::cout << face_stack.size() << " faces" << std::endl;
#endif

  // Introduce new vertices, convert the patch in vertex patches
  std::vector<Vertex_face> patch_with_vertices;
  patch_with_vertices.reserve(patch.size());

  std::map<Point, vertex_descriptor> point_to_vs;

  // first, add those for which the vertex will not change
  for(const vertex_descriptor v : border_vertices)
    point_to_vs[get(vpm, v)] = v;

  // Do a check to see if we are not introducing any new edge that might be incompatible
  // with a face graph data structure
  for(const Point_face& pf : patch)
  {
    CGAL_assertion(!pf.empty());

    typename Point_face::const_iterator pit = pf.begin(),
                                        pend = pf.end(),
                                        last = std::prev(pf.end());
    for(; pit!=pend; ++pit)
    {
      const Point& p = *pit;
      typename std::map<Point, vertex_descriptor>::iterator vit = point_to_vs.find(p);
      if(vit == point_to_vs.end())
        break;

      // @todo can avoid one find() since it's consecutive edges within the same face
      typename Point_face::const_iterator npit = (pit == last) ? pf.begin() : std::next(pit);
      const Point& np = *npit;
      typename std::map<Point, vertex_descriptor>::iterator nvit = point_to_vs.find(np);
      if(nvit == point_to_vs.end())
        break;

      const std::pair<edge_descriptor, bool> eb = edge(vit->second, nvit->second, pmesh);
      if(!eb.second)
        continue;

      if(is_border(eb.first, pmesh))
        continue;

      // edge is not a border, but one of the faces incident might be up for removal
      const halfedge_descriptor h = halfedge(eb.first, pmesh);
      if(face_range.count(face(h, pmesh)) == 0 && face_range.count(face(opposite(h, pmesh), pmesh)) == 0)
      {
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
        std::cout << "Cannot deal with replacement patch (nm edge)" << std::endl;
#endif
        return false;
      }
    }
  }

  // now build a correspondence map and the faces with vertices
  const vertex_descriptor null_v = boost::graph_traits<PolygonMesh>::null_vertex();
  for(const Point_face& pf : patch)
  {
    Vertex_face vf;
    vf.reserve(pf.size());

    for(const Point& p : pf)
    {
      bool success;
      typename std::map<Point, vertex_descriptor>::iterator it;
      std::tie(it, success) = point_to_vs.insert(std::make_pair(p, null_v));
      vertex_descriptor& v = it->second;

      if(success) // first time we meet that point
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

      vf.push_back(v);
    }

    patch_with_vertices.push_back(vf);
  }

  typedef std::pair<vertex_descriptor, vertex_descriptor>                        Vertex_pair;
  typedef std::map<Vertex_pair, halfedge_descriptor>                             Vertex_pair_halfedge_map;

  Vertex_pair_halfedge_map halfedge_map;

  // register border halfedges
  for(halfedge_descriptor h : border_hedges)
  {
    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    halfedge_map.insert(std::make_pair(std::make_pair(vs, vt), h));
  }

  face_descriptor f = boost::graph_traits<PolygonMesh>::null_face();
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::vector<face_descriptor> new_faces;
#endif

  for(const Vertex_face& vf : patch_with_vertices)
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

    std::vector<halfedge_descriptor> hedges;
    hedges.reserve(vf.size());

    for(std::size_t i=0, n=vf.size(); i<n; ++i)
    {
      const vertex_descriptor vi = vf[i];
      const vertex_descriptor vj = vf[(i+1)%n];

      // get the corresponding halfedge (either a new one or an already created)
      bool success;
      typename Vertex_pair_halfedge_map::iterator it;
      std::tie(it, success) = halfedge_map.insert(std::make_pair(std::make_pair(vi, vj),
                                                  boost::graph_traits<PolygonMesh>::null_halfedge()));
      halfedge_descriptor& h = it->second;

      if(success) // this halfedge is an interior halfedge, not seen through its opposite before
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

        CGAL_assertion(halfedge_map.count(std::make_pair(vj, vi)) == 0);
        halfedge_map[std::make_pair(vj, vi)] = opposite(h, pmesh);

        set_face(opposite(h, pmesh), boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      }

      hedges.push_back(h);
    }

    CGAL_assertion(vf.size() == hedges.size());

    // update halfedge connections + face pointers
    for(std::size_t i=0, n=vf.size(); i<n; ++i)
    {
      set_next(hedges[i], hedges[(i+1)%n], pmesh);
      set_face(hedges[i], f, pmesh);

      set_target(hedges[i], vf[(i+1)%n], pmesh);
      set_halfedge(vf[(i+1)%n], hedges[i], pmesh);
      set_target(opposite(hedges[i], pmesh), vf[i], pmesh);
      set_halfedge(vf[i], opposite(hedges[i], pmesh), pmesh);
    }

    set_halfedge(f, hedges[0], pmesh);

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
    new_faces.push_back(f);
#endif
  }

  // Remove the remaining superfluous vertices, edges, faces
#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << vertex_stack.size() << " vertices to remove" << std::endl;
  std::cout << edge_stack.size() << " edges to remove" << std::endl;
  std::cout << face_stack.size() << " faces to remove" << std::endl;
#endif

  for(vertex_descriptor v : vertex_stack)
    remove_vertex(v, pmesh);
  for(edge_descriptor e : edge_stack)
    remove_edge(e, pmesh);
  for(face_descriptor f : face_stack)
    remove_face(f, pmesh);

  // If 'border_hedges' describes a topological circle, which is possible if we don't want
  // to conform on true border halfedges (in the sense of halfedges that are on the border
  // of not only the patch, but also of the mesh), we must link these new true border halfedges
  // by turning around the vertex in the interior of the mesh
  const halfedge_descriptor null_h = boost::graph_traits<PolygonMesh>::null_halfedge();
  for(const auto& p : halfedge_map)
  {
    const halfedge_descriptor h = p.second;
    if(is_border(h, pmesh))
      set_next(h, null_h, pmesh);
  }

  for(const auto& p : halfedge_map)
  {
    halfedge_descriptor h = p.second;
    if(!is_border(h, pmesh))
      continue;

    // forced to go both way because we need to fix the 'next' also for one halfedge
    // which is not part of the patch
    halfedge_descriptor prev_h = opposite(h, pmesh);
    do
    {
      prev_h = opposite(next(prev_h, pmesh), pmesh);
    }
    while(!is_border(prev_h, pmesh));
    set_next(prev_h, h, pmesh);

    if(next(h, pmesh) != null_h) // already done the 'next' as part of some halfedge's 'prev'
      continue;

    prev_h = h;
    h = opposite(prev_h, pmesh);
    do
    {
      h = opposite(prev(h, pmesh), pmesh);
    }
    while(!is_border(h, pmesh));
    set_next(prev_h, h, pmesh);
  }

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
  std::ofstream res_out("results/last_patch_replacement.off");
  res_out << std::setprecision(17) << pmesh;
  res_out.close();
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  std::cout << "  DEBUG: Replacing range with patch: ";
  std::cout << face_range.size() << " triangles removed, " << new_faces.size() << " created\n";
  CGAL_assertion(new_faces.size() == patch.size());
#endif

#ifdef CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
  if(Polygon_mesh_processing::does_self_intersect(new_faces, pmesh))
    std::cout << "!! NEW FACES SELF INTERSECT !!" << std::endl;
#endif

  return true;
}

template <typename PolygonMesh, typename VPM, typename Point, typename FaceOutputIterator>
bool replace_faces_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& face_range,
                              const std::vector<std::vector<Point> >& patch,
                              PolygonMesh& pmesh,
                              VPM& vpm,
                              FaceOutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor       face_descriptor;

  CGAL_assertion(!face_range.empty());

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

      // The patch has to conform if there is a face on the opposite side, but not necessarily to the border
      if(!is_border(opp_h, pmesh) && face_range.count(opp_f) == 0)
      {
        border_hedges.push_back(h);

        interior_vertices.erase(source(h, pmesh));
        border_vertices.insert(source(h, pmesh));
        interior_vertices.erase(target(h, pmesh));
        border_vertices.insert(target(h, pmesh));
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
bool replace_faces_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& face_range,
                              const std::vector<std::vector<Point> >& patch,
                              PolygonMesh& pmesh,
                              VPM& vpm)
{
  CGAL::Emptyset_iterator out;
  return replace_faces_with_patch(face_range, patch, pmesh, vpm, out);
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
                                     cc_border_hedges, cc_faces, patch, tmesh, vpm, gt))
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
  bool success = replace_faces_with_patch(cc_border_vertices, cc_interior_vertices,
                                          cc_border_hedges, cc_interior_edges,
                                          cc_faces, patch, tmesh, vpm,
                                          std::inserter(working_face_range, working_face_range.end()));
  if(!success)
    return false;

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

// Same as above, but we don't care about maintaining the working face range
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool fill_hole(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& cc_faces,
               TriangleMesh& tmesh,
               VertexPointMap vpm,
               const GeomTraits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  std::set<face_descriptor> unused_working_face_range;
  return fill_hole(cc_faces, unused_working_face_range, tmesh, vpm, gt);
}

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_HELPER_H
