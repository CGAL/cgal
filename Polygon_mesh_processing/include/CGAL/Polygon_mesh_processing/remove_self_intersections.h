// Copyright (c) 2015-2019 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMOVE_SELF_INTERSECTIONS_H
#define CGAL_POLYGON_MESH_PROCESSING_REMOVE_SELF_INTERSECTIONS_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/remove_degeneracies.h> // only for the preconditions
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/utility.h>

#include <array>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename PolygonMesh, typename VPM, typename Point, typename FaceOutputIterator>
FaceOutputIterator replace_face_range_with_patch(const std::vector<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& border_vertices,
                                                 const std::set<typename boost::graph_traits<PolygonMesh>::vertex_descriptor>& interior_vertices,
                                                 const std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor>& border_hedges,
                                                 const std::set<typename boost::graph_traits<PolygonMesh>::edge_descriptor>& interior_edges,
                                                 const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& pol_faces,
                                                 const std::vector<std::vector<Point> >& patch,
                                                 PolygonMesh& pmesh,
                                                 VPM& vpm,
                                                 FaceOutputIterator out,
                                                 const bool verbose)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor        edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor        face_descriptor;

  CGAL_static_assertion((std::is_same<typename boost::property_traits<VPM>::value_type, Point>::value));

  typedef std::vector<Point>                                                Point_face;
  typedef std::vector<vertex_descriptor>                                    Vertex_face;

  // To be used to create new elements
  std::vector<vertex_descriptor> vertex_stack(interior_vertices.begin(), interior_vertices.end());
  std::vector<edge_descriptor> edge_stack(interior_edges.begin(), interior_edges.end());
  std::vector<face_descriptor> face_stack(pol_faces.begin(), pol_faces.end());

  // Introduce new vertices, convert the patch in vertex patches
  std::vector<Vertex_face> vertex_patch;
  vertex_patch.reserve(patch.size());

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

      if(success) // first time we meet that point, means it's an interior point and we need to make a new vertex
      {
        if(vertex_stack.empty())
        {
          v = add_vertex(pmesh);;
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

    vertex_patch.push_back(vface);
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

  for(const Vertex_face& vface : vertex_patch)
  {
    if(face_stack.empty())
    {
      f = add_face(pmesh);
      *out++ = f;
    }
    else
    {
      f = face_stack.back();
      face_stack.pop_back();
      *out++ = f;
    }

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
    for(int i=0, n=vface.size(); i<n; ++i)
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

  if(verbose)
  {
    std::cout << "  DEBUG: Replacing range with patch: ";
    std::cout << pol_faces.size() << " triangles removed, " << patch.size() << " created\n";
  }

  CGAL_postcondition(pmesh.is_valid());
  CGAL_postcondition(is_valid_polygon_mesh(pmesh));

  return out;
}

template <typename PolygonMesh, typename VPM, typename Point, typename FaceOutputIterator>
FaceOutputIterator replace_face_range_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                                                 const std::vector<std::vector<Point> >& patch,
                                                 PolygonMesh& pmesh,
                                                 VPM& vpm,
                                                 FaceOutputIterator out,
                                                 const bool verbose)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor       face_descriptor;

  std::vector<vertex_descriptor> border_vertices;
  std::set<vertex_descriptor> interior_vertices;
  std::vector<halfedge_descriptor> border_hedges;
  std::set<edge_descriptor> interior_edges;

  for(face_descriptor fh : faces)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(fh, pmesh), pmesh))
    {
      if(halfedge(target(h, pmesh), pmesh) == h) // limit the number of insertions
        interior_vertices.insert(target(h, pmesh));
    }
  }

  for(face_descriptor fh : faces)
  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(fh, pmesh), pmesh))
    {
      CGAL_assertion(!is_border(h, pmesh));

      const edge_descriptor e = edge(h, pmesh);
      const halfedge_descriptor opp_h = opposite(h, pmesh);
      const face_descriptor opp_f = face(opp_h, pmesh);

      if(is_border(opp_h, pmesh) || faces.count(opp_f) == 0)
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

  return replace_face_range_with_patch(border_vertices, interior_vertices,
                                       border_hedges, interior_edges, faces, patch,
                                       pmesh, vpm, out, verbose);
}

template <typename PolygonMesh, typename VPM, typename Point>
void replace_face_range_with_patch(const std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& faces,
                                   const std::vector<std::vector<Point> >& patch,
                                   PolygonMesh& pmesh,
                                   VPM& vpm,
                                   const bool verbose = false)
{
  CGAL::Emptyset_iterator out;
  replace_face_range_with_patch(faces, patch, pmesh, vpm, out, verbose);
}

template <typename TriangleMesh>
bool is_new_patch_close_enough_to_initial_patch(const TriangleMesh& old_patch,
                                                const TriangleMesh& new_patch)
{
#if defined(CGAL_LINKED_WITH_TBB)
  typedef CGAL::Parallel_tag                                        Tag
#else
  typedef CGAL::Sequential_tag                                      Tag;
#endif

  double d = Polygon_mesh_processing::approximate_Hausdorff_distance<Tag>( // @todo should it be symmetric?
               new_patch, old_patch, parameters::number_of_points_on_edges(10)
                                                .number_of_points_on_faces(100));

  std::cout << "d: " << d << std::endl;
  return (d <= 0.1); // @fixme hardcoded, must depend on local edge length
}

template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
bool remove_self_intersections_with_smoothing(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& faces,
                                              const bool constrain_sharp_edges,
                                              TriangleMesh& tmesh,
                                              VertexPointMap& vpmap,
                                              const GeomTraits& gt,
                                              const bool verbose)
{
  if(verbose)
  {
    std::cout << "  DEBUG: repair with smoothing... (constraining sharp edges: ";
    std::cout << std::boolalpha << constrain_sharp_edges << ")" << std::endl;
  }

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;

  typedef typename GeomTraits::FT                                           FT;
  typedef typename boost::property_traits<VertexPointMap>::value_type       Point;
  typedef typename GeomTraits::Vector_3                                     Vector;

  typedef CGAL::Face_filtered_graph<TriangleMesh>                           Filtered_graph;

  CGAL_precondition(does_self_intersect(faces, tmesh));

  // keep in memory the face patch so that we can restore if smoothing fails
  std::vector<std::vector<Point> > patch;
  patch.reserve(faces.size());
  for(const face_descriptor f : faces)
  {
    std::vector<Point> face_points;
    for(const halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
      face_points.push_back(get(vpmap, target(h, tmesh)));

    patch.push_back(face_points);
  }

  // also extract the full mesh because we need to compute the Hausdorff distance
  Filtered_graph ffg(tmesh, faces);
  TriangleMesh patch_mesh;
  CGAL::copy_face_graph(ffg, patch_mesh, parameters::vertex_point_map(vpmap));

  typedef typename boost::property_map<TriangleMesh, CGAL::edge_is_feature_t>::type  EIFMap;
  EIFMap eif = get(CGAL::edge_is_feature, tmesh);

  std::map<edge_descriptor, bool> is_border_of_selection;
  for(face_descriptor f : faces)
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      // meet it once --> switch to 'true'; meet it twice --> switch back to 'false'
      const edge_descriptor e = edge(h, tmesh);
      is_border_of_selection[e] = !(is_border_of_selection[e]);
    }
  }

  const double bound = 60.; // @fixme hardcoded
  const double cos_angle = std::cos(bound * CGAL_PI / 180.);

  for(const auto& ep : is_border_of_selection)
  {
    bool flag = ep.second;
    if(constrain_sharp_edges && !flag)
    {
      const halfedge_descriptor h = halfedge(ep.first, tmesh);
      const face_descriptor f1 = face(h, tmesh);
      const face_descriptor f2 = face(opposite(h, tmesh), tmesh);

      const Vector n1 = compute_face_normal(f1, tmesh, parameters::vertex_point_map(vpmap).geom_traits(gt));
      const Vector n2 = compute_face_normal(f2, tmesh, parameters::vertex_point_map(vpmap).geom_traits(gt));
      const FT cos = n1 * n2;

      // do not mark as sharp edges with a dihedral angle that is almost 'pi' because this is likely
      // due to a foldness on the mesh rather than a sharp edge that we wish to preserve
      // (Ideally this would be pre-treated as part of the flatness treatment)
      flag = (cos <= cos_angle && cos >= -cos_angle);
    }

    put(eif, ep.first, flag);
  }

  // @todo choice of number of iterations? Till convergence && max of 100?
  smooth_mesh(faces, tmesh, CGAL::parameters::edge_is_constrained_map(eif)
                                             .number_of_iterations(10)
                                             .use_safety_constraints(false));

  Filtered_graph ffg_post(tmesh, faces);
  TriangleMesh patch_mesh_post;
  CGAL::copy_face_graph(ffg_post, patch_mesh_post, parameters::vertex_point_map(vpmap));

  bool is_acceptable = (!does_self_intersect(faces, tmesh) &&
                        is_new_patch_close_enough_to_initial_patch(patch_mesh, patch_mesh_post));
  if(!is_acceptable)
    replace_face_range_with_patch(faces, patch, tmesh, vpmap, verbose);

  return is_acceptable;
}

// the parameter 'step' controls how many extra layers of faces we take around the range 'faces_to_remove'
template <typename TriangleMesh, typename VertexPointMap, typename GeomTraits>
std::pair<bool, bool>
remove_self_intersections_one_step(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& faces_to_remove,
                                   std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor>& working_face_range,
                                   TriangleMesh& tmesh,
                                   VertexPointMap& vpmap,
                                   const GeomTraits& gt,
                                   const int step,
                                   const bool preserve_genus,
                                   const bool only_treat_self_intersections_locally,
                                   const bool verbose)
{
  typedef boost::graph_traits<TriangleMesh>                               graph_traits;
  typedef typename graph_traits::vertex_descriptor                        vertex_descriptor;
  typedef typename graph_traits::halfedge_descriptor                      halfedge_descriptor;
  typedef typename graph_traits::edge_descriptor                          edge_descriptor;
  typedef typename graph_traits::face_descriptor                          face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type     Point;

  std::set<face_descriptor> faces_to_remove_copy = faces_to_remove;

  if(verbose)
  {
    std::cout << "DEBUG: running remove_self_intersections_one_step, step " << step
              << " with " << faces_to_remove.size() << " intersecting faces\n";
  }

  CGAL_assertion(tmesh.is_valid());

  bool something_was_done = false; // indicates if a region was successfully remeshed
  bool all_fixed = true; // indicates if all removal went well
  // indicates if a removal was not possible because the region handle has
  // some boundary cycle of halfedges
  bool topology_issue = false;

  if(verbose)
  {
    std::cout << "  DEBUG: is_valid in one_step(tmesh)? ";
    std::cout.flush();
    std::cout << is_valid_polygon_mesh(tmesh) << "\n";
    std::cout << "  DEBUG: " << faces_to_remove.empty() << " faces to remove" << std::endl;
  }

  if(!faces_to_remove.empty())
  {
    while(!faces_to_remove.empty())
    {
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

      if(verbose)
        std::cout << "  DEBUG: " << cc_faces.size() << " faces in CC\n";

      // expand the region to be filled
      if(step > 0)
      {
        expand_face_selection(cc_faces, tmesh, step,
                              make_boolean_property_map(cc_faces),
                              Emptyset_iterator());
      }

      // try to compactify the selection region by also selecting all the faces included
      // in the bounding box of the initial selection
      std::vector<halfedge_descriptor> stack_for_expension;
      Bbox_3 bb;
      for(face_descriptor fd : cc_faces)
      {
        for(halfedge_descriptor h : halfedges_around_face(halfedge(fd, tmesh), tmesh))
        {
          bb += get(vpmap, target(h, tmesh)).bbox();
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

        if(do_overlap(bb, get(vpmap, target(next(h, tmesh), tmesh)).bbox()))
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
        if(!does_self_intersect(cc_faces, tmesh, parameters::vertex_point_map(vpmap).geom_traits(gt)))
        {
          if(verbose)
            std::cout << "  DEBUG: No self-intersection in CC\n";

          for(face_descriptor f : cc_faces)
            faces_to_remove.erase(f);

          continue;
        }
      }

      if(verbose)
        std::cout << "  DEBUG: " << cc_faces.size() << " faces in expanded CC\n";

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
        if(verbose)
          std::cout << "  DEBUG: CC not handled due to the presence at least one non-manifold vertex\n";

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
        for(int i=0; i<3;++i)
        {
          if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh)) == 0)
            cc_border_hedges.push_back(h);

          h = next(h, tmesh);
        }
      }

      if(cc_faces.size() == 1) // it is a triangle nothing better can be done
        continue;

      working_face_range.insert(cc_faces.begin(), cc_faces.end());

      // Now, we have a proper selection that we can work on.

      // First, try to smooth if we only care about local self-intersections
      // Two different approaches:
      // - First, try to constrain edges that are in the zone to smooth and whose dihedral angle is large,
      //   but not too large (we don't want to constrain edges that are created by foldings)
      // - If that fails, try to smooth without any constraints, but make sure that the deviation from
      //   the first zone is small
      //
      // If smoothing fails, geometry / combinatoris are restored to pre-smoothing state.
      //
      // Note that there is no need to update the working range because smoothing doesn't change
      // the number of faces (and old faces are re-used).
      bool fixed_by_smoothing = false;
      if(only_treat_self_intersections_locally)
      {
        fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, true, tmesh, vpmap, gt, verbose);
        if(!fixed_by_smoothing) // try again, but without constraining sharp edges
          fixed_by_smoothing = remove_self_intersections_with_smoothing(cc_faces, false, tmesh, vpmap, gt, verbose);
      }

      std::ofstream out2("results/post_fixup.off");
      out2 << std::setprecision(17) << tmesh;
      out2.close();

      // remove faces from the set to process
      for(face_descriptor f : cc_faces)
        faces_to_remove.erase(f);

      if(fixed_by_smoothing)
      {
        something_was_done = true;
        continue;
      }

      if(verbose)
        std::cout << "  DEBUG: Could not be solved via smoothing, trying hole-filling based approach\n";

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

        int nb_cycles=0;
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
          if(verbose)
            std::cout << "  DEBUG: CC not handled due to the presence of  "
                      << nb_cycles << " of boundary edges\n";

          topology_issue = true;
          continue;
        }
        else
        {
          if(preserve_genus)
          {
            if(verbose)
              std::cout << "  DEBUG: CC not handled because it is not a topological disk (preserve_genus=true)\n";

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

          if(nbc!=1)
          {
            if(verbose)
              std::cout << "  DEBUG: CC not handled because it is not a topological disk("
                        << nbc << " boundary cycles)\n";

            all_fixed = false;
            continue;
          }
          else
          {
            if(verbose)
              std::cout << "  DEBUG: CC that is not a topological disk but has only one boundary cycle(preserve_genus=false)\n";
          }
        }
      }

      // sort halfedges so that they describe the sequence
      // of halfedges of the hole to be made
      CGAL_assertion(cc_border_hedges.size() > 2);
      for(std::size_t i=0; i < cc_border_hedges.size()-2; ++i)
      {
        vertex_descriptor tgt = target(cc_border_hedges[i], tmesh);
        for(std::size_t j=i+1; j<cc_border_hedges.size(); ++j)
        {
          if(tgt == source(cc_border_hedges[j], tmesh))
          {
            std::swap(cc_border_hedges[i+1], cc_border_hedges[j]);
            break;
          }

          CGAL_assertion(j!=cc_border_hedges.size()-1);
        }
      }

      CGAL_assertion(source(cc_border_hedges.front(), tmesh) == target(cc_border_hedges.back(), tmesh));

      // collect vertices and edges inside the current selection cc
      std::set<vertex_descriptor> cc_interior_vertices;
      std::set<edge_descriptor>  cc_interior_edges;

      // first collect all vertices and edges incident to the faces to remove
      for(face_descriptor fh : cc_faces)
      {
        for(halfedge_descriptor h : halfedges_around_face(halfedge(fh, tmesh), tmesh))
        {
          if(halfedge(target(h, tmesh), tmesh) == h) // limit the number of insertions
            cc_interior_vertices.insert(target(h, tmesh));

          cc_interior_edges.insert(edge(h, tmesh));
        }
      }

      // and then remove those on the boundary
      for(halfedge_descriptor h : cc_border_hedges)
      {
        cc_interior_vertices.erase(target(h, tmesh));
        cc_interior_edges.erase(edge(h, tmesh));
      }

      if(verbose)
      {
        std::cout << "  DEBUG: is_valid(tmesh) in one_step, before mesh changes? ";
        std::cout << is_valid_polygon_mesh(tmesh) << std::endl;
      }

      //try hole_filling.
      typedef CGAL::Triple<int, int, int>                                       Face_indices;
      typedef typename boost::property_traits<VertexPointMap>::value_type       Point;

      std::vector<Point> hole_points, third_points;
      hole_points.reserve(cc_border_hedges.size());
      third_points.reserve(cc_border_hedges.size());
      std::vector<vertex_descriptor> cc_border_vertices;

      for(halfedge_descriptor h : cc_border_hedges)
      {
        vertex_descriptor v = source(h, tmesh);
        hole_points.push_back(get(vpmap, v));
        cc_border_vertices.push_back(v);
        third_points.push_back(get(vpmap, target(next(opposite(h, tmesh), tmesh), tmesh))); // TODO fix me for mesh border edges
      }

      CGAL_assertion(hole_points.size() >= 3);

      // try to triangulate the hole using default parameters
      //(using Delaunay search space if CGAL_HOLE_FILLING_DO_NOT_USE_DT3 is not defined)
      std::vector<Face_indices> patch;
      if(hole_points.size()>3)
        triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch));
      else
        patch.push_back(Face_indices(0,1,2)); // trivial hole filling

      if(patch.empty())
      {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
        if(verbose)
          std::cout << "  DEBUG: Failed to fill a hole using Delaunay search space.\n";

        triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch),
                                  parameters::use_delaunay_triangulation(false));
#endif
        if(patch.empty())
        {
          if(verbose)
            std::cout << "  DEBUG: Failed to fill a hole using the whole search space.\n";
          all_fixed = false;
          continue;
        }
      }

      // make sure that the hole filling is valid, we check that no
      // edge already in the mesh is present in patch.
      bool non_manifold_edge_found = false;
      for(const Face_indices& triangle : patch)
      {
        std::array<int, 6> edges = make_array(triangle.first, triangle.second,
                                              triangle.second, triangle.third,
                                              triangle.third, triangle.first);
        for(int k=0; k<3; ++k)
        {
          int vi=edges[2*k], vj=edges[2*k+1];
          // ignore boundary edges
          if(vi+1==vj || (vj==0 && static_cast<std::size_t>(vi) == cc_border_vertices.size()-1))
            continue;

          halfedge_descriptor h = halfedge(cc_border_vertices[vi], cc_border_vertices[vj], tmesh).first;
          if(h != boost::graph_traits<TriangleMesh>::null_halfedge() &&
              cc_interior_edges.count(edge(h, tmesh)) == 0)
          {
            non_manifold_edge_found=true;
            break;
          }
        }

        if(non_manifold_edge_found)
          break;
      }

      if(non_manifold_edge_found)
      {
        if(verbose)
          std::cout << "  DEBUG: Triangulation produced is non-manifold when plugged into the mesh.\n";

        all_fixed = false;
        continue;
      }

      if(verbose)
        std::cout << "  DEBUG: Plugging new triangles in...\n";

      // plug the new triangles in the mesh, reusing previous edges and faces
      std::vector<std::vector<Point> > point_patch;
      point_patch.reserve(patch.size());
      for(const Face_indices& face : patch)
      {
        std::vector<Point> point_face = { hole_points[face.first],
                                          hole_points[face.second],
                                          hole_points[face.third] };
        point_patch.push_back(point_face);
      }

      // Could renew the range directly within the patch replacement function
      // to avoid erasing and re-adding the same face
      for(const face_descriptor f : cc_faces)
        working_face_range.erase(f);

      replace_face_range_with_patch(cc_border_vertices, cc_interior_vertices,
                                    cc_border_hedges, cc_interior_edges,
                                    cc_faces, point_patch, tmesh, vpmap,
                                    std::inserter(working_face_range, working_face_range.end()),
                                    verbose);

      CGAL_postcondition(is_valid_polygon_mesh(tmesh));

      something_was_done = true;
    }
  }

  if(!something_was_done)
  {
    faces_to_remove.swap(faces_to_remove_copy);
    if(verbose)
      std::cout << "  DEBUG: Nothing was changed during this step, self-intersections won't be recomputed." << std::endl;
  }

  return std::make_pair(all_fixed, topology_issue);
}

} // namespace internal

/// \cond SKIP_IN_MANUAL

template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_self_intersections(const FaceRange& face_range,
                               TriangleMesh& tmesh,
                               const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh>                                 graph_traits;
  typedef typename graph_traits::face_descriptor                            face_descriptor;

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                           get_property_map(vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       GeomTraits;
  GeomTraits gt = boost::choose_param(boost::get_param(np, internal_np::geom_traits), GeomTraits());

  const int max_steps = boost::choose_param(boost::get_param(np, internal_np::number_of_iterations), 7);
  bool verbose = boost::choose_param(boost::get_param(np, internal_np::verbosity_level), 0) > 0;
  bool preserve_genus = boost::choose_param(boost::get_param(np, internal_np::preserve_genus), true);
  verbose = true;
  const bool only_treat_self_intersections_locally = true;

  std::set<face_descriptor> working_face_range(face_range.begin(), face_range.end());

  if(verbose)
    std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";

  CGAL_precondition_code(std::set<face_descriptor> degenerate_face_set;)
  CGAL_precondition_code(degenerate_faces(working_face_range, tmesh, std::inserter(degenerate_face_set, degenerate_face_set.begin()), np);)
  CGAL_precondition(degenerate_face_set.empty());

  if(!preserve_genus)
    duplicate_non_manifold_vertices(tmesh, np);

  // Look for self-intersections in the mesh and remove them
  int step = -1;
  bool all_fixed = true; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> faces_to_remove;

  while(++step < max_steps)
  {
    if(faces_to_remove.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;
      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(working_face_range, tmesh, std::back_inserter(self_inter));

      for(const Face_pair& fp : self_inter)
      {
        faces_to_remove.insert(fp.first);
        faces_to_remove.insert(fp.second);
      }
    }

    if(faces_to_remove.empty() && all_fixed)
    {
      if(verbose)
        std::cout << "DEBUG: There is no more face to remove." << std::endl;
      break;
    }

    std::tie(all_fixed, topology_issue) =
      internal::remove_self_intersections_one_step(
          faces_to_remove, working_face_range, tmesh, vpm, gt,
          step, preserve_genus, only_treat_self_intersections_locally, verbose);

    if(all_fixed && topology_issue)
    {
      if(verbose)
        std::cout << "DEBUG: boundary cycles of boundary edges involved in self-intersections.\n";
    }
  }

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

/// \endcond

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REMOVE_SELF_INTERSECTIONS_H
