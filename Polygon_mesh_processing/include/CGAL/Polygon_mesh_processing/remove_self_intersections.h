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

#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/remove_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/utility.h>

#include <array>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

/// \cond SKIP_IN_MANUAL
template <typename TriangleMesh, typename face_descriptor, typename VertexPointMap>
std::pair<bool, bool>
remove_self_intersections_one_step(TriangleMesh& tmesh,
                                   std::set<face_descriptor>& faces_to_remove,
                                   VertexPointMap& vpmap,
                                   int step,
                                   bool preserve_genus,
                                   bool verbose)
{
  std::set<face_descriptor> faces_to_remove_copy = faces_to_remove;

  if(verbose)
  {
    std::cout << "DEBUG: running remove_self_intersections_one_step, step " << step
              << " with " << faces_to_remove.size() << " intersecting faces\n";
  }

  if(!does_self_intersect(faces_to_remove, tmesh, parameters::vertex_point_map(vpmap)))
  {
    if(verbose)
      std::cout << "DEBUG: Range has no self-intersections\n";

    return std::make_pair(true, false);
  }

  CGAL_assertion(tmesh.is_valid());

  typedef boost::graph_traits<TriangleMesh>                       graph_traits;
  typedef typename graph_traits::vertex_descriptor                vertex_descriptor;
  typedef typename graph_traits::halfedge_descriptor              halfedge_descriptor;
  typedef typename graph_traits::edge_descriptor                  edge_descriptor;

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

      // remove faces from the set to process
      for(face_descriptor f : cc_faces)
        faces_to_remove.erase(f);

      if(cc_faces.size() == 1)
        continue; // it is a triangle nothing better can be done

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
          if(is_border(opposite(h, tmesh), tmesh) || cc_faces.count(face(opposite(h, tmesh), tmesh))== 0)
            cc_border_hedges.push_back(h);

          h = next(h, tmesh);
        }
      }

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
              typename std::map<vertex_descriptor, vertex_descriptor>::iterator
                  it_find = bhs.find(top.second);
              if(it_find == bhs.end()) break;
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
      std::vector<vertex_descriptor> border_vertices;

      for(halfedge_descriptor h : cc_border_hedges)
      {
        vertex_descriptor v = source(h, tmesh);
        hole_points.push_back(get(vpmap, v));
        border_vertices.push_back(v);
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
          if(vi+1==vj || (vj==0 && static_cast<std::size_t>(vi)==border_vertices.size()-1))
            continue;

          halfedge_descriptor h = halfedge(border_vertices[vi], border_vertices[vj], tmesh).first;
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

      // plug the new triangles in the mesh, reusing previous edges and faces
      std::vector<edge_descriptor> edge_stack(cc_interior_edges.begin(), cc_interior_edges.end());
      std::vector<face_descriptor> face_stack(cc_faces.begin(), cc_faces.end());

      std::map<std::pair<int, int>, halfedge_descriptor> halfedge_map;
      int i = 0;

      // register border halfedges
      for(halfedge_descriptor h : cc_border_hedges)
      {
        int j = static_cast<int>(std::size_t(i+1)%cc_border_hedges.size());
        halfedge_map.insert(std::make_pair(std::make_pair(i, j), h));
        set_halfedge(target(h, tmesh), h, tmesh); // update vertex halfedge pointer
        CGAL_assertion(border_vertices[i] == source(h, tmesh) && border_vertices[j] == target(h, tmesh));
        ++i;
      }

      std::vector<halfedge_descriptor> hedges;
      hedges.reserve(4);
      face_descriptor f = boost::graph_traits<TriangleMesh>::null_face();
      for(const Face_indices& triangle : patch)
      {
        // get the new face
        if(face_stack.empty())
        {
          f = add_face(tmesh);
        }
        else
        {
          f = face_stack.back();
          face_stack.pop_back();
        }

        std::array<int, 4> indices = make_array(triangle.first, triangle.second, triangle.third, triangle.first);
        for(int i=0; i<3; ++i)
        {
          // get the corresponding halfedge (either a new one or an already created)
          typename std::map<std::pair<int, int>, halfedge_descriptor >::iterator insert_res =
              halfedge_map.insert(std::make_pair(std::make_pair(indices[i], indices[i+1]),
                                  boost::graph_traits<TriangleMesh>::null_halfedge())).first;

          if(insert_res->second == boost::graph_traits<TriangleMesh>::null_halfedge())
          {
            if(edge_stack.empty())
            {
              insert_res->second = halfedge(add_edge(tmesh), tmesh);
            }
            else
            {
              insert_res->second = halfedge(edge_stack.back(), tmesh);
              edge_stack.pop_back();
            }

            halfedge_map[std::make_pair(indices[i+1], indices[i])] = opposite(insert_res->second, tmesh);
          }

          hedges.push_back(insert_res->second);
        }

        hedges.push_back(hedges.front());

        // update halfedge connections + face pointers
        for(int i=0; i<3;++i)
        {
          set_next(hedges[i], hedges[i+1], tmesh);
          set_face(hedges[i], f, tmesh);
          set_target(hedges[i], border_vertices[indices[i+1]], tmesh);
        }

        set_halfedge(f, hedges[0], tmesh);
        hedges.clear();
      }

      // now remove remaining edges,
      for(edge_descriptor e : edge_stack)
        remove_edge(e, tmesh);

      // vertices,
      for(vertex_descriptor vh : cc_interior_vertices)
        remove_vertex(vh, tmesh);

      // and remaning faces
      for(face_descriptor f : face_stack)
        remove_face(f, tmesh);

      if(verbose)
        std::cout << "  DEBUG: " << cc_faces.size() << " triangles removed, " << patch.size() << " created\n";

      CGAL_assertion(is_valid_polygon_mesh(tmesh));

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

template <typename TriangleMesh, typename NamedParameters>
bool remove_self_intersections(TriangleMesh& tmesh,
                               const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh>                                 graph_traits;
  typedef typename graph_traits::face_descriptor                            face_descriptor;

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                           get_property_map(vertex_point, tmesh));

  const int max_steps = boost::choose_param(boost::get_param(np, internal_np::number_of_iterations), 7);
  bool verbose = boost::choose_param(boost::get_param(np, internal_np::verbosity_level), 0) > 0;
  bool preserve_genus = boost::choose_param(boost::get_param(np, internal_np::preserve_genus), true);

  if(verbose)
    std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";

  CGAL_precondition_code(std::set<face_descriptor> degenerate_face_set;)
  CGAL_precondition_code(degenerate_faces(tmesh, std::inserter(degenerate_face_set, degenerate_face_set.begin()), np);)
  CGAL_precondition(degenerate_face_set.empty());

  if(!preserve_genus)
    duplicate_non_manifold_vertices(tmesh, np);

  // Look for self-intersections in the mesh and remove them
  int step = -1;
  bool all_fixed = true; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> faces_to_remove;

  while(++step<max_steps)
  {
    if(faces_to_remove.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;
      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(tmesh, std::back_inserter(self_inter));

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
        remove_self_intersections_one_step(tmesh, faces_to_remove, vpm, step, preserve_genus, verbose);

    if(all_fixed && topology_issue)
    {
      if(verbose)
        std::cout<< "DEBUG: Process stopped because of boundary cycles"
                    " of boundary edges involved in self-intersections.\n";
      return false;
    }
  }

  return step < max_steps;
}

template <typename TriangleMesh>
bool remove_self_intersections(TriangleMesh& tmesh)
{
  return remove_self_intersections(tmesh, parameters::all_default());
}

template <typename FaceContainer, typename TriangleMesh, typename SICCContainer>
int identify_self_intersecting_ccs(const FaceContainer& all_intersecting_faces,
                                   TriangleMesh& mesh,
                                   SICCContainer& self_intersecting_ccs)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor            face_descriptor;

  typedef CGAL::dynamic_face_property_t<bool>                                    Face_bool_tag;
  typedef typename boost::property_map<TriangleMesh, Face_bool_tag>::type        Is_self_intersecting_face_map;

  typedef CGAL::dynamic_face_property_t<int>                                     Face_color_tag;
  typedef typename boost::property_map<TriangleMesh, Face_color_tag>::type       SI_CC_face_map;

  std::list<face_descriptor> si_faces(all_intersecting_faces.begin(), all_intersecting_faces.end());

  Is_self_intersecting_face_map is_si_fmap = get(Face_bool_tag(), mesh);
  SI_CC_face_map colors = get(Face_color_tag(), mesh);

  // multiple iterations means we must reset properties
  for(face_descriptor f : faces(mesh))
  {
    put(is_si_fmap, f, false);
    put(colors, f, 0);
  }

  for(face_descriptor f : si_faces)
    put(is_si_fmap, f, true);

  int color = 0;
  while(!si_faces.empty())
  {
    face_descriptor seed_f = si_faces.front();
    si_faces.pop_front();

    if(get(colors, seed_f) != 0) // 'seed_f' is already colored, nothing to do
      continue;
    put(colors, seed_f, ++color);

    std::stack<face_descriptor> to_color;
    to_color.push(seed_f);

    while(!to_color.empty())
    {
      face_descriptor f = to_color.top();
      to_color.pop();

      // add neighbors to stack @todo should this be vertex-neighbors and not edge neighbors
      for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, mesh), mesh))
      {
        halfedge_descriptor opp_h = opposite(h, mesh);
        face_descriptor opp_f = face(opp_h, mesh);

        if(is_border(opp_h, mesh) || !get(is_si_fmap, opp_f))
          continue;

        // 'f' is already colored; nothing to do
        if(get(colors, opp_f) != 0) {
          CGAL_assertion(get(colors, opp_f) == color);
        } else {
          put(colors, opp_f, color);
          to_color.push(opp_f);
        }
      }
    }
  }

  self_intersecting_ccs.resize(color);
  for(face_descriptor f : all_intersecting_faces)
  {
    CGAL_assertion(get(is_si_fmap, f));
    CGAL_assertion(0 < get(colors, f) && get(colors, f) <= color);
    self_intersecting_ccs[get(colors, f) - 1].insert(f); // -1 because colors start at '1'
  }

  return color;
}

template <typename TriangleMesh, typename NamedParameters>
bool remove_self_intersections_locally(TriangleMesh& tmesh,
                                       const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh>                                 graph_traits;
  typedef typename graph_traits::face_descriptor                            face_descriptor;

  // named parameter extraction
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                           get_property_map(vertex_point, tmesh));

  const int max_steps = boost::choose_param(boost::get_param(np, internal_np::number_of_iterations), 7);
  bool verbose = boost::choose_param(boost::get_param(np, internal_np::verbosity_level), 0) > 0;
  bool preserve_genus = boost::choose_param(boost::get_param(np, internal_np::preserve_genus), true);

  if(verbose)
    std::cout << "DEBUG: Starting remove_self_intersections, is_valid(tmesh)? " << is_valid_polygon_mesh(tmesh) << "\n";

  CGAL_precondition_code(std::set<face_descriptor> degenerate_face_set;)
  CGAL_precondition_code(degenerate_faces(tmesh, std::inserter(degenerate_face_set, degenerate_face_set.begin()), np);)
  CGAL_precondition(degenerate_face_set.empty());

  if(!preserve_genus)
    duplicate_non_manifold_vertices(tmesh, np);

  // Look for self-intersections in the mesh and remove them
  int step = -1;
  bool all_fixed = false; // indicates if the filling of all created holes went fine
  bool topology_issue = false; // indicates if some boundary cycles of edges are blocking the fixing
  std::set<face_descriptor> all_intersecting_faces;

  while(++step < max_steps)
  {
    if(all_fixed)
    {
      if(verbose)
        std::cout << "DEBUG: Fixed all local self-intersections!" << std::endl;
      break;
    }

    if(all_intersecting_faces.empty()) // the previous round might have been blocked due to topological constraints
    {
      typedef std::pair<face_descriptor, face_descriptor> Face_pair;
      std::vector<Face_pair> self_inter;
      // TODO : possible optimization to reduce the range to check with the bbox
      // of the previous patches or something.
      self_intersections(tmesh, std::back_inserter(self_inter));

      for(const Face_pair& fp : self_inter)
      {
        all_intersecting_faces.insert(fp.first);
        all_intersecting_faces.insert(fp.second);
      }
    }

    if(verbose)
      std::cout << all_intersecting_faces.size() << " faces part of self-intersections" << std::endl;

    if(all_intersecting_faces.empty())
    {
      if(verbose)
        std::cout << "DEBUG: There is no more face to remove." << std::endl;
      break;
    }

    std::vector<std::set<face_descriptor> > faces_to_remove_ccs;
    identify_self_intersecting_ccs(all_intersecting_faces, tmesh, faces_to_remove_ccs);
    std::cout << faces_to_remove_ccs.size() << " CCs to handle" << std::endl;

    all_fixed = true;

    for(std::set<face_descriptor>& cc : faces_to_remove_ccs)
    {
      bool local_fixed;
      std::tie(local_fixed, topology_issue) =
          remove_self_intersections_one_step(tmesh, cc, vpm, step, preserve_genus, verbose);

      if(topology_issue)
      {
        if(verbose)
          std::cout << "DEBUG: Process stopped because of boundary cycles"
                       " of boundary edges involved in self-intersections.\n";

        return false;
      }

      all_fixed = (all_fixed && local_fixed);
    }
  }

  return step < max_steps;
}

template <typename TriangleMesh>
bool remove_self_intersections_locally(TriangleMesh& tmesh)
{
  return remove_self_intersections_locally(tmesh, parameters::all_default());
}

/// \endcond

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REMOVE_SELF_INTERSECTIONS_H
