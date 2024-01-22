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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SURFACE_SNAPPING_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SURFACE_SNAPPING_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <map>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace experimental {

template <class K, class AABB_tree_, class TriangleMesh, class VPM, class MatchingMap>
void surface_snapping_impl_tm1_on_tm2(const AABB_tree_& tm1_tree,
                                      TriangleMesh& tm1, TriangleMesh& tm2,
                                      VPM vpm1, VPM vpm2,
                                      MatchingMap matching_vertices_1, MatchingMap matching_vertices_2,
                                      double eps2)
{
  using BGT = boost::graph_traits<TriangleMesh>;
  using vertex_descriptor = typename BGT::vertex_descriptor;
  using face_descriptor = typename BGT::face_descriptor;
  using halfedge_descriptor = typename BGT::halfedge_descriptor;

  auto split_edge_and_incident_faces =
  [](halfedge_descriptor h, const typename K::Point_3& pt, TriangleMesh& tm, VPM vpm)
  {
    halfedge_descriptor res = Euler::split_edge(h, tm);
    put(vpm, target(res, tm), pt);

    if(!is_border(res, tm))
      Euler::split_face(res, next(h, tm), tm);

    halfedge_descriptor opp_h = opposite(h, tm);
    if(!is_border(opp_h, tm))
      Euler::split_face(opp_h, next(next(opp_h, tm), tm), tm);

    return res;
  };

  // TODO PARALLEL FOR
  for (vertex_descriptor vd : vertices(tm2))
  {
    if (get(matching_vertices_2, vd)) continue; // already on tm1
    typename K::Sphere_3 epsilon_sphere(get(vpm2, vd), eps2);
    std::vector<face_descriptor> closed_faces;
    tm1_tree.all_intersected_primitives(epsilon_sphere, std::back_inserter(closed_faces));

    switch(closed_faces.size())
    {
      case 0:
        break;
      case 1:
      {
        // check if close to border edges
        face_descriptor f1 = closed_faces[0];
        halfedge_descriptor h1 = halfedge(f1, tm1);
        std::array<halfedge_descriptor, 3> f1_hedges;
        std::vector<int> closed_hedges;
        for (int i=0; i<3; ++i)
        {
          f1_hedges[i]=h1;
          if (face(opposite(h1, tm1), tm1) == BGT::null_face())
          {
            if (do_intersect(typename K::Segment_3(get(vpm1,source(h1,tm1)),
                                                   get(vpm1,target(h1,tm1))),
                             epsilon_sphere))
            {
              closed_hedges.push_back(i);
            }
          }
          h1 = next(h1, tm1);
        }
        switch(closed_hedges.size())
        {
          case 0:
          {
            halfedge_descriptor hnew = Euler::add_center_vertex(h1,tm1);
            put(vpm1, target(hnew, tm1), get(vpm2, vd));
            put(matching_vertices_1, target(hnew, tm1), true);
            put(matching_vertices_2, vd, true);
            break;
          }
          case 1:
          {
            h1 = f1_hedges[closed_hedges[0]];
            // try to snap on h1 edge points
            if (compare_squared_distance(get(vpm1, target(h1, tm1)), get(vpm2, vd), eps2)!=LARGER)
            {
              put(vpm2, vd, get(vpm1,  target(h1, tm1)));
              put(matching_vertices_1, target(h1, tm1), true);
              put(matching_vertices_2, vd, true);
              break;
            }
            if (compare_squared_distance(get(vpm1, source(h1, tm1)), get(vpm2, vd), eps2)!=LARGER)
            {
              put(vpm2, vd, get(vpm1,  source(h1, tm1)));
              put(matching_vertices_1, source(h1, tm1), true);
              put(matching_vertices_2, vd, true);
              break;
            }
            // split and snap onto the edge
            halfedge_descriptor hnew =
              split_edge_and_incident_faces(h1, get(vpm2, vd), tm1, vpm1);
            put(matching_vertices_1, target(hnew, tm1), true);
            put(matching_vertices_2, vd, true);
            break;
          }
          case 2:
          {
            // snap onto the common vertex
            h1=next(f1_hedges[closed_hedges[0]+closed_hedges[1]-3], tm1);
            put(vpm2, vd, get(vpm1,  target(h1, tm1)));
            put(matching_vertices_1, target(h1, tm1), true);
            put(matching_vertices_2, vd, true);
            break;
          }
          default:
            //snap to any of the face vertex
            put(vpm2, vd, get(vpm1,target(h1, tm1)));
            put(matching_vertices_1, target(h1, tm1), true);
            put(matching_vertices_2, vd, true);
        }
        break;
      }
      case 2:
      {
        // TODO: in case some edges are boundary edges, we might snap to a vertex
        // look for a common edge
        face_descriptor f1 = closed_faces[0], f2 = closed_faces[1];
        halfedge_descriptor h1 = halfedge(f1, tm1);
        halfedge_descriptor h_to_split = BGT::null_halfedge();
        for (int i=0; i<3; ++i)
        {
          if (face(opposite(h1, tm1), tm1) == f2)
          {
            h_to_split=h1;
            break;
          }
          h1 = next(h1, tm1);
        }
        if (h_to_split != BGT::null_halfedge())
        {
          halfedge_descriptor hnew = split_edge_and_incident_faces(h_to_split, get(vpm2, vd), tm1, vpm1);
          put(matching_vertices_1, target(hnew, tm1), true);
          put(matching_vertices_2, vd, true);
          break;
        }
        break; // TODO: decide what to do here. snap on the closest face?
      }
      default:
      {
        // look for a vertex common to all faces
        std::set<vertex_descriptor> vset;
        halfedge_descriptor h1 = halfedge(closed_faces[0], tm1);
        for (int i=0;i<3; ++i)
        {
          vset.insert(target(h1, tm1));
          h1 = next(h1, tm1);
        }
        for (auto fit = std::next(closed_faces.begin()); fit!=closed_faces.end(); ++fit)
        {
          h1 = halfedge(*fit, tm1);
          std::set<vertex_descriptor> common;
          for (int i=0;i<3; ++i)
          {
            if (vset.count(target(h1, tm1))==1)
              common.insert(target(h1,tm1));
            h1 = next(h1, tm1);
          }
          common.swap(vset);
          if (vset.empty()) break;
        }
        if (vset.size()!=1) break; // TODO snap on the closest?
        //snap on the vertex
        put(vpm2, vd, get(vpm1,  *vset.begin()));
        put(matching_vertices_1, *vset.begin(), true);
        put(matching_vertices_2, vd, true);
      }
    }
  }
}

// TODO test on Thingi10K/raw_meshes/162100.stl after splitting it
template <class TriangleMesh>
void surface_snapping(TriangleMesh& tm1, TriangleMesh& tm2, double epsilon)
{
  using VPM = typename property_map_selector<TriangleMesh, boost::vertex_point_t>::type;
  using K = typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;

  using AABB_primitive_ = AABB_face_graph_triangle_primitive<TriangleMesh, VPM>;
  using AABB_traits_ = AABB_traits<K, AABB_primitive_>;
  using AABB_tree_ = AABB_tree<AABB_traits_>;

  using BGT = boost::graph_traits<TriangleMesh>;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  VPM vpm1 = get(boost::vertex_point, tm1);
  VPM vpm2 = get(boost::vertex_point, tm2);

// check for identical points in tm1 and tm2
  // Note does not work well in parallel: bool --> int
  auto matching_vertices_1 = get(CGAL::dynamic_vertex_property_t<bool>(), tm1);
  auto matching_vertices_2 = get(CGAL::dynamic_vertex_property_t<bool>(), tm2);
  std::map<typename K::Point_3, vertex_descriptor> pts_map;

  for (vertex_descriptor vd : vertices(tm1))
  {
    auto insert_res = pts_map.emplace(get(vpm1, vd), vd);
    if (!insert_res.second)
      std::cerr << "Warning duplicated vertex in tm1, better solve that first!\n";
  }

  for (vertex_descriptor vd : vertices(tm2))
  {
    auto insert_res = pts_map.emplace(get(vpm2, vd), BGT::null_vertex());
    if (!insert_res.second)
    {
      if (insert_res.first->second == BGT::null_vertex())
        std::cerr << "Warning duplicated vertex in tm2, better solve that first!\n";
      else
      {
        if (get(matching_vertices_1, insert_res.first->second))
        {
          std::cerr << "Two vertices of tm2 closed to the same vertex of tm1\n";
        }
        put(matching_vertices_1, insert_res.first->second, true);
        put(matching_vertices_2, vd, true);
      }
    }
  }

  AABB_tree_ tm1_tree(faces(tm1).begin(), faces(tm1).end(), tm1);
  AABB_tree_ tm2_tree(faces(tm2).begin(), faces(tm2).end(), tm2);

  double eps2 = epsilon * epsilon;
  surface_snapping_impl_tm1_on_tm2<K>(tm1_tree, tm1, tm2, vpm1, vpm2,
                                      matching_vertices_1, matching_vertices_2,
                                      eps2);
  surface_snapping_impl_tm1_on_tm2<K>(tm2_tree, tm2, tm1, vpm2, vpm1,
                                      matching_vertices_2, matching_vertices_1,
                                      eps2);



  // TODO: apply the same thing using tm2_tree (

}

} } } /// end of CGAL::Polygon_mesh_processing::experimental

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SURFACE_SNAPPING_H