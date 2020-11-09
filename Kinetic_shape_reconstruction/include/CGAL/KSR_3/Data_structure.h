// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_DATA_STRUCTURE_H
#define CGAL_KSR_3_DATA_STRUCTURE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <queue>

#include <boost/function_output_iterator.hpp>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/verbosity.h>
#include <CGAL/KSR/debug.h>

#include <CGAL/KSR_3/Support_plane.h>
#include <CGAL/KSR_3/Intersection_graph.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/centroid.h>

#include <CGAL/Delaunay_triangulation_2.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Data_structure
{
public:

  typedef GeomTraits Kernel;

private:

  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Direction_2 Direction_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Triangle_2 Triangle_2;

  typedef KSR_3::Support_plane<Kernel> Support_plane;
  typedef typename Support_plane::Mesh Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Face_index Face_index;
  typedef typename Mesh::Edge_index Edge_index;
  typedef typename Mesh::Halfedge_index Halfedge_index;

  typedef KSR_3::Intersection_graph<Kernel> Intersection_graph;

  typedef KSR::vector<Support_plane> Support_planes;

public:

  typedef std::pair<KSR::size_t, Vertex_index> PVertex;
  typedef std::pair<KSR::size_t, Edge_index> PEdge;
  typedef std::pair<KSR::size_t, Face_index> PFace;

  template <typename PSimplex>
  struct Make_PSimplex
  {
    typedef typename PSimplex::second_type argument_type;
    typedef PSimplex result_type;
    KSR::size_t support_plane_idx;

    Make_PSimplex (KSR::size_t support_plane_idx) : support_plane_idx (support_plane_idx) { }

    result_type operator() (const argument_type& arg) const
    {
      return result_type(support_plane_idx, arg);
    }
  };

  typedef boost::transform_iterator<Make_PSimplex<PVertex>,
                                    typename Mesh::Vertex_range::iterator> PVertex_iterator;
  typedef boost::transform_iterator<Make_PSimplex<PEdge>,
                                    typename Mesh::Edge_range::iterator> PEdge_iterator;
  typedef boost::transform_iterator<Make_PSimplex<PFace>,
                                    typename Mesh::Face_range::iterator> PFace_iterator;

  typedef Iterator_range<PVertex_iterator> PVertices;
  typedef Iterator_range<PEdge_iterator> PEdges;
  typedef Iterator_range<PFace_iterator> PFaces;

  struct Halfedge_to_pvertex
  {
    typedef Halfedge_index argument_type;
    typedef PVertex result_type;
    KSR::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pvertex (KSR::size_t support_plane_idx, const Mesh& mesh)
      : support_plane_idx (support_plane_idx), mesh (mesh) { }

    result_type operator() (const argument_type& arg) const
    {
      return result_type(support_plane_idx, mesh.target(arg));
    }
  };

  typedef boost::transform_iterator<Halfedge_to_pvertex,
                                    Halfedge_around_face_iterator<Mesh> > PVertex_of_pface_iterator;
  typedef Iterator_range<PVertex_of_pface_iterator> PVertices_of_pface;

  struct Halfedge_to_pedge
  {
    typedef Halfedge_index argument_type;
    typedef PEdge result_type;
    KSR::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pedge (KSR::size_t support_plane_idx, const Mesh& mesh)
      : support_plane_idx (support_plane_idx), mesh (mesh) { }

    result_type operator() (const argument_type& arg) const
    {
      return result_type(support_plane_idx, mesh.edge(arg));
    }
  };

  typedef boost::transform_iterator<Halfedge_to_pedge,
                                    Halfedge_around_target_iterator<Mesh> > PEdge_around_pvertex_iterator;
  typedef Iterator_range<PEdge_around_pvertex_iterator> PEdges_around_pvertex;

  typedef boost::transform_iterator<Halfedge_to_pedge,
                                    Halfedge_around_face_iterator<Mesh> > PEdge_of_pface_iterator;
  typedef Iterator_range<PEdge_of_pface_iterator> PEdges_of_pface;

  typedef typename Intersection_graph::Vertex_descriptor IVertex;
  typedef typename Intersection_graph::Edge_descriptor IEdge;
  typedef typename Intersection_graph::Vertices IVertices;
  typedef typename Intersection_graph::Edges IEdges;
  typedef typename Intersection_graph::Incident_edges Incident_iedges;

  struct Volume_cell {
    std::vector<PFace> pfaces;
    std::vector<int> neighbors;
    std::set<PVertex> pvertices;

    void add_pface(const PFace& pface, const int neighbor) {
      pfaces.push_back(pface);
      neighbors.push_back(neighbor);
    }
  };

private:

  // Main data structure
  Support_planes m_support_planes;
  Intersection_graph m_intersection_graph;
  std::vector<Volume_cell> m_volumes;

  // Helping data structures
  std::map<Point_3, KSR::size_t> m_meta_map;

  FT m_current_time;
  FT m_previous_time;

public:

  Data_structure()
    : m_current_time(0)
  { }

  void print() const
  {
    CGAL_assertion_msg(false, "TODO: ADD PRINTING!");
  }

  void init (std::size_t number_of_polygons)
  {
    m_support_planes.reserve (number_of_polygons + 6);
  }

  const FT& current_time() const { return m_current_time; }

  void set_last_event_time(const PVertex& pvertex, const FT time) {
    support_plane(pvertex).set_last_event_time(pvertex.second, time);
  }

  const FT last_event_time(const PVertex& pvertex) {
    return support_plane(pvertex).last_event_time(pvertex.second);
  }

  const std::vector<Volume_cell>& polyhedrons() const {
    return m_volumes;
  }

  /*******************************
   * Support planes
   *******************************/

  KSR::size_t number_of_support_planes() const { return m_support_planes.size(); }

  bool is_bbox_support_plane (KSR::size_t support_plane_idx) const
  { return (support_plane_idx < 6); }

  bool mesh_is_valid (KSR::size_t support_plane_idx) const
  {
    bool is_valid = mesh(support_plane_idx).is_valid();
    if (!is_valid)
      return false;

    for (PFace pface : pfaces(support_plane_idx))
    {
      std::function<Point_2(PVertex)> unary_f = [&](const PVertex& pvertex) -> Point_2 { return point_2(pvertex); };
      CGAL::Polygon_2<Kernel> polygon
        (boost::make_transform_iterator
         (pvertices_of_pface(pface).begin(), unary_f),
         boost::make_transform_iterator
         (pvertices_of_pface(pface).end(), unary_f));

      // if (!polygon.is_simple())
      // {
      //   std::cerr << "PFace(" << pface.first << ":" << pface.second << ") is not simple" << std::endl;
      //   for (const Point_2& p : polygon)
      //     std::cerr << to_3d(support_plane_idx,p) << " ";
      //   std::cerr << to_3d(support_plane_idx,polygon[0]) << " ";
      //   std::cerr << std::endl;
      //   return false;
      // }

      // if (!polygon.is_convex())
      // {
      //   std::cerr << "PFace(" << pface.first << ":" << pface.second << ") is not convex" << std::endl;
      //   for (const Point_2& p : polygon)
      //     std::cerr << to_3d(support_plane_idx,p) << " ";
      //   std::cerr << to_3d(support_plane_idx,polygon[0]) << " ";
      //   std::cerr << std::endl;
      //   return false;
      // }

      PVertex prev = null_pvertex();

      for (const PVertex pvertex : pvertices_of_pface (pface))
      {
        if (prev == null_pvertex())
        {
          prev = pvertex;
          continue;
        }

        if (point_2(prev) == point_2(pvertex)
            && direction(prev) == direction(pvertex))

        {
          std::cerr << "PFace(" << pface.first << ":" << pface.second << ") has two consequent identical vertices "
                    << str(prev) << " and " << str(pvertex) << std::endl;
          return false;
        }

        prev = pvertex;
      }
    }

    return true;
  }

  KSR::size_t add_support_plane (const Support_plane& new_support_plane)
  {
    KSR::size_t support_plane_idx = KSR::no_element();
    for (KSR::size_t i = 0; i < number_of_support_planes(); ++ i)
      if (new_support_plane == support_plane(i))
      {
        support_plane_idx = i;
        break;
      }

    if (support_plane_idx == KSR::no_element())
    {
      support_plane_idx = number_of_support_planes();
      m_support_planes.push_back (new_support_plane);
    }

    if (support_plane_idx >= 6) // Intersect planes with bbox...
    {
      std::vector<std::pair<IEdge, Point_3> > intersections;

      Point_3 centroid = CGAL::ORIGIN;
      for (const IEdge edge : m_intersection_graph.edges())
      {
        Point_3 point;
        if (!KSR::intersection_3 (support_plane(support_plane_idx).plane(),
                                  m_intersection_graph.segment_3 (edge), point))
          continue;

        centroid = CGAL::barycenter (centroid, intersections.size(), point, 1);
        intersections.push_back (std::make_pair (edge, point));
      }

      Point_2 centroid_2 = support_plane(support_plane_idx).to_2d (centroid);
      std::sort (intersections.begin(), intersections.end(),
                 [&] (const std::pair<IEdge, Point_3>& a,
                      const std::pair<IEdge, Point_3>& b) -> bool
                 {
                   return (Direction_2 (Segment_2 (centroid_2, support_plane(support_plane_idx).to_2d (a.second)))
                           < Direction_2 (Segment_2 (centroid_2, support_plane(support_plane_idx).to_2d (b.second))));
                 });

      KSR::vector<KSR::size_t> common_planes_idx;
      std::map<KSR::size_t, KSR::size_t> map_lines_idx;
      KSR::vector<IVertex> vertices;
      vertices.reserve (intersections.size());
      for (std::size_t i = 0; i < intersections.size(); ++ i)
      {
        const IEdge& e0 = intersections[i].first;
        const IEdge& e1 = intersections[(i+1)%intersections.size()].first;

        KSR::size_t common_plane_idx = KSR::no_element();
        std::set_intersection (m_intersection_graph.intersected_planes(e0).begin(),
                               m_intersection_graph.intersected_planes(e0).end(),
                               m_intersection_graph.intersected_planes(e1).begin(),
                               m_intersection_graph.intersected_planes(e1).end(),
                               boost::make_function_output_iterator
                               ([&](const KSR::size_t& idx) -> void
                                {
                                  if (idx < 6)
                                  {
                                    CGAL_assertion (common_plane_idx == KSR::no_element());
                                    common_plane_idx = idx;
                                  }
                                }));
        CGAL_assertion (common_plane_idx != KSR::no_element());
        common_planes_idx.push_back (common_plane_idx);

        typename std::map<KSR::size_t, KSR::size_t>::iterator iter;
        bool inserted;
        std::tie (iter, inserted)
          = map_lines_idx.insert (std::make_pair (common_plane_idx, KSR::no_element()));
        if (inserted)
          iter->second = m_intersection_graph.add_line();

        vertices.push_back (m_intersection_graph.add_vertex (intersections[i].second).first);
      }

      for (std::size_t i = 0; i < intersections.size(); ++ i)
      {
        for (KSR::size_t sp_idx : m_intersection_graph.intersected_planes(intersections[i].first))
          support_plane(sp_idx).iedges().erase (intersections[i].first);
        const auto edges = m_intersection_graph.split_edge (intersections[i].first, vertices[i]);
        // for (const IEdge& edge : { edge_0, edge_1 })
        //   for (KSR::size_t sp_idx : m_intersection_graph.intersected_planes(edge))
        //     support_plane(sp_idx).iedges().insert (edge); // bugs!

        const auto& iplanes_1 = m_intersection_graph.intersected_planes(edges.first);
        for (const KSR::size_t sp_idx : iplanes_1) {
          support_plane(sp_idx).iedges().insert(edges.first);
        }

        const auto& iplanes_2 = m_intersection_graph.intersected_planes(edges.second);
        for (const KSR::size_t sp_idx : iplanes_2) {
          support_plane(sp_idx).iedges().insert(edges.second);
        }

        IEdge new_edge =
          m_intersection_graph.add_edge (vertices[i], vertices[(i+1)%vertices.size()], support_plane_idx).first;
        m_intersection_graph.intersected_planes(new_edge).insert (common_planes_idx[i]);
        m_intersection_graph.set_line (new_edge, map_lines_idx[common_planes_idx[i]]);

        support_plane(support_plane_idx).iedges().insert (new_edge);
        support_plane(common_planes_idx[i]).iedges().insert (new_edge);
      }
    }

    return support_plane_idx;
  }

  void add_bbox_polygon (const std::array<Point_3, 4>& polygon)
  {
    KSR::size_t support_plane_idx = add_support_plane (Support_plane (polygon));

    std::array<IVertex, 4> ivertices;
    std::array<Point_2, 4> points;
    for (std::size_t i = 0; i < 4; ++ i)
    {
      points[i] = support_plane(support_plane_idx).to_2d(polygon[i]);
      ivertices[i] = m_intersection_graph.add_vertex(polygon[i]).first;
    }

    std::array<Vertex_index, 4> vertices
      = support_plane(support_plane_idx).add_bbox_polygon (points, ivertices);

    for (std::size_t i = 0; i < 4; ++ i)
    {
      IEdge iedge;
      bool inserted;
      std::tie (iedge, inserted)
        = m_intersection_graph.add_edge (ivertices[i], ivertices[(i+1)%4], support_plane_idx);
      if (inserted)
        m_intersection_graph.set_line (iedge, m_intersection_graph.add_line());

      support_plane(support_plane_idx).set_iedge
        (vertices[i], vertices[(i+1)%4], iedge);

      support_plane(support_plane_idx).iedges().insert (iedge);
    }
  }

  template <typename PointRange>
  void add_polygon (const PointRange& polygon, KSR::size_t input_idx)
  {
    KSR::size_t support_plane_idx = add_support_plane (Support_plane (polygon));

    // Create ordered polygon
    std::vector<Point_2> points;
    points.reserve (polygon.size());
    for (const Point_3& p : polygon)
      points.push_back (support_plane(support_plane_idx).to_2d(p));

    // const auto centroid = CGAL::centroid(points.begin(), points.end());

    using TRI = CGAL::Delaunay_triangulation_2<Kernel>;
    TRI tri(points.begin(), points.end());
    std::vector<Triangle_2> triangles;
    triangles.reserve(tri.number_of_faces());
    for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
      triangles.push_back(Triangle_2(
          fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point()));
    }
    const auto centroid = CGAL::centroid(triangles.begin(), triangles.end());

    std::sort (points.begin(), points.end(),
               [&](const Point_2& a, const Point_2& b) -> bool
               {
                 return (Direction_2 (Segment_2 (centroid, a))
                         < Direction_2 (Segment_2 (centroid, b)));
               });

    support_plane(support_plane_idx).add_polygon (points, centroid, input_idx);
  }

  /*******************************
   * PSimplices
   *******************************/

  static PVertex null_pvertex() { return PVertex(KSR::no_element(), Vertex_index()); }
  static PEdge null_pedge() { return PEdge(KSR::no_element(), Edge_index()); }
  static PFace null_pface() { return PFace(KSR::no_element(), Face_index()); }

  PVertices pvertices (KSR::size_t support_plane_idx) const
  {
    return PVertices (boost::make_transform_iterator
                      (mesh(support_plane_idx).vertices().begin(),
                       Make_PSimplex<PVertex>(support_plane_idx)),
                      boost::make_transform_iterator
                      (mesh(support_plane_idx).vertices().end(),
                       Make_PSimplex<PVertex>(support_plane_idx)));

  }

  PEdges pedges (KSR::size_t support_plane_idx) const
  {
    return PEdges (boost::make_transform_iterator
                   (mesh(support_plane_idx).edges().begin(),
                    Make_PSimplex<PEdge>(support_plane_idx)),
                   boost::make_transform_iterator
                   (mesh(support_plane_idx).edges().end(),
                    Make_PSimplex<PEdge>(support_plane_idx)));

  }

  PFaces pfaces (KSR::size_t support_plane_idx) const
  {
    return PFaces (boost::make_transform_iterator
                   (mesh(support_plane_idx).faces().begin(),
                    Make_PSimplex<PFace>(support_plane_idx)),
                   boost::make_transform_iterator
                   (mesh(support_plane_idx).faces().end(),
                    Make_PSimplex<PFace>(support_plane_idx)));

  }

  // Get prev and next of free vertex
  PVertex prev (const PVertex& pvertex) const
  {
    return PVertex (pvertex.first, support_plane(pvertex).prev(pvertex.second));
  }
  PVertex next (const PVertex& pvertex) const
  {
    return PVertex (pvertex.first, support_plane(pvertex).next(pvertex.second));
  }

  // Get prev and next of constrained vertex
  std::pair<PVertex, PVertex> prev_and_next (const PVertex& pvertex) const
  {
    std::pair<PVertex, PVertex> out (null_pvertex(), null_pvertex());

    for (Halfedge_index hi : halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)))
    {
      IEdge iedge = support_plane(pvertex).iedge (mesh(pvertex).edge(hi));
      if (iedge == this->iedge(pvertex))
        continue;
      if (out.first == null_pvertex())
        out.first = PVertex (pvertex.first, mesh(pvertex).source(hi));
      else
      {
        out.second = PVertex (pvertex.first, mesh(pvertex).source(hi));
        return out;
      }
    }

    return out;
  }

  const std::pair<PVertex, PVertex> border_prev_and_next(const PVertex& pvertex) const {

    // std::cout << point_3(pvertex) << std::endl;
    Halfedge_index he = mesh(pvertex).halfedge(pvertex.second);
    const auto end = he;

    // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).source(he))) << std::endl;
    // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).target(he))) << std::endl;

    // If the assertion below fails, it probably means that we need to circulate
    // longer until we hit the border edge!

    std::size_t count = 0;
    while (true) {
      if (mesh(pvertex).face(he) != Face_index()) {
        he = mesh(pvertex).prev(mesh(pvertex).opposite(he));

        // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).source(he))) << std::endl;
        // std::cout << point_3(PVertex(pvertex.first, mesh(pvertex).target(he))) << std::endl;

        ++count;
      } else { break; }

      // std::cout << "count: " << count << std::endl;
      CGAL_assertion(count <= 2);
      if (he == end) {
        CGAL_assertion_msg(false, "ERROR: BORDER HALFEDGE IS NOT FOUND, FULL CIRCLE!");
        break;
      }
      if (count == 100) {
        CGAL_assertion_msg(false, "ERROR: BORDER HALFEDGE IS NOT FOUND, LIMIT ITERATIONS!");
        break;
      }
    }

    CGAL_assertion(mesh(pvertex).face(he) == Face_index());
    return std::make_pair(
      PVertex(pvertex.first, mesh(pvertex).source(he)),
      PVertex(pvertex.first, mesh(pvertex).target(mesh(pvertex).next(he))));
  }

  PVertex add_pvertex (KSR::size_t support_plane_idx, const Point_2& point)
  {
    CGAL_assertion(support_plane_idx != KSR::uninitialized());
    CGAL_assertion(support_plane_idx != KSR::no_element());

    auto& m = mesh(support_plane_idx);
    const auto vertex_index = m.add_vertex(point);
    CGAL_assertion(vertex_index != typename Support_plane::Mesh::Vertex_index());
    return PVertex(support_plane_idx, vertex_index);
  }

  template <typename VertexRange>
  PFace add_pface (const VertexRange& pvertices)
  {
    auto support_plane_idx = pvertices.front().first;
    CGAL_assertion(support_plane_idx != KSR::uninitialized());
    CGAL_assertion(support_plane_idx != KSR::no_element());

    auto& m = mesh(support_plane_idx);
    const auto range = CGAL::make_range(
      boost::make_transform_iterator(pvertices.begin(),
      CGAL::Property_map_to_unary_function<CGAL::Second_of_pair_property_map<PVertex> >()),
      boost::make_transform_iterator(pvertices.end(),
      CGAL::Property_map_to_unary_function<CGAL::Second_of_pair_property_map<PVertex> >()));
    const auto face_index = m.add_face(range);
    CGAL_assertion(face_index != Support_plane::Mesh::null_face());
    return PFace(support_plane_idx, face_index);
  }

  void clear_polygon_faces (KSR::size_t support_plane_idx)
  {
    Mesh& m = mesh(support_plane_idx);
    for (Face_index fi : m.faces())
      m.remove_face(fi);
    for (Edge_index ei : m.edges())
      m.remove_edge(ei);
    for (Vertex_index vi : m.vertices())
      m.set_halfedge(vi, Halfedge_index());
  }

  PVertex source (const PEdge& pedge) const
  { return PVertex (pedge.first, mesh(pedge).source(mesh(pedge).halfedge(pedge.second))); }
  PVertex target (const PEdge& pedge) const
  { return PVertex (pedge.first, mesh(pedge).target(mesh(pedge).halfedge(pedge.second))); }
  PVertex opposite (const PEdge& pedge, const PVertex& pvertex) const
  {
    if (mesh(pedge).target(mesh(pedge).halfedge(pedge.second)) == pvertex.second)
      return PVertex (pedge.first, mesh(pedge).source(mesh(pedge).halfedge(pedge.second)));

    CGAL_assertion (mesh(pedge).source(mesh(pedge).halfedge(pedge.second)) == pvertex.second);
    return PVertex (pedge.first, mesh(pedge).target(mesh(pedge).halfedge(pedge.second)));
  }

  PFace pface_of_pvertex (const PVertex& pvertex) const
  {
    return PFace (pvertex.first, support_plane(pvertex).face (pvertex.second));
  }

  std::pair<PFace, PFace> pfaces_of_pvertex (const PVertex& pvertex) const
  {
    std::pair<PFace, PFace> out (null_pface(), null_pface());

    std::tie (out.first.second, out.second.second)
      = support_plane(pvertex).faces (pvertex.second);

    if (out.first.second != Face_index())
      out.first.first = pvertex.first;
    if (out.second.second != Face_index())
      out.second.first = pvertex.first;
    return out;
  }

  PVertices_of_pface pvertices_of_pface (const PFace& pface) const
  {
    return PVertices_of_pface (boost::make_transform_iterator
                               (halfedges_around_face(halfedge(pface.second, mesh(pface)),
                                                      mesh(pface)).begin(),
                                Halfedge_to_pvertex(pface.first, mesh(pface))),
                               boost::make_transform_iterator
                               (halfedges_around_face(halfedge(pface.second, mesh(pface)),
                                                      mesh(pface)).end(),
                                Halfedge_to_pvertex(pface.first, mesh(pface))));
  }

  PEdges_of_pface pedges_of_pface(const PFace& pface) const {

    return PEdges_of_pface(
      boost::make_transform_iterator(
      halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).begin(),
      Halfedge_to_pedge(pface.first, mesh(pface))),
      boost::make_transform_iterator(
      halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).end(),
      Halfedge_to_pedge(pface.first, mesh(pface))));
  }

  PEdges_around_pvertex pedges_around_pvertex (const PVertex& pvertex) const
  {
    return PEdges_around_pvertex (boost::make_transform_iterator
                                  (halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)),
                                                           mesh(pvertex)).begin(),
                                   Halfedge_to_pedge(pvertex.first, mesh(pvertex))),
                                  boost::make_transform_iterator
                                  (halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)),
                                                           mesh(pvertex)).end(),
                                   Halfedge_to_pedge(pvertex.first, mesh(pvertex))));
  }

  const KSR::size_t& input (const PFace& pface) const
  { return support_plane(pface).input(pface.second); }
  KSR::size_t& input (const PFace& pface)
  { return support_plane(pface).input(pface.second); }

  const unsigned int& k (const PFace& pface) const
  { return support_plane(pface).k(pface.second); }
  unsigned int& k (const PFace& pface)
  { return support_plane(pface).k(pface.second); }

  bool is_frozen (const PVertex& pvertex) const
  { return support_plane(pvertex).is_frozen (pvertex.second); }
  const Vector_2& direction (const PVertex& pvertex) const
  { return support_plane(pvertex).direction (pvertex.second); }
  Vector_2& direction (const PVertex& pvertex)
  { return support_plane(pvertex).direction (pvertex.second); }
  FT speed (const PVertex& pvertex)
  { return support_plane(pvertex).speed (pvertex.second); }

  void freeze (PVertex& pvertex)
  {
    Point_2 p = point_2 (pvertex, m_current_time);
    support_plane(pvertex).direction (pvertex.second) = CGAL::NULL_VECTOR;
    support_plane(pvertex).set_point (pvertex.second, p);
  }

  bool is_active (const PVertex& pvertex) const
  { return support_plane(pvertex).is_active (pvertex.second); }

  void deactivate (const PVertex& pvertex)
  {
    support_plane(pvertex).set_active (pvertex.second, false);
    if (iedge(pvertex) != null_iedge())
      m_intersection_graph.is_active(iedge(pvertex)) = false;
    if (ivertex(pvertex) != null_ivertex())
      m_intersection_graph.is_active(ivertex(pvertex)) = false;
  }
  void activate (const PVertex& pvertex)
  {
    support_plane(pvertex).set_active (pvertex.second, true);
    if (iedge(pvertex) != null_iedge())
      m_intersection_graph.is_active(iedge(pvertex)) = true;
    if (ivertex(pvertex) != null_ivertex())
      m_intersection_graph.is_active(ivertex(pvertex)) = true;
  }

#ifdef CGAL_KSR_DEBUG
  template <typename PSimplex>
  const Mesh& dbg_mesh (const PSimplex& psimplex) const { return dbg_mesh(psimplex.first); }
  const Mesh& dbg_mesh (KSR::size_t support_plane_idx) const { return support_plane(support_plane_idx).dbg_mesh(); }

  PVertices dbg_pvertices (KSR::size_t support_plane_idx) const
  {
    return PVertices (boost::make_transform_iterator
                      (dbg_mesh(support_plane_idx).vertices().begin(),
                       Make_PSimplex<PVertex>(support_plane_idx)),
                      boost::make_transform_iterator
                      (dbg_mesh(support_plane_idx).vertices().end(),
                       Make_PSimplex<PVertex>(support_plane_idx)));
  }
  PFaces dbg_pfaces (KSR::size_t support_plane_idx) const
  {
    return PFaces (boost::make_transform_iterator
                   (dbg_mesh(support_plane_idx).faces().begin(),
                    Make_PSimplex<PFace>(support_plane_idx)),
                   boost::make_transform_iterator
                   (dbg_mesh(support_plane_idx).faces().end(),
                    Make_PSimplex<PFace>(support_plane_idx)));

  }
  PVertices_of_pface dbg_pvertices_of_pface (const PFace& pface) const
  {
    return PVertices_of_pface (boost::make_transform_iterator
                               (halfedges_around_face(halfedge(pface.second, dbg_mesh(pface)),
                                                      dbg_mesh(pface)).begin(),
                                Halfedge_to_pvertex(pface.first, dbg_mesh(pface))),
                               boost::make_transform_iterator
                               (halfedges_around_face(halfedge(pface.second, dbg_mesh(pface)),
                                                      dbg_mesh(pface)).end(),
                                Halfedge_to_pvertex(pface.first, dbg_mesh(pface))));
  }
  Point_3 dbg_point_3 (const PVertex& pvertex) const
  { return support_plane(pvertex).dbg_point_3 (pvertex.second, m_current_time); }
#endif

  /*******************************
   * ISimplices
   *******************************/

  static IVertex null_ivertex() { return Intersection_graph::null_ivertex(); }
  static IEdge null_iedge() { return Intersection_graph::null_iedge(); }

  IVertices ivertices() const { return m_intersection_graph.vertices(); }
  IEdges iedges() const { return m_intersection_graph.edges(); }

  KSR::size_t nb_intersection_lines() const { return m_intersection_graph.nb_lines(); }
  KSR::size_t line_idx (const IEdge& iedge) const { return m_intersection_graph.line (iedge); }
  KSR::size_t line_idx (const PVertex& pvertex) const { return line_idx (iedge(pvertex)); }

  IVertex add_ivertex (const Point_3& point, const KSR::Idx_set& support_planes_idx)
  {
    KSR::Idx_vector vec_planes;
    std::copy (support_planes_idx.begin(), support_planes_idx.end(),
               std::back_inserter (vec_planes));

    IVertex vertex;
    bool inserted;
    std::tie (vertex, inserted) = m_intersection_graph.add_vertex (point, vec_planes);
    return vertex;
  }

  void add_iedge (const KSR::Idx_set& support_planes_idx,
                  KSR::vector<IVertex>& vertices)
  {
    Point_3 source = m_intersection_graph.point_3 (vertices.front());

    std::sort (vertices.begin(), vertices.end(),
               [&](const IVertex& a, const IVertex& b) -> bool
               {
                 return (CGAL::squared_distance (source, m_intersection_graph.point_3(a))
                         < CGAL::squared_distance (source, m_intersection_graph.point_3(b)));
               });

    KSR::size_t line_idx = m_intersection_graph.add_line();

    for (KSR::size_t i = 0; i < vertices.size() - 1; ++ i)
    {
      IEdge iedge;
      bool inserted;
      std::tie (iedge, inserted)
        = m_intersection_graph.add_edge (vertices[i],
                                         vertices[i+1],
                                         support_planes_idx);
      CGAL_assertion (inserted);
      m_intersection_graph.set_line (iedge, line_idx);

      for (KSR::size_t support_plane_idx : support_planes_idx)
        support_plane(support_plane_idx).iedges().insert (iedge);
    }
  }

  IVertex source (const IEdge& edge) const
  { return m_intersection_graph.source (edge); }
  IVertex target (const IEdge& edge) const
  { return m_intersection_graph.target (edge); }
  IVertex opposite (const IEdge& edge, const IVertex& ivertex) const
  {
    IVertex out = source (edge);
    if (out == ivertex)
      return target (edge);
    CGAL_assertion (target(edge) == ivertex);
    return out;
  }

  Incident_iedges incident_iedges (const IVertex& ivertex) const
  { return m_intersection_graph.incident_edges(ivertex); }

  const std::set<IEdge>& iedges (KSR::size_t support_plane_idx) const
  { return support_plane(support_plane_idx).iedges(); }

  const KSR::Idx_set& intersected_planes (const IEdge& iedge) const
  { return m_intersection_graph.intersected_planes(iedge); }

  KSR::Idx_set intersected_planes (const IVertex& ivertex, bool keep_bbox = true) const
  {
    KSR::Idx_set out;
    for (const IEdge incident_iedge : incident_iedges (ivertex))
      for (KSR::size_t support_plane_idx : intersected_planes (incident_iedge))
      {
        if (!keep_bbox && support_plane_idx < 6)
          continue;
        out.insert (support_plane_idx);
      }
    return out;
  }

  bool is_iedge (const IVertex& source, const IVertex& target) const
  { return m_intersection_graph.is_edge (source, target); }

  bool is_active (const IEdge& iedge) const
  { return m_intersection_graph.is_active (iedge); }
  bool is_active (const IVertex& ivertex) const
  { return m_intersection_graph.is_active (ivertex); }

  bool is_bbox_iedge (const IEdge& edge) const
  {
    for (KSR::size_t support_plane_idx : m_intersection_graph.intersected_planes(edge))
      if (support_plane_idx < 6)
        return true;
    return false;
  }

  /*******************************
   * Connectivity
   *******************************/

  bool has_ivertex (const PVertex& pvertex) const
  { return support_plane(pvertex).has_ivertex(pvertex.second); }
  IVertex ivertex (const PVertex& pvertex) const
  { return support_plane(pvertex).ivertex(pvertex.second); }

  bool has_iedge (const PVertex& pvertex) const
  { return support_plane(pvertex).has_iedge(pvertex.second); }
  IEdge iedge (const PVertex& pvertex) const
  { return support_plane(pvertex).iedge(pvertex.second); }

  bool has_iedge (const PEdge& pedge) const
  { return support_plane(pedge).has_iedge(pedge.second); }
  IEdge iedge (const PEdge& pedge) const
  { return support_plane(pedge).iedge(pedge.second); }


  void connect (const PVertex& pvertex, const IVertex& ivertex)
  { support_plane(pvertex).set_ivertex (pvertex.second, ivertex); }
  void connect (const PVertex& pvertex, const IEdge& iedge)
  { support_plane(pvertex).set_iedge (pvertex.second, iedge); }
  void connect (const PVertex& a, const PVertex& b, const IEdge& iedge)
  { support_plane(a).set_iedge (a.second, b.second, iedge); }
  void connect (const PEdge& pedge, const IEdge& iedge)
  { support_plane(pedge).set_iedge (pedge.second, iedge); }

  IVertex disconnect_ivertex (const PVertex& pvertex)
  {
    IVertex out = ivertex (pvertex);
    support_plane(pvertex).set_ivertex (pvertex.second, null_ivertex());
    return out;
  }
  IEdge disconnect_iedge (const PVertex& pvertex)
  {
    IEdge out = iedge (pvertex);
    support_plane(pvertex).set_iedge (pvertex.second, null_iedge());
    return out;
  }

  struct Queue_element
  {
    PVertex previous;
    PVertex pvertex;
    bool front;
    bool previous_was_free;

    Queue_element (const PVertex& previous, const PVertex& pvertex, bool front,
                   bool previous_was_free)
      : previous (previous), pvertex (pvertex), front (front),
        previous_was_free(previous_was_free) { }
  };

  std::vector<PVertex> pvertices_around_ivertex (const PVertex& pvertex, const IVertex& ivertex) const
  {
    // std::cout.precision(20);
    std::deque<PVertex> vertices;
    vertices.push_back (pvertex);

    std::cout << "came from: " <<
    str(iedge(pvertex)) << " " << segment_3(iedge(pvertex)) << std::endl;

    std::queue<Queue_element> todo;
    PVertex prev, next;
    std::tie (prev, next) = border_prev_and_next (pvertex);
    // std::cout << "prev in: " << str(prev) << " " << point_3(prev) << std::endl;
    // std::cout << "next in: " << str(next) << " " << point_3(next) << std::endl;
    // std::cout << "curr in: " << str(pvertex) << " " << point_3(pvertex) << std::endl;

    todo.push (Queue_element (pvertex, prev, true, false));
    todo.push (Queue_element (pvertex, next, false, false));

    while (!todo.empty())
    {
      // std::cout << std::endl;
      PVertex previous = todo.front().previous;
      PVertex current = todo.front().pvertex;
      bool front = todo.front().front;
      bool previous_was_free = todo.front().previous_was_free;
      todo.pop();

      IEdge iedge = this->iedge (current);
      bool is_free = (iedge == null_iedge());
      // std::cout << "is free 1: " << is_free << std::endl;

      // std::cout << "iedge: " << segment_3(iedge) << std::endl;
      if (!is_free && source(iedge) != ivertex && target(iedge) != ivertex) {
        // std::cout << "is free 2: " << is_free << std::endl;
        is_free = true;
      }

      if (!is_free)
      {
        IVertex other = source(iedge);
        if (other == ivertex)
          other = target(iedge);
        else
          CGAL_assertion (target(iedge) == ivertex);

        // Filter backwards vertex.
        const Vector_2 dir1 = direction(current);
        // std::cout << "dir1: " << dir1 << std::endl;
        const Vector_2 dir2(
          point_2(current.first, other), point_2(current.first, ivertex));
        // std::cout << "dir2: " << dir2 << std::endl;
        const FT dot_product = dir1 * dir2;
        // std::cout << "dot: " << dot_product << std::endl;

        if (dot_product < FT(0))
        {
          std::cerr << str(current) << " is backwards" << std::endl;
          // std::cout << point_3(current) << std::endl;
          is_free = true;
        }
        if (is_frozen(current)) {
          std::cerr << str(current) << " is frozen" << std::endl;
          // std::cout << point_3(current) << std::endl;
          is_free = true;
        }
        // std::cout << "is free 3: " << is_free << std::endl;
      }

      if (previous_was_free && is_free)
      {
        std::cerr << str(current) << " has no iedge, stopping there" << std::endl;
        // std::cout << point_3(current) << std::endl;
        continue;
      }

      if (is_free)
      {
        std::cerr << str(current) << " has no iedge" << std::endl;
        // std::cout << point_3(current) << std::endl;
      }
      else
      {
        std::cerr << str(current) << " has iedge " << str(iedge)
                  << " from " << str(source(iedge)) << " to " << str(target(iedge)) << std::endl;
        // std::cout << segment_3(iedge) << std::endl;
        // std::cout << point_3(current) << std::endl;
      }

      if (front) {
        vertices.push_front (current);
        // std::cout << "pushed front" << std::endl;
      }
      else {
        vertices.push_back (current);
        // std::cout << "pushed back" << std::endl;
      }

      std::tie (prev, next) = border_prev_and_next (current);

      if (prev == previous)
      {
        CGAL_assertion (next != previous);
        todo.push (Queue_element (current, next, front, is_free));
        // std::cout << "pushed next" << std::endl;
      }
      else {
        todo.push (Queue_element (current, prev, front, is_free));
        // std::cout << "pushed prev" << std::endl;
      }
    }

    std::vector<PVertex> out;
    out.reserve (vertices.size());
    std::copy (vertices.begin(), vertices.end(),
               std::back_inserter (out));

    std::cout << "*** Found vertices:";
    for (const PVertex& pv : out)
      std::cout << " " << str(pv);
    std::cout << std::endl;
    return out;
  }

  /*******************************
   * Conversions
   *******************************/

  Point_2 to_2d (KSR::size_t support_plane_idx, const IVertex& ivertex) const
  { return support_plane(support_plane_idx).to_2d (point_3(ivertex)); }
  Segment_2 to_2d (KSR::size_t support_plane_idx, const Segment_3& segment_3) const
  { return support_plane(support_plane_idx).to_2d (segment_3); }

  Point_2 point_2 (const PVertex& pvertex, FT time) const
  { return support_plane(pvertex).point_2 (pvertex.second, time); }
  Point_2 point_2 (const PVertex& pvertex) const
  { return point_2 (pvertex, m_current_time); }
  Point_2 point_2 (KSR::size_t support_plane_idx, const IVertex& ivertex) const
  { return support_plane(support_plane_idx).to_2d(point_3 (ivertex)); }

  Segment_2 segment_2 (KSR::size_t support_plane_idx, const IEdge& iedge) const
  { return support_plane(support_plane_idx).to_2d(segment_3(iedge)); }

  Point_3 to_3d (KSR::size_t support_plane_idx, const Point_2& point_2) const
  { return support_plane(support_plane_idx).to_3d (point_2); }

  Point_3 point_3 (const PVertex& pvertex, FT time) const
  { return support_plane(pvertex).point_3 (pvertex.second, time); }
  Point_3 point_3 (const PVertex& pvertex) const
  { return point_3 (pvertex, m_current_time); }
  Point_3 point_3 (const IVertex& vertex) const
  { return m_intersection_graph.point_3 (vertex); }

  Segment_3 segment_3 (const PEdge& pedge, FT time) const
  { return support_plane(pedge).segment_3 (pedge.second, time); }
  Segment_3 segment_3 (const PEdge& pedge) const
  { return segment_3 (pedge, m_current_time); }
  Segment_3 segment_3 (const IEdge& edge) const
  { return m_intersection_graph.segment_3 (edge); }

  /*******************************
   * Predicates
   *******************************/

  // Check if there is a collision with another polygon.
  std::pair<bool, bool> collision_occured (
    const PVertex& pvertex, const IEdge& iedge) const {

    // const FT tol = FT(1) / FT(100000);
    bool collision = false;
    for (const auto support_plane_idx : intersected_planes(iedge)) {
      if (support_plane_idx < 6)
        return std::make_pair(true, true);

      for (const auto pedge : pedges(support_plane_idx)) {
        if (this->iedge(pedge) == iedge) {
          const auto pedge_segment = Segment_3(point_3(source(pedge)), point_3(target(pedge)));

          const Segment_3 source_to_pvertex(pedge_segment.source(), point_3(pvertex));
          // if (CGAL::sqrt(source_to_pvertex.squared_length()) < tol) {
          //   // std::cout << "WARNING: POINTS ARE ALMOST EQUAL!" << std::endl;
          //   collision = true;
          //   break;
          // }

          // std::cout << point_3(source(pedge)) << std::endl;
          // std::cout << point_3(target(pedge)) << std::endl;
          // std::cout << point_3(pvertex) << std::endl;

          const FT dot_product = pedge_segment.to_vector() * source_to_pvertex.to_vector();
          if (dot_product < FT(0)) {
            continue;
          }
          // std::cout << source_to_pvertex.squared_length() << std::endl;
          // std::cout << pedge_segment.squared_length() << std::endl;

          if (pedge_segment.squared_length() == FT(0))
            std::cout << "ERROR: SOURCE_TO_PVERTEX/PEDGE SEGMENT SQ LENGTH = "
            << source_to_pvertex.squared_length() << std::endl;
          CGAL_assertion(pedge_segment.squared_length() != FT(0));

          // if (pedge_segment.squared_length() == FT(0)) {
          //   if (pedge_segment.source() == point_3(pvertex)) {
          //     collision = false;
          //     break;
          //   }
          // }
          // CGAL_assertion(pedge_segment.squared_length() != FT(0));

          if (source_to_pvertex.squared_length() <= pedge_segment.squared_length()) {
            collision = true;
            break;
          }
        }
      }
    }
    return std::make_pair(collision, false);
  }

  std::pair<bool, bool> is_occupied(
    const PVertex& pvertex,
    const IEdge& query_iedge) {

    // std::cout << str(query_iedge) << " " << segment_3(query_iedge) << std::endl;
    KSR::size_t num_adjacent_faces = 0;
    for (const auto plane_idx : intersected_planes(query_iedge)) {
      if (plane_idx == pvertex.first) continue; // current plane
      if (plane_idx < 6) return std::make_pair(true, true); // bbox plane

      for (const auto pedge : pedges(plane_idx)) {
        // std::cout << str(iedge(pedge)) << std::endl;
        if (iedge(pedge) == query_iedge) {
          const auto& m = mesh(plane_idx);
          const auto he = m.halfedge(pedge.second);
          const auto op = m.opposite(he);
          const auto face1 = m.face(he);
          const auto face2 = m.face(op);
          if (face1 != Support_plane::Mesh::null_face()) ++num_adjacent_faces;
          if (face2 != Support_plane::Mesh::null_face()) ++num_adjacent_faces;
        }
      }
    }

    std::cout << "num adjacent faces: " << num_adjacent_faces << std::endl;
    if (num_adjacent_faces <= 1)
      return std::make_pair(false, false);
    return std::make_pair(true, false);
  }

  /*******************************
   * Operations on polygons
   *******************************/

  PVertex crop_polygon (const PVertex& pvertex, const IEdge& iedge)
  {
    std::cout << "*** Cropping " << str(pvertex) << " along " << str(iedge) << std::endl;

    Point_2 future_point_a, future_point_b;
    Vector_2 direction_a, direction_b;

    compute_future_points_and_directions (pvertex, iedge,
                                          future_point_a, future_point_b,
                                          direction_a, direction_b);

    PEdge pedge (pvertex.first, support_plane(pvertex).split_vertex(pvertex.second));
    CGAL_assertion (source(pedge) == pvertex || target(pedge) == pvertex);

    PVertex other = opposite(pedge, pvertex);

    std::cout << "*** New edge " << str(pedge) << " between " << str(pvertex)
                     << " and " << str(other) << std::endl;

    connect (pedge, iedge);
    connect (pvertex, iedge);
    connect (other, iedge);

    support_plane(pvertex).set_point (pvertex.second, future_point_a);
    support_plane(other).set_point (other.second, future_point_b);

    direction(pvertex) = direction_a;
    direction(other) = direction_b;

    // std::cout << "pvertex: " << point_3(pvertex) << std::endl;
    // std::cout << "pvertex dir: " << direction_a << std::endl;
    // std::cout << "other: " << point_3(other) << std::endl;
    // std::cout << "other dir: " << direction_b << std::endl;

    std::cout << "New vertices: " << str(other) << std::endl;
    return other;
  }

  std::array<PVertex, 3> propagate_polygon (
    const unsigned int last_k,
    const PVertex& pvertex, const IEdge& iedge)
  {
    std::cout << "*** Propagating " << str(pvertex) << " along " << str(iedge) << std::endl;

    Point_2 original_point = point_2 (pvertex, 0);
    Vector_2 original_direction = direction(pvertex);

    PVertex other = crop_polygon (pvertex, iedge);

    PVertex propagated = add_pvertex (pvertex.first, original_point);
    direction(propagated) = original_direction;

    std::array<PVertex, 3> pvertices;

    pvertices[0] = pvertex;
    pvertices[1] = other;
    pvertices[2] = propagated;

    PFace new_pface = add_pface (pvertices);
    this->k(new_pface) = last_k;
    CGAL_assertion (new_pface.second != Face_index());

    std::cout << "*** New face " << lstr(new_pface) << std::endl;

    return pvertices;
  }

  void crop_polygon (const PVertex& pv0, const PVertex& pv1, const IEdge& iedge)
  {
    std::cout << "*** Cropping " << str(pv0) << "/" << str(pv1) << " along " << str(iedge) << std::endl;

    std::cout.precision(20);
    // std::cout << "pv0: " << point_3(pv0) << std::endl;
    // std::cout << "pv1: " << point_3(pv1) << std::endl;

    Point_2 future_point;
    Vector_2 future_direction;
    // const Line_2 iedge_line = segment_2(pv0.first, iedge).supporting_line();
    CGAL_assertion(pv0.first == pv1.first);

    {
      // const Point_2 pinit = iedge_line.projection(point_2(pv0, m_current_time));
      // const Point_2 future_point = iedge_line.projection(point_2(pv0, m_current_time + FT(1)));

      const PVertex prev(pv0.first, support_plane(pv0).prev(pv0.second));
      const PVertex next(pv0.first, support_plane(pv0).next(pv0.second));

      if (prev == pv1)
        compute_future_point_and_direction(0, pv0, next, iedge, future_point, future_direction);
      else {
        CGAL_assertion(next == pv1);
        compute_future_point_and_direction(0, pv0, prev, iedge, future_point, future_direction);
      }

      direction(pv0) = future_direction;
      std::cout << "pv0 dir: " << direction(pv0) << std::endl;
      support_plane(pv0).set_point(pv0.second, future_point);
      connect(pv0, iedge);
    }

    {
      // const Point_2 pinit = iedge_line.projection(point_2(pv1, m_current_time));
      // const Point_2 future_point = iedge_line.projection(point_2(pv1, m_current_time + FT(1)));

      const PVertex prev(pv1.first, support_plane(pv1).prev(pv1.second));
      const PVertex next(pv1.first, support_plane(pv1).next(pv1.second));

      if (prev == pv0)
        compute_future_point_and_direction(0, pv1, next, iedge, future_point, future_direction);
      else {
        CGAL_assertion(next == pv0);
        compute_future_point_and_direction(0, pv1, prev, iedge, future_point, future_direction);
      }

      direction(pv1) = future_direction;
      std::cout << "pv1 dir: " << direction(pv1) << std::endl;
      support_plane(pv1).set_point(pv1.second, future_point);
      connect(pv1, iedge);
    }

    const PEdge pedge(pv0.first, support_plane(pv0).edge(pv0.second, pv1.second));
    connect(pedge, iedge);
  }

  std::pair<PVertex, PVertex> propagate_polygon(
    const unsigned int, // last_k,
    const PVertex& pvertex, const PVertex& pother, const IEdge& iedge)
  {
    std::cout << "*** Propagating " << str(pvertex) << "/" << str(pother) << " along " << str(iedge) << std::endl;
    CGAL_assertion_msg(false, "TODO: PROPAGATE POLYGON VIA THE EDGE!");
    return std::make_pair(null_pvertex(), null_pvertex());
  }

  bool transfer_vertex (const PVertex& pvertex, const PVertex& pother)
  {
    std::cout << "*** Transfering " << str(pother) << " through " << str(pvertex) << std::endl;

    // If pvertex is adjacent to one or two
    PFace source_face, target_face;
    std::tie (source_face, target_face) = pfaces_of_pvertex (pvertex);

    PFace common_pface = pface_of_pvertex (pother);

    if (common_pface == target_face)
      std::swap (source_face, target_face);
    CGAL_assertion (common_pface == source_face);

    std::cout << "*** Initial faces: " << lstr(source_face)
                     << " and " << lstr(target_face) << std::endl;

    PVertex pthird = next(pother);
    if (pthird == pvertex)
      pthird = prev(pother);

    if (target_face == null_pface())
    {
      Vector_2 new_direction;

      Line_2 iedge_line = segment_2(pother.first, iedge(pvertex)).supporting_line();
      Point_2 pinit = iedge_line.projection(point_2 (pother, m_current_time));

      Line_2 future_line (point_2 (pother, m_current_time + 1),
                          point_2 (pthird, m_current_time + 1));

      Point_2 future_point = KSR::intersection_2<Point_2> (future_line, iedge_line);

      direction(pvertex) = Vector_2 (pinit, future_point);
      support_plane(pvertex).set_point (pvertex.second,
                                        pinit - direction(pvertex) * m_current_time);

      Halfedge_index hi = mesh(pvertex).halfedge(pother.second, pvertex.second);
      CGAL::Euler::join_vertex(hi, mesh(pvertex));
    }
    else
    {
      IEdge iedge = disconnect_iedge (pvertex);
      // std::cerr << "Disconnect " << str(pvertex) << " from " << str(iedge) << std::endl;

      PEdge pedge = null_pedge();
      for (PEdge pe : pedges_around_pvertex (pvertex))
        if (this->iedge(pe) == iedge)
        {
          pedge = pe;
          break;
        }
      CGAL_assertion (pedge != null_pedge());

      Halfedge_index hi = mesh(pedge).halfedge(pedge.second);
      if (mesh(pedge).face(hi) != common_pface.second)
        hi = mesh(pedge).opposite(hi);
      CGAL_assertion (mesh(pedge).face(hi) == common_pface.second);

      if (mesh(pedge).target(hi) == pvertex.second)
      {
        // std::cerr << "Shift target" << std::endl;
        CGAL::Euler::shift_target (hi, mesh(pedge));
      }
      else
      {
        CGAL_assertion (mesh(pedge).source(hi) == pvertex.second);
        // std::cerr << "Shift source" << std::endl;
        CGAL::Euler::shift_source (hi, mesh(pedge));
      }

      Vector_2 new_direction;

      Line_2 iedge_line = segment_2(pother.first, iedge).supporting_line();
      Point_2 pinit = iedge_line.projection(point_2 (pother, m_current_time));

      direction(pvertex) = direction(pother);
      support_plane(pother).set_point (pvertex.second,
                                       pinit - direction(pvertex) * m_current_time);

      Line_2 future_line (point_2 (pvertex, m_current_time + 1),
                          point_2 (pthird, m_current_time + 1));

      Point_2 future_point = KSR::intersection_2<Point_2> (future_line, iedge_line);

      direction(pother) = Vector_2 (pinit, future_point);
      support_plane(pother).set_point (pother.second,
                                       pinit - direction(pother) * m_current_time);

      // std::cerr << "Connect " << str(pother) << " to " << str(iedge) << std::endl;
      connect (pother, iedge);
    }

    std::cout << "*** New faces: " << lstr(source_face)
                     << " and " << lstr(target_face) << std::endl;

    return (target_face != null_pface());
  }

  void merge_pvertices (const PVertex& pvertex, const PVertex& pother)
  {
    std::cout << "*** Merging " << str(pvertex) << " with " << str(pother) << std::endl;

    Halfedge_index hi = mesh(pvertex).halfedge(pother.second, pvertex.second);
    disconnect_ivertex (pother);
    CGAL::Euler::join_vertex(hi, mesh(pvertex));
  }

  std::vector<PVertex> merge_pvertices_on_ivertex (const FT min_time,
                                                   const FT max_time,
                                                   std::vector<PVertex>& pvertices,
                                                   const IVertex& ivertex,
                                                   std::vector<IEdge>& crossed)
  {
    crossed.clear();
    KSR::size_t support_plane_idx = pvertices.front().first;

    PVertex prev = pvertices.front();
    PVertex next = pvertices.back();

    IEdge prev_iedge = null_iedge(), next_iedge = null_iedge();

    // std::ofstream("came_from.polylines.txt")
    // << "2 " << segment_3(iedge(pvertices[1])) << std::endl;

    std::cout << "came from: " <<
    str(iedge(pvertices[1])) << " " << segment_3(iedge(pvertices[1])) << std::endl;

    // Copy front/back to remember position/direction.
    PVertex front, back;
    if (pvertices.size() < 3) {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    } else if (pvertices.size() == 3 || pvertices.size() == 4) {

      // BUG: In this case, the point that is duplicated twice is not always copied.
      // To fix it, we copy the second point not from the original vertex but from the first
      // copy of that vertex.

      const auto& initial = pvertices[1];
      front = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial.second));
      support_plane(support_plane_idx).set_point(
        front.second, support_plane(support_plane_idx).get_point(initial.second));
      back  = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(front.second));
      support_plane(support_plane_idx).set_point(
        back.second, support_plane(support_plane_idx).get_point(front.second));

    } else if (pvertices.size() >= 5) {

      const auto& initial1 = pvertices[1];
      front = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial1.second));
      support_plane(support_plane_idx).set_point(
        front.second, support_plane(support_plane_idx).get_point(initial1.second));

      const auto& initial2 = pvertices[pvertices.size() - 2];
      back  = PVertex(support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(initial2.second));
      support_plane(support_plane_idx).set_point(
        back.second, support_plane(support_plane_idx).get_point(initial2.second));

    } else {
      CGAL_assertion_msg(false, "ERROR: INVALID CASE!");
    }

    // auto pvertex_to_point =
    //   [&](const PVertex& a) -> Point_2 {
    //     return point_2(a);
    //   };

    // PFace fprev = pface_of_pvertex(prev);
    // Point_2 pprev = CGAL::centroid
    //   (boost::make_transform_iterator (pvertices_of_pface(fprev).begin(), pvertex_to_point),
    //    boost::make_transform_iterator (pvertices_of_pface(fprev).end(), pvertex_to_point));
    // PFace fnext = pface_of_pvertex(next);
    // Point_2 pnext = CGAL::centroid
    //   (boost::make_transform_iterator (pvertices_of_pface(fnext).begin(), pvertex_to_point),
    //    boost::make_transform_iterator (pvertices_of_pface(fnext).end(), pvertex_to_point));

    bool was_swapped = false;
    // if (CGAL::orientation(pprev, point_2(support_plane_idx, ivertex), pnext) == CGAL::LEFT_TURN) {
    //   std::cout << "Swapped!" << std::endl;
    //   was_swapped = true;
    //   std::swap(prev, next);
    //   std::swap(front, back);
    // }

    // Freeze vertices.
    for (std::size_t i = 1; i < pvertices.size() - 1; ++i) {
      PVertex& pvertex = pvertices[i];
      Point_2 point = point_2(support_plane_idx, ivertex);
      support_plane(pvertex).direction(pvertex.second) = CGAL::NULL_VECTOR;
      support_plane(pvertex).set_point(pvertex.second, point);
    }

    PVertex pvertex = pvertices[1];
    connect (pvertex, ivertex);

    std::cout << "*** Frozen vertex: " << str(pvertex) << std::endl;
    // std::cout << point_3(pvertex) << std::endl;
    // std::cout << "*** Removed vertices:";

    // Join vertices
    for (std::size_t i = 2; i < pvertices.size() - 1; ++ i)
    {
      // std::cout << " " << str(pvertices[i]) << std::endl;
      // std::cout << point_3(pvertices[i]) << std::endl;

      Halfedge_index hi = mesh(support_plane_idx).halfedge(pvertices[i].second, pvertex.second);
      disconnect_ivertex (pvertices[i]);
      CGAL::Euler::join_vertex(hi, mesh(support_plane_idx));
    }
    // std::cout << std::endl;

    Incident_iedges i_iedges = incident_iedges (ivertex);
    std::vector<std::pair<IEdge, Direction_2> > iedges;
    std::copy (i_iedges.begin(), i_iedges.end(),
               boost::make_function_output_iterator
               ([&](const IEdge& ie) -> void
                {
                  if (intersected_planes(ie).find (support_plane_idx)
                      == intersected_planes(ie).end())
                    return;

                  Direction_2 dir (point_2 (support_plane_idx, opposite (ie, ivertex))
                                   - point_2 (support_plane_idx, ivertex));
                  iedges.push_back (std::make_pair (ie, dir));
                }));

    std::sort (iedges.begin(), iedges.end(),
               [&](const std::pair<IEdge, Direction_2>& a,
                   const std::pair<IEdge, Direction_2>& b) -> bool
               {
                 return a.second < b.second;
               });
    CGAL_assertion(iedges.size() != 0);

    std::cout.precision(20);
    std::cout << "Prev = "  << point_3 (prev)  << " / " << direction (prev)  << std::endl
              << "Front = " << point_3 (front) << " / " << direction (front) << std::endl
              << "Back = "  << point_3 (back)  << " / " << direction (back)  << std::endl
              << "Next = "  << point_3 (next)  << " / " << direction (next)  << std::endl;

    // std::cout << (iedge(next) != null_iedge()) << std::endl;
    // std::cout << "source: " << point_3(source(iedge(next))) << std::endl;
    // std::cout << "target: " << point_3(target(iedge(next))) << std::endl;
    // std::cout << "ivertex: " << point_3(ivertex) << std::endl;

    bool back_constrained = false;
    if ((iedge(next) != null_iedge()
         && (source(iedge(next)) == ivertex || target(iedge(next)) == ivertex))
        || (this->ivertex(next) != null_ivertex()
            && is_iedge (this->ivertex(next), ivertex)))
      back_constrained = true;

    bool front_constrained = false;
    if ((iedge(prev) != null_iedge()
         && (source(iedge(prev)) == ivertex || target(iedge(prev)) == ivertex))
        || (this->ivertex(prev) != null_ivertex()
            && is_iedge (this->ivertex(prev), ivertex)))
      front_constrained = true;

    std::vector<PVertex> new_vertices;
    if (back_constrained && front_constrained) // Closing case
    {
      std::cout << "*** Closing case" << std::endl;
    }
    else if (back_constrained) // Border case
    {
      std::cout << "*** Back border case" << std::endl;

      CGAL_assertion(has_iedge(pvertex));
      // std::ofstream("limit.polylines.txt")
      // << "2 " << segment_3(iedge(pvertex)) << std::endl;
      const KSR::size_t other_side_limit = line_idx(pvertex);

      // const Direction_2 dir(point_2(prev) - point_2(pvertex));

      const FT prev_time = last_event_time(prev);
      CGAL_assertion(prev_time < m_current_time);
      CGAL_assertion(prev_time >= FT(0));

      const auto pp_last = point_2(prev, prev_time);
      const auto pp_curr = point_2(prev, m_current_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      const auto tmp_prev = pp_curr - dirp / FT(10);

      const Direction_2 tmp_dir(tmp_prev - point_2(pvertex.first, ivertex));
      // std::cout << to_3d(prev.first, tmp_prev) << std::endl;

      std::reverse(iedges.begin(), iedges.end());

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++i) {
        if (tmp_dir.counterclockwise_in_between(
          iedges[(i + 1) % iedges.size()].second, iedges[i].second)) {

          first_idx = (i + 1) % iedges.size();
          break;
        }
      }

      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;

        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = collision_occured(pvertex, iedge);
        bool limit_reached = (line_idx(iedge) == other_side_limit);
        std::cout << "limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;

        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);

        if (limit_reached || bbox_reached) {
          break;
        }

        iedge_idx = (iedge_idx + 1) % iedges.size();

        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: BACK WHY SO MANY ITERATIONS?");
        }
        ++iter;
      }

      CGAL_assertion(crossed.size() != 0);
      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;
      for (const auto& iedge : crossed)
        std::cout << segment_3(iedge) << std::endl;

      std::vector<Point_2> future_points(crossed.size());
      std::vector<Vector_2> future_directions(crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        const bool is_parallel = compute_future_point_and_direction(
          i, back, prev, crossed[i], future_points[i], future_directions[i]);
        if (is_parallel) {
          if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
            prev_iedge = crossed[i];
          }
        }
      }

      PVertex previous = null_pvertex();
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        if (i == 0) // crop
        {
          std::cout << "Cropping" << std::endl;
          PVertex cropped; Point_2 future_point; Vector_2 future_direction;
          if (prev_iedge != null_iedge() && prev_iedge == crossed[i]) {
            std::cout << "prev parallel case" << std::endl;

            cropped = prev;
            const auto pair = this->border_prev_and_next(prev);
            const auto pprev = pair.first;
            compute_future_point_and_direction(
              i, prev, pprev, prev_iedge, future_point, future_direction);

          } else {
            std::cout << "standard case" << std::endl;
            cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, prev.second));
            future_point = future_points[i];
            future_direction = future_directions[i];
          }

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
          new_vertices.push_back(cropped);

          connect(pedge, crossed[i]);
          connect(cropped, crossed[i]);

          support_plane(cropped).set_point(cropped.second, future_point);
          direction(cropped) = future_direction;
          previous = cropped;
          // std::cerr << "cropped point -> direction: " << point_2 (cropped) << " -> " << direction(cropped) << std::endl;
          std::cout << "cropped: " << point_3(cropped) << std::endl;
        }
        else // create triangle face
        {
          bool is_occupied_edge, bbox_reached;
          std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, crossed[0]);
          // std::tie(is_occupied_edge, bbox_reached) = collision_occured(pvertex, crossed[0]);
          std::cout << "is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;

          // Stop.
          const auto pface = pface_of_pvertex(pvertex);
          std::cout << "k intersections: " << this->k(pface) << std::endl;
          if (bbox_reached) {
            std::cout << "stop bbox" << std::endl;
            CGAL_assertion_msg(false, "ERROR: THIS CASE CANNOT HAPPEN!");
            break;
          } else if (is_occupied_edge && this->k(pface) == 1) {
            std::cout << "stop k" << std::endl;
            break;
          }

          // Create a new face.
          std::cout << "adding new face!" << std::endl;
          if (is_occupied_edge && this->k(pface) > 1) {
            std::cout << "continue k > 1" << std::endl;
            this->k(pface)--;
          } else {
            std::cout << "continue k = 1" << std::endl;
          }
          CGAL_assertion(this->k(pface) >= 1);

          PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back(propagated);

          std::cout << "propagated: " << point_3(propagated) << std::endl;

          PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, propagated, previous});
          // this->k(new_pface) = k;
          this->k(new_pface) = this->k(pface);
          previous = propagated;

          PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
          connect(pedge, crossed[i]);
          connect(propagated, crossed[i]);
        }
      }
    }
    else if (front_constrained) // Border case
    {
      std::cout << "*** Front border case" << std::endl;

      CGAL_assertion(has_iedge(pvertex));
      // std::ofstream("limit.polylines.txt")
      // << "2 " << segment_3(iedge(pvertex)) << std::endl;
      const KSR::size_t other_side_limit = line_idx(pvertex);

      // const Direction_2 dir(point_2(next) - point_2(pvertex));

      const FT next_time = last_event_time(next);
      CGAL_assertion(next_time < m_current_time);
      CGAL_assertion(next_time >= FT(0));

      const auto pn_last = point_2(next, next_time);
      const auto pn_curr = point_2(next, m_current_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      const auto tmp_next = pn_curr - dirn / FT(10);

      const Direction_2 tmp_dir(tmp_next - point_2(pvertex.first, ivertex));
      // std::cout << to_3d(next.first, tmp_next) << std::endl;

      if (was_swapped) {
        std::reverse(iedges.begin(), iedges.end());
      }

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++ i)
      {
        if (!was_swapped) {
          if (tmp_dir.counterclockwise_in_between(
            iedges[i].second, iedges[(i + 1) % iedges.size()].second)) {
            first_idx = (i + 1) % iedges.size();
            break;
          }
        } else {
          if (tmp_dir.counterclockwise_in_between(
            iedges[(i + 1) % iedges.size()].second, iedges[i].second)) {
            first_idx = (i + 1) % iedges.size();
            break;
          }
        }
      }

      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true) {
        const IEdge& iedge = iedges[iedge_idx].first;

        bool collision, bbox_reached;
        std::tie(collision, bbox_reached) = collision_occured(pvertex, iedge);
        bool limit_reached = (line_idx(iedge) == other_side_limit);
        std::cout << "limit/bbox: " << limit_reached << "/" << bbox_reached << std::endl;

        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);

        if (limit_reached || bbox_reached) {
          break;
        }

        iedge_idx = (iedge_idx + 1) % iedges.size();

        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: FRONT WHY SO MANY ITERATIONS?");
        }
        ++iter;
      }

      CGAL_assertion(crossed.size() != 0);
      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;
      for (const auto& iedge : crossed)
        std::cout << segment_3(iedge) << std::endl;

      std::vector<Point_2> future_points(crossed.size());
      std::vector<Vector_2> future_directions(crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        const bool is_parallel = compute_future_point_and_direction(
          i, front, next, crossed[i], future_points[i], future_directions[i]);

        if (is_parallel) {
          if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
            next_iedge = crossed[i];
          }
        }
      }

      PVertex previous = null_pvertex();
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        if (i == 0) // crop
        {
          std::cout << "Cropping" << std::endl;
          PVertex cropped; Point_2 future_point; Vector_2 future_direction;
          if (next_iedge != null_iedge() && next_iedge == crossed[i]) {
            std::cout << "next parallel case" << std::endl;

            cropped = next;
            const auto pair = this->border_prev_and_next(next);
            const auto nnext = pair.second;
            compute_future_point_and_direction(
              i, next, nnext, next_iedge, future_point, future_direction);

          } else {
            std::cout << "standard case" << std::endl;
            cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, next.second));
            future_point = future_points[i];
            future_direction = future_directions[i];
          }

          const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
          CGAL_assertion(cropped != pvertex);
          new_vertices.push_back(cropped);

          connect(pedge, crossed[i]);
          connect(cropped, crossed[i]);

          support_plane(cropped).set_point(cropped.second, future_point);
          direction(cropped) = future_direction;
          previous = cropped;
          // std::cerr << point_2 (cropped) << " -> " << direction(cropped) << std::endl;
        }
        else // create triangle face
        {
          bool is_occupied_edge, bbox_reached;
          std::tie(is_occupied_edge, bbox_reached) = is_occupied(pvertex, crossed[0]);
          // std::tie(is_occupied_edge, bbox_reached) = collision_occured(pvertex, crossed[0]);
          std::cout << "is already occupied / bbox: " << is_occupied_edge << "/" << bbox_reached << std::endl;

          // Stop.
          const auto pface = pface_of_pvertex(pvertex);
          std::cout << "k intersections: " << this->k(pface) << std::endl;
          if (bbox_reached) {
            std::cout << "stop bbox" << std::endl;
            CGAL_assertion_msg(false, "ERROR: THIS CASE CANNOT HAPPEN!");
            break;
          } else if (is_occupied_edge && this->k(pface) == 1) {
            std::cout << "stop k" << std::endl;
            break;
          }

          // Create a new face.
          std::cout << "adding new face!" << std::endl;
          if (is_occupied_edge && this->k(pface) > 1) {
            std::cout << "continue k > 1" << std::endl;
            this->k(pface)--;
          } else {
            std::cout << "continue k = 1" << std::endl;
          }
          CGAL_assertion(this->k(pface) >= 1);

          PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back(propagated);

          std::cout << "propagated: " << point_3(propagated) << std::endl;

          PFace new_pface = add_pface(std::array<PVertex, 3>{pvertex, previous, propagated});
          // this->k(new_pface) = k;
          this->k(new_pface) = this->k(pface);
          previous = propagated;

          PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, propagated.second));
          connect(pedge, crossed[i]);
          connect(propagated, crossed[i]);
        }
      }
    }
    else // Open case
    {
      std::cout << "*** Open case" << std::endl;

      // const Direction_2 dir_prev(point_2(prev) - point_2(pvertex));
      // const Direction_2 dir_next(point_2(next) - point_2(pvertex));

      const FT prev_time = last_event_time(prev);
      const FT next_time = last_event_time(next);
      CGAL_assertion(prev_time < m_current_time);
      CGAL_assertion(next_time < m_current_time);
      CGAL_assertion(prev_time >= FT(0));
      CGAL_assertion(next_time >= FT(0));

      const auto pp_last = point_2(prev, prev_time);
      const auto pp_curr = point_2(prev, m_current_time);
      const auto dirp = Vector_2(pp_last, pp_curr);
      const auto tmp_prev = pp_curr - dirp / FT(10);

      const auto pn_last = point_2(next, next_time);
      const auto pn_curr = point_2(next, m_current_time);
      const auto dirn = Vector_2(pn_last, pn_curr);
      const auto tmp_next = pn_curr - dirn / FT(10);

      const Direction_2 dir_prev(tmp_prev - point_2(pvertex.first, ivertex));
      const Direction_2 dir_next(tmp_next - point_2(pvertex.first, ivertex));

      // std::cout << to_3d(prev.first, tmp_prev) << std::endl;
      // std::cout << to_3d(next.first, tmp_next) << std::endl;

      // std::cout << "initial iedges: " << std::endl;
      // for (const auto& iedge : iedges) {
      //   std::cout << segment_3(iedge.first) << std::endl;
      // }

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++i) {
        if (dir_next.counterclockwise_in_between(
          iedges[i].second, iedges[(i + 1) % iedges.size()].second)) {

          first_idx = (i + 1) % iedges.size();
          break;
        }
      }

      CGAL_assertion(first_idx != KSR::no_element());
      crossed.clear();

      // std::ofstream("first.polylines.txt")
      // << "2 " << segment_3(iedges[first_idx].first) << std::endl;

      KSR::size_t iedge_idx = first_idx;
      std::size_t iter = 0;
      while (true)
      {
        const IEdge& iedge = iedges[iedge_idx].first;
        const Direction_2& dir = iedges[iedge_idx].second;

        if (!dir.counterclockwise_in_between (dir_next, dir_prev))
          break;

        // std::ofstream("next" + std::to_string(iter) + ".polylines.txt")
        // << "2 " << segment_3(iedge) << std::endl;
        crossed.push_back(iedge);

        iedge_idx = (iedge_idx + 1) % iedges.size();

        if (iter == 100) {
          CGAL_assertion_msg(false, "ERROR: OPEN WHY SO MANY ITERATIONS?");
        }
        ++iter;
      }

      CGAL_assertion(crossed.size() != 0);
      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;
      for (const auto& iedge : crossed)
        std::cout << segment_3(iedge) << std::endl;

      std::vector<Point_2> future_points(crossed.size());
      std::vector<Vector_2> future_directions(crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++i) {
        const bool is_parallel = compute_future_point_and_direction(
          pvertex, prev, next, crossed[i], future_points[i], future_directions[i]);

        if (is_parallel) {
          if (is_intersecting_iedge(min_time, max_time, prev, crossed[i])) {
            prev_iedge = crossed[i];
          }
          if (is_intersecting_iedge(min_time, max_time, next, crossed[i])) {
            next_iedge = crossed[i];
          }
        }
      }

      {
        PVertex cropped; Point_2 future_point; Vector_2 future_direction;
        if (next_iedge != null_iedge() && next_iedge == crossed.front()) {
          std::cout << "next parallel case" << std::endl;

          cropped = next;
          const auto pair = this->border_prev_and_next(next);
          const auto nnext = pair.second;
          compute_future_point_and_direction(
            0, next, nnext, next_iedge, future_point, future_direction);

        } else {
          std::cout << "standard case" << std::endl;
          cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, next.second));
          future_point = future_points.front();
          future_direction = future_directions.front();
        }

        const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
        new_vertices.push_back(cropped);

        connect(pedge, crossed.front());
        connect(cropped, crossed.front());

        support_plane(cropped).set_point(cropped.second, future_point);
        direction(cropped) = future_direction;
        std::cout << direction(cropped) << std::endl;
        std::cout << "cropped 1: " << point_3(cropped) << std::endl;
      }

      for (std::size_t i = 1; i < crossed.size() - 1; ++i)
      {
        const PVertex propagated = add_pvertex(pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        connect(propagated, crossed[i]);
        new_vertices.push_back(propagated);
        std::cout << "propagated " << std::to_string(i) << ": " << point_3(propagated) << std::endl;
      }

      {
        PVertex cropped; Point_2 future_point; Vector_2 future_direction;
        if (prev_iedge != null_iedge() && prev_iedge == crossed.back()) {
          std::cout << "prev parallel case" << std::endl;

          cropped = prev;
          const auto pair = this->border_prev_and_next(prev);
          const auto pprev = pair.first;
          compute_future_point_and_direction(
            0, prev, pprev, prev_iedge, future_point, future_direction);

        } else {
          std::cout << "standard case" << std::endl;
          cropped = PVertex(support_plane_idx, support_plane(pvertex).split_edge(pvertex.second, prev.second));
          future_point = future_points.back();
          future_direction = future_directions.back();
        }

        const PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, cropped.second));
        new_vertices.push_back(cropped);

        connect(pedge, crossed.back());
        connect(cropped, crossed.back());

        support_plane(cropped).set_point(cropped.second, future_point);
        direction(cropped) = future_direction;
        std::cout << direction(cropped) << std::endl;
        std::cout << "cropped 2: " << point_3(cropped) << std::endl;
      }

      std::cerr << new_vertices.size() << " new vertice(s)" << std::endl;

      bool is_occupied_edge_back, bbox_reached_back;
      std::tie(is_occupied_edge_back, bbox_reached_back) = is_occupied(pvertex, crossed.back());
      // std::tie(is_occupied_edge_back, bbox_reached_back) = collision_occured(pvertex, crossed.back());
      std::cout << "is already occupied back / bbox: " << is_occupied_edge_back << "/" << bbox_reached_back << std::endl;

      bool is_occupied_edge_front, bbox_reached_front;
      std::tie(is_occupied_edge_front, bbox_reached_front) = is_occupied(pvertex, crossed.front());
      // std::tie(is_occupied_edge_front, bbox_reached_front) = collision_occured(pvertex, crossed.front());
      std::cout << "is already occupied front / bbox: " << is_occupied_edge_front << "/" << bbox_reached_front << std::endl;

      const auto pface = pface_of_pvertex(pvertex);
      std::cout << "k intersections: " << this->k(pface) << std::endl;
      if (bbox_reached_back) {

        CGAL_assertion(bbox_reached_front);
        std::cout << "stop bbox back" << std::endl;

      } else if (bbox_reached_front) {

        CGAL_assertion(bbox_reached_back);
        std::cout << "stop bbox front" << std::endl;

      } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) == 1) {

        add_new_faces(this->k(pface), pvertex, new_vertices, pface);
        std::cout << "back && front k = 1" << std::endl;

      } else if ((is_occupied_edge_back && is_occupied_edge_front) && this->k(pface) > 1) {

        // this->k(pface)--;
        // CGAL_assertion(this->k(pface) >= 1);
        add_new_faces(this->k(pface), pvertex, new_vertices, pface);
        std::cout << "back && front k > 1" << std::endl;

      } else if ((!is_occupied_edge_back && !is_occupied_edge_front)) {

        add_new_faces(this->k(pface), pvertex, new_vertices, pface);
        std::cout << "!back && !front" << std::endl;

      } else if (is_occupied_edge_back && !is_occupied_edge_front) {

        add_new_faces(this->k(pface), pvertex, new_vertices, pface);
        std::cout << "back && !front" << std::endl;

      } else if (!is_occupied_edge_back && is_occupied_edge_front) {

        add_new_faces(this->k(pface), pvertex, new_vertices, pface);
        std::cout << "!back && front" << std::endl;

        // if (this->k(pface) > 1) {
        //   this->k(pface)--;
        //   CGAL_assertion(this->k(pface) >= 1);
        //   add_new_faces(this->k(pface), pvertex, new_vertices, pface);
        // }

      } else {
        CGAL_assertion_msg(false, "TODO: ADD NEW OPEN CASE! DO NOT FORGET TO UPDATE K!");
      }

      for (std::size_t i = 1; i < crossed.size() - 1; ++i) {
        PEdge pedge(support_plane_idx, support_plane(pvertex).edge(pvertex.second, new_vertices[i].second));
        connect(pedge, crossed[i]);
        connect(new_vertices[i], crossed[i]);
      }
    }

    support_plane(support_plane_idx).remove_vertex(front.second);
    support_plane(support_plane_idx).remove_vertex(back.second);

    // push also remaining vertex so that its events are recomputed
    // std::cout << "pushing new pv: " << str(pvertex) << std::endl;
    // std::cout << "pv direction: " << direction(pvertex) << std::endl;
    new_vertices.push_back(pvertex);
    crossed.push_back(iedge(pvertex));

    std::cout << "*** New vertices:";
    for (const PVertex& pv : new_vertices)
      std::cout << " " << str(pv);
    std::cout << std::endl;

    // for (const PVertex& pv : new_vertices)
    //   std::cout << point_3(pv) << std::endl;

    // if (has_iedge(prev) && !is_frozen(prev)) {
    // // if (iedge(prev) != iedge(pvertex)) {
    //   std::cout << "pushing new prev: " << str(prev) << std::endl;
    //   new_vertices.push_back (prev);
    // }

    // if (has_iedge(next) && !is_frozen(next)) {
    // // if (back_constrained) {
    //   std::cout << "pushing new next: " << str(next) << std::endl;
    //   new_vertices.push_back (next);
    // }

    return new_vertices;
  }

  void add_new_faces(
    const unsigned int k,
    const PVertex& pvertex,
    const std::vector<PVertex>& new_vertices,
    const PFace& pface) {

    CGAL_assertion(new_vertices.size() >= 2);
    for (std::size_t i = 0; i < new_vertices.size() - 1; ++i) {
      std::cout << "adding a new face" << std::endl;
      const PFace new_pface = add_pface(std::array<PVertex, 3>{new_vertices[i], new_vertices[i + 1], pvertex});
      this->k(new_pface) = k;
    }
  }

  void create_polyhedrons() {
    // for (std::size_t i = 0; i < number_of_support_planes(); ++i)
    //   std::cout << "num faces sp " << i << ": " << pfaces(i).size() << std::endl;

    // Check vertices.
    for (const auto vertex : m_intersection_graph.vertices()) {
      const auto nedges = m_intersection_graph.incident_edges(vertex);
      if (nedges.size() <= 2)
        std::cerr << "current num edges = " << nedges.size() << std::endl;
      CGAL_assertion_msg(nedges.size() > 2,
      "ERROR: VERTEX MUST HAVE AT LEAST 3 NEIGHBORS!");
    }

    // Check edges.
    std::vector<PFace> nfaces;
    for (const auto edge : m_intersection_graph.edges()) {
      incident_faces(edge, nfaces);
      if (nfaces.size() <= 1)
        std::cerr << "current num faces = " << nfaces.size() << std::endl;
      CGAL_assertion_msg(nfaces.size() > 1,
      "ERROR: EDGE MUST HAVE AT LEAST 2 NEIGHBORS!");
    }

    // Check faces.
    create_volumes();
    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        const auto nvolumes = incident_volumes(pface);
        if (nvolumes.size() == 0 || nvolumes.size() > 2)
          std::cerr << "current num volumes = " << nvolumes.size() << std::endl;
        CGAL_assertion_msg(nvolumes.size() == 1 || nvolumes.size() == 2,
        "ERROR: FACE MUST HAVE 1 OR 2 NEIGHBORS!");
      }
    }
  }

  void create_volumes() {

    m_volumes.clear();
    std::map<PFace, std::pair<int, int> > map_volumes;
    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces)
        map_volumes[pface] = std::make_pair(-1, -1);
    }

    int volume_index = 0;
    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        const bool success = traverse_pface(pface, volume_index, map_volumes);
        if (success) ++volume_index;
      }
    }
    std::cout << "Found " << volume_index << " polyhedrons!" << std::endl;

    for (const auto& item : map_volumes) {
      const auto& pair = item.second;

      CGAL_assertion(pair.first != -1);
      if (m_volumes.size() <= pair.first)
        m_volumes.resize(pair.first + 1);
      m_volumes[pair.first].add_pface(item.first, pair.second);
      if (pair.second == -1) continue;

      CGAL_assertion(pair.second != -1);
      if (m_volumes.size() <= pair.second)
        m_volumes.resize(pair.second + 1);
      m_volumes[pair.second].add_pface(item.first, pair.first);
    }
    for (auto& volume : m_volumes)
      create_cell_pvertices(volume);
    std::cout << "Created " << m_volumes.size() << " polyhedrons!" << std::endl;

    for (const auto& volume : m_volumes) {
      std::cout <<
      " pvertices: " << volume.pvertices.size() <<
      " pfaces: " << volume.pfaces.size() << std::endl;
    }

    dump_polyhedrons(*this, "iter_1000");
  }

  const bool traverse_pface(
    const PFace& pface,
    const int volume_index,
    std::map<PFace, std::pair<int, int> >& map_volumes) {

    auto& pair = map_volumes.at(pface);
    if (pair.first != -1 && pair.second != -1) return false;
    CGAL_assertion(pair.second == -1);
    if (pair.first == volume_index) return false;
    CGAL_assertion(pair.first != volume_index);

    if (pair.first != -1) pair.second = volume_index;
    else pair.first = volume_index;

    const auto pedges = pedges_of_pface(pface);
    const std::size_t n = pedges.size();

    std::vector<PEdge> edges;
    edges.reserve(n);
    for (const auto pedge : pedges) {
      edges.push_back(pedge);
    }
    CGAL_assertion(edges.size() == n);

    std::vector<PFace> nfaces;
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t im = (i + n - 1) % n;
      const std::size_t ip = (i + 1) % n;

      const auto& edgem = edges[im];
      const auto& edgei = edges[i];
      const auto& edgep = edges[ip];

      CGAL_assertion(has_iedge(edgei));
      incident_faces(this->iedge(edgei), nfaces);
      for (const auto& nface : nfaces) {
        if (belongs_to_this_volume(pface, nface, edgem, edgep)) {
          traverse_pface(nface, volume_index, map_volumes);
        }
      }
    }
    return true;
  }

  const bool belongs_to_this_volume(
    const PFace& pface, const PFace& nface,
    const PEdge& edgem, const PEdge& edgep) const {

    if (pface.first == nface.first) return false;

    std::vector<PFace> facesm, facesp;
    CGAL_assertion(has_iedge(edgem));
    incident_faces(this->iedge(edgem), facesm);
    CGAL_assertion(has_iedge(edgep));
    incident_faces(this->iedge(edgep), facesp);

    const bool found_prev = check_neighbor_faces(pface, nface, facesm);
    const bool found_next = check_neighbor_faces(pface, nface, facesp);
    return (found_prev && found_next);
  }

  const bool check_neighbor_faces(
    const PFace& pface,
    const PFace& nface,
    const std::vector<PFace>& faces) const {

    CGAL_assertion(pface != nface);
    for (const auto& face : faces) {
      CGAL_assertion(face != nface);
      if (face.first == nface.first) continue;
      if (face == pface) return true;
      if (has_equal_edge(pface, face)) {
        return true;
      }
    }
    return false;
  }

  const bool has_equal_edge(const PFace& pface, const PFace& nface) const {

    CGAL_assertion(pface != nface);
    for (const auto pedge : pedges_of_pface(pface)) {
      for (const auto nedge : pedges_of_pface(nface)) {
        if (pedge == nedge) return true;
      }
    }
    return false;
  }

  void create_cell_pvertices(Volume_cell& cell) {
    cell.pvertices.clear();
    for (const auto& pface : cell.pfaces) {
      for (const auto pvertex : pvertices_of_pface(pface)) {
        cell.pvertices.insert(pvertex);
      }
    }
  }

  const std::vector<Volume_cell> incident_volumes(const PFace& query_pface) const {
    std::vector<Volume_cell> nvolumes;
    for (const auto& volume : m_volumes) {
      for (const auto& pface : volume.pfaces) {
        if (pface == query_pface) nvolumes.push_back(volume);
      }
    }
    return nvolumes;
  }

  void incident_faces(const IEdge& query_iedge, std::vector<PFace>& nfaces) const {

    nfaces.clear();
    for (const auto plane_idx : intersected_planes(query_iedge)) {
      for (const auto pedge : pedges(plane_idx)) {
        if (iedge(pedge) == query_iedge) {
          const auto& m = mesh(plane_idx);
          const auto he = m.halfedge(pedge.second);
          const auto op = m.opposite(he);
          const auto face1 = m.face(he);
          const auto face2 = m.face(op);
          if (face1 != Support_plane::Mesh::null_face()) {
            nfaces.push_back(PFace(plane_idx, face1));
          }
          if (face2 != Support_plane::Mesh::null_face()) {
            nfaces.push_back(PFace(plane_idx, face2));
          }
        }
      }
    }
  }

  void update_positions (FT time)
  {
    m_previous_time = m_current_time;
    m_current_time = time;
  }

  inline std::string str (const PVertex& pvertex) const
  { return "PVertex(" + std::to_string(pvertex.first) + ":v" + std::to_string(pvertex.second) + ")"; }
  inline std::string str (const PEdge& pedge) const
  { return "PEdge(" + std::to_string(pedge.first) + ":e" + std::to_string(pedge.second) + ")"; }
  inline std::string str (const PFace& pface) const
  { return "PFace(" + std::to_string(pface.first) + ":f" + std::to_string(pface.second) + ")"; }
  inline std::string str (const IVertex& ivertex) const
  { return "IVertex(" + std::to_string(ivertex) + ")"; }
  inline std::string str (const IEdge& iedge) const
  { std::ostringstream oss; oss << "IEdge" << iedge; return oss.str(); }

  inline std::string lstr (const PFace& pface) const
  {
    if (pface == null_pface())
      return "PFace(null)";
    std::string out = "PFace(" + std::to_string(pface.first) + ":f" + std::to_string(pface.second) + ")[";
    for (PVertex pvertex : pvertices_of_pface (pface))
      out += "v" + std::to_string(pvertex.second);
    out += "]";
    return out;
  }
  inline std::string lstr (const PEdge& pedge) const
  { return "PEdge(" + std::to_string(pedge.first) + ":e" + std::to_string(pedge.second)
      + ")[v" + std::to_string(source(pedge).second) + "->v" + std::to_string(target(pedge).second) + "]"; }

  template <typename PSimplex>
  const Support_plane& support_plane (const PSimplex& psimplex) const { return support_plane(psimplex.first); }
  const Support_plane& support_plane (KSR::size_t idx) const { return m_support_planes[idx]; }
  template <typename PSimplex>
  Support_plane& support_plane (const PSimplex& psimplex) { return support_plane(psimplex.first); }
  Support_plane& support_plane (KSR::size_t idx) { return m_support_planes[idx]; }

private:

  template <typename PSimplex>
  const Mesh& mesh (const PSimplex& psimplex) const { return mesh(psimplex.first); }
  const Mesh& mesh (KSR::size_t support_plane_idx) const { return support_plane(support_plane_idx).mesh(); }
  template <typename PSimplex>
  Mesh& mesh (const PSimplex& psimplex) { return mesh(psimplex.first); }
  Mesh& mesh (KSR::size_t support_plane_idx) { return support_plane(support_plane_idx).mesh(); }

  void compute_future_points_and_directions (const PVertex& pvertex, const IEdge& iedge,
                                             Point_2& future_point_a, Point_2& future_point_b,
                                             Vector_2& direction_a, Vector_2& direction_b) const
  {
    const PVertex prev(pvertex.first, support_plane(pvertex).prev(pvertex.second));
    const PVertex next(pvertex.first, support_plane(pvertex).next(pvertex.second));
    const auto& curr = pvertex;

    const Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    const Point_2 pinit = iedge_line.projection(point_2(pvertex));

    // std::cout << "iedge segment: " << segment_3(iedge) << std::endl;

    const auto prev_p = point_2(prev);
    const auto next_p = point_2(next);
    const auto curr_p = point_2(curr);

    // std::cout << "prev: " << point_3(prev) << std::endl;
    // std::cout << "next: " << point_3(next) << std::endl;
    // std::cout << "curr: " << point_3(curr) << std::endl;

    const Line_2 future_line_prev(
      point_2(prev, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));
    const Line_2 future_line_next(
      point_2(next, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));

    // std::cout << "future line prev: " <<
    // Segment_3(
    //   to_3d(pvertex.first, point_2(prev, m_current_time + FT(1))),
    //   to_3d(pvertex.first, point_2(curr, m_current_time + FT(1)))) << std::endl;
    // std::cout << "future line next: " <<
    // Segment_3(
    //   to_3d(pvertex.first, point_2(next, m_current_time + FT(1))),
    //   to_3d(pvertex.first, point_2(curr, m_current_time + FT(1)))) << std::endl;

    const Vector_2 current_vec_prev(prev_p, curr_p);
    const Vector_2 current_vec_next(next_p, curr_p);

    const auto source_p = point_2(pvertex.first, source(iedge));
    const auto target_p = point_2(pvertex.first, target(iedge));
    const Vector_2 iedge_vec(source_p, target_p);

    const FT tol = FT(1) / FT(100000);
    FT m1 = FT(100000), m2 = FT(100000), m3 = FT(100000);

    const FT prev_d = (curr_p.x() - prev_p.x());
    const FT next_d = (curr_p.x() - next_p.x());
    const FT edge_d = (target_p.x() - source_p.x());

    if (CGAL::abs(prev_d) > tol)
      m1 = (curr_p.y() - prev_p.y()) / prev_d;
    if (CGAL::abs(next_d) > tol)
      m2 = (curr_p.y() - next_p.y()) / next_d;
    if (CGAL::abs(edge_d) > tol)
      m3 = (target_p.y() - source_p.y()) / edge_d;

    // std::cout << "prev slope: "  << m1 << std::endl;
    // std::cout << "next slope: "  << m2 << std::endl;
    // std::cout << "iedge slope: " << m3 << std::endl;

    if (CGAL::abs(m1 - m3) < tol) {
      std::cout << "prev parallel lines" << std::endl;
      const FT prev_dot = current_vec_prev * iedge_vec;
      if (prev_dot < FT(0)) {
        std::cout << "prev moves backwards" << std::endl;
        future_point_a = target_p;
      } else {
        std::cout << "prev moves forwards" << std::endl;
        future_point_a = source_p;
      }
    } else {

      std::cout << "prev intersected lines" << std::endl;
      const bool a_found = KSR::intersection_2(future_line_prev, iedge_line, future_point_a);
      if (!a_found)
      {
        std::cerr << "Warning: a not found" << std::endl;
        future_point_b = pinit + (pinit - future_point_a);
      }
    }

    direction_a = Vector_2(pinit, future_point_a);
    future_point_a = pinit - m_current_time * direction_a;
    std::cout << "future point a: " << to_3d(pvertex.first, future_point_a + m_current_time * direction_a) << std::endl;
    std::cout << "dir a: " << direction_a << std::endl;

    if (CGAL::abs(m2 - m3) < tol) {
      std::cout << "next parallel lines" << std::endl;
      const FT next_dot = current_vec_next * iedge_vec;
      if (next_dot < FT(0)) {
        std::cout << "next moves backwards" << std::endl;
        future_point_b = target_p;
      } else {
        std::cout << "next moves forwards" << std::endl;
        future_point_b = source_p;
      }

    } else {

      std::cout << "next intersected lines" << std::endl;
      const bool b_found = KSR::intersection_2(future_line_next, iedge_line, future_point_b);
      if (!b_found)
      {
        std::cerr << "Warning: b not found" << std::endl;
        future_point_a = pinit + (pinit - future_point_b);
      }
    }

    direction_b = Vector_2(pinit, future_point_b);
    future_point_b = pinit - m_current_time * direction_b;
    std::cout << "future point b: " << to_3d(pvertex.first, future_point_b + m_current_time * direction_b) << std::endl;
    std::cout << "dir b: " << direction_b << std::endl;
  }

  bool compute_future_point_and_direction (const std::size_t idx,
                                           const PVertex& pvertex, const PVertex& next, // back prev
                                           const IEdge& iedge,
                                           Point_2& future_point, Vector_2& direction) const
  {
    bool is_parallel = false;
    // if (this->iedge(pvertex) != null_iedge()
    //     && line_idx(pvertex) == line_idx(iedge))
    // {
    //   std::cerr << "found limit" << std::endl;
    //   future_point = point_2(pvertex, FT(0));
    //   direction = this->direction(pvertex);
    //   return is_parallel;
    // }

    const Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    const Point_2 pinit = iedge_line.projection(point_2(pvertex));

    const auto& curr = pvertex;
    const auto next_p = point_2(next);
    const auto curr_p = point_2(curr);

    const Line_2 future_line_next(
      point_2(next, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));
    const Vector_2 current_vec_next(next_p, curr_p);

    const auto source_p = point_2(pvertex.first, source(iedge));
    const auto target_p = point_2(pvertex.first, target(iedge));
    const Vector_2 iedge_vec(source_p, target_p);

    const FT tol = FT(1) / FT(100000);
    FT m2 = FT(100000), m3 = FT(100000);

    const FT next_d = (curr_p.x() - next_p.x());
    const FT edge_d = (target_p.x() - source_p.x());

    if (CGAL::abs(next_d) > tol)
      m2 = (curr_p.y() - next_p.y()) / next_d;
    if (CGAL::abs(edge_d) > tol)
      m3 = (target_p.y() - source_p.y()) / edge_d;

    // std::cout << "m2: " << m2 << std::endl;
    // std::cout << "m3: " << m3 << std::endl;

    if (CGAL::abs(m2 - m3) < tol) {
      std::cout << "back/front parallel lines" << std::endl;

      is_parallel = true;
      const FT next_dot = current_vec_next * iedge_vec;
      if (next_dot < FT(0)) {
        std::cout << "back/front moves backwards" << std::endl;
        future_point = target_p;
        // std::cout << point_3(target(iedge)) << std::endl;
      } else {
        std::cout << "back/front moves forwards" << std::endl;
        future_point = source_p;
        // std::cout << point_3(source(iedge)) << std::endl;
      }

    } else {
      std::cout << "back/front intersected lines" << std::endl;
      future_point = KSR::intersection_2<Point_2>(future_line_next, iedge_line);
    }

    direction = Vector_2(pinit, future_point);
    future_point = pinit - m_current_time * direction;
    return is_parallel;
  }

  bool compute_future_point_and_direction (const PVertex& pvertex,
                                           const PVertex& prev, const PVertex& next,
                                           const IEdge& iedge,
                                           Point_2& future_point, Vector_2& direction) const
  {
    const Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    const auto pv_point = point_2(pvertex);
    const Point_2 pinit = iedge_line.projection(pv_point);

    const auto& curr = prev;
    const auto next_p = point_2(next);
    const auto curr_p = point_2(curr);

    const Line_2 future_line_next(
      point_2(next, m_current_time + FT(1)),
      point_2(curr, m_current_time + FT(1)));

    const auto source_p = point_2(pvertex.first, source(iedge));
    const auto target_p = point_2(pvertex.first, target(iedge));
    const Vector_2 iedge_vec(source_p, target_p);

    const FT tol = FT(1) / FT(100000);
    FT m2 = FT(100000), m3 = FT(100000);

    const FT next_d = (curr_p.x() - next_p.x());
    const FT edge_d = (target_p.x() - source_p.x());

    if (CGAL::abs(next_d) > tol)
      m2 = (curr_p.y() - next_p.y()) / next_d;
    if (CGAL::abs(edge_d) > tol)
      m3 = (target_p.y() - source_p.y()) / edge_d;

    // std::cout << "m2: " << m2 << std::endl;
    // std::cout << "m3: " << m3 << std::endl;
    // std::cout << "mm: " << m2 - m3 << std::endl;

    bool is_parallel = false;
    if (CGAL::abs(m2 - m3) < tol) {
      std::cout << "open parallel lines" << std::endl;

      is_parallel = true;
      if (source_p == pv_point)
        future_point = target_p;
      else
        future_point = source_p;

    } else {
      std::cout << "open intersected lines" << std::endl;
      future_point = KSR::intersection_2<Point_2>(future_line_next, iedge_line);
    }

    direction = Vector_2(pinit, future_point);
    future_point = pinit - m_current_time * direction;
    return is_parallel;
  }

  bool is_intersecting_iedge(
    const FT min_time, const FT max_time,
    const PVertex& pvertex, const IEdge& iedge) {

    const FT time_step = (max_time - min_time) / FT(100);
    const FT time_1 = m_current_time - time_step;
    const FT time_2 = m_current_time + time_step;
    CGAL_assertion(time_1 != time_2);

    const Segment_2 pv_seg(
      point_2(pvertex, time_1), point_2(pvertex, time_2));
    const auto pv_bbox = pv_seg.bbox();

    const auto iedge_seg  = segment_2(pvertex.first, iedge);
    const auto iedge_bbox = iedge_seg.bbox();

    if (has_iedge(pvertex)) {
      std::cout << "constrained pvertex case" << std::endl;
      return false;
    }

    if (!is_active(pvertex)) {
      std::cout << "pvertex no active case" << std::endl;
      return false;
    }

    if (!is_active(iedge)) {
      std::cout << "iedge no active case" << std::endl;
      return false;
    }

    if (!CGAL::do_overlap(pv_bbox, iedge_bbox)) {
      std::cout << "no overlap case" << std::endl;
      return false;
    }

    Point_2 point;
    if (!KSR::intersection_2(pv_seg, iedge_seg, point)) {
      std::cout << "no intersection case" << std::endl;
      return false;
    }

    std::cout << "found intersection" << std::endl;
    return true;
  }
};

}} // namespace CGAL::KSR_3

#endif // CGAL_KSR_3_DATA_STRUCTURE_H
