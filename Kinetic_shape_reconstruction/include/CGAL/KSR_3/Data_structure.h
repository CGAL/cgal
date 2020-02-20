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

  typedef typename Intersection_graph::Vertex_descriptor IVertex;
  typedef typename Intersection_graph::Edge_descriptor IEdge;
  typedef typename Intersection_graph::Vertices IVertices;
  typedef typename Intersection_graph::Edges IEdges;
  typedef typename Intersection_graph::Incident_edges Incident_iedges;

private:

  // Main data structure
  Support_planes m_support_planes;
  Intersection_graph m_intersection_graph;

  // Helping data structures
  std::map<Point_3, KSR::size_t> m_meta_map;
  
  FT m_current_time;
  
public:

  Data_structure()
    : m_current_time(0)
  { }

  void print() const
  {
    // TODO
  }

  void init (std::size_t number_of_polygons)
  {
    m_support_planes.reserve (number_of_polygons + 6);
  }
  
  const FT& current_time() const { return m_current_time; }

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

      for (const PVertex& pvertex : pvertices_of_pface (pface))
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
      for (const IEdge& edge : m_intersection_graph.edges())
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
        IEdge edge_0, edge_1;
        std::tie (edge_0, edge_1)
          = m_intersection_graph.split_edge (intersections[i].first, vertices[i]);
        for (const IEdge& edge : { edge_0, edge_1 })
          for (KSR::size_t sp_idx : m_intersection_graph.intersected_planes(edge))
            support_plane(sp_idx).iedges().insert (edge);

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
    
    Point_2 centroid = CGAL::centroid (points.begin(), points.end());

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

  std::pair<PVertex, PVertex> border_prev_and_next (const PVertex& pvertex) const
  {
    Halfedge_index hi = mesh(pvertex).halfedge(pvertex.second);
    if (mesh(pvertex).face(hi) != Face_index())
      hi = mesh(pvertex).prev (mesh(pvertex).opposite(hi));
    
    CGAL_assertion (mesh(pvertex).face(hi) == Face_index());
    return std::make_pair (PVertex (pvertex.first, mesh(pvertex).source (hi)),
                           PVertex (pvertex.first, mesh(pvertex).target (mesh(pvertex).next(hi))));
  }

  PVertex add_pvertex (KSR::size_t support_plane_idx, const Point_2& point)
  {
    return PVertex (support_plane_idx, mesh(support_plane_idx).add_vertex(point));
  }

  template <typename VertexRange>
  PFace add_pface (const VertexRange& pvertices)
  {
    return PFace (pvertices.front().first,
                  mesh(pvertices.front()).add_face
                  (CGAL::make_range
                   (boost::make_transform_iterator
                    (pvertices.begin(),
                     CGAL::Property_map_to_unary_function
                     <CGAL::Second_of_pair_property_map<PVertex> >()),
                    boost::make_transform_iterator
                    (pvertices.end(),
                     CGAL::Property_map_to_unary_function
                     <CGAL::Second_of_pair_property_map<PVertex> >()))));
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
    for (const IEdge& incident_iedge : incident_iedges (ivertex))
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
    
    std::deque<PVertex> vertices;
    vertices.push_back (pvertex);
    
    std::queue<Queue_element> todo;
    PVertex prev, next;
    std::tie (prev, next) = border_prev_and_next (pvertex);
    todo.push (Queue_element (pvertex, prev, true, false));
    todo.push (Queue_element (pvertex, next, false, false));

    while (!todo.empty())
    {
      PVertex previous = todo.front().previous;
      PVertex current = todo.front().pvertex;
      bool front = todo.front().front;
      bool previous_was_free = todo.front().previous_was_free;
      todo.pop();

      IEdge iedge = this->iedge (current);
      bool is_free = (iedge == null_iedge());

      if (!is_free && source(iedge) != ivertex && target(iedge) != ivertex)
        is_free = true;

      if (!is_free)
      {
        IVertex other = source(iedge);
        if (other == ivertex)
          other = target(iedge);
        else
          CGAL_assertion (target(iedge) == ivertex);

        // Filter backwards vertex
        if (direction (current) * Vector_2 (point_2 (current.first, other),
                                            point_2 (current.first, ivertex))
            < 0)
        {
          std::cerr << str(current) << " is backwards" << std::endl;
          is_free = true;
        }
      }
      
      if (previous_was_free && is_free)
      {
        std::cerr << str(current) << " has no iedge, stopping there" << std::endl;
        continue;
      }
      
      if (is_free)
      {
        std::cerr << str(current) << " has no iedge" << std::endl;

      }
      else
      {
        std::cerr << str(current) << " has iedge " << str(iedge)
                  << " from " << str(source(iedge)) << " to " << str(target(iedge)) << std::endl;

      }

      if (front)
        vertices.push_front (current);
      else
        vertices.push_back (current);
      
      std::tie (prev, next) = border_prev_and_next (current);

      if (prev == previous)
      {
        CGAL_assertion (next != previous);
        todo.push (Queue_element (current, next, front, is_free));
      }
      else
        todo.push (Queue_element (current, prev, front, is_free));
    }

    std::vector<PVertex> out;
    out.reserve (vertices.size());
    std::copy (vertices.begin(), vertices.end(),
               std::back_inserter (out));
    
    CGAL_KSR_CERR(3) << "*** Found vertices:";
    for (const PVertex& pv : out)
      CGAL_KSR_CERR(3) << " " << str(pv);
    CGAL_KSR_CERR(3) << std::endl;
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

  std::pair<bool, bool> collision_occured (const PVertex& pvertex, const IEdge& iedge) const
  {
    bool collision = false;
    
    for (KSR::size_t support_plane_idx : intersected_planes(iedge))
    {
      if (support_plane_idx < 6)
        return std::make_pair (true, true);
      
      for (const PEdge& pedge : pedges(support_plane_idx))
        if (this->iedge(pedge) == iedge)
        {
          Segment_2 iedge_segment = segment_2 (support_plane_idx, iedge);
          
          Vector_2 source_2_vertex (iedge_segment.source(), point_2(pvertex));

          FT prod = iedge_segment.to_vector() * source_2_vertex;
          if (prod < 0)
            continue;

          if (source_2_vertex.squared_length() <= iedge_segment.squared_length())
          {
            collision = true;
            break;
          }
        }
    }
    
    return std::make_pair (collision, false);
  }
  
  /*******************************
   * Operations on polygons
   *******************************/

  PVertex crop_polygon (const PVertex& pvertex, const IEdge& iedge)
  {
    CGAL_KSR_CERR(3) << "*** Cropping " << str(pvertex) << " along " << str(iedge) << std::endl;
    
    Point_2 future_point_a, future_point_b;
    Vector_2 direction_a, direction_b;

    compute_future_points_and_directions (pvertex, iedge,
                                          future_point_a, future_point_b,
                                          direction_a, direction_b);
    
    PEdge pedge (pvertex.first, support_plane(pvertex).split_vertex(pvertex.second));
    CGAL_assertion (source(pedge) == pvertex || target(pedge) == pvertex);

    PVertex other = opposite(pedge, pvertex);
    
    CGAL_KSR_CERR(3) << "*** New edge " << str(pedge) << " between " << str(pvertex)
                     << " and " << str(other) << std::endl;

    connect (pedge, iedge);
    connect (pvertex, iedge);
    connect (other, iedge);

    support_plane(pvertex).set_point (pvertex.second, future_point_a);
    support_plane(other).set_point (other.second, future_point_b);
    
    direction(pvertex) = direction_a;
    direction(other) = direction_b;

    return other;
  }

  std::array<PVertex, 3> propagate_polygon (const PVertex& pvertex, const IEdge& iedge)
  {
    CGAL_KSR_CERR(3) << "*** Propagating " << str(pvertex) << " along " << str(iedge) << std::endl;
    
    Point_2 original_point = point_2 (pvertex, 0);
    Vector_2 original_direction = direction(pvertex);

    PVertex other = crop_polygon (pvertex, iedge);

    PVertex propagated = add_pvertex (pvertex.first, original_point);
    direction(propagated) = original_direction;
    
    std::array<PVertex, 3> pvertices;

    pvertices[0] = pvertex;
    pvertices[1] = other;
    pvertices[2] = propagated;

    PFace pface = add_pface (pvertices);
    CGAL_assertion (pface.second != Face_index());
    
    CGAL_KSR_CERR(3) << "*** New face " << lstr(pface) << std::endl;

    return pvertices;
  }

  void crop_polygon (const PVertex& pv0, const PVertex& pv1, const IEdge& iedge)
  {
    CGAL_KSR_CERR(3) << "*** Cropping " << str(pv0) << "/" << str(pv1) << " along " << str(iedge) << std::endl;

    Line_2 iedge_line = segment_2(pv0.first, iedge).supporting_line();
    
    for (const PVertex& pvertex : { pv0, pv1 })
    {
      Point_2 init = iedge_line.projection (point_2 (pvertex, m_current_time));
      Point_2 future = iedge_line.projection (point_2 (pvertex, m_current_time + 1));

      direction(pvertex) = (future - init);
      support_plane(pvertex).set_point (pvertex.second, init - direction(pvertex) * m_current_time);

      connect (pvertex, iedge);
    }

    PEdge pedge (pv0.first, support_plane(pv0).edge (pv0.second, pv1.second));
    connect (pedge, iedge);
  }

  std::pair<PVertex, PVertex> propagate_polygon
  (const PVertex& pvertex, const PVertex& pother, const IEdge& iedge)
  {
    CGAL_KSR_CERR(3) << "*** Propagating " << str(pvertex) << "/" << str(pother) << " along " << str(iedge) << std::endl;

    CGAL_assertion_msg (false, "TODO: propagate polygon via edge");

    return std::make_pair (null_pvertex(), null_pvertex());
  }

  bool transfer_vertex (const PVertex& pvertex, const PVertex& pother)
  {
    CGAL_KSR_CERR(3) << "*** Transfering " << str(pother) << " through " << str(pvertex) << std::endl;
    
    // If pvertex is adjacent to one or two
    PFace source_face, target_face;
    std::tie (source_face, target_face) = pfaces_of_pvertex (pvertex);

    PFace common_pface = pface_of_pvertex (pother);

    if (common_pface == target_face)
      std::swap (source_face, target_face);
    CGAL_assertion (common_pface == source_face);

    CGAL_KSR_CERR(3) << "*** Initial faces: " << lstr(source_face)
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
//      std::cerr << "Disconnect " << str(pvertex) << " from " << str(iedge) << std::endl;

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
//        std::cerr << "Shift target" << std::endl;
        CGAL::Euler::shift_target (hi, mesh(pedge));
      }
      else
      {
        CGAL_assertion (mesh(pedge).source(hi) == pvertex.second);
//        std::cerr << "Shift source" << std::endl;
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

//      std::cerr << "Connect " << str(pother) << " to " << str(iedge) << std::endl;
      connect (pother, iedge);
    }
    
    CGAL_KSR_CERR(3) << "*** New faces: " << lstr(source_face)
                     << " and " << lstr(target_face) << std::endl;

    return (target_face != null_pface());
  }
  
  void merge_pvertices (const PVertex& pvertex, const PVertex& pother)
  {
    CGAL_KSR_CERR(3) << "*** Merging " << str(pvertex) << " with " << str(pother) << std::endl;
    
    Halfedge_index hi = mesh(pvertex).halfedge(pother.second, pvertex.second);
    disconnect_ivertex (pother);
    CGAL::Euler::join_vertex(hi, mesh(pvertex));
  }

  std::vector<PVertex> merge_pvertices_on_ivertex (std::vector<PVertex>& pvertices,
                                                   const IVertex& ivertex)
  {
    KSR::size_t support_plane_idx = pvertices.front().first;

    PVertex prev = pvertices.front();
    PVertex next = pvertices.back();
    
    // Copy front/back to remember position/direction
    PVertex front (support_plane_idx, support_plane(support_plane_idx).duplicate_vertex(pvertices[1].second));
    PVertex back (support_plane_idx,support_plane(support_plane_idx).duplicate_vertex(pvertices[pvertices.size() - 2].second));

    auto pvertex_to_point =
      [&](const PVertex& a) -> Point_2
      {
        return point_2(a);
      };

    PFace fprev = pface_of_pvertex(prev);
    Point_2 pprev = CGAL::centroid
      (boost::make_transform_iterator (pvertices_of_pface(fprev).begin(), pvertex_to_point),
       boost::make_transform_iterator (pvertices_of_pface(fprev).end(), pvertex_to_point));
    PFace fnext = pface_of_pvertex(next);
    Point_2 pnext = CGAL::centroid
      (boost::make_transform_iterator (pvertices_of_pface(fnext).begin(), pvertex_to_point),
       boost::make_transform_iterator (pvertices_of_pface(fnext).end(), pvertex_to_point));

    if (CGAL::orientation (pprev, point_2 (support_plane_idx, ivertex), pnext) == CGAL::LEFT_TURN)
    {
      std::swap (prev, next);
      std::swap (front, back);
    }

    // Freeze vertices
    for (std::size_t i = 1; i < pvertices.size() - 1; ++ i)
    {
      PVertex& pvertex = pvertices[i];
      Point_2 p = point_2 (support_plane_idx, ivertex);
      support_plane(pvertex).direction (pvertex.second) = CGAL::NULL_VECTOR;
      support_plane(pvertex).set_point (pvertex.second, p);
    }

    PVertex pvertex = pvertices[1];
    connect (pvertex, ivertex);
    
    CGAL_KSR_CERR(3) << "*** Frozen vertex: " << str(pvertex) << std::endl;
    CGAL_KSR_CERR(3) << "*** Removed vertices:";

    // Join vertices
    for (std::size_t i = 2; i < pvertices.size() - 1; ++ i)
    {
      CGAL_KSR_CERR(3) << " " << str(pvertices[i]);

      Halfedge_index hi = mesh(support_plane_idx).halfedge(pvertices[i].second, pvertex.second);
      disconnect_ivertex (pvertices[i]);
      CGAL::Euler::join_vertex(hi, mesh(support_plane_idx));
    }
    CGAL_KSR_CERR(3) << std::endl;


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

    std::cerr << "Prev = " << point_2 (prev) << " / " << direction (prev) << std::endl
              << "Front = " << point_2 (front) << " / " << direction (front) << std::endl
              << "Back = " << point_2 (back) << " / " << direction (back) << std::endl
              << "Next = " << point_2 (next) << " / " << direction (next) << std::endl;

    std::vector<PVertex> new_vertices;

    if (back_constrained && front_constrained) // Closing case
    {
      CGAL_assertion_msg (false, "TODO: closing");
    }
    else if (back_constrained) // Border case
    {
      CGAL_KSR_CERR(3) << "*** Back border case" << std::endl;
      KSR::size_t other_side_limit = line_idx(next);

      Direction_2 dir (point_2(prev) - point_2 (pvertex));

      std::reverse (iedges.begin(), iedges.end());
      
      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++ i)
      {
        if (dir.counterclockwise_in_between(iedges[(i+1)%iedges.size()].second,
                                            iedges[i].second))
        {
          first_idx = (i+1)%iedges.size();
          break;
        }
      }

      std::ofstream ("first.polylines.txt")
        << "2 " << segment_3 (iedges[first_idx].first) << std::endl;

      CGAL_assertion (first_idx != KSR::no_element());

      std::vector<IEdge> crossed;

      KSR::size_t iedge_idx = first_idx;
      while (true)
      {
        const IEdge& iedge = iedges[iedge_idx].first;

        bool collision, bbox_reached;
        std::tie (collision, bbox_reached) = collision_occured (pvertex, iedge);
        bool limit_reached = (line_idx(iedge) == other_side_limit);

        crossed.push_back (iedge);

        std::ofstream ("next.polylines.txt")
          << "2 " << segment_3 (iedge) << std::endl;
        if (limit_reached || bbox_reached)
          break;
        
        iedge_idx = (iedge_idx + 1) % iedges.size();
      }

      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;

      std::vector<Point_2> future_points (crossed.size());
      std::vector<Vector_2> future_directions (crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++ i)
        compute_future_point_and_direction (back, prev, crossed[i], future_points[i], future_directions[i]);

      PVertex previous = null_pvertex();
      
      for (std::size_t i = 0; i < crossed.size(); ++ i)
      {
        if (i == 0) // crop
        {
          PVertex cropped (support_plane_idx, support_plane(pvertex).split_edge (pvertex.second, prev.second));

          PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, cropped.second));

          new_vertices.push_back (cropped);
    
          connect (pedge, crossed[i]);
          connect (cropped, crossed[i]);

          support_plane(cropped).set_point (cropped.second, future_points[i]);
          direction(cropped) = future_directions[i];
          previous = cropped;
          std::cerr << point_2 (cropped) << " -> " << direction(cropped) << std::endl;
        }
        else // create triangle face
        {
          PVertex propagated = add_pvertex (pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back (propagated);
          
          add_pface (std::array<PVertex, 3>{ pvertex, previous, propagated });
          previous = propagated;
          
          PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, propagated.second));
          connect (pedge, crossed[i]);
          connect (propagated, crossed[i]);
        }
      }
    }
    else if (front_constrained) // Border case
    {
      CGAL_KSR_CERR(3) << "*** Front border case" << std::endl;
      
      KSR::size_t other_side_limit = line_idx(prev);

      Direction_2 dir (point_2(next) - point_2 (pvertex));

      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++ i)
      {
        if (dir.counterclockwise_in_between(iedges[i].second,
                                            iedges[(i+1)%iedges.size()].second))
                                            
        {
          first_idx = (i+1)%iedges.size();
          break;
        }
      }

      std::ofstream ("first.polylines.txt")
        << "2 " << segment_3 (iedges[first_idx].first) << std::endl;

      CGAL_assertion (first_idx != KSR::no_element());

      std::vector<IEdge> crossed;

      KSR::size_t iedge_idx = first_idx;
      while (true)
      {
        const IEdge& iedge = iedges[iedge_idx].first;

        bool collision, bbox_reached;
        std::tie (collision, bbox_reached) = collision_occured (pvertex, iedge);
        bool limit_reached = (line_idx(iedge) == other_side_limit);

        std::ofstream ("next.polylines.txt")
          << "2 " << segment_3 (iedge) << std::endl;
        crossed.push_back (iedge);

        if (limit_reached || bbox_reached)
          break;
        
        iedge_idx = (iedge_idx + 1) % iedges.size();
      }

      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;

      std::vector<Point_2> future_points (crossed.size());
      std::vector<Vector_2> future_directions (crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++ i)
        compute_future_point_and_direction (front, next, crossed[i], future_points[i], future_directions[i]);

      PVertex previous = null_pvertex();
      
      for (std::size_t i = 0; i < crossed.size(); ++ i)
      {
        if (i == 0) // crop
        {
          PVertex cropped (support_plane_idx, support_plane(pvertex).split_edge (pvertex.second, next.second));

          PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, cropped.second));

          CGAL_assertion (cropped != pvertex);

          new_vertices.push_back (cropped);
    
          connect (pedge, crossed[i]);
          connect (cropped, crossed[i]);

          support_plane(cropped).set_point (cropped.second, future_points[i]);
          direction(cropped) = future_directions[i];
          previous = cropped;
          std::cerr << point_2 (cropped) << " -> " << direction(cropped) << std::endl;
        }
        else // create triangle face
        {
          PVertex propagated = add_pvertex (pvertex.first, future_points[i]);
          direction(propagated) = future_directions[i];
          new_vertices.push_back (propagated);
          
          add_pface (std::array<PVertex, 3>{ pvertex, previous, propagated });
          previous = propagated;
          
          PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, propagated.second));
          connect (pedge, crossed[i]);
          connect (propagated, crossed[i]);
        }
      }
    }
    else // Open case
    {
      CGAL_KSR_CERR(3) << "*** Open case" << std::endl;

      Direction_2 dir_next (point_2(next) - point_2 (pvertex));
      Direction_2 dir_prev (point_2(prev) - point_2 (pvertex));
      KSR::size_t first_idx = KSR::no_element();
      for (std::size_t i = 0; i < iedges.size(); ++ i)
      {
        if (dir_next.counterclockwise_in_between(iedges[i].second,
                                                 iedges[(i+1)%iedges.size()].second))
                                            
        {
          first_idx = (i+1)%iedges.size();
          break;
        }
      }

      std::vector<IEdge> crossed;

      KSR::size_t iedge_idx = first_idx;
      while (true)
      {
        const IEdge& iedge = iedges[iedge_idx].first;
        const Direction_2& dir = iedges[iedge_idx].second;

        if (!dir.counterclockwise_in_between (dir_next, dir_prev))
          break;

        crossed.push_back (iedge);

        iedge_idx = (iedge_idx + 1) % iedges.size();
      }

      std::cerr << "IEdges crossed = " << crossed.size() << std::endl;

      std::vector<Point_2> future_points (crossed.size());
      std::vector<Vector_2> future_directions (crossed.size());
      for (std::size_t i = 0; i < crossed.size(); ++ i)
        compute_future_point_and_direction (pvertex, prev, next, crossed[i], future_points[i], future_directions[i]);

      {
        PVertex cropped (support_plane_idx, support_plane(pvertex).split_edge (pvertex.second, next.second));

        PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, cropped.second));

        new_vertices.push_back (cropped);
    
        connect (pedge, crossed.front());
        connect (cropped, crossed.front());

        support_plane(cropped).set_point (cropped.second, future_points.front());
        direction(cropped) = future_directions.front();
      }
      
      for (std::size_t i = 1; i < crossed.size() - 1; ++ i)
      {
        PVertex propagated = add_pvertex (pvertex.first, future_points[i]);
        direction(propagated) = future_directions[i];
        connect (propagated, crossed[i]);
        new_vertices.push_back (propagated);
      }

      {
        PVertex cropped (support_plane_idx, support_plane(pvertex).split_edge (pvertex.second, prev.second));

        PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, cropped.second));

        new_vertices.push_back (cropped);
    
        connect (pedge, crossed.back());
        connect (cropped, crossed.back());

        support_plane(cropped).set_point (cropped.second, future_points.back());
        direction(cropped) = future_directions.back();
      }
      std::cerr << new_vertices.size() << " new vertice(s)" << std::endl;
      for (std::size_t i = 0; i < new_vertices.size() - 1; ++ i)
        add_pface (std::array<PVertex, 3>{ new_vertices[i], new_vertices[i+1], pvertex });
      
      for (std::size_t i = 1; i < crossed.size() - 1; ++ i)
      {
        PEdge pedge (support_plane_idx, support_plane(pvertex).edge (pvertex.second, new_vertices[i].second));
        connect (pedge, crossed[i]);
        connect (new_vertices[i], crossed[i]);
      }
    }

    support_plane(support_plane_idx).remove_vertex(front.second);
    support_plane(support_plane_idx).remove_vertex(back.second);

    CGAL_KSR_CERR(3) << "*** New vertices:";
    for (const PVertex& pv : new_vertices)
      CGAL_KSR_CERR(3) << " " << str(pv);
    CGAL_KSR_CERR(3) << std::endl;

    // push also remaining vertex so that its events are recomputed
    new_vertices.push_back (pvertex);
    
    return new_vertices;
  }

  
  void update_positions (FT time)
  {
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
private:

  template <typename PSimplex>
  const Support_plane& support_plane (const PSimplex& psimplex) const { return support_plane(psimplex.first); }
  const Support_plane& support_plane (KSR::size_t idx) const { return m_support_planes[idx]; }
  template <typename PSimplex>
  Support_plane& support_plane (const PSimplex& psimplex) { return support_plane(psimplex.first); }
  Support_plane& support_plane (KSR::size_t idx) { return m_support_planes[idx]; }
  
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
    PVertex prev (pvertex.first, support_plane(pvertex).prev(pvertex.second));
    PVertex next (pvertex.first, support_plane(pvertex).next(pvertex.second));
    
    Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    Point_2 pinit = iedge_line.projection(point_2 (pvertex, m_current_time));
    
    Line_2 future_line_prev (point_2 (prev, m_current_time + 1),
                             point_2 (pvertex, m_current_time + 1));
    Line_2 future_line_next (point_2 (next, m_current_time + 1),
                             point_2 (pvertex, m_current_time + 1));

    bool a_found = KSR::intersection_2 (future_line_prev, iedge_line, future_point_a);
    bool b_found = KSR::intersection_2 (future_line_next, iedge_line, future_point_b);

    if (!a_found)
    {
      std::cerr << "Warning: a not found" << std::endl;
      CGAL_assertion (b_found);
      future_point_b = pinit + (pinit - future_point_a);
    }
    if (!b_found)
    {
      std::cerr << "Warning: b not found" << std::endl;
      CGAL_assertion (a_found);
      future_point_a = pinit + (pinit - future_point_b);
    }
        
    direction_a = Vector_2 (pinit, future_point_a);
    direction_b = Vector_2 (pinit, future_point_b);
    future_point_a = pinit - m_current_time * direction_a;
    future_point_b = pinit - m_current_time * direction_b;
  }

  void compute_future_point_and_direction (const PVertex& pvertex, const PVertex& next,
                                           const IEdge& iedge,
                                           Point_2& future_point, Vector_2& direction) const
  {
    if (this->iedge(pvertex) != null_iedge()
        && line_idx(pvertex) == line_idx(iedge))
    {
      std::cerr << "Found limit" << std::endl;
      future_point = point_2 (pvertex, 0);
      direction = this->direction (pvertex);
      return;
    }
    Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    Point_2 pinit = iedge_line.projection(point_2 (pvertex, m_current_time));
    
    Line_2 future_line_next (point_2 (next, m_current_time + 1),
                             point_2 (pvertex, m_current_time + 1));

    future_point = KSR::intersection_2<Point_2> (future_line_next, iedge_line);
    direction = Vector_2 (pinit, future_point);
    future_point = pinit - m_current_time * direction;
  }

  void compute_future_point_and_direction (const PVertex& pvertex,
                                           const PVertex& prev, const PVertex& next,
                                           const IEdge& iedge,
                                           Point_2& future_point, Vector_2& direction) const
  {
    Line_2 iedge_line = segment_2(pvertex.first, iedge).supporting_line();
    Point_2 pinit = iedge_line.projection(point_2 (pvertex, m_current_time));
    
    Line_2 future_line (point_2 (next, m_current_time + 1),
                        point_2 (prev, m_current_time + 1));

    future_point = KSR::intersection_2<Point_2> (future_line, iedge_line);
    direction = Vector_2 (pinit, future_point);
    future_point = pinit - m_current_time * direction;
  }

  
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
