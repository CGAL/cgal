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

#include <boost/function_output_iterator.hpp>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/verbosity.h>
#include <CGAL/KSR/debug.h>

#include <CGAL/KSR_3/Support_plane.h>
#include <CGAL/KSR_3/Intersection_line.h>
#include <CGAL/KSR_3/Meta_vertex.h>

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
  { return mesh(support_plane_idx).is_valid(); }

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

      Point_3 centroid;
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
      IEdge iedge
        = m_intersection_graph.add_edge (ivertices[i], ivertices[(i+1)%4], support_plane_idx).first;
      
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

  PVertex add_pvertex (KSR::size_t support_plane_idx, const Point_2& point)
  {
    return PVertex (support_plane_idx, mesh(support_plane_idx).add_vertex(point));
  }
  PFace add_pface (KSR::size_t support_plane_idx, const std::vector<PVertex>& pvertices)
  {
    return PFace (support_plane_idx,
                  mesh(support_plane_idx).add_face
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

  bool is_frozen (const PVertex& pvertex) const
  { return support_plane(pvertex).is_frozen (pvertex.second); }
  const Vector_2& direction (const PVertex& pvertex) const
  { return support_plane(pvertex).direction (pvertex.second); }
  Vector_2& direction (const PVertex& pvertex)
  { return support_plane(pvertex).direction (pvertex.second); }
  FT speed (const PVertex& pvertex)
  { return support_plane(pvertex).speed (pvertex.second); }

  /*******************************
   * ISimplices
   *******************************/

  static IVertex null_ivertex() { return Intersection_graph::null_ivertex(); }
  static IEdge null_iedge() { return Intersection_graph::null_iedge(); }
  
  IVertices ivertices() const { return m_intersection_graph.vertices(); }
  IEdges iedges() const { return m_intersection_graph.edges(); }

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

    for (KSR::size_t i = 0; i < vertices.size() - 1; ++ i)
    {
      IEdge iedge;
      bool inserted;
      std::tie (iedge, inserted)
        = m_intersection_graph.add_edge (vertices[i],
                                         vertices[i+1],
                                         support_planes_idx);
      CGAL_assertion (inserted);
      
      for (KSR::size_t support_plane_idx : support_planes_idx)
        support_plane(support_plane_idx).iedges().insert (iedge);
    }
  }

  IVertex source (const IEdge& edge) const
  { return m_intersection_graph.source (edge); }
  IVertex target (const IEdge& edge) const
  { return m_intersection_graph.target (edge); }

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

  



#if 0
  inline bool point_is_inside_bbox_section_of_intersection_line
  (const Point_3& point, KSR::size_t intersection_line_idx) const
  {
    Vector_3 ref (meta_vertex(intersection_line(intersection_line_idx).meta_vertices_idx()[0]).point(),
                  meta_vertex(intersection_line(intersection_line_idx).meta_vertices_idx()[1]).point());
    Vector_3 position (meta_vertex(intersection_line(intersection_line_idx).meta_vertices_idx()[0]).point(),
                       point);

    if (ref * position < 0)
      return false;

    return (position * position) < (ref * ref);
  }

  KSR::size_t add_intersection_line (KSR::size_t support_plane_idx,
                                     KSR::size_t meta_vertex_idx_0,
                                     KSR::size_t meta_vertex_idx_1)
  {
    KSR::size_t intersection_line_idx
      = add_intersection_line (Intersection_line (Line_3 (meta_vertex(meta_vertex_idx_0).point(),
                                                          meta_vertex(meta_vertex_idx_1).point())));
    
    KSR::size_t common_support_plane_idx = KSR::no_element();
    for (KSR::size_t support_plane_idx_0 : meta_vertex(meta_vertex_idx_0).support_planes_idx())
      for (KSR::size_t support_plane_idx_1 : meta_vertex(meta_vertex_idx_1).support_planes_idx())
        if (support_plane_idx_0 != support_plane_idx
            && support_plane_idx_0 == support_plane_idx_1)
        {
          common_support_plane_idx = support_plane_idx_0;
          break;
        }
    CGAL_assertion (common_support_plane_idx != KSR::no_element());
    CGAL_assertion (common_support_plane_idx < 6);

    intersection_line(intersection_line_idx).meta_vertices_idx().push_back (meta_vertex_idx_0);
    intersection_line(intersection_line_idx).meta_vertices_idx().push_back (meta_vertex_idx_1);

    intersection_line(intersection_line_idx).support_planes_idx().push_back (support_plane_idx);
    intersection_line(intersection_line_idx).support_planes_idx().push_back (common_support_plane_idx);
    support_plane(support_plane_idx).intersection_lines_idx().push_back (intersection_line_idx);
    support_plane(common_support_plane_idx).intersection_lines_idx().push_back (intersection_line_idx);

    KSR::Idx_vector polygons_idx = support_plane(common_support_plane_idx).polygons_idx();
    for (KSR::size_t polygon_idx : polygons_idx)
      if (do_intersect (polygon_idx, line_on_support_plane (intersection_line_idx, common_support_plane_idx)))
        cut_polygon (polygon_idx, intersection_line_idx);

    return intersection_line_idx;
  }

  // Add segment on full intersection line, using 2 extrem meta vertices
  KSR::size_t add_segment (KSR::size_t intersection_line_idx, KSR::size_t source_idx, KSR::size_t target_idx,
                           KSR::size_t other_source_idx = KSR::no_element(),
                           KSR::size_t other_target_idx = KSR::no_element())
  {
    KSR::size_t segment_idx = add_segment (Segment (intersection_line_idx, source_idx, target_idx,
                                                    other_source_idx, other_target_idx));
    intersection_line(intersection_line_idx).segments_idx().push_back (segment_idx);
    return segment_idx;
  }
  
  void partition (KSR::size_t polygon_idx, const Line_2& line,
                  KSR::Idx_vector& positive_side, KSR::Idx_vector& negative_side) const
  {
    std::vector<bool> has_on_positive_side;
    has_on_positive_side.reserve(polygon(polygon_idx).vertices_idx().size());
    
    for (KSR::size_t vertex_idx : polygon(polygon_idx).vertices_idx())
      has_on_positive_side.push_back (line.has_on_positive_side(vertex(vertex_idx).point(m_current_time)));

    KSR::size_t first_positive = KSR::no_element();
        
    for (std::size_t i = 0; i <= has_on_positive_side.size(); ++ i)
      if (!has_on_positive_side[i] && has_on_positive_side[(i+1) % has_on_positive_side.size()])
      {
        first_positive = KSR::size_t((i+1) % has_on_positive_side.size());
        break;
      }

    KSR::size_t current_position = first_positive;
    do
    {
      if (has_on_positive_side[std::size_t(current_position)])
        positive_side.push_back (polygon(polygon_idx).vertices_idx()[current_position]);
      else
        negative_side.push_back (polygon(polygon_idx).vertices_idx()[current_position]);
      
      current_position = (current_position + 1) % has_on_positive_side.size();
    }
    while (current_position != first_positive);

    CGAL_assertion (!positive_side.empty() && !negative_side.empty());
    CGAL_assertion (positive_side.size() + negative_side.size() >= 3);
  }

  std::tuple<Point_2, Vector_2, Point_2, Vector_2>
  compute_constrained_points_along_line (const Line_2& line_2,
                                         const KSR::Idx_vector& positive_side,
                                         const KSR::Idx_vector& negative_side) const
  {
    Line_2 line_0 (vertex(positive_side.back()).point(m_current_time),
                   vertex(negative_side.front()).point(m_current_time));
    Line_2 line_1 (vertex(negative_side.back()).point(m_current_time),
                   vertex(positive_side.front()).point(m_current_time));
    
    Point_2 inter_0 = KSR::intersection_2<Point_2> (line_0, line_2);
    Point_2 inter_1 = KSR::intersection_2<Point_2> (line_1, line_2);

    // Compute speeds
    Line_2 future_line_0 (vertex(positive_side.back()).point(m_current_time + 1),
                          vertex(negative_side.front()).point(m_current_time + 1));
    Line_2 future_line_1 (vertex(negative_side.back()).point(m_current_time + 1),
                          vertex(positive_side.front()).point(m_current_time + 1));

    Point_2 future_inter_0 = KSR::intersection_2<Point_2> (future_line_0, line_2);
    Point_2 future_inter_1 = KSR::intersection_2<Point_2> (future_line_1, line_2);
    
    Vector_2 direction_0 (inter_0, future_inter_0);
    Vector_2 direction_1 (inter_1, future_inter_1);

    return std::make_tuple (inter_0, direction_0, inter_1, direction_1);
  }

  void cut_polygon (KSR::size_t polygon_idx, KSR::size_t intersection_line_idx)
  {
    CGAL_KSR_CERR(3) << "** Cutting " << polygon_str(polygon_idx) << std::endl;

    Line_2 line_2 = line_on_support_plane (intersection_line_idx, polygon(polygon_idx).support_plane_idx());
    // Partition
    KSR::Idx_vector positive_side, negative_side;
    partition (polygon_idx, line_2, positive_side, negative_side);

    // Position + Direction of V1 + V2
    Point_2 point_0, point_1;
    Vector_2 direction_0, direction_1;
    std::tie (point_0, direction_0, point_1, direction_1)
      = compute_constrained_points_along_line (line_2, positive_side, negative_side);

    KSR::size_t new_polygon_idx = add_polygon (Polygon(polygon(polygon_idx).input_idx(), polygon(polygon_idx).support_plane_idx()));
    support_plane(polygon(polygon_idx).support_plane_idx()).polygons_idx().push_back (new_polygon_idx);

    for (KSR::size_t vertex_idx : positive_side)
      vertex(vertex_idx).polygon_idx() = polygon_idx;
    for (KSR::size_t vertex_idx : negative_side)
      vertex(vertex_idx).polygon_idx() = new_polygon_idx;

    KSR::size_t segment_idx = add_segment (intersection_line_idx,
                                           number_of_vertices(), number_of_vertices() + 1,
                                           number_of_vertices() + 3, number_of_vertices() + 2);
    
    positive_side.push_back (add_vertex (Vertex (point_0 - m_current_time * direction_0, polygon_idx)));
    m_vertices.back().segment_idx() = segment_idx;
    m_vertices.back().direction() = direction_0;
    
    positive_side.push_back (add_vertex (Vertex (point_1 - m_current_time * direction_1, polygon_idx)));
    m_vertices.back().segment_idx() = segment_idx;
    m_vertices.back().direction() = direction_1;

    negative_side.push_back (add_vertex (Vertex (point_1 - m_current_time * direction_1, new_polygon_idx)));
    m_vertices.back().segment_idx() = segment_idx;
    m_vertices.back().direction() = direction_1;
    
    negative_side.push_back (add_vertex (Vertex (point_0 - m_current_time * direction_0, new_polygon_idx)));
    m_vertices.back().segment_idx() = segment_idx;
    m_vertices.back().direction() = direction_0;

    if (direction_0 == CGAL::NULL_VECTOR)
    {
      KSR::size_t meta_vertex_idx = add_meta_vertex (support_plane_of_polygon(polygon_idx).to_3d (point_0),
                                                     polygon(polygon_idx).support_plane_idx());
      attach_vertex_to_meta_vertex (m_vertices.size() - 4, meta_vertex_idx);
      attach_vertex_to_meta_vertex (m_vertices.size() - 1, meta_vertex_idx);
    }
    
    if (direction_1 == CGAL::NULL_VECTOR)
    {
      KSR::size_t meta_vertex_idx = add_meta_vertex (support_plane_of_polygon(polygon_idx).to_3d (point_1),
                                                     polygon(polygon_idx).support_plane_idx());
      attach_vertex_to_meta_vertex (m_vertices.size() - 3, meta_vertex_idx);
      attach_vertex_to_meta_vertex (m_vertices.size() - 2, meta_vertex_idx);
    }
    
    polygon(polygon_idx).vertices_idx().swap (positive_side);
    polygon(new_polygon_idx).vertices_idx().swap (negative_side);
    
    CGAL_KSR_CERR(3) << "*** new polygons:";
    for (KSR::size_t i : { polygon_idx, new_polygon_idx })
      CGAL_KSR_CERR(3) << " " << polygon_str(i);
    CGAL_KSR_CERR(3) << std::endl;
  }

  KSR::size_t crop_polygon (KSR::size_t polygon_idx, KSR::size_t intersection_line_idx,
                            KSR::Idx_vector& vertices,
                            const Point_2& point_0, const Vector_2& direction_0,
                            const Point_2& point_1, const Vector_2& direction_1)
  {
    m_vertices[vertices.back()] = Vertex (point_1 - m_current_time * direction_1, polygon_idx);
    m_vertices[vertices.back()].segment_idx() = number_of_segments();
    m_vertices[vertices.back()].direction() = direction_1;
    
    vertices.push_back (add_vertex (Vertex (point_0 - m_current_time * direction_0, polygon_idx)));
    m_vertices.back().segment_idx() = number_of_segments();
    m_vertices.back().direction() = direction_0;

    KSR::size_t segment_idx = add_segment (intersection_line_idx,
                                           vertices[vertices.size() - 2],
                                           vertices[vertices.size() - 1]);

    polygon(polygon_idx).vertices_idx().swap (vertices);

    return segment_idx;
  }

  KSR::size_t propagate_polygon (KSR::size_t segment_idx,
                                 const Point_2& point,
                                 const Vector_2& direction)
  {
    KSR::size_t source_idx = segment(segment_idx).source_idx();
    KSR::size_t target_idx = segment(segment_idx).target_idx();
    KSR::size_t polygon_idx = vertex(source_idx).polygon_idx();
    
    KSR::size_t new_polygon_idx = add_polygon (Polygon(polygon(polygon_idx).input_idx(), polygon(polygon_idx).support_plane_idx()));
    support_plane(polygon(polygon_idx).support_plane_idx()).polygons_idx().push_back (new_polygon_idx);

    // Copy segment vertices
    segment(segment_idx).other_source_idx() = add_vertex (vertex(source_idx));
    segment(segment_idx).other_target_idx() = add_vertex (vertex(target_idx));

    vertex(segment(segment_idx).other_source_idx()).segment_idx() = segment_idx;
    vertex(segment(segment_idx).other_target_idx()).segment_idx() = segment_idx;

    vertex(segment(segment_idx).other_source_idx()).polygon_idx() = new_polygon_idx;
    vertex(segment(segment_idx).other_target_idx()).polygon_idx() = new_polygon_idx;

    KSR::size_t new_vertex_idx = add_vertex (Vertex (point - m_current_time * direction, new_polygon_idx));

    polygon(new_polygon_idx).vertices_idx().push_back (segment(segment_idx).other_target_idx());
    polygon(new_polygon_idx).vertices_idx().push_back (segment(segment_idx).other_source_idx());
    polygon(new_polygon_idx).vertices_idx().push_back (new_vertex_idx);

    return new_vertex_idx;
  }

  std::pair<KSR::size_t, KSR::size_t>
  transfer_vertex (KSR::size_t vertex_idx, KSR::size_t segment_border_idx,
                   const Point_2& point, const Vector_2& direction)
  {
    KSR::size_t segment_idx = vertex(segment_border_idx).segment_idx();
    KSR::size_t mirror_segment_border_idx = segment(segment_idx).mirror_vertex (segment_border_idx);
    KSR::size_t polygon_idx = vertex(vertex_idx).polygon_idx();

    if (mirror_segment_border_idx == KSR::no_element()) // Border reached
    {
      vertex(segment_border_idx) = Vertex (point - m_current_time * direction, polygon_idx);
      vertex(segment_border_idx).segment_idx() = segment_idx;
      vertex(segment_border_idx).direction() = CGAL::NULL_VECTOR;
      
      vertex(vertex_idx) = Vertex (point - m_current_time * direction, polygon_idx);
      vertex(vertex_idx).direction() = direction;
      
      return std::make_pair(vertex_idx, KSR::no_element());
    }
    else
    {
      KSR::size_t other_polygon_idx = vertex(mirror_segment_border_idx).polygon_idx();

      polygon(polygon_idx).vertices_idx().erase
        (std::find (polygon(polygon_idx).vertices_idx().begin(), polygon(polygon_idx).vertices_idx().end(), vertex_idx));

      vertex(vertex_idx).polygon_idx() = other_polygon_idx;

      vertex(segment_border_idx) = Vertex (point - m_current_time * direction, polygon_idx);
      vertex(segment_border_idx).segment_idx() = segment_idx;
      vertex(segment_border_idx).direction() = direction;
    
      vertex(mirror_segment_border_idx) = Vertex (point - m_current_time * direction, other_polygon_idx);
      vertex(mirror_segment_border_idx).segment_idx() = segment_idx;
      vertex(mirror_segment_border_idx).direction() = direction;

      for (KSR::size_t i = 0; i < polygon(other_polygon_idx).vertices_idx().size(); ++ i)
      {
        KSR::size_t vertex_idx_0 = polygon(other_polygon_idx).vertices_idx()[i];
        KSR::size_t vertex_idx_1 = polygon(other_polygon_idx).vertices_idx()[(i+1) % polygon(other_polygon_idx).vertices_idx().size()];

        if ((vertex_idx_0 == mirror_segment_border_idx && vertex(vertex_idx_1).segment_idx() != segment_idx)
            || (vertex_idx_1 == mirror_segment_border_idx && vertex(vertex_idx_0).segment_idx() != segment_idx))
        {
          polygon(other_polygon_idx).vertices_idx().insert
            (polygon(other_polygon_idx).vertices_idx().begin() + i + 1, vertex_idx);
          return std::make_pair (segment_border_idx, mirror_segment_border_idx);
        }
      }
    }
    
    // This should never be reached
    CGAL_assertion(false);
    return std::make_pair (KSR::no_element(), KSR::no_element());
  }

  void cut_segment_of_vertex (KSR::size_t vertex_idx)
  {
    KSR::size_t segment_idx = vertex(vertex_idx).segment_idx();
    
  }
#endif
  
  void update_positions (FT time)
  {
    m_current_time = time;
  }

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
  

  
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
