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

#ifndef CGAL_KSR_3_INTERSECTION_GRAPH_H
#define CGAL_KSR_3_INTERSECTION_GRAPH_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

#include <boost/graph/adjacency_list.hpp>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Intersection_graph
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Plane_3 Plane_3;

  struct Vertex_property
  {
    Point_3 point;

    Vertex_property () { }
    Vertex_property (const Point_3& point) : point (point) { }
  };
  
  struct Edge_property
  {
    KSR::Idx_set planes;
  };

  typedef boost::adjacency_list <boost::setS,
                                 boost::vecS,
                                 boost::undirectedS,
                                 Vertex_property,
                                 Edge_property> Graph;

  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::vertex_iterator Vertex_iterator;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_iterator Edge_iterator;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator Incident_edge_iterator;
  typedef CGAL::Iterator_range<Vertex_iterator> Vertices;
  typedef CGAL::Iterator_range<Edge_iterator> Edges;
  typedef CGAL::Iterator_range<Incident_edge_iterator> Incident_edges;
  
private:

  Graph m_graph;
  std::map<Point_3, Vertex_descriptor> m_map_points;
  std::map<KSR::Idx_vector, Vertex_descriptor> m_map_vertices;
  
public:

  Intersection_graph() { }

  static Vertex_descriptor null_ivertex()
  { return boost::graph_traits<Graph>::null_vertex(); }
    
  static Edge_descriptor null_iedge()
  { return Edge_descriptor(null_ivertex(), null_ivertex(), nullptr); }
  
  std::pair<Vertex_descriptor, bool> add_vertex (const Point_3& point)
  {
    typename std::map<Point_3, Vertex_descriptor>::iterator iter;
    bool inserted;
    std::tie (iter, inserted) = m_map_points.insert (std::make_pair (point, Vertex_descriptor()));
    if (inserted)
    {
      iter->second = boost::add_vertex(m_graph);
      m_graph[iter->second].point = point;
    }

    return std::make_pair (iter->second, inserted);
  }

  std::pair<Vertex_descriptor, bool> add_vertex (const Point_3& point,
                                                 const KSR::Idx_vector& intersected_planes)
  {
    typename std::map<KSR::Idx_vector, Vertex_descriptor>::iterator iter;
    bool inserted;
    std::tie (iter, inserted) = m_map_vertices.insert (std::make_pair (intersected_planes, Vertex_descriptor()));
    if (inserted)
    {
      iter->second = boost::add_vertex(m_graph);
      m_graph[iter->second].point = point;
    }

    return std::make_pair (iter->second, inserted);
  }

  std::pair<Edge_descriptor, bool> add_edge (const Vertex_descriptor& source, const Vertex_descriptor& target,
                                             KSR::size_t support_plane_idx)
  {
    std::pair<Edge_descriptor, bool> out = boost::add_edge (source, target, m_graph);
    m_graph[out.first].planes.insert (support_plane_idx);
    return out;
  }

  template <typename IndexContainer>
  std::pair<Edge_descriptor, bool> add_edge (const Vertex_descriptor& source, const Vertex_descriptor& target,
                                             const IndexContainer& support_planes_idx)
  {
    std::pair<Edge_descriptor, bool> out = boost::add_edge (source, target, m_graph);
    for (KSR::size_t support_plane_idx : support_planes_idx)
      m_graph[out.first].planes.insert (support_plane_idx);
    return out;
  }
  
  std::pair<Edge_descriptor, bool> add_edge (const Point_3& source, const Point_3& target)
  {
    return add_edge (add_vertex (source).first, add_vertex (target).first);
  }

  std::pair<Edge_descriptor, Edge_descriptor>
  split_edge (const Edge_descriptor& edge, const Vertex_descriptor& vertex)
  {
    Vertex_descriptor source = boost::source (edge, m_graph);
    Vertex_descriptor target = boost::target (edge, m_graph);
    Edge_property prop = m_graph[edge];
    
    boost::remove_edge (edge, m_graph);

    bool inserted;
    Edge_descriptor sedge;
    std::tie (sedge, inserted) = boost::add_edge (source, vertex, m_graph);

    if (!inserted)
    {
      std::cerr << segment_3 (edge) << " " << point_3 (vertex) << std::endl;
    }
    CGAL_assertion (inserted);
    
    m_graph[sedge] = prop;
    
    Edge_descriptor tedge;
    std::tie (tedge, inserted) = boost::add_edge (vertex, target, m_graph);
    if (!inserted)
    {
      std::cerr << segment_3 (edge) << " " << point_3 (vertex) << std::endl;
    }
    CGAL_assertion (inserted);
    
    m_graph[tedge] = prop;

    return std::make_pair (sedge, tedge);
  }

  Vertices vertices() const { return CGAL::make_range(boost::vertices (m_graph)); }
  Edges edges() const { return CGAL::make_range(boost::edges (m_graph)); }

  Vertex_descriptor source (const Edge_descriptor& edge) const
  { return boost::source (edge, m_graph); }
  Vertex_descriptor target (const Edge_descriptor& edge) const
  { return boost::target (edge, m_graph); }

  Incident_edges incident_edges (const Vertex_descriptor& vertex) const
  { return CGAL::make_range (boost::out_edges(vertex, m_graph)); }

  const KSR::Idx_set& intersected_planes (const Edge_descriptor& edge) const
  { return m_graph[edge].planes; }
  KSR::Idx_set& intersected_planes (const Edge_descriptor& edge)
  { return m_graph[edge].planes; }

  const Point_3& point_3 (const Vertex_descriptor& vertex) const
  {
    return m_graph[vertex].point;
  }

  Segment_3 segment_3 (const Edge_descriptor& edge) const
  {
    return Segment_3 (m_graph[boost::source (edge, m_graph)].point,
                      m_graph[boost::target (edge, m_graph)].point);
  }
  
  Line_3 line_3 (const Edge_descriptor& edge) const
  {
    return Line_3 (m_graph[source (edge, m_graph)].point,
                   m_graph[target (edge, m_graph)].point);
  }
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_INTERSECTION_GRAPH_H
