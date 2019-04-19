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
#include <CGAL/KSR_3/Support_plane.h>
#include <CGAL/KSR_3/Support_line.h>
#include <CGAL/KSR_3/Polygon.h>
#include <CGAL/KSR_3/Vertex.h>

#include <CGAL/KSR/Event.h>
#include <CGAL/KSR/Event_queue.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Data_structure
{
public:
  
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Line_3 Line_3;

  typedef KSR_3::Support_plane<Kernel> Support_plane;
  typedef KSR_3::Support_line<Kernel> Support_line;
  typedef KSR_3::Vertex<Kernel> Vertex;
  typedef KSR_3::Polygon<Kernel> Polygon;

  typedef std::vector<Support_plane> Support_planes;
  typedef std::vector<Support_line> Support_lines;
  typedef std::vector<Vertex> Vertices;
  typedef std::vector<Polygon> Polygons;

  typedef KSR::Event<Kernel> Event;
  typedef KSR::Event_queue<Kernel> Event_queue;

private:

  // Utilies / Property maps / Unary functions
  struct Point_3_to_2
  {
    typedef Point_3 argument_type;
    typedef Point_2 return_type;
    const Plane_3& plane;
    Point_3_to_2 (const Plane_3& plane) : plane (plane) { }
    Point_2 operator() (const Point_3& point) const { return plane.to_2d(point); }
  };

  struct Vertex_index_to_point_2
  {
    typedef KSR::size_t key_type;
    typedef Point_2 return_type;
    const Data_structure& ds;
    Vertex_index_to_point_2 (const Data_structure& ds) : ds (ds) { }
    const return_type& operator() (const key_type& k) const { return ds.vertices()[std::size_t(k)].point(); }
  };
  Vertex_index_to_point_2 vertex_index_to_point_2() const
  {
    return Vertex_index_to_point_2(*this);
  }

  // Main data structure
  Support_planes m_support_planes;
  Support_lines m_support_lines;
  Vertices m_vertices;
  Polygons m_polygons;

  Event_queue m_queue;

public:

  Data_structure() { }

  const Support_planes& support_planes() const { return m_support_planes; }
  const Support_lines& support_lines() const { return m_support_lines; }
  const Vertices& vertices() const { return m_vertices; }
  const Polygons& polygons() const { return m_polygons; }

  Support_planes& support_planes() { return m_support_planes; }
  Support_lines& support_lines() { return m_support_lines; }
  Vertices& vertices() { return m_vertices; }
  Polygons& polygons() { return m_polygons; }

  std::size_t number_of_vertices() const { return m_vertices.size(); }
  std::size_t number_of_polygons() const { return m_polygons.size(); }

  Point_3 point (std::size_t vertex_idx) const
  {
    const Vertex& vertex = m_vertices[vertex_idx];
    const Polygon& polygon = m_polygons[std::size_t(vertex.polygon())];
    const Support_plane& plane = m_support_planes[std::size_t(polygon.support_plane())];
    return plane.to_3d(vertex.point());
  }

  std::vector<std::size_t> polygon (std::size_t polygon_idx) const
  {
    std::vector<std::size_t> out;
    out.reserve (m_polygons.size());
    std::transform (m_polygons[polygon_idx].vertices().begin(), m_polygons[polygon_idx].vertices().end(),
                    std::back_inserter(out),
                    [&](const KSR::size_t& idx) -> std::size_t { return std::size_t(idx); });
    return out;
  }

  void add_bbox_as_polygons (const CGAL::Bbox_3& bbox)
  {
    FT xmed = (bbox.xmin() + bbox.xmax()) / 2.;
    FT ymed = (bbox.ymin() + bbox.ymax()) / 2.;
    FT zmed = (bbox.zmin() + bbox.zmax()) / 2.;
    FT dx = (bbox.xmax() - bbox.xmin()) / 2.;
    FT dy = (bbox.ymax() - bbox.ymin()) / 2.;
    FT dz = (bbox.zmax() - bbox.zmin()) / 2.;

    FT ratio = 1.1;
    FT xmin = xmed - ratio * dx;
    FT xmax = xmed + ratio * dx;
    FT ymin = ymed - ratio * dy;
    FT ymax = ymed + ratio * dy;
    FT zmin = zmed - ratio * dz;
    FT zmax = zmed + ratio * dz;
    
    std::array<Point_3, 8> bbox_points
      = { Point_3 (xmin, ymin, zmin),
          Point_3 (xmin, ymin, zmax),
          Point_3 (xmin, ymax, zmin),
          Point_3 (xmin, ymax, zmax),
          Point_3 (xmax, ymin, zmin),
          Point_3 (xmax, ymin, zmax),
          Point_3 (xmax, ymax, zmin),
          Point_3 (xmax, ymax, zmax) };

    std::array<Point_3, 4> facet_points = { bbox_points[0], bbox_points[1], bbox_points[3], bbox_points[2] };
    add_polygon (facet_points);
    
    facet_points = { bbox_points[4], bbox_points[5], bbox_points[7], bbox_points[6] };
    add_polygon (facet_points);
    
    facet_points = { bbox_points[0], bbox_points[1], bbox_points[5], bbox_points[4] };
    add_polygon (facet_points);

    facet_points = { bbox_points[2], bbox_points[3], bbox_points[7], bbox_points[6] };
    add_polygon (facet_points);
    
    facet_points = { bbox_points[1], bbox_points[5], bbox_points[7], bbox_points[3] };
    add_polygon (facet_points);
    
    facet_points = { bbox_points[0], bbox_points[4], bbox_points[6], bbox_points[2] };
    add_polygon (facet_points);
  }
  
  template <typename PolygonRange, typename PolygonMap>
  void add_polygons (const PolygonRange& polygons, PolygonMap polygon_map)
  {
    for (const typename PolygonRange::const_iterator::value_type& vt : polygons)
    {
      add_polygon (get (polygon_map, vt));
      initialize_vertices_directions(m_polygons.size() - 1);
    }
  }

  void compute_support_lines()
  {
    for (Support_plane& p : m_support_planes)
      p.support_lines().resize (m_support_planes.size(), KSR::no_element());

    for (std::size_t i = 0; i < m_support_planes.size() - 1; ++ i)
    {
      const Plane_3& pi = m_support_planes[i].plane();
      
      for (std::size_t j = i+1; j < m_support_planes.size(); ++ j)
      {
        const Plane_3& pj = m_support_planes[j].plane();

        Line_3 line_inter;
        if (!KSR::intersection_3 (pi, pj, line_inter))
          continue;
        
        // TODO check if line already exist and do not duplicate

        m_support_planes[i].support_lines()[j] = KSR::size_t(m_support_lines.size());
        m_support_planes[j].support_lines()[i] = KSR::size_t(m_support_lines.size());
        m_support_lines.push_back (Support_line(line_inter));
      }
    }

    // Compute intersections
    for (std::size_t p = 0; p < m_support_planes.size(); ++ p)
    {
      const std::vector<KSR::size_t>& support_lines = m_support_planes[p].support_lines();

      for (std::size_t i = 0; i < support_lines.size()-1; ++ i)
      {
        Line_2 li = m_support_planes[p].to_2d (m_support_lines[support_lines[i]].line());
        for (std::size_t j = i+1; j < support_lines.size(); ++ j)
        {
          Line_2 lj = m_support_planes[p].to_2d (m_support_lines[support_lines[j]].line());

          Point_2 point_inter;
          if (!KSR::intersection_2 (li, lj, point_inter))
            continue;

          m_vertices.push_back (Vertex(point_inter, KSR::no_element(), p));
        }
      }
    }

  }

  void make_polygons_intersection_free()
  {
    // TODO
  }

  void initialize_queue()
  {
    // Loop over vertices and schedule events
    for (std::size_t i = 0; i < m_vertices.size(); ++ i)
    {
      Vertex& vertex = m_vertices[i];
      if (vertex.is_frozen())
        continue;

      Polygon& polygon = m_polygons[vertex.polygon()];
      Support_plane& support_plane = m_support_planes[polygon.support_plane()];

      Ray_2 ray = vertex.ray();

      for (KSR::size_t sl : support_plane.support_lines())
      {
        if (sl == KSR::no_element())
          continue;
        Support_line& support_line = m_support_lines[sl];
        Line_2 line = support_plane.to_2d (support_line.line());

        Point_2 point;
        if (!KSR::intersection_2 (ray, line, point))
          continue;

        FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (vertex.point(), point));
        FT time = dist / vertex.speed();
          
        m_queue.push (Event (i, sl, time));
      }
    }

    m_queue.print();
  }

  void run()
  {
    FT latest_time = FT(0);

    std::size_t iterations = 0;
    while (!m_queue.empty())
    {
      Event ev = m_queue.pop();
      std::cerr << " * Applying " << ev << std::endl;

      FT ellapsed_time = ev.time() - latest_time;
      latest_time = ev.time();

      advance_time (ellapsed_time);

      Vertex& vertex = m_vertices[ev.vertex()];
      if (!vertex.is_constrained())
        split_vertex(ev.vertex(), ev.intersection_line());

      ++ iterations;
      if (iterations == 5)
        break;
    }
  }
  

private:

  template <typename PointRange>
  void add_polygon (const PointRange& points)
  {
    // Compute support plane
    Vector_3 normal = CGAL::NULL_VECTOR;

    //Newell's method
    for (std::size_t i = 0; i < points.size(); ++ i)
    {
      const Point_3& pa = points[i];
      const Point_3& pb = points[(i+1) % points.size()];
      FT x = normal.x() + (pa.y()-pb.y())*(pa.z()+pb.z());
      FT y = normal.y() + (pa.z()-pb.z())*(pa.x()+pb.x());
      FT z = normal.z() + (pa.x()-pb.x())*(pa.y()+pb.y());
      normal = Vector_3 (x,y,z);
    }
    CGAL_assertion_msg (normal != CGAL::NULL_VECTOR, "Polygon is flat");

    m_support_planes.push_back (Support_plane (CGAL::centroid (points.begin(), points.end()), normal));

    m_polygons.push_back (Polygon(m_support_planes.size() - 1));
    
    CGAL::convex_hull_2 (boost::make_transform_iterator
                         (points.begin(), Point_3_to_2(m_support_planes.back().plane())),
                         boost::make_transform_iterator
                         (points.end(), Point_3_to_2(m_support_planes.back().plane())),
                         boost::make_function_output_iterator
                         ([&](const Point_2& point) -> void
                          {
                            m_polygons.back().add_vertex(m_vertices.size());
                            m_vertices.push_back(Vertex(point, m_polygons.size() - 1, m_support_planes.size() - 1));
                          }));
  }

  void initialize_vertices_directions (std::size_t polygon_idx)
  {
    Polygon& polygon = m_polygons[polygon_idx];
    
    Point_2 centroid = CGAL::centroid(boost::make_transform_iterator
                                      (polygon.vertices().begin(),
                                       vertex_index_to_point_2()),
                                      boost::make_transform_iterator
                                      (polygon.vertices().end(),
                                       vertex_index_to_point_2()));

    for (KSR::size_t vidx : polygon.vertices())
    {
      Vector_2 d (centroid, m_vertices[vidx].point());
      m_vertices[vidx].direction() = KSR::normalize(d);
    }
  }

  void advance_time (FT time)
  {
    for (Vertex& v : m_vertices)
    {
      if (v.is_frozen())
        continue;

      v.point() = v.point() + time * v.direction();

    }
  }

  void split_vertex (std::size_t vertex_idx, std::size_t line_idx)
  {
    std::ofstream file ("split.xyz");
    file << point(vertex_idx) << std::endl;
    
    KSR::size_t polygon_idx = m_vertices[vertex_idx].polygon();
    Polygon& polygon = m_polygons[polygon_idx];
    
    KSR::size_t previous_vertex_idx, next_vertex_idx;
    std::tie (previous_vertex_idx, next_vertex_idx)
      = polygon.previous_and_next_vertex(vertex_idx);
    
    std::size_t new_vertex_idx = m_vertices.size();
    m_vertices.push_back (Vertex (m_vertices[vertex_idx].point(), polygon_idx));
    Vertex& v0 = m_vertices[vertex_idx];
    Vertex& v1 = m_vertices.back();

    std::cerr << "Polygon p" << polygon_idx << " before:";
    for (KSR::size_t v : polygon.vertices())
      std::cerr << " v" << v;
    std::cerr << std::endl;

    std::cerr << "Splitting v" << vertex_idx << " between v" << previous_vertex_idx
              << " and v" << next_vertex_idx << ": new vertex v" << new_vertex_idx << std::endl;
    
    typename std::vector<KSR::size_t>::iterator iter = polygon.vertices().begin();
    while (*iter != vertex_idx)
      ++ iter;
    polygon.vertices().insert(iter, new_vertex_idx);

    std::cerr << "Polygon p" << polygon_idx << " after:";
    for (KSR::size_t v : polygon.vertices())
      std::cerr << " v" << v;
    std::cerr << std::endl;
    
    const Support_line& support_line = m_support_lines[line_idx];
    const Support_plane& support_plane = m_support_planes[polygon.support_plane()];

    Line_2 line = support_plane.to_2d (support_line.line());

    Point_2 point = line.projection (v0.point());
    Vector_2 direction = v0.direction();
    v0.point() = point;
    v1.point() = point;

    const Point_2& previous_point = m_vertices[previous_vertex_idx].point();
    const Vector_2&  previous_direction = m_vertices[previous_vertex_idx].direction();
    
    const Point_2& next_point = m_vertices[next_vertex_idx].point();
    const Vector_2& next_direction = m_vertices[next_vertex_idx].direction();

    Point_2 moved_point = point + direction;
    Point_2 moved_previous_point = previous_point + previous_direction;
    Point_2 moved_next_point = next_point + next_direction;

    Point_2 target_previous = KSR::intersection_2<Point_2> (Line_2 (moved_previous_point, moved_point), line);
    Point_2 target_next = KSR::intersection_2<Point_2> (Line_2 (moved_next_point, moved_point), line);

    v1.direction() = Vector_2 (point, target_previous);
    v0.direction() = Vector_2 (point, target_next);

  }

  
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
