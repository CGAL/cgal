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

#ifndef CGAL_KSR_2_DATA_STRUCTURE_H
#define CGAL_KSR_2_DATA_STRUCTURE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <boost/function_output_iterator.hpp>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_2/Support_line.h>
#include <CGAL/KSR_2/Segment.h>
#include <CGAL/KSR_2/Vertex.h>
#include <CGAL/KSR/Event.h>
#include <CGAL/KSR/Event_queue.h>


namespace CGAL
{

namespace KSR_2
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
  typedef typename Kernel::Segment_2 Segment_2;

  typedef KSR_2::Support_line<Kernel> Support_line;
  typedef KSR_2::Segment<Kernel> Segment;
  typedef KSR_2::Vertex<Kernel> Vertex;

  typedef std::vector<Support_line> Support_lines;
  typedef std::vector<Segment> Segments;
  typedef std::vector<Vertex> Vertices;

  typedef KSR::Event<Kernel> Event;
  typedef KSR::Event_queue<Kernel> Event_queue;

private:

  // Main data structure
  Support_lines m_support_lines;
  Segments m_segments;
  Vertices m_vertices;
  Event_queue m_queue;

public:

  Data_structure() { }

  const Support_lines& support_lines() const { return m_support_lines; }
  const Vertices& vertices() const { return m_vertices; }

  Support_lines& support_lines() { return m_support_lines; }
  Vertices& vertices() { return m_vertices; }

  std::size_t number_of_vertices() const { return m_vertices.size(); }
  std::size_t number_of_segments() const { return m_segments.size(); }

  Segment_2 segment (std::size_t segment_idx) const
  {
    const Segment& segment = m_segments[segment_idx];
    const Support_line& support_line = m_support_lines[segment.support_line()];
    const Vertex& source = m_vertices[segment.source()];
    const Vertex& target = m_vertices[segment.target()];
    
    return Segment_2 (support_line.to_2d(source.point()), support_line.to_2d(target.point()));
  }

  bool is_bbox_segment (std::size_t segment_idx) const
  {
    return segment_idx < 4;
  }

  void add_bbox_as_segments (const CGAL::Bbox_2& bbox)
  {
    FT xmed = (bbox.xmin() + bbox.xmax()) / 2.;
    FT ymed = (bbox.ymin() + bbox.ymax()) / 2.;
    FT dx = (bbox.xmax() - bbox.xmin()) / 2.;
    FT dy = (bbox.ymax() - bbox.ymin()) / 2.;

    FT ratio = 1.1;
    FT xmin = xmed - ratio * dx;
    FT xmax = xmed + ratio * dx;
    FT ymin = ymed - ratio * dy;
    FT ymax = ymed + ratio * dy;
    
    std::array<Point_2, 4> bbox_points
      = { Point_2 (xmin, ymin),
          Point_2 (xmin, ymax),
          Point_2 (xmax, ymin),
          Point_2 (xmax, ymax) };
    
    add_segment (Segment_2 (bbox_points[0], bbox_points[1]));
    add_segment (Segment_2 (bbox_points[1], bbox_points[3]));
    add_segment (Segment_2 (bbox_points[3], bbox_points[2]));
    add_segment (Segment_2 (bbox_points[2], bbox_points[0]));
  }
  
  template <typename SegmentRange, typename SegmentMap>
  void add_segments (const SegmentRange& segments, SegmentMap segment_map)
  {
    for (const typename SegmentRange::const_iterator::value_type& vt : segments)
    {
      add_segment (get (segment_map, vt));
      initialize_vertices_directions(m_segments.size() - 1);
    }
  }

  void compute_support_lines()
  {
    // TODO
  }

  void make_segments_intersection_free()
  {
    // TODO
  }

  void initialize_queue(unsigned int k)
  {
    // Loop over vertices and schedule events
    for (std::size_t i = 0; i < m_vertices.size(); ++ i)
    {
      Vertex& vertex = m_vertices[i];
      if (vertex.is_frozen())
        continue;

      vertex.remaining_intersections() = k;

      Segment& segment = m_segments[vertex.segment()];
      Support_line& sli = m_support_lines[segment.support_line()];

      Ray_2 ray = sli.to_ray (vertex);

      for (std::size_t j = 0; j < m_support_lines.size(); ++ j)
      {
        if (j == vertex.support_line())
          continue;

        Support_line& slj_line = m_support_lines[j];
        Line_2 line = slj_line.line();

        Point_2 point;
        if (!KSR::intersection_2 (ray, line, point))
          continue;

        FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (sli.to_2d(vertex.point()), point));
        FT time = dist / CGAL::abs(vertex.speed());
        
        m_queue.push (Event (i, j, time));
      }
    }

//    m_queue.print();
  }

  void run()
  {
    FT latest_time = FT(0);

    std::size_t iterations = 0;
    while (!m_queue.empty())
    {
      Event ev = m_queue.pop();
//      std::cerr << " * Applying " << ev << std::endl;

      FT ellapsed_time = ev.time() - latest_time;
      latest_time = ev.time();

      advance_time (ellapsed_time);

      if (stop_vertex_if_intersection(ev.vertex(), ev.intersection_line()))
      {
//        std::cerr << "  -> Intersection happened" << std::endl;

        m_vertices[ev.vertex()].remaining_intersections() --;
        if (is_bbox_segment (ev.intersection_line()))
          m_vertices[ev.vertex()].remaining_intersections() = 0;
        
        // If there are still intersections to be made
        if (m_vertices[ev.vertex()].remaining_intersections() != 0)
        {
          // Create a new segment and transfer events
          m_segments.push_back (m_vertices[ev.vertex()].support_line());
          m_support_lines[m_vertices[ev.vertex()].support_line()].segments().push_back (m_segments.size() - 1);
          
          m_vertices.push_back (Vertex (m_vertices[ev.vertex()]));
          m_vertices.push_back (Vertex (m_vertices[ev.vertex()]));
          
          m_segments.back().source() = m_vertices.size() - 2;
          m_segments.back().target() = m_vertices.size() - 1;
          
          m_vertices[m_vertices.size() - 2].speed() = 0.; // Freeze other end
          
          m_queue.transfer_vertex_events (ev.vertex(), m_vertices.size() - 1);
        }
        else
          m_queue.remove_vertex_events (ev.vertex());

        m_vertices[ev.vertex()].speed() = 0.;
        
      }
      else
      {
//        std::cerr << "  -> Nothing happened" << std::endl;
      }
      
      ++ iterations;
      // if (iterations == 6)
      //   break;
    }
  }
  

private:

  void add_segment (const Segment_2 segment)
  {
    m_support_lines.push_back (Support_line(segment));
    m_segments.push_back (m_support_lines.size() - 1);
    m_support_lines.back().segments().push_back (m_segments.size() - 1);

    m_vertices.push_back (Vertex (m_support_lines.back().to_1d (segment.source()),
                                  m_segments.size() - 1, m_support_lines.size() - 1));
    m_vertices.push_back (Vertex (m_support_lines.back().to_1d (segment.target()),
                                  m_segments.size() - 1, m_support_lines.size() - 1));

    m_segments.back().source() = m_vertices.size() - 2;
    m_segments.back().target() = m_vertices.size() - 1;
  }

  void initialize_vertices_directions (std::size_t segment_idx)
  {
    Segment& segment = m_segments[segment_idx];
    Support_line& support_line = m_support_lines[segment.support_line()];

    Vertex& source = m_vertices[segment.source()];
    Vertex& target = m_vertices[segment.target()];

    Point_2 psource = support_line.to_2d (source.point());
    Point_2 ptarget = support_line.to_2d (target.point());

    if (Vector_2 (psource, ptarget) * support_line.line().to_vector() > 0.)
    {
      source.speed() = -1;
      target.speed() = 1;
    }
    else
    {
      source.speed() = 1;
      target.speed() = -1;
    }
    
  }

  void advance_time (FT time)
  {
    for (Vertex& v : m_vertices)
    {
      if (v.is_frozen())
        continue;

      v.point() = v.point() + time * v.speed();

    }
  }

  bool stop_vertex_if_intersection (std::size_t vertex_idx, std::size_t line_idx)
  {
    const Vertex& vertex = m_vertices[vertex_idx];
    const Support_line& intersecting = m_support_lines[vertex.support_line()];
    const Support_line& intersected = m_support_lines[line_idx];

    Point_2 point_inter = KSR::intersection_2<Point_2> (intersecting.line(), intersected.line());

    KSR::size_t intersected_segment = KSR::no_element();
    
    for (KSR::size_t sg : intersected.segments())
    {
      const Segment& segment = m_segments[sg];

      FT source = m_vertices[segment.source()].point();
      FT target = m_vertices[segment.target()].point();
      
      FT point = intersected.to_1d (point_inter);

      if ((source <= point && point <= target) ||
          (target <= point && point <= source))
      {
        intersected_segment = sg;
        break;
      }
    }

    if (intersected_segment == KSR::no_element()) // No intersection happened
      return false;

    m_vertices[vertex_idx].point() = intersecting.to_1d(point_inter);
    return true;
  }

  
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
