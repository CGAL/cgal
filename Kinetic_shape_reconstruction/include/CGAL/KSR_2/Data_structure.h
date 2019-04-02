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
#include <CGAL/KSR/verbosity.h>
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
  const Vertex& vertex (std::size_t idx) const { return m_vertices[idx]; }
  Vertex& vertex (std::size_t idx) { return m_vertices[idx]; }
  
  std::size_t number_of_segments() const { return m_segments.size(); }
  const Segment& segment (std::size_t idx) const { return m_segments[idx]; }
  Segment& segment (std::size_t idx) { return m_segments[idx]; }

  std::size_t number_of_support_lines() const { return m_support_lines.size(); }
  const Support_line& support_line (std::size_t idx) const { return m_support_lines[idx]; }
  Support_line& support_line (std::size_t idx) { return m_support_lines[idx]; }

  // Vertex/idx -> Point_2
  inline Point_2 point_of_vertex (const Vertex& vertex) const
  { return support_line_of_vertex(vertex).to_2d(vertex.point()); }
  inline Point_2 point_of_vertex (std::size_t vertex_idx) const
  { return point_of_vertex (m_vertices[vertex_idx]); }

  // Vertex/idx -> Segment
  inline const Segment& segment_of_vertex (const Vertex& vertex) const
  { return m_segments[vertex.segment_idx()]; }
  inline Segment& segment_of_vertex(const Vertex& vertex)
  { return m_segments[vertex.segment_idx()]; }
  inline const Segment& segment_of_vertex (std::size_t vertex_idx) const
  { return segment_of_vertex(m_vertices[vertex_idx]); }
  inline Segment& segment_of_vertex (std::size_t vertex_idx)
  { return segment_of_vertex(m_vertices[vertex_idx]); }

  // Segment -> source Vertex
  inline const Vertex& source_of_segment (const Segment& segment) const
  { return m_vertices[segment.source_idx()]; }
  inline Vertex& source_of_segment (const Segment& segment)
  { return m_vertices[segment.source_idx()]; }

  // Segment -> target Vertex
  inline const Vertex& target_of_segment (const Segment& segment) const
  { return m_vertices[segment.target_idx()]; }
  inline Vertex& target_of_segment (const Segment& segment)
  { return m_vertices[segment.target_idx()]; }

  // Segment/idx -> Support_line
  inline const Support_line& support_line_of_segment (const Segment& segment) const
  { return m_support_lines[segment.support_line_idx()]; }
  inline Support_line& support_line_of_segment (const Segment& segment)
  { return m_support_lines[segment.support_line_idx()]; }
  inline const Support_line& support_line_of_segment (std::size_t segment_idx) const
  { return support_line_of_segment(m_segments[segment_idx]); }
  inline Support_line& support_line_of_segment (std::size_t segment_idx)
  { return support_line_of_segment(m_segments[segment_idx]); }
  
  // Vertex/idx -> Support_line
  inline const Support_line& support_line_of_vertex (const Vertex& vertex) const
  { return support_line_of_segment(vertex.segment_idx()); }
  inline Support_line& support_line_of_vertex (const Vertex& vertex)
  { return support_line_of_segment(vertex.segment_idx()); }
  inline const Support_line& support_line_of_vertex (std::size_t vertex_idx) const
  { return support_line_of_vertex(m_vertices[vertex_idx]); }
  inline Support_line& support_line_of_vertex (std::size_t vertex_idx)
  { return support_line_of_vertex(m_vertices[vertex_idx]); }

  // Event -> Vertex
  const Vertex& vertex_of_event (const Event& ev) const
  { return m_vertices[ev.vertex_idx()]; }
  Vertex& vertex_of_event (const Event& ev)
  { return m_vertices[ev.vertex_idx()]; }

  // idx -> Segment_2
  Segment_2 segment_2 (std::size_t segment_idx) const
  {
    const Segment& segment = m_segments[segment_idx];
    const Support_line& support_line = m_support_lines[segment.support_line_idx()];
    const Vertex& source = m_vertices[segment.source_idx()];
    const Vertex& target = m_vertices[segment.target_idx()];
    
    return Segment_2 (support_line.to_2d(source.point()), support_line.to_2d(target.point()));
  }

  bool is_bbox_segment (std::size_t segment_idx) const
  {
    return segment_idx < 4;
  }

  KSR::size_t add_support_line (const Segment_2& segment)
  {
    m_support_lines.push_back (Support_line(segment));
    return KSR::size_t(m_support_lines.size() - 1);
  }

  Segment& add_segment (const Segment_2 segment)
  {
    m_support_lines.push_back (Support_line(segment));
    m_segments.push_back (m_support_lines.size() - 1);
    m_support_lines.back().segments_idx().push_back (m_segments.size() - 1);

    m_vertices.push_back (Vertex (m_support_lines.back().to_1d (segment.source()),
                                  m_segments.size() - 1, m_support_lines.size() - 1));
    m_vertices.push_back (Vertex (m_support_lines.back().to_1d (segment.target()),
                                  m_segments.size() - 1, m_support_lines.size() - 1));

    m_segments.back().source_idx() = m_vertices.size() - 2;
    m_segments.back().target_idx() = m_vertices.size() - 1;
    return m_segments.back();
  }

  Segment& propagate_segment (const Vertex& vertex)
  {
    // Create a new segment and transfer events
    m_segments.push_back (Segment(segment_of_vertex(vertex).support_line_idx()));
    support_line_of_vertex(vertex).segments_idx().push_back (m_segments.size() - 1);
    
    m_vertices.push_back (Vertex (vertex));
    m_vertices.push_back (Vertex (vertex));
          
    m_segments.back().source_idx() = m_vertices.size() - 2;
    m_segments.back().target_idx() = m_vertices.size() - 1;

    // Freeze one end
    m_vertices[m_vertices.size() - 2].remaining_intersections() = 0;
    m_vertices[m_vertices.size() - 2].direction() = 0.;
    
    return m_segments.back();
  }

  void push_to_queue (const Event& ev) { m_queue.push(ev); }
  bool queue_is_empty() const { return m_queue.empty(); }
  Event queue_pop() { return m_queue.pop(); }
  void transfer_events (std::size_t old_vertex, std::size_t new_vertex)
  { m_queue.transfer_vertex_events (old_vertex, new_vertex); }
  void remove_events (std::size_t vertex_idx) { m_queue.remove_vertex_events (vertex_idx); };

  void advance_time (FT time)
  {
    for (Vertex& v : m_vertices)
    {
      if (v.is_frozen())
        continue;

      v.point() = v.point() + time * v.direction();

    }
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
