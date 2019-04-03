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

#include <CGAL/KSR_2/Meta_vertex.h>

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
  
  typedef KSR_2::Meta_vertex<Kernel> Meta_vertex;

  typedef std::vector<Support_line> Support_lines;
  typedef std::vector<Segment> Segments;
  typedef std::vector<Vertex> Vertices;
  
  typedef std::vector<Meta_vertex> Meta_vertices;

  typedef KSR::Event<Kernel> Event;
  typedef KSR::Event_queue<Kernel> Event_queue;

private:

  // Main data structure
  Support_lines m_support_lines;
  Segments m_segments;
  Vertices m_vertices;
  
  Meta_vertices m_meta_vertices;
  
  Event_queue m_queue;

  // Helping data structures
  std::map<Point_2, KSR::size_t> m_meta_map;
  
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

  std::size_t number_of_meta_vertices() const { return m_meta_vertices.size(); }
  const Meta_vertex& meta_vertex (std::size_t idx) const { return m_meta_vertices[idx]; }
  Meta_vertex& meta_vertex (std::size_t idx) { return m_meta_vertices[idx]; }
  
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

  // Segment/idx -> source Vertex
  inline const Vertex& source_of_segment (const Segment& segment) const
  { return m_vertices[segment.source_idx()]; }
  inline Vertex& source_of_segment (const Segment& segment)
  { return m_vertices[segment.source_idx()]; }
  inline const Vertex& source_of_segment (std::size_t segment_idx) const
  { return source_of_segment(m_segments[segment_idx]); }
  inline Vertex& source_of_segment (std::size_t segment_idx)
  { return source_of_segment(m_segments[segment_idx]); }

  // Segment/idx -> target Vertex
  inline const Vertex& target_of_segment (const Segment& segment) const
  { return m_vertices[segment.target_idx()]; }
  inline Vertex& target_of_segment (const Segment& segment)
  { return m_vertices[segment.target_idx()]; }
  inline const Vertex& target_of_segment (std::size_t segment_idx) const
  { return target_of_segment(m_segments[segment_idx]); }
  inline Vertex& target_of_segment (std::size_t segment_idx)
  { return target_of_segment(m_segments[segment_idx]); }


  // idx -> opposite Vertex
  inline const Vertex& opposite_vertex (std::size_t vertex_idx) const
  {
    const Segment& segment = segment_of_vertex(vertex_idx);

    CGAL_assertion (segment.source_idx() == vertex_idx
                    || segment.target_idx() == vertex_idx);
    
    return (segment.source_idx() == vertex_idx ?
            m_vertices[segment.target_idx()] :
            m_vertices[segment.source_idx()]);
  }
    

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

  // Vertex/idx -> Meta_vertex
  inline const Meta_vertex& meta_vertex_of_vertex (const Vertex& vertex) const
  { return m_meta_vertices[vertex.meta_vertex_idx()]; }
  inline Meta_vertex& meta_vertex_of_vertex (const Vertex& vertex)
  { return m_meta_vertices[vertex.meta_vertex_idx()]; }
  inline const Meta_vertex& meta_vertex_of_vertex (std::size_t vertex_idx) const
  { return meta_vertex_of_vertex(m_vertices[vertex_idx]); }
  inline Meta_vertex& meta_vertex_of_vertex (std::size_t vertex_idx)
  { return meta_vertex_of_vertex(m_vertices[vertex_idx]); }

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

  void add_meta_vertex (const Point_2& p, std::size_t vertex_idx)
  {
    CGAL_assertion (m_vertices[vertex_idx].meta_vertex_idx() == KSR::no_element());
    
    typename std::map<Point_2, KSR::size_t>::iterator iter;
    bool inserted = false;
    std::tie (iter, inserted) = m_meta_map.insert (std::make_pair (p, KSR::size_t(m_meta_vertices.size())));
    if (inserted)
      m_meta_vertices.push_back (Meta_vertex(p));
    
    m_meta_vertices[iter->second].vertices_idx().push_back (vertex_idx);

    m_vertices[vertex_idx].meta_vertex_idx() = iter->second;
  }

  void add_meta_vertex (std::size_t vertex_idx)
  {
    add_meta_vertex (point_of_vertex(vertex_idx), vertex_idx);
  }

  void cut_segment (KSR::size_t segment_idx, const Point_2& point)
  {
    std::vector<Point_2> vec (1, point);
    cut_segment (segment_idx, vec);
  }
    
  void cut_segment (KSR::size_t segment_idx, std::vector<Point_2>& points)
  {
    Segment& segment = m_segments[segment_idx];
    KSR::size_t source_idx = segment.source_idx();
    KSR::size_t target_idx = segment.target_idx();
    
    Support_line& support_line = support_line_of_segment(segment_idx);

    points.push_back (point_of_vertex(source_idx));
    points.push_back (point_of_vertex(target_idx));

    std::sort (points.begin(), points.end(),
               [&](const Point_2& a, const Point_2& b) -> bool
               {
                 return support_line.to_1d(a) < support_line.to_1d(b);
               });

    for (std::size_t i = 0; i < points.size() - 1; ++ i)
    {
      KSR::size_t sidx = segment_idx;
      if (i != 0)
      {
        sidx = m_segments.size();
        m_segments.push_back (Segment (segment.support_line_idx()));
        support_line.segments_idx().push_back (m_segments.size() - 1);
      }

      for (std::size_t j = 0; j < 2; ++ j)
      {
        KSR::size_t vertex_idx = KSR::no_element();
        if (points[i+j] == point_of_vertex(source_idx))
          vertex_idx = source_idx;
        else if (points[i+j] == point_of_vertex(target_idx))
          vertex_idx = target_idx;
        else
        {
          vertex_idx = m_vertices.size();
          m_vertices.push_back (Vertex (support_line.to_1d(points[i+j])));
          add_meta_vertex(points[i+j], m_vertices.size() - 1);
        }

        m_vertices[vertex_idx].segment_idx() = sidx;

        if (j == 0)
          m_segments[sidx].source_idx() = vertex_idx;
        else
          m_segments[sidx].target_idx() = vertex_idx;
      }
    }

  }

  Segment& propagate_segment (std::size_t vertex_idx)
  {
    // Create a new segment
    m_segments.push_back (Segment(segment_of_vertex(vertex_idx).support_line_idx()));
    support_line_of_vertex(vertex_idx).segments_idx().push_back (m_segments.size() - 1);

    // Create new vertices
    m_vertices.push_back (Vertex (m_vertices[vertex_idx]));
    m_vertices.push_back (Vertex (m_vertices[vertex_idx]));

    // Connect segments and vertices
    m_segments.back().source_idx() = m_vertices.size() - 2;
    m_segments.back().target_idx() = m_vertices.size() - 1;
    m_vertices[m_vertices.size() - 2].segment_idx() = m_segments.size() - 1;
    m_vertices[m_vertices.size() - 1].segment_idx() = m_segments.size() - 1;

    // Freeze one end
    meta_vertex_of_vertex(vertex_idx).vertices_idx().push_back (m_vertices.size() - 2);
    m_vertices[m_vertices.size() - 2].remaining_intersections() = 0;
    m_vertices[m_vertices.size() - 2].direction() = 0.;
    
    // Release other end
    m_vertices[m_vertices.size() - 1].meta_vertex_idx() = KSR::no_element();
    
    return m_segments.back();
  }

  void push_to_queue (const Event& ev) { m_queue.push(ev); }
  bool queue_is_empty() const { return m_queue.empty(); }
  Event queue_pop() { return m_queue.pop(); }
  void transfer_events (std::size_t old_vertex, std::size_t new_vertex)
  { m_queue.transfer_vertex_events (old_vertex, new_vertex); }
  void remove_events (std::size_t vertex_idx) { m_queue.remove_vertex_events (vertex_idx); };

  void update_positions (FT time)
  {
    for (Vertex& v : m_vertices)
      v.update_position(time);
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
