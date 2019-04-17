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
  typedef KSR_2::Segment Segment;
  typedef KSR_2::Vertex<FT> Vertex;
  
  typedef KSR_2::Meta_vertex<Point_2> Meta_vertex;

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
  
  FT m_current_time;
  
public:

  Data_structure()
    : m_current_time(0)
  { }

  void print() const
  {
    for (std::size_t i = 0; i < m_support_lines.size(); ++ i)
    {
      std::cerr << "* Support_line[" << i << "]" << std::endl;
      
      for (KSR::size_t segment_idx : m_support_lines[i].segments_idx())
      {
        std::cerr << "** Segment[" << segment_idx << "]" << std::endl;
        std::cerr << "*** Vertex[" << segment(segment_idx).source_idx() << "]" << std::endl;
        std::cerr << "*** Vertex[" << segment(segment_idx).target_idx() << "]" << std::endl;
      }
    }
  }

  Event_queue& queue() { return m_queue; }
  
  const FT& current_time() const { return m_current_time; }

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
  
  std::string segment_str (std::size_t segment_idx) const
  {
    return "Segment[" + std::to_string(segment_idx)
      + " from " + (segment(segment_idx).input_idx() == KSR::no_element() ?
                    "bbox" : std::to_string(segment(segment_idx).input_idx()))
      + "](v" + std::to_string(segment(segment_idx).source_idx())
      + "->v" + std::to_string(segment(segment_idx).target_idx())
      + ")";
  }
  std::string vertex_str (std::size_t vertex_idx) const
  {
    return "Vertex[" + std::to_string(vertex_idx) + "]";
  }

  // Vertex/idx -> Point_2
  inline Point_2 point_of_vertex (const Vertex& vertex, FT time) const
  { return support_line_of_vertex(vertex).to_2d(vertex.point(time)); }
  inline Point_2 point_of_vertex (std::size_t vertex_idx, FT time) const
  { return point_of_vertex (m_vertices[vertex_idx], time); }
  inline Point_2 point_of_vertex (const Vertex& vertex) const
  { return point_of_vertex (vertex, m_current_time); }
  inline Point_2 point_of_vertex (std::size_t vertex_idx) const
  { return point_of_vertex (vertex_idx, m_current_time); }

  // Vertex/idx -> Vector_2
  inline Vector_2 direction_of_vertex (const Vertex& vertex) const
  { return Vector_2 (support_line_of_vertex(vertex).to_2d(vertex.point(m_current_time)),
                     support_line_of_vertex(vertex).to_2d(vertex.point(m_current_time) + vertex.direction())); }
  inline Vector_2 direction_of_vertex (std::size_t vertex_idx) const
  { return direction_of_vertex (m_vertices[vertex_idx]); }

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
  
  bool has_meta_vertex (const Vertex& vertex) const
  { return vertex.meta_vertex_idx() != KSR::no_element(); }
  bool has_meta_vertex (KSR::size_t vertex_idx) const
  { return has_meta_vertex (m_vertices[vertex_idx]); }

  FT position_of_meta_vertex_on_support_line (KSR::size_t meta_vertex_idx, KSR::size_t support_line_idx) const
  { return support_line(support_line_idx).to_1d(meta_vertex(meta_vertex_idx).point()); }
   
  inline bool meta_vertex_exists (const Point_2& point) const
  { return m_meta_map.find(point) != m_meta_map.end(); }

  void get_vertices_of_meta_vertex (KSR::size_t meta_vertex_idx,
                                    std::vector<KSR::size_t>& vertices_idx) const
  {
    const Meta_vertex& meta_vertex = m_meta_vertices[meta_vertex_idx];
    for (KSR::size_t support_line_idx : meta_vertex.support_lines_idx())
    {
      const Support_line& support_line = m_support_lines[support_line_idx];
      for (KSR::size_t segment_idx : support_line.segments_idx())
      {
        const Segment& segment = m_segments[segment_idx];
        for (KSR::size_t vertex_idx : { segment.source_idx() , segment.target_idx() })
          if (m_vertices[vertex_idx].meta_vertex_idx() == meta_vertex_idx)
            vertices_idx.push_back (vertex_idx);
      }
    }
  }

  // Event -> Vertex
  const Vertex& vertex_of_event (const Event& ev) const
  { return m_vertices[ev.vertex_idx()]; }
  Vertex& vertex_of_event (const Event& ev)
  { return m_vertices[ev.vertex_idx()]; }

  inline CGAL::Bbox_2 bbox (const Vertex& vertex) const
  {
    return point_of_vertex(vertex).bbox();
  }
  inline CGAL::Bbox_2 bbox (const Support_line& support_line) const
  {
    return std::accumulate (support_line.segments_idx().begin(), support_line.segments_idx().end(),
                            CGAL::Bbox_2(),
                            [&](const CGAL::Bbox_2& bbox_2, const KSR::size_t& segment_idx) -> CGAL::Bbox_2
                            {
                              return bbox_2
                                + bbox(source_of_segment(segment_idx))
                                + bbox(target_of_segment(segment_idx));
                            });
  }

  bool is_segment_frozen (std::size_t segment_idx) const
  { return (source_of_segment(segment_idx).is_frozen() && target_of_segment(segment_idx).is_frozen()); }
  
  // idx -> Segment_2
  Segment_2 segment_2 (std::size_t segment_idx) const
  {
    const Segment& segment = m_segments[segment_idx];
    const Support_line& support_line = m_support_lines[segment.support_line_idx()];
    const Vertex& source = m_vertices[segment.source_idx()];
    const Vertex& target = m_vertices[segment.target_idx()];
    
    return Segment_2 (support_line.to_2d(source.point(m_current_time)), support_line.to_2d(target.point(m_current_time)));
  }

  bool is_bbox_support_line (std::size_t support_line_idx) const
  {
    return support_line_idx < 4;
  }

  bool is_bbox_segment (std::size_t segment_idx) const
  {
    return is_bbox_support_line(segment(segment_idx).support_line_idx());
  }

  bool is_bbox_meta_vertex (KSR::size_t meta_vertex_idx) const
  {
    for (KSR::size_t support_line_idx : meta_vertex(meta_vertex_idx).support_lines_idx())
      if (is_bbox_support_line(support_line_idx))
        return true;
    return false;
  }

  bool is_bbox_meta_edge (KSR::size_t source_idx, KSR::size_t target_idx) const
  {
    KSR::size_t common_line_idx = KSR::no_element();
    
    for (KSR::size_t support_line_idx : meta_vertex(source_idx).support_lines_idx())
      if (m_meta_vertices[target_idx].support_lines_idx().find(support_line_idx)
          != m_meta_vertices[target_idx].support_lines_idx().end())
      {
        common_line_idx = support_line_idx;
        break;
      }

    CGAL_assertion (common_line_idx != KSR::no_element());

    return is_bbox_support_line (common_line_idx);
  }

  bool is_meta_vertex_active (KSR::size_t meta_vertex_idx) const
  {
    for (KSR::size_t support_line_idx : meta_vertex(meta_vertex_idx).support_lines_idx())
      for (KSR::size_t segment_idx : support_line(support_line_idx).segments_idx())
        for (KSR::size_t vertex_idx : { segment(segment_idx).source_idx(), segment(segment_idx).target_idx() })
          if (vertex(vertex_idx).meta_vertex_idx() == meta_vertex_idx)
            return true;
    return false;
  }

  bool is_meta_vertex_intersection (KSR::size_t meta_vertex_idx) const
  {
    bool found_one = false;
    
    for (KSR::size_t support_line_idx : meta_vertex(meta_vertex_idx).support_lines_idx())
    {
      bool broken = false;
      for (KSR::size_t segment_idx : support_line(support_line_idx).segments_idx())
      {
        for (KSR::size_t vertex_idx : { segment(segment_idx).source_idx(), segment(segment_idx).target_idx() })
        {
          if (vertex(vertex_idx).meta_vertex_idx() == meta_vertex_idx)
          {
            if (found_one)
              return true;
            found_one = true;
            broken = true;
            break;
          }
        }
        if (broken)
          break;
      }
    }
    return false;
  }

  bool is_meta_vertex_deadend_of_vertex (KSR::size_t meta_vertex_idx, KSR::size_t vertex_idx) const
  {
    return meta_vertex(meta_vertex_idx).is_deadend_of (segment_of_vertex(vertex_idx).support_line_idx());
  }

  void make_meta_vertex_deadend_of_vertex (KSR::size_t meta_vertex_idx, KSR::size_t vertex_idx)
  {
    meta_vertex(meta_vertex_idx).make_deadend_of (segment_of_vertex(vertex_idx).support_line_idx());
  }

  void make_meta_vertex_no_longer_deadend_of_vertex (KSR::size_t meta_vertex_idx, KSR::size_t vertex_idx)
  {
    meta_vertex(meta_vertex_idx).make_no_longer_deadend_of (segment_of_vertex(vertex_idx).support_line_idx());
  }

  bool are_support_lines_connected (KSR::size_t support_line_0,
                                    KSR::size_t support_line_1) const
  {
    for (KSR::size_t meta_vertex_idx : support_line(support_line_0).meta_vertices_idx())
      if (m_meta_vertices[meta_vertex_idx].support_lines_idx().find(support_line_1)
          != m_meta_vertices[meta_vertex_idx].support_lines_idx().end())
        return true;
    return false;
  }

  KSR::size_t add_support_line (const Segment_2& segment)
  {
    m_support_lines.push_back (Support_line(segment));
    return KSR::size_t(m_support_lines.size() - 1);
  }

  Segment& add_segment (const Segment_2 segment, KSR::size_t input_idx = KSR::no_element())
  {
    // Check if support line exists first
    Support_line new_support_line (segment);
    KSR::size_t support_line_idx = KSR::no_element();
    for (std::size_t i = 0; i < m_support_lines.size(); ++ i)
      if (new_support_line == m_support_lines[i])
      {
        support_line_idx = i;
        break;
      }

    if (support_line_idx == KSR::no_element())
    {
      support_line_idx = m_support_lines.size();
      m_support_lines.push_back (new_support_line);

      if (input_idx == KSR::no_element())
      {
        m_support_lines.back().minimum() = m_support_lines.back().to_1d (segment.source());
        m_support_lines.back().maximum() = m_support_lines.back().to_1d (segment.target());
      }
      else
      {
        FT max_negative = -std::numeric_limits<FT>::max();
        FT min_positive = std::numeric_limits<FT>::max();

        for (std::size_t i = 0; i < 4; ++ i)
        {
          Point_2 point;
          if (!KSR::intersection_2 (m_support_lines[i].line(), m_support_lines.back().line(), point))
            continue;

          FT position = m_support_lines.back().to_1d (point);
          if (position < 0 && position > max_negative)
            max_negative = position;
          if (position > 0 && position < min_positive)
            min_positive = position;
        }

        CGAL_assertion (max_negative != -std::numeric_limits<FT>::max()
                        && min_positive != -std::numeric_limits<FT>::min());

        m_support_lines.back().minimum() = max_negative;
        m_support_lines.back().maximum() = min_positive;
      }
    }
    else
      m_support_lines[support_line_idx].connected_components() ++;

    
    KSR::size_t segment_idx = m_segments.size();
    m_segments.push_back (Segment(input_idx, support_line_idx));
    m_support_lines[support_line_idx].segments_idx().push_back (segment_idx);
    
    KSR::size_t source_idx = m_vertices.size();
    m_vertices.push_back (Vertex (m_support_lines[support_line_idx].to_1d (segment.source()),
                                  segment_idx, m_support_lines.size() - 1));
    KSR::size_t target_idx = m_vertices.size();
    m_vertices.push_back (Vertex (m_support_lines[support_line_idx].to_1d (segment.target()),
                                  segment_idx, m_support_lines.size() - 1));

    // Keep segment ordered
    if (m_vertices[source_idx].point(0) > m_vertices[target_idx].point(0))
      std::swap (source_idx, target_idx);
    
    m_segments[segment_idx].source_idx() = source_idx;
    m_segments[segment_idx].target_idx() = target_idx;
    return m_segments.back();
  }

  KSR::size_t add_meta_vertex (const Point_2& point,
                               KSR::size_t support_line_idx_0,
                               KSR::size_t support_line_idx_1 = KSR::no_element())
  {
    // Avoid several points almost equal
    Point_2 p (CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.x()) / CGAL_KSR_SAME_POINT_TOLERANCE),
               CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.y()) / CGAL_KSR_SAME_POINT_TOLERANCE));
      
    typename std::map<Point_2, KSR::size_t>::iterator iter;
    bool inserted = false;
    std::tie (iter, inserted) = m_meta_map.insert (std::make_pair (p, KSR::size_t(m_meta_vertices.size())));
    if (inserted)
      m_meta_vertices.push_back (Meta_vertex(p));

    KSR::size_t meta_vertex_idx = iter->second;

    CGAL_KSR_CERR_3 << "** Adding meta vertex " << meta_vertex_idx << " between "
                    << support_line_idx_0 << " and " << support_line_idx_1
                    << " at point " << p << std::endl;

    for (KSR::size_t support_line_idx : { support_line_idx_0, support_line_idx_1 })
    {
      if (support_line_idx != KSR::no_element())
      {
        m_meta_vertices[meta_vertex_idx].support_lines_idx().insert (support_line_idx);

        if (std::find(m_support_lines[support_line_idx].meta_vertices_idx().begin(),
                                 m_support_lines[support_line_idx].meta_vertices_idx().end(),
                      meta_vertex_idx) == m_support_lines[support_line_idx].meta_vertices_idx().end())
          m_support_lines[support_line_idx].meta_vertices_idx().push_back (meta_vertex_idx);
      }
    }

    // Special case = meta vertex is deadend of one line
    if (support_line_idx_1 == KSR::no_element())
    {
      meta_vertex(meta_vertex_idx).make_deadend_of (support_line_idx_0);
      CGAL_KSR_CERR_3 << "*** Meta_vertex[" << meta_vertex_idx
                      << "] is deadend of Support_line[" << support_line_idx_0 << "]" << std::endl;
    }

    return meta_vertex_idx;
  }

  KSR::size_t add_meta_vertex (KSR::size_t support_line_idx_0, KSR::size_t support_line_idx_1)
  {
    Point_2 inter_point = KSR::intersection_2<Point_2> (support_line(support_line_idx_0).line(),
                                                        support_line(support_line_idx_1).line());

    return add_meta_vertex (inter_point, support_line_idx_0, support_line_idx_1);
  }

  void attach_vertex_to_meta_vertex (KSR::size_t vertex_idx, KSR::size_t meta_vertex_idx)
  {
    CGAL_assertion (!has_meta_vertex(vertex_idx));
    CGAL_assertion_msg (m_meta_vertices[meta_vertex_idx].support_lines_idx().find
                        (segment_of_vertex(vertex_idx).support_line_idx())
                        != m_meta_vertices[meta_vertex_idx].support_lines_idx().end(),
                        "Trying to attach a vertex to a meta vertex not on its support line");
    m_vertices[vertex_idx].meta_vertex_idx() = meta_vertex_idx;
  }

  void cut_segment (KSR::size_t segment_idx, KSR::size_t meta_vertex_idx)
  {
    std::vector<KSR::size_t> vec (1, meta_vertex_idx);
    cut_segment (segment_idx, vec);
  }

  void cut_segment (KSR::size_t segment_idx, std::vector<KSR::size_t>& meta_vertices_idx)
  {
    CGAL_KSR_CERR_3 << "** Cutting " << segment_str(segment_idx) << std::endl;

    Segment& segment = m_segments[segment_idx];
    KSR::size_t source_idx = segment.source_idx();
    KSR::size_t target_idx = segment.target_idx();
    
    Support_line& support_line = support_line_of_segment(segment_idx);

    std::sort (meta_vertices_idx.begin(), meta_vertices_idx.end(),
               [&](const KSR::size_t& a,
                   const KSR::size_t& b) -> bool
               {
                 return (position_of_meta_vertex_on_support_line(a, segment.support_line_idx())
                         < position_of_meta_vertex_on_support_line(b, segment.support_line_idx()));
               });

    std::size_t nb_segments_before = m_segments.size();
    std::size_t nb_vertices_before = m_vertices.size();

    // Attach to existing endpoint
    KSR::size_t new_target_idx = m_vertices.size();
    m_vertices.push_back (Vertex (position_of_meta_vertex_on_support_line(meta_vertices_idx.front(),
                                                                          segment.support_line_idx())));
    m_vertices[new_target_idx].segment_idx() = segment_idx;
    segment.target_idx() = new_target_idx;
    attach_vertex_to_meta_vertex (new_target_idx, meta_vertices_idx.front());

    // Create new segments
    for (std::size_t i = 0; i < meta_vertices_idx.size() - 1; ++ i)
    {
      KSR::size_t sidx = m_segments.size();
      m_segments.push_back (Segment (segment.input_idx(), segment.support_line_idx()));
      support_line.segments_idx().push_back (sidx);

      KSR::size_t source_idx = m_vertices.size();
      m_vertices.push_back (Vertex (position_of_meta_vertex_on_support_line(meta_vertices_idx[i],
                                                                            segment.support_line_idx())));
      m_vertices[source_idx].segment_idx() = sidx;
      m_segments[sidx].source_idx() = source_idx;
      attach_vertex_to_meta_vertex (source_idx, meta_vertices_idx[i]);
      
      KSR::size_t target_idx = m_vertices.size();
      m_vertices.push_back (Vertex (position_of_meta_vertex_on_support_line(meta_vertices_idx[i+1],
                                                                            segment.support_line_idx())));
      m_vertices[target_idx].segment_idx() = sidx;
      m_segments[sidx].target_idx() = target_idx;
      attach_vertex_to_meta_vertex (source_idx, meta_vertices_idx[i+1]);
    }

    // Create final segment and attach to existing endpoint
    KSR::size_t sidx = m_segments.size();
    m_segments.push_back (Segment (segment.input_idx(), segment.support_line_idx()));
    support_line.segments_idx().push_back (sidx);

    KSR::size_t new_source_idx = m_vertices.size();
    m_vertices.push_back (Vertex (position_of_meta_vertex_on_support_line(meta_vertices_idx.back(),
                                                                            segment.support_line_idx())));
    m_vertices[new_source_idx].segment_idx() = sidx;
    m_segments[sidx].source_idx() = new_source_idx;
    attach_vertex_to_meta_vertex (new_source_idx, meta_vertices_idx.back());
      
    m_vertices[target_idx].segment_idx() = sidx;
    m_segments[sidx].target_idx() = target_idx;

    CGAL_KSR_CERR_3 << "*** new vertices:";
    for (std::size_t i = nb_vertices_before; i < m_vertices.size(); ++ i)
      CGAL_KSR_CERR_3 << " " << vertex_str(i);
    CGAL_KSR_CERR_3 << std::endl;
    
    CGAL_KSR_CERR_3 << "*** new segments: " << segment_str(segment_idx);
    for (std::size_t i = nb_segments_before; i < m_segments.size(); ++ i)
      CGAL_KSR_CERR_3 << " " << segment_str(i);
    CGAL_KSR_CERR_3 << std::endl;
  }

  void connect_vertices (KSR::size_t vertex_1_idx, KSR::size_t vertex_2_idx, const Point_2& point)
  {
    Vertex& vertex_1 = vertex(vertex_1_idx);
    Vertex& vertex_2 = vertex(vertex_2_idx);

    KSR::size_t meta_vertex_idx;
    if (vertex_1.meta_vertex_idx() != KSR::no_element())
    {
      CGAL_assertion (vertex_2.meta_vertex_idx() == KSR::no_element());
      vertex_2.meta_vertex_idx() = vertex_1.meta_vertex_idx();
      m_meta_vertices[vertex_1.meta_vertex_idx()].vertices_idx().push_back(vertex_2_idx);
      m_meta_map[point] = vertex_1.meta_vertex_idx();
    }
    else if (vertex_2.meta_vertex_idx() != KSR::no_element())
    {
      CGAL_assertion (vertex_1.meta_vertex_idx() == KSR::no_element());
      vertex_1.meta_vertex_idx() = vertex_2.meta_vertex_idx();
      m_meta_vertices[vertex_2.meta_vertex_idx()].vertices_idx().push_back(vertex_1_idx);
      m_meta_map[point] = vertex_2.meta_vertex_idx();
    }
    else
    {
      KSR::size_t meta_vertex_idx = add_meta_vertex(point, vertex_1_idx);
      vertex_2.meta_vertex_idx() = vertex_1.meta_vertex_idx();
      m_meta_vertices[vertex_1.meta_vertex_idx()].vertices_idx().push_back(vertex_2_idx);
      m_meta_map[point] = vertex_1.meta_vertex_idx();
    }
  }

  KSR::size_t propagate_segment (std::size_t vertex_idx)
  {
    CGAL_KSR_CERR_3 << "** Propagating " << vertex_str(vertex_idx) << std::endl;

    // Create a new segment
    KSR::size_t segment_idx = m_segments.size();
    m_segments.push_back (Segment(segment_of_vertex(vertex_idx).input_idx(),
                                  segment_of_vertex(vertex_idx).support_line_idx()));
    support_line_of_vertex(vertex_idx).segments_idx().push_back (segment_idx);

    // Create new vertices
    KSR::size_t source_idx = m_vertices.size();
    m_vertices.push_back (Vertex (m_vertices[vertex_idx]));
    KSR::size_t target_idx = m_vertices.size();
    m_vertices.push_back (Vertex (m_vertices[vertex_idx]));

    // Connect segments and vertices
    m_segments[segment_idx].source_idx() = source_idx;
    m_segments[segment_idx].target_idx() = target_idx;
    m_vertices[source_idx].segment_idx() = segment_idx;
    m_vertices[target_idx].segment_idx() = segment_idx;

    CGAL_assertion (m_vertices[vertex_idx].direction() != 0);

    // Keep vertices ordered on the segment
    if (m_vertices[vertex_idx].direction() < 0)
      std::swap (source_idx, target_idx);

    // Freeze one end
    m_vertices[source_idx].freeze(m_current_time);

    // Release other end
    m_vertices[target_idx].meta_vertex_idx() = KSR::no_element();

    CGAL_KSR_CERR_3 << "*** new vertices: " << vertex_str (source_idx)
                  << " " << vertex_str (target_idx) << std::endl;
    CGAL_KSR_CERR_3 << "*** new segment: " << segment_str(segment_idx) << std::endl;

    return target_idx;
  }

  void update_positions (FT time)
  {
    m_current_time = time;
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
