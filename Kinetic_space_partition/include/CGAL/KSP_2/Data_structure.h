// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSP_2_DATA_STRUCTURE_H
#define CGAL_KSP_2_DATA_STRUCTURE_H

#include <CGAL/license/Kinetic_space_partition.h>

#include <CGAL/KSP/utils.h>
#include <CGAL/KSP_2/Support_line.h>
#include <CGAL/KSP_2/Segment.h>
#include <CGAL/KSP_2/Vertex.h>

#include <CGAL/KSP_2/Meta_vertex.h>

namespace CGAL {
namespace KSP_2 {
namespace internal {

template <typename GeomTraits>
class Data_structure {
public:

  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Segment_2 Segment_2;

  typedef Support_line<Kernel> Support_line;
  typedef Segment Segment;
  typedef Vertex<FT> Vertex;

  typedef Meta_vertex<Point_2> Meta_vertex;

  typedef std::vector<Support_line> Support_lines;
  typedef std::vector<Segment> Segments;
  typedef std::vector<Vertex> Vertices;

  typedef std::vector<Meta_vertex> Meta_vertices;

private:

  // Main data structure
  Support_lines m_support_lines;
  Segments m_segments;
  Vertices m_vertices;

  Meta_vertices m_meta_vertices;

  // Helping data structures
  std::map<Point_2, std::size_t> m_meta_map;

  FT m_current_time;

public:

  Data_structure()
    : m_current_time(0)
  { }

  void print() const
  {
    for (std::size_t i = 0; i < m_support_lines.size(); ++i)
    {
      std::cerr << "* Support_line[" << i << "]" << std::endl;

      for (std::size_t segment_idx : m_support_lines[i].segments_idx())
      {
        std::cerr << "** Segment[" << segment_idx << "]" << std::endl;
        std::cerr << "*** Vertex[" << segment(segment_idx).source_idx() << "]" << std::endl;
        std::cerr << "*** Vertex[" << segment(segment_idx).target_idx() << "]" << std::endl;
      }
    }
  }

  const FT& current_time() const { return m_current_time; }

  std::size_t number_of_vertices() const { return m_vertices.size(); }
  const Vertex& vertex(std::size_t idx) const { return m_vertices[idx]; }
  Vertex& vertex(std::size_t idx) { return m_vertices[idx]; }

  std::size_t number_of_segments() const { return m_segments.size(); }
  const Segment& segment(std::size_t idx) const { return m_segments[idx]; }
  Segment& segment(std::size_t idx) { return m_segments[idx]; }

  std::size_t number_of_support_lines() const { return m_support_lines.size(); }
  const Support_line& support_line(std::size_t idx) const { return m_support_lines[idx]; }
  Support_line& support_line(std::size_t idx) { return m_support_lines[idx]; }

  std::size_t number_of_meta_vertices() const { return m_meta_vertices.size(); }
  const Meta_vertex& meta_vertex(std::size_t idx) const { return m_meta_vertices[idx]; }
  Meta_vertex& meta_vertex(std::size_t idx) { return m_meta_vertices[idx]; }

  std::string segment_str(std::size_t segment_idx) const
  {
    return "Segment[" + std::to_string(segment_idx)
      + " from " + (segment(segment_idx).input_idx() == std::size_t(-1) ?
        "bbox" : std::to_string(segment(segment_idx).input_idx()))
      + "](v" + std::to_string(segment(segment_idx).source_idx())
      + "->v" + std::to_string(segment(segment_idx).target_idx())
      + ")";
  }
  std::string vertex_str(std::size_t vertex_idx) const
  {
    return "Vertex[" + std::to_string(vertex_idx) + "]";
  }

  // Vertex/idx -> Point_2
  inline Point_2 point_of_vertex(const Vertex& vertex, FT time) const
  {
    return support_line_of_vertex(vertex).to_2d(vertex.point(time));
  }
  inline Point_2 point_of_vertex(std::size_t vertex_idx, FT time) const
  {
    return point_of_vertex(m_vertices[vertex_idx], time);
  }
  inline Point_2 point_of_vertex(const Vertex& vertex) const
  {
    return point_of_vertex(vertex, m_current_time);
  }
  inline Point_2 point_of_vertex(std::size_t vertex_idx) const
  {
    return point_of_vertex(vertex_idx, m_current_time);
  }

  // Vertex/idx -> Vector_2
  inline Vector_2 direction_of_vertex(const Vertex& vertex) const
  {
    return Vector_2(support_line_of_vertex(vertex).to_2d(vertex.point(m_current_time)),
      support_line_of_vertex(vertex).to_2d(vertex.point(m_current_time) + vertex.direction()));
  }
  inline Vector_2 direction_of_vertex(std::size_t vertex_idx) const
  {
    return direction_of_vertex(m_vertices[vertex_idx]);
  }

  // Vertex/idx -> Segment
  inline const Segment& segment_of_vertex(const Vertex& vertex) const
  {
    return m_segments[vertex.segment_idx()];
  }
  inline Segment& segment_of_vertex(const Vertex& vertex)
  {
    return m_segments[vertex.segment_idx()];
  }
  inline const Segment& segment_of_vertex(std::size_t vertex_idx) const
  {
    return segment_of_vertex(m_vertices[vertex_idx]);
  }
  inline Segment& segment_of_vertex(std::size_t vertex_idx)
  {
    return segment_of_vertex(m_vertices[vertex_idx]);
  }

  // Segment/idx -> source Vertex
  inline const Vertex& source_of_segment(const Segment& segment) const
  {
    return m_vertices[segment.source_idx()];
  }
  inline Vertex& source_of_segment(const Segment& segment)
  {
    return m_vertices[segment.source_idx()];
  }
  inline const Vertex& source_of_segment(std::size_t segment_idx) const
  {
    return source_of_segment(m_segments[segment_idx]);
  }
  inline Vertex& source_of_segment(std::size_t segment_idx)
  {
    return source_of_segment(m_segments[segment_idx]);
  }

  // Segment/idx -> target Vertex
  inline const Vertex& target_of_segment(const Segment& segment) const
  {
    return m_vertices[segment.target_idx()];
  }
  inline Vertex& target_of_segment(const Segment& segment)
  {
    return m_vertices[segment.target_idx()];
  }
  inline const Vertex& target_of_segment(std::size_t segment_idx) const
  {
    return target_of_segment(m_segments[segment_idx]);
  }
  inline Vertex& target_of_segment(std::size_t segment_idx)
  {
    return target_of_segment(m_segments[segment_idx]);
  }


  // idx -> opposite Vertex
  inline const Vertex& opposite_vertex(std::size_t vertex_idx) const
  {
    const Segment& segment = segment_of_vertex(vertex_idx);

    CGAL_assertion(segment.source_idx() == vertex_idx
      || segment.target_idx() == vertex_idx);

    return (segment.source_idx() == vertex_idx ?
      m_vertices[segment.target_idx()] :
      m_vertices[segment.source_idx()]);
  }


  // Segment/idx -> Support_line
  inline const Support_line& support_line_of_segment(const Segment& segment) const
  {
    return m_support_lines[segment.support_line_idx()];
  }
  inline Support_line& support_line_of_segment(const Segment& segment)
  {
    return m_support_lines[segment.support_line_idx()];
  }
  inline const Support_line& support_line_of_segment(std::size_t segment_idx) const
  {
    return support_line_of_segment(m_segments[segment_idx]);
  }
  inline Support_line& support_line_of_segment(std::size_t segment_idx)
  {
    return support_line_of_segment(m_segments[segment_idx]);
  }

  // Vertex/idx -> Support_line
  inline const Support_line& support_line_of_vertex(const Vertex& vertex) const
  {
    return support_line_of_segment(vertex.segment_idx());
  }
  inline Support_line& support_line_of_vertex(const Vertex& vertex)
  {
    return support_line_of_segment(vertex.segment_idx());
  }
  inline const Support_line& support_line_of_vertex(std::size_t vertex_idx) const
  {
    return support_line_of_vertex(m_vertices[vertex_idx]);
  }
  inline Support_line& support_line_of_vertex(std::size_t vertex_idx)
  {
    return support_line_of_vertex(m_vertices[vertex_idx]);
  }

  // Vertex/idx -> Meta_vertex
  inline const Meta_vertex& meta_vertex_of_vertex(const Vertex& vertex) const
  {
    return m_meta_vertices[vertex.meta_vertex_idx()];
  }
  inline Meta_vertex& meta_vertex_of_vertex(const Vertex& vertex)
  {
    return m_meta_vertices[vertex.meta_vertex_idx()];
  }
  inline const Meta_vertex& meta_vertex_of_vertex(std::size_t vertex_idx) const
  {
    return meta_vertex_of_vertex(m_vertices[vertex_idx]);
  }
  inline Meta_vertex& meta_vertex_of_vertex(std::size_t vertex_idx)
  {
    return meta_vertex_of_vertex(m_vertices[vertex_idx]);
  }

  bool has_meta_vertex(const Vertex& vertex) const
  {
    return vertex.meta_vertex_idx() != std::size_t(-1);
  }
  bool has_meta_vertex(std::size_t vertex_idx) const
  {
    return has_meta_vertex(m_vertices[vertex_idx]);
  }

  FT position_of_meta_vertex_on_support_line(std::size_t meta_vertex_idx, std::size_t support_line_idx) const
  {
    return support_line(support_line_idx).to_1d(meta_vertex(meta_vertex_idx).point());
  }

  inline bool meta_vertex_exists(const Point_2& point) const
  {
    return m_meta_map.find(point) != m_meta_map.end();
  }

  void get_vertices_of_meta_vertex(std::size_t meta_vertex_idx,
    std::vector<std::size_t>& vertices_idx) const
  {
    const Meta_vertex& meta_vertex = m_meta_vertices[meta_vertex_idx];
    for (std::size_t support_line_idx : meta_vertex.support_lines_idx())
    {
      const Support_line& support_line = m_support_lines[support_line_idx];
      for (std::size_t segment_idx : support_line.segments_idx())
      {
        const Segment& segment = m_segments[segment_idx];
        for (std::size_t vertex_idx : { segment.source_idx(), segment.target_idx() })
          if (m_vertices[vertex_idx].meta_vertex_idx() == meta_vertex_idx)
            vertices_idx.push_back(vertex_idx);
      }
    }
  }

  inline CGAL::Bbox_2 bbox(const Vertex& vertex) const
  {
    return point_of_vertex(vertex).bbox();
  }
  inline CGAL::Bbox_2 bbox(const Support_line& support_line) const
  {
    return std::accumulate(support_line.segments_idx().begin(), support_line.segments_idx().end(),
      CGAL::Bbox_2(),
      [&](const CGAL::Bbox_2& bbox_2, const std::size_t& segment_idx) -> CGAL::Bbox_2
      {
        return bbox_2
          + bbox(source_of_segment(segment_idx))
          + bbox(target_of_segment(segment_idx));
      });
  }

  bool is_segment_frozen(std::size_t segment_idx) const
  {
    return (source_of_segment(segment_idx).is_frozen() && target_of_segment(segment_idx).is_frozen());
  }

  // idx -> Segment_2
  Segment_2 segment_2(std::size_t segment_idx) const
  {
    const Segment& segment = m_segments[segment_idx];
    const Support_line& support_line = m_support_lines[segment.support_line_idx()];
    const Vertex& source = m_vertices[segment.source_idx()];
    const Vertex& target = m_vertices[segment.target_idx()];

    return Segment_2(support_line.to_2d(source.point(m_current_time)), support_line.to_2d(target.point(m_current_time)));
  }

  bool is_bbox_support_line(std::size_t support_line_idx) const
  {
    return support_line_idx < 4;
  }

  bool is_bbox_segment(std::size_t segment_idx) const
  {
    return is_bbox_support_line(segment(segment_idx).support_line_idx());
  }

  bool is_bbox_meta_vertex(std::size_t meta_vertex_idx) const
  {
    for (std::size_t support_line_idx : meta_vertex(meta_vertex_idx).support_lines_idx())
      if (is_bbox_support_line(support_line_idx))
        return true;
    return false;
  }

  bool is_bbox_meta_edge(std::size_t source_idx, std::size_t target_idx) const
  {
    std::size_t common_line_idx = std::size_t(-1);

    for (std::size_t support_line_idx : meta_vertex(source_idx).support_lines_idx())
      if (m_meta_vertices[target_idx].support_lines_idx().find(support_line_idx)
        != m_meta_vertices[target_idx].support_lines_idx().end())
      {
        common_line_idx = support_line_idx;
        break;
      }

    CGAL_assertion(common_line_idx != std::size_t(-1));

    return is_bbox_support_line(common_line_idx);
  }

  bool is_meta_vertex_active(std::size_t meta_vertex_idx) const
  {
    for (std::size_t support_line_idx : meta_vertex(meta_vertex_idx).support_lines_idx())
      for (std::size_t segment_idx : support_line(support_line_idx).segments_idx())
        for (std::size_t vertex_idx : { segment(segment_idx).source_idx(), segment(segment_idx).target_idx() })
          if (vertex(vertex_idx).meta_vertex_idx() == meta_vertex_idx)
            return true;
    return false;
  }

  bool is_meta_vertex_intersection(std::size_t meta_vertex_idx) const
  {
    bool found_one = false;

    for (std::size_t support_line_idx : meta_vertex(meta_vertex_idx).support_lines_idx())
    {
      bool broken = false;
      for (std::size_t segment_idx : support_line(support_line_idx).segments_idx())
      {
        for (std::size_t vertex_idx : { segment(segment_idx).source_idx(), segment(segment_idx).target_idx() })
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

  bool is_meta_vertex_deadend_of_vertex(std::size_t meta_vertex_idx, std::size_t vertex_idx) const
  {
    return meta_vertex(meta_vertex_idx).is_deadend_of(segment_of_vertex(vertex_idx).support_line_idx());
  }

  void make_meta_vertex_deadend_of_vertex(std::size_t meta_vertex_idx, std::size_t vertex_idx)
  {
    meta_vertex(meta_vertex_idx).make_deadend_of(segment_of_vertex(vertex_idx).support_line_idx());
  }

  void make_meta_vertex_no_longer_deadend_of_vertex(std::size_t meta_vertex_idx, std::size_t vertex_idx)
  {
    meta_vertex(meta_vertex_idx).make_no_longer_deadend_of(segment_of_vertex(vertex_idx).support_line_idx());
  }

  std::size_t add_support_line(const Segment_2& segment)
  {
    m_support_lines.push_back(Support_line(segment));
    return std::size_t(m_support_lines.size() - 1);
  }

  Segment& add_segment(const Segment_2 segment, std::size_t input_idx = std::size_t(-1))
  {
    // Check if support line exists first
    Support_line new_support_line(segment);
    std::size_t support_line_idx = std::size_t(-1);
    for (std::size_t i = 0; i < number_of_support_lines(); ++i)
      if (new_support_line == support_line(i))
      {
        support_line_idx = i;
        break;
      }

    if (support_line_idx == std::size_t(-1))
    {
      support_line_idx = number_of_support_lines();
      m_support_lines.push_back(new_support_line);

      if (input_idx == std::size_t(-1))
      {
        m_support_lines.back().minimum() = m_support_lines.back().to_1d(segment.source());
        m_support_lines.back().maximum() = m_support_lines.back().to_1d(segment.target());
      }
      else
      {
        FT max_negative = -(std::numeric_limits<FT>::max)();
        FT min_positive = (std::numeric_limits<FT>::max)();

        for (std::size_t i = 0; i < 4; ++i)
        {
          Point_2 point;
          if (!KSP::internal::intersection(m_support_lines[i].line(), m_support_lines.back().line(), point))
            continue;

          FT position = m_support_lines.back().to_1d(point);
          if (position < 0 && position > max_negative)
            max_negative = position;
          if (position > 0 && position < min_positive)
            min_positive = position;
        }

        CGAL_assertion(max_negative != -(std::numeric_limits<FT>::max)()
          && min_positive != -(std::numeric_limits<FT>::min)());

        m_support_lines.back().minimum() = max_negative;
        m_support_lines.back().maximum() = min_positive;
      }
    }
    else
      support_line(support_line_idx).connected_components()++;


    std::size_t segment_idx = m_segments.size();
    m_segments.push_back(Segment(input_idx, support_line_idx));
    m_support_lines[support_line_idx].segments_idx().push_back(segment_idx);

    std::size_t source_idx = m_vertices.size();
    m_vertices.push_back(Vertex(m_support_lines[support_line_idx].to_1d(segment.source()),
      segment_idx));
    std::size_t target_idx = m_vertices.size();
    m_vertices.push_back(Vertex(m_support_lines[support_line_idx].to_1d(segment.target()),
      segment_idx));

    // Keep segment ordered
    if (m_vertices[source_idx].point(0) > m_vertices[target_idx].point(0))
      std::swap(source_idx, target_idx);

    m_segments[segment_idx].source_idx() = source_idx;
    m_segments[segment_idx].target_idx() = target_idx;
    return m_segments.back();
  }

  std::size_t add_meta_vertex(const Point_2& point,
    std::size_t support_line_idx_0,
    std::size_t support_line_idx_1 = std::size_t(-1))
  {
    // Avoid several points almost equal
    Point_2 p(1e-10 * std::floor(CGAL::to_double(point.x()) / 1e-10),
      1e-10 * std::floor(CGAL::to_double(point.y()) / 1e-10));

    typename std::map<Point_2, std::size_t>::iterator iter;
    bool inserted = false;
    std::tie(iter, inserted) = m_meta_map.insert(std::make_pair(p, number_of_meta_vertices()));
    if (inserted)
      m_meta_vertices.push_back(Meta_vertex(p));

    std::size_t meta_vertex_idx = iter->second;

    for (std::size_t support_line_idx : { support_line_idx_0, support_line_idx_1 })
    {
      if (support_line_idx != std::size_t(-1))
      {
        meta_vertex(meta_vertex_idx).support_lines_idx().insert(support_line_idx);

        if (std::find(support_line(support_line_idx).meta_vertices_idx().begin(),
          support_line(support_line_idx).meta_vertices_idx().end(),
          meta_vertex_idx) == support_line(support_line_idx).meta_vertices_idx().end())
          support_line(support_line_idx).meta_vertices_idx().push_back(meta_vertex_idx);
      }
    }

    // Special case = meta vertex is deadend of one line
    if (support_line_idx_1 == std::size_t(-1))
    {
      meta_vertex(meta_vertex_idx).make_deadend_of(support_line_idx_0);
    }

    return meta_vertex_idx;
  }

  void attach_vertex_to_meta_vertex(std::size_t vertex_idx, std::size_t meta_vertex_idx)
  {
    CGAL_assertion(!has_meta_vertex(vertex_idx));
    CGAL_assertion_msg(meta_vertex(meta_vertex_idx).support_lines_idx().find
    (segment_of_vertex(vertex_idx).support_line_idx())
      != meta_vertex(meta_vertex_idx).support_lines_idx().end(),
      "Trying to attach a vertex to a meta vertex not on its support line");
    vertex(vertex_idx).meta_vertex_idx() = meta_vertex_idx;
  }

  void cut_segment(std::size_t segment_idx, std::size_t meta_vertex_idx)
  {
    std::vector<std::size_t> vec(1, meta_vertex_idx);
    cut_segment(segment_idx, vec);
  }

  void cut_segment(std::size_t segment_idx, std::vector<std::size_t>& meta_vertices_idx)
  {
    Segment& segment = m_segments[segment_idx];
    std::size_t input_idx = segment.input_idx();
    std::size_t support_line_idx = segment.support_line_idx();
    // std::size_t source_idx = segment.source_idx();
    std::size_t target_idx = segment.target_idx();

    Support_line& support_line = support_line_of_segment(segment_idx);

    std::sort(meta_vertices_idx.begin(), meta_vertices_idx.end(),
      [&](const std::size_t& a,
        const std::size_t& b) -> bool
      {
        return (position_of_meta_vertex_on_support_line(a, support_line_idx)
          < position_of_meta_vertex_on_support_line(b, support_line_idx));
      });

    std::size_t nb_segments_before = m_segments.size();
    std::size_t nb_vertices_before = m_vertices.size();

    // Attach to existing endpoint
    std::size_t new_target_idx = m_vertices.size();
    m_vertices.push_back(Vertex(position_of_meta_vertex_on_support_line(meta_vertices_idx.front(),
      support_line_idx)));
    m_vertices[new_target_idx].segment_idx() = segment_idx;
    segment.target_idx() = new_target_idx;
    attach_vertex_to_meta_vertex(new_target_idx, meta_vertices_idx.front());

    // Create new segments
    for (std::size_t i = 0; i < meta_vertices_idx.size() - 1; ++i)
    {
      std::size_t sidx = m_segments.size();
      m_segments.push_back(Segment(input_idx, support_line_idx));
      support_line.segments_idx().push_back(sidx);

      std::size_t source_idx = m_vertices.size();
      m_vertices.push_back(Vertex(position_of_meta_vertex_on_support_line(meta_vertices_idx[i],
        support_line_idx)));
      m_vertices[source_idx].segment_idx() = sidx;
      m_segments[sidx].source_idx() = source_idx;
      attach_vertex_to_meta_vertex(source_idx, meta_vertices_idx[i]);

      std::size_t target_idx = m_vertices.size();
      m_vertices.push_back(Vertex(position_of_meta_vertex_on_support_line(meta_vertices_idx[i + 1],
        support_line_idx)));
      m_vertices[target_idx].segment_idx() = sidx;
      m_segments[sidx].target_idx() = target_idx;
      attach_vertex_to_meta_vertex(source_idx, meta_vertices_idx[i + 1]);
    }

    // Create final segment and attach to existing endpoint
    std::size_t sidx = m_segments.size();
    m_segments.push_back(Segment(input_idx, support_line_idx));
    support_line.segments_idx().push_back(sidx);

    std::size_t new_source_idx = m_vertices.size();
    m_vertices.push_back(Vertex(position_of_meta_vertex_on_support_line(meta_vertices_idx.back(),
      support_line_idx)));
    m_vertices[new_source_idx].segment_idx() = sidx;
    m_segments[sidx].source_idx() = new_source_idx;
    attach_vertex_to_meta_vertex(new_source_idx, meta_vertices_idx.back());

    m_vertices[target_idx].segment_idx() = sidx;
    m_segments[sidx].target_idx() = target_idx;
  }

  std::size_t propagate_segment(std::size_t vertex_idx)
  {
    // Create a new segment
    std::size_t segment_idx = m_segments.size();
    m_segments.push_back(Segment(segment_of_vertex(vertex_idx).input_idx(),
      segment_of_vertex(vertex_idx).support_line_idx()));
    support_line_of_vertex(vertex_idx).segments_idx().push_back(segment_idx);

    // Create new vertices
    std::size_t source_idx = m_vertices.size();
    m_vertices.push_back(Vertex(m_vertices[vertex_idx]));
    std::size_t target_idx = m_vertices.size();
    m_vertices.push_back(Vertex(m_vertices[vertex_idx]));

    // Connect segments and vertices
    m_segments[segment_idx].source_idx() = source_idx;
    m_segments[segment_idx].target_idx() = target_idx;
    m_vertices[source_idx].segment_idx() = segment_idx;
    m_vertices[target_idx].segment_idx() = segment_idx;

    CGAL_assertion(m_vertices[vertex_idx].direction() != 0);

    // Keep vertices ordered on the segment
    if (m_vertices[vertex_idx].direction() < 0)
      std::swap(source_idx, target_idx);

    // Freeze one end
    m_vertices[source_idx].freeze(m_current_time);

    // Release other end
    m_vertices[target_idx].meta_vertex_idx() = std::size_t(-1);

    return target_idx;
  }

  void update_positions(FT time)
  {
    m_current_time = time;
  }

};

} // namespace internal
} // namespace KSP_2
} // namespace CGAL


#endif // CGAL_KSP_2_DATA_STRUCTURE_H
