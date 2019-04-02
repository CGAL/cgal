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

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_2_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_2_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR_2/Data_structure.h>

namespace CGAL
{

template <typename GeomTraits>
class Kinetic_shape_reconstruction_2
{
public:

  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Vector_2 Vector_2;

  typedef KSR_2::Data_structure<Kernel> Data;
  typedef typename Data::Support_line Support_line;
  typedef typename Data::Segment Segment;
  typedef typename Data::Vertex Vertex;
  typedef typename Data::Event Event;

private:

  Data m_data;

public:

  Kinetic_shape_reconstruction_2()
  {

  }


  template <typename SegmentRange, typename SegmentMap>
  void partition (const SegmentRange& segments, SegmentMap segment_map, unsigned int k = 2)
  {
    CGAL::Bbox_2 bbox;
    for (const auto& vt : segments)
    {
      const Segment_2& segment = get (segment_map, vt);
      bbox += segment.bbox();
    }

    add_bbox_as_segments (bbox);

    // Add input as segments
    for (const typename SegmentRange::const_iterator::value_type& vt : segments)
    {
      Segment& segment = m_data.add_segment (get (segment_map, vt));
      initialize_vertices_directions (segment);
    }

    make_segments_intersection_free();

    initialize_queue(k);

    run();
  }

  
  template <typename PointRange, typename PointMap, typename VectorMap>
  void reconstruct (const PointRange& points, PointMap point_map, VectorMap normal_map)
  {

  }

  template <typename OutputIterator>
  OutputIterator output_partition_edges_to_segment_soup (OutputIterator output) const
  {
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
      *(output ++) = m_data.segment_2(i);
    return output;
  }

  template <typename OutputIterator>
  OutputIterator output_partition_cells_to_surface_meshes (OutputIterator output) const
  {

  }

private:
  
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
    
    m_data.add_segment (Segment_2 (bbox_points[0], bbox_points[1]));
    m_data.add_segment (Segment_2 (bbox_points[1], bbox_points[3]));
    m_data.add_segment (Segment_2 (bbox_points[3], bbox_points[2]));
    m_data.add_segment (Segment_2 (bbox_points[2], bbox_points[0]));
  }

  void initialize_vertices_directions (Segment& segment)
  {
    const Support_line& support_line = m_data.support_line_of_segment (segment);

    Vertex& source = m_data.source_of_segment (segment);
    Vertex& target = m_data.target_of_segment (segment);

    Point_2 psource = m_data.point_of_vertex(source);
    Point_2 ptarget = m_data.point_of_vertex(target);

    if (Vector_2 (psource, ptarget) * support_line.line().to_vector() > 0.)
    {
      source.direction() = -1;
      target.direction() = 1;
    }
    else
    {
      source.direction() = 1;
      target.direction() = -1;
    }
  }

  void make_segments_intersection_free()
  {
    // TODO
  }

  void initialize_queue(unsigned int k)
  {
    // Loop over vertices and schedule events
    for (std::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
    {
      Vertex& vertex = m_data.vertex(i);
      if (vertex.is_frozen())
        continue;

      vertex.remaining_intersections() = k;

      Support_line& sli = m_data.support_line_of_vertex(vertex);

      Ray_2 ray = sli.to_ray (vertex);

      for (std::size_t j = 0; j < m_data.number_of_support_lines(); ++ j)
      {
        if (j == m_data.segment_of_vertex(vertex).support_line_idx())
          continue;

        Support_line& slj_line = m_data.support_line(j);
        Line_2 line = slj_line.line();

        Point_2 point;
        if (!KSR::intersection_2 (ray, line, point))
          continue;

        FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (sli.to_2d(vertex.point()), point));
        FT time = dist / vertex.speed();
        
        m_data.push_to_queue (Event (i, j, time));
      }
    }
  }

  void run()
  {
    FT latest_time = FT(0);

    std::size_t iterations = 0;
    while (!m_data.queue_is_empty())
    {
      Event ev = m_data.queue_pop();

      CGAL_KSR_CERR << " * Applying " << ev << std::endl;

      FT ellapsed_time = ev.time() - latest_time;
      latest_time = ev.time();

      m_data.advance_time (ellapsed_time);

      if (stop_vertex_if_intersection(ev.vertex_idx(), ev.intersection_line_idx()))
      {
        CGAL_KSR_CERR << "  -> Intersection happened" << std::endl;

        m_data.vertex_of_event(ev).remaining_intersections() --;
        if (m_data.is_bbox_segment (ev.intersection_line_idx()))
          m_data.vertex_of_event(ev).remaining_intersections() = 0;
        CGAL_KSR_CERR << " -> Remaining intersections = " << m_data.vertex_of_event(ev).remaining_intersections() << std::endl;
        
        // If there are still intersections to be made
        if (m_data.vertex_of_event(ev).remaining_intersections() != 0)
        {
          // Create a new segment
          Segment& segment = m_data.propagate_segment (m_data.vertex_of_event(ev));

          // Transfer events to new moving vertex
          m_data.transfer_events (ev.vertex_idx(), segment.target_idx());
        }
        else
          m_data.remove_events (ev.vertex_idx());

        m_data.vertex_of_event(ev).direction() = 0.;
      }
      else
      {
        CGAL_KSR_CERR << "  -> Nothing happened" << std::endl;
      }
      
      ++ iterations;
      // if (iterations == 6)
      //   break;
    }
  }
  
  bool stop_vertex_if_intersection (std::size_t vertex_idx, std::size_t line_idx)
  {
    Vertex& vertex = m_data.vertex(vertex_idx);
    const Support_line& intersecting = m_data.support_line_of_vertex(vertex);
    const Support_line& intersected = m_data.support_line(line_idx);

    Point_2 point_inter = KSR::intersection_2<Point_2> (intersecting.line(), intersected.line());

    KSR::size_t intersected_segment = KSR::no_element();
    
    for (KSR::size_t sg : intersected.segments_idx())
    {
      const Segment& segment = m_data.segment(sg);

      FT source = m_data.source_of_segment(segment).point();
      FT target = m_data.target_of_segment(segment).point();
      
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

    vertex.point() = intersecting.to_1d(point_inter);
    return true;
  }

  
};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_2_H
