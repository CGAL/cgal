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

#ifndef CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
#define CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/box_intersection_d.h>

#include <CGAL/KSR_3/Data_structure.h>
#include <CGAL/KSR_3/Polygon_splitter.h>

#include <CGAL/KSR_3/Event.h>
#include <CGAL/KSR_3/Event_queue.h>

#include <CGAL/KSR/debug.h>

#include <unordered_set>

namespace CGAL
{

template <typename GeomTraits>
class Kinetic_shape_reconstruction_3
{
public:

  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Direction_2 Direction_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Direction_3 Direction_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Vector_3 Vector_3;

  typedef KSR_3::Data_structure<Kernel> Data;
  typedef typename Data::Support_plane Support_plane;
  typedef typename Support_plane::Mesh Mesh;
  typedef typename Data::Intersection_graph Intersection_graph;
  typedef typename Intersection_graph::Vertex_descriptor Intersection_vertex;
  typedef typename Intersection_graph::Edge_descriptor Intersection_edge;

  typedef KSR_3::Event<Kernel> Event;
  typedef KSR_3::Event_queue<Kernel> Event_queue;


private:

  Data m_data;
  Event_queue m_queue;
  FT m_min_time;
  FT m_max_time;

public:

  Kinetic_shape_reconstruction_3()
  {

  }


  template <typename PolygonRange, typename PolygonMap>
  void partition (const PolygonRange& polygons, PolygonMap polygon_map,
                  unsigned int k = 2, FT enlarge_bbox_ratio = 1.1)
  {
    CGAL::Bbox_3 bbox;
    for (const auto& poly : polygons)
    {
      const std::vector<Point_3>& polygon = get (polygon_map, poly);
      bbox += CGAL::bbox_3 (polygon.begin(), polygon.end());
    }

    m_data.init (polygons.size());

    CGAL_KSR_CERR(1) << "Adding bbox as polygons" << std::endl;
    add_bbox_as_polygons (bbox, enlarge_bbox_ratio);

    CGAL_KSR_CERR(1) << "Adding input as polygons" << std::endl;

    KSR::size_t polygon_idx = 0;
    for (const typename PolygonRange::const_iterator::value_type& poly : polygons)
    {
      m_data.add_polygon (get (polygon_map, poly), polygon_idx);
      ++ polygon_idx;
    }

    FT time_step = CGAL::approximate_sqrt(CGAL::squared_distance(Point_3 (bbox.xmin(), bbox.ymin(), bbox.zmin()),
                                                                 Point_3 (bbox.xmax(), bbox.ymax(), bbox.zmax())));
    
    time_step /= 50;
    
    CGAL_KSR_CERR(1) << "Making input polygons intersection free" << std::endl;
    
    KSR_3::dump (m_data, "init");
    
    CGAL_assertion(check_integrity(true));
    make_polygons_intersection_free();
    CGAL_assertion(check_integrity(true));
    
    KSR_3::dump (m_data, "intersected");

    std::size_t iter = 0;
    m_min_time = 0;
    m_max_time = time_step;
    while (initialize_queue())
    {
      run();
      m_min_time = m_max_time;
      m_max_time += time_step;

      CGAL_assertion(check_integrity(true));
      ++ iter;
    }
    CGAL_assertion(check_integrity(true));
  }

  
  template <typename PointRange, typename PointMap, typename VectorMap>
  void reconstruct (const PointRange& points, PointMap point_map, VectorMap normal_map)
  {
    // TODO
  }

  bool check_integrity(bool verbose = false) const
  {
    // TODO
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++ i)
    {
      if (!m_data.mesh(i).is_valid())
      {
        if (verbose)
          std::cerr << "ERROR: Mesh " << i << "is invalid" << std::endl;
        return false;
      }

      for (const Intersection_edge& intersection_edge : m_data.support_plane(i).intersection_edges())
        if (m_data.intersected_planes(intersection_edge).find (i)
            == m_data.intersected_planes(intersection_edge).end())
        {
          if (verbose)
            std::cerr << "ERROR: Support_plane[" << i
                      << "] is intersected by Intersection_edge[" << intersection_edge
                      << "] which claims it does not intersect it" << std::endl;
          for (KSR::size_t spi : m_data.intersected_planes(intersection_edge))
            std::cerr << spi << " ";
          std::cerr << std::endl;
          return false;
        }
    }

    for (const Intersection_edge intersection_edge : m_data.intersection_edges())
    {
      for (KSR::size_t support_plane_idx : m_data.intersected_planes (intersection_edge))
      {
        if (m_data.support_plane(support_plane_idx).intersection_edges().find (intersection_edge)
            == m_data.support_plane(support_plane_idx).intersection_edges().end())
        {
          if (verbose)
            std::cerr << "ERROR: Intersection_edge[" << intersection_edge
                      << "] intersects Support_plane[" << support_plane_idx
                      << "] which claims it's not intersected by it" << std::endl;
          return false;
        }
      }
    }
    
    return true;
  }

  template <typename OutputIterator>
  OutputIterator output_partition_edges_to_segment_soup (OutputIterator output) const
  {
  }

  template <typename VertexOutputIterator, typename FacetOutputIterator>
  void output_partition_facets_to_polygon_soup (VertexOutputIterator vertices,
                                                FacetOutputIterator facets) //const
  {
  }


private:
  
  void add_bbox_as_polygons (const CGAL::Bbox_3& bbox, FT ratio)
  {
    FT xmed = (bbox.xmin() + bbox.xmax()) / 2.;
    FT ymed = (bbox.ymin() + bbox.ymax()) / 2.;
    FT zmed = (bbox.zmin() + bbox.zmax()) / 2.;
    FT dx = (bbox.xmax() - bbox.xmin()) / 2.;
    FT dy = (bbox.ymax() - bbox.ymin()) / 2.;
    FT dz = (bbox.zmax() - bbox.zmin()) / 2.;

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

    std::array<Point_3, 4> facet_points;

    facet_points = { bbox_points[0], bbox_points[1], bbox_points[3], bbox_points[2] };
    m_data.add_bbox_polygon (facet_points);
    
    facet_points = { bbox_points[4], bbox_points[5], bbox_points[7], bbox_points[6] };
    m_data.add_bbox_polygon (facet_points);

    facet_points = { bbox_points[0], bbox_points[1], bbox_points[5], bbox_points[4] };
    m_data.add_bbox_polygon (facet_points);

    facet_points = { bbox_points[2], bbox_points[3], bbox_points[7], bbox_points[6] };
    m_data.add_bbox_polygon (facet_points);
    
    facet_points = { bbox_points[1], bbox_points[5], bbox_points[7], bbox_points[3] };
    m_data.add_bbox_polygon (facet_points);
    
    facet_points = { bbox_points[0], bbox_points[4], bbox_points[6], bbox_points[2] };
    m_data.add_bbox_polygon (facet_points);

    CGAL_assertion (m_data.intersection_vertices().size() == 8);
    CGAL_assertion (m_data.intersection_edges().size() == 12);
  }

  struct Box_with_idx : public CGAL::Box_intersection_d::Box_d<FT,3>
  {
    typedef CGAL::Box_intersection_d::Box_d<FT,3> Base;
    KSR::size_t idx;

    Box_with_idx () { }

    Box_with_idx (const Bbox_3& bbox, KSR::size_t idx)
      : Base(bbox), idx(idx)
    { }
  };

  struct Intersection
  {
    Line_3 line;
    KSR::size_t support_plane_idx_0;
    Intersection_edge source_0;
    Intersection_edge target_0;
    KSR::size_t support_plane_idx_1;
    Intersection_edge source_1;
    Intersection_edge target_1;
  };

  void make_polygons_intersection_free()
  {
    // First, generate all transverse intersection lines
    typedef std::map<KSR::Idx_set, std::pair<Intersection_vertex, Intersection_vertex> > Map;
    Map map_p2vv;

    for (const Intersection_vertex& intersection_vertex : m_data.intersection_vertices())
    {
      KSR::Idx_set key = m_data.intersected_planes (intersection_vertex, false);
      if (key.size() < 2)
        continue;

      typename Map::iterator iter;
      bool inserted;
      std::tie (iter, inserted) = map_p2vv.insert (std::make_pair (key,
                                                                   std::make_pair (intersection_vertex,
                                                                                   Intersection_vertex())));
      if (!inserted)
        iter->second.second = intersection_vertex;
    }


    // Then, intersect these lines to find internal intersection vertices
    KSR::vector<std::pair<KSR::Idx_set, KSR::vector<Intersection_vertex> > > todo;
    for (typename Map::iterator it_a = map_p2vv.begin(); it_a != map_p2vv.end(); ++ it_a)
    {
      const KSR::Idx_set& set_a = it_a->first;

      todo.push_back (std::make_pair (set_a, KSR::vector<Intersection_vertex>()));

      KSR::vector<Intersection_vertex>& crossed_vertices = todo.back().second;
      crossed_vertices.push_back (it_a->second.first);

      std::set<KSR::Idx_set> done;
      
      for (typename Map::iterator it_b = map_p2vv.begin() ; it_b != map_p2vv.end(); ++ it_b)
      {
        const KSR::Idx_set& set_b = it_b->first;
        KSR::size_t common_plane_idx = KSR::no_element();
        std::set_intersection (set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                               boost::make_function_output_iterator
                               ([&](const KSR::size_t& idx) -> void { common_plane_idx = idx; }));

        if (common_plane_idx != KSR::no_element())
        {
          KSR::Idx_set union_set = set_a;
          union_set.insert (set_b.begin(), set_b.end());
          if (!done.insert (union_set).second)
            continue;
          
          Point_2 inter;
          if (!KSR::intersection_2 (m_data.support_plane(common_plane_idx).to_2d
                                    (Segment_3 (m_data.point_3 (it_a->second.first),
                                                m_data.point_3 (it_a->second.second))),
                                    m_data.support_plane(common_plane_idx).to_2d
                                    (Segment_3 (m_data.point_3 (it_b->second.first),
                                                m_data.point_3 (it_b->second.second))),
                                    inter))
            continue;

          crossed_vertices.push_back (m_data.add_vertex
                                      (m_data.support_plane(common_plane_idx).to_3d(inter), union_set));
        }
      }
      crossed_vertices.push_back (it_a->second.second);
    }

    for (auto& t : todo)
      m_data.add_intersection_edge (t.first, t.second);

    // Refine polygons
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++ i)
    {
      KSR_3::Polygon_splitter<GeomTraits> splitter (m_data);
      splitter.split_support_plane (i);
    }
  }

  bool initialize_queue()
  {
#if 0
    CGAL_KSR_CERR(1) << "Initializing queue for events in [" << m_min_time << ";" << m_max_time << "]" << std::endl;

    m_data.update_positions(m_max_time);

    bool still_running = false;

    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++ i)
    {
      const Support_plane& support_plane = m_data.support_plane(i);
      
      // Precompute segments and bboxes
      KSR::vector<Segment_2> segments_2;
      segments_2.reserve (support_plane.intersection_lines_idx().size());
      KSR::vector<CGAL::Bbox_2> segment_bboxes;
      segment_bboxes.reserve (support_plane.intersection_lines_idx().size());
      for (KSR::size_t intersection_line_idx : support_plane.intersection_lines_idx())
      {
        segments_2.push_back (m_data.segment_of_intersection_line_on_support_plane (intersection_line_idx, i));
        segment_bboxes.push_back (segments_2.back().bbox());
      }

      for (KSR::size_t polygon_idx : support_plane.polygons_idx())
      {
        const Polygon& polygon = m_data.polygon (polygon_idx);
        for (KSR::size_t v = 0; v < polygon.vertices_idx().size(); ++ v)
        {
          KSR::size_t vertex_idx = polygon.vertices_idx()[v];
          const Vertex& vertex = m_data.vertex (vertex_idx);
          if (vertex.is_frozen())
            continue;

          still_running = true;

          Segment_2 si (vertex.point (m_min_time), vertex.point (m_max_time));
          CGAL::Bbox_2 si_bbox = si.bbox();

          for (std::size_t j = 0; j < segments_2.size(); ++ j)
          {
            KSR::size_t intersection_line_idx = support_plane.intersection_lines_idx()[j];

            if (m_data.intersection_line_idx_of_vertex(vertex_idx) == intersection_line_idx)
              continue;

            if (!CGAL::do_overlap (si_bbox, segment_bboxes[j]))
              continue;

            Point_2 point;
            if (!KSR::intersection_2 (si, segments_2[j], point))
              continue;

            FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (vertex.point (m_min_time), point));
            FT time = dist / vertex.speed();

            m_queue.push (Event (vertex_idx, intersection_line_idx, m_min_time + time));
          }
        }
      }
    }

    m_data.update_positions(m_min_time);

    return still_running;
#endif
    return false;
  }

  void run()
  {
#if 0
    CGAL_KSR_CERR(1) << "Unstacking queue" << std::endl;
    
    KSR::size_t iterations = 0;
    
    static int iter = 0;

    while (!m_queue.empty())
    {
      Event ev = m_queue.pop();

      FT current_time = ev.time();

      m_data.update_positions (current_time);

      CGAL_KSR_CERR(2) << "* Applying " << iter << ": " << ev << std::endl;

      if (iter < 10)
      {
        dump (m_data, "iter_0" + std::to_string(iter));
        dump_event (m_data, ev, "iter_0" + std::to_string(iter));
      }
      else
      {
        dump (m_data, "iter_" + std::to_string(iter));
        dump_event (m_data, ev, "iter_" + std::to_string(iter));
      }
      ++ iter;
      if (iter == 30)
        exit(0);
      
      apply(ev);
              
      ++ iterations;
    }
#endif
  }

  void apply (const Event& ev)
  {
#if 0
    const Vertex& vertex = m_data.vertex (ev.vertex_idx());
    const Intersection_line& intersection_line = m_data.intersection_line (ev.intersection_line_idx());
    
    bool is_vertex_along_line = (vertex.segment_idx() != KSR::no_element());
    
    bool is_intersection_occupied = false;
    bool is_segment_bbox = false;

    const Point_2& point = vertex.point(m_data.current_time());
    for (KSR::size_t segment_idx : intersection_line.segments_idx())
    {
      Point_2 psource = m_data.point_on_plane (m_data.segment (segment_idx).source_idx(),
                                               m_data.polygon_of_vertex(vertex).support_plane_idx());
      Point_2 ptarget = m_data.point_on_plane (m_data.segment (segment_idx).target_idx(),
                                               m_data.polygon_of_vertex(vertex).support_plane_idx());

      Vector_2 ref (psource, ptarget);
      Vector_2 vec (psource, point);

      if (ref * vec < 0)
        continue;
      if (vec * vec < ref * ref)
      {
        is_intersection_occupied = true;

        if (m_data.is_bbox_segment (segment_idx))
          is_segment_bbox = true;
        
        break;
      }
    }
    

    CGAL_KSR_CERR(3) << "** Vertex " << ev.vertex_idx()
                     << (is_vertex_along_line ? " (constrained on " + std::to_string(vertex.segment_idx()) + ")" : " (not constrained)")
                     << (is_intersection_occupied ? " reaching intersection " : " reaching free line ")
                     << ev.intersection_line_idx() << std::endl;


    KSR::size_t polygon_idx = m_data.vertex(ev.vertex_idx()).polygon_idx();
    
    KSR::Idx_vector positive_side (1, ev.vertex_idx()), negative_side;
    negative_side.reserve (m_data.polygon(polygon_idx).vertices_idx().size() - 1);

    KSR::size_t idx = 0;
    bool inside = false;
    while (negative_side.size() < m_data.polygon(polygon_idx).vertices_idx().size() - 1)
    {
      KSR::size_t current_vertex_idx = m_data.polygon(polygon_idx).vertices_idx()[idx];

      if (inside)
        negative_side.push_back (current_vertex_idx);
      else if (current_vertex_idx == ev.vertex_idx())
        inside = true;

      idx = (idx + 1) % m_data.polygon(polygon_idx).vertices_idx().size();
    }

    Line_2 line_2 = m_data.line_on_support_plane (ev.intersection_line_idx(), m_data.polygon(polygon_idx).support_plane_idx());
      
    Point_2 new_point_0, new_point_1;
    Vector_2 new_direction_0, new_direction_1;
    std::tie (new_point_0, new_direction_0, new_point_1, new_direction_1)
      = m_data.compute_constrained_points_along_line (line_2, positive_side, negative_side);
    
    // If one of the neighbor vertices is already on line, then we
    // have to transfer the vertex to the other polygon
    if (m_data.intersection_line_idx_of_vertex(negative_side.front()) == ev.intersection_line_idx()
        || m_data.intersection_line_idx_of_vertex(negative_side.back()) == ev.intersection_line_idx())
    {
      CGAL_KSR_CERR(3) << "** Transfering to the other side" << std::endl;

      KSR::size_t changed_vertex_idx_0, changed_vertex_idx_1;
      if (m_data.intersection_line_idx_of_vertex(negative_side.front()) == ev.intersection_line_idx())
        std::tie (changed_vertex_idx_0, changed_vertex_idx_1)
          = m_data.transfer_vertex (ev.vertex_idx(), negative_side.front(), new_point_1, new_direction_1);
      else
        std::tie (changed_vertex_idx_0, changed_vertex_idx_1)
          = m_data.transfer_vertex (ev.vertex_idx(), negative_side.back(), new_point_0, new_direction_0);

      for (KSR::size_t changed_vertex_idx : { changed_vertex_idx_0, changed_vertex_idx_1 })
        if (changed_vertex_idx != KSR::no_element())
          update_events (changed_vertex_idx);
    }
    else
    {
      CGAL_KSR_CERR(3) << "** Cropping" << std::endl;

      // Remember position/Direction of point
      Point_2 point = m_data.vertex(ev.vertex_idx()).point(m_data.current_time());
      Vector_2 direction = m_data.vertex(ev.vertex_idx()).direction();
      
      negative_side.push_back (positive_side.front());
      KSR::size_t segment_idx
        = m_data.crop_polygon (polygon_idx, ev.intersection_line_idx(),
                               negative_side,
                               new_point_0, new_direction_0, new_point_1, new_direction_1);
      
      if (!is_intersection_occupied)
      {
        CGAL_KSR_CERR(3) << "** Propagating" << std::endl;

        KSR::size_t new_vertex_idx
          = m_data.propagate_polygon (segment_idx, point, direction);

        transfer_events (ev.vertex_idx(), new_vertex_idx);
      }

      for (KSR::size_t vertex_idx : { m_data.segment(segment_idx).source_idx(),
                                      m_data.segment(segment_idx).target_idx(),
                                      m_data.segment(segment_idx).other_source_idx(),
                                      m_data.segment(segment_idx).other_target_idx() })
        if (vertex_idx != KSR::no_element())
          update_events (vertex_idx);
    }
#endif
  }

#if 0
  void transfer_events (KSR::size_t old_vertex, KSR::size_t new_vertex)
  {
    CGAL_KSR_CERR(3) << "** Transfering events of vertex " << old_vertex << " to " << new_vertex << std::endl;

    KSR::vector<Event> events;
    m_queue.erase_vertex_events (old_vertex, events);

    for (Event& ev : events)
    {
      ev.vertex_idx() = new_vertex;
      CGAL_KSR_CERR(4) << "****   - Pushing " << ev << std::endl;
      m_queue.push (ev);
    }
  }

  void update_events (KSR::size_t vertex_idx)
  {
    remove_events (vertex_idx);
    compute_new_events (vertex_idx);
  }

  void remove_events (KSR::size_t vertex_idx)
  {
    CGAL_KSR_CERR(3) << "** Removing events of vertex " << vertex_idx << std::endl;
    m_queue.erase_vertex_events (vertex_idx);
  }

  void compute_new_events (KSR::size_t vertex_idx)
  {
    const Vertex& vertex = m_data.vertex (vertex_idx);
    if (vertex.is_frozen())
      return;
    
    CGAL_KSR_CERR(3) << "** Computing new events of vertex " << vertex_idx << std::endl;

    FT current_time = m_data.current_time();
    
    m_data.update_positions(m_max_time);

    KSR::size_t support_plane_idx = m_data.polygon_of_vertex(vertex_idx).support_plane_idx();
    const Support_plane& support_plane = m_data.support_plane(support_plane_idx);
      
    // Precompute segments and bboxes
    KSR::vector<Segment_2> segments_2;
    segments_2.reserve (support_plane.intersection_lines_idx().size());
    KSR::vector<CGAL::Bbox_2> segment_bboxes;
    segment_bboxes.reserve (support_plane.intersection_lines_idx().size());
    for (KSR::size_t intersection_line_idx : support_plane.intersection_lines_idx())
    {
      segments_2.push_back (m_data.segment_of_intersection_line_on_support_plane (intersection_line_idx, support_plane_idx));
      segment_bboxes.push_back (segments_2.back().bbox());
    }
    
    Segment_2 si (vertex.point (current_time), vertex.point (m_max_time));
    CGAL::Bbox_2 si_bbox = si.bbox();

    for (std::size_t j = 0; j < segments_2.size(); ++ j)
    {
      KSR::size_t intersection_line_idx = support_plane.intersection_lines_idx()[j];

      if (m_data.intersection_line_idx_of_vertex(vertex_idx) == intersection_line_idx)
        continue;

      if (!CGAL::do_overlap (si_bbox, segment_bboxes[j]))
        continue;

      Point_2 point;
      if (!KSR::intersection_2 (si, segments_2[j], point))
        continue;

      FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (vertex.point (current_time), point));
      FT time = dist / vertex.speed();

      m_queue.push (Event (vertex_idx, intersection_line_idx, current_time + time));
    }

    m_data.update_positions(current_time);
  }

  void get_meta_neighbors (KSR::vector<KSR::vector<KSR::size_t> >& neighbors) const
  {
  }
#endif

};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
