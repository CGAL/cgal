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

#include <unordered_set>

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
  typedef typename Kernel::Direction_2 Direction_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Vector_2 Vector_2;

  typedef KSR_2::Data_structure<Kernel> Data;
  typedef typename Data::Support_line Support_line;
  typedef typename Data::Segment Segment;
  typedef typename Data::Vertex Vertex;
  
  typedef typename Data::Meta_vertex Meta_vertex;
  
  typedef typename Data::Event Event;
  typedef typename Data::Event_queue Event_queue;

private:

  Data m_data;

public:

  Kinetic_shape_reconstruction_2()
  {

  }


  template <typename SegmentRange, typename SegmentMap>
  void partition (const SegmentRange& segments, SegmentMap segment_map,
                  unsigned int k = 2, FT enlarge_bbox_ratio = 1.1)
  {
    CGAL::Bbox_2 bbox;
    for (const auto& vt : segments)
    {
      const Segment_2& segment = get (segment_map, vt);
      bbox += segment.bbox();
    }

    CGAL_KSR_CERR_1 << "Adding bbox as segments" << std::endl;
    add_bbox_as_segments (bbox, enlarge_bbox_ratio);

    CGAL_KSR_CERR_1 << "Adding input as segments" << std::endl;
    // Add input as segments
    KSR::size_t segment_idx = 0;
    for (const typename SegmentRange::const_iterator::value_type& vt : segments)
    {
      Segment& segment = m_data.add_segment (get (segment_map, vt), segment_idx);
      initialize_vertices_directions (segment, k);
      ++ segment_idx;
    }

    CGAL_KSR_CERR_1 << "Making input segments intersection free" << std::endl;
    make_segments_intersection_free();

    FT time_step = CGAL::approximate_sqrt(CGAL::squared_distance(Point_2 (bbox.xmin(), bbox.ymin()),
                                                                 Point_2 (bbox.xmax(), bbox.ymax())));
    time_step /= 50;

    m_data.print();
    
    FT min_time = 0;
    while (initialize_queue(min_time, min_time + time_step))
    {
      run();
      min_time += time_step;
    }
  }

  
  template <typename PointRange, typename PointMap, typename VectorMap>
  void reconstruct (const PointRange& points, PointMap point_map, VectorMap normal_map)
  {

  }

  bool check_integrity(bool verbose = false) const
  {
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
    {
      const Support_line& support_line = m_data.support_line(i);
      for (KSR::size_t s : support_line.segments_idx())
      {
        if (s == KSR::no_element())
        {
          if (verbose)
            std::cerr << "ERROR: Support_line[" << i
                      << "] supports Segment[-1]" << std::endl;
          return false;
        }
        const Segment& segment = m_data.segment(s);
        if (segment.support_line_idx() != i)
        {
          if (verbose)
            std::cerr << "ERROR: Support_line[" << i
                      << "] supports Segment[" << s
                      << "] which claims to be supported by Support_line[" << segment.support_line_idx()
                      << "]" << std::endl;
          return false;
        }
      }
    }

    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      const Segment& segment = m_data.segment(i);

      if (segment.source_idx() == KSR::no_element())
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has source Vertex[-1]" << std::endl;
        return false;
      }
      if (segment.target_idx() == KSR::no_element())
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has source Vertex[-1]" << std::endl;
        return false;
      }
      if (segment.support_line_idx() == KSR::no_element())
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has support line Support_line[-1]" << std::endl;
        return false;
      }
      if (m_data.source_of_segment(segment).segment_idx() != i)
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has source Vertex[" << segment.source_idx()
                    << "] which claims to belong to Segment[" << m_data.source_of_segment(segment).segment_idx()
                    << "]" << std::endl;
        return false;
      }
      if (m_data.target_of_segment(segment).segment_idx() != i)
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has target Vertex[" << segment.target_idx()
                    << "] which claims to belong to Segment[" << m_data.target_of_segment(segment).segment_idx()
                    << "]" << std::endl;
        return false;
      }
      if (segment.source_idx() == segment.target_idx())
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has Vertex[" << segment.source_idx()
                    << "] acting both as source and target" << std::endl;
        return false;
      }
      if (m_data.source_of_segment(segment).meta_vertex_idx() == m_data.target_of_segment(segment).meta_vertex_idx())
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] joins Vertex[" << segment.source_idx()
                    << "] to Vertex[" << segment.target_idx()
                    << "] which have the same meta vertex Meta_vertex["
                    << m_data.source_of_segment(segment).meta_vertex_idx() << "]" << std::endl;
        return false;
      }
      if (m_data.source_of_segment(segment).point(m_data.current_time())
          == m_data.target_of_segment(segment).point(m_data.current_time()))
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] joins Vertex[" << segment.source_idx()
                    << "] to Vertex[" << segment.target_idx()
                    << "] which represent the same point" << std::endl;
        return false;
      }
      if (m_data.source_of_segment(segment).point(m_data.current_time())
          > m_data.target_of_segment(segment).point(m_data.current_time()))
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] joins Vertex[" << segment.source_idx()
                    << "] to Vertex[" << segment.target_idx()
                    << "] which are wrongly ordered" << std::endl;
        return false;
      }

      if (std::find(m_data.support_line_of_segment(segment).segments_idx().begin(),
                    m_data.support_line_of_segment(segment).segments_idx().end(),
                    i) == m_data.support_line_of_segment(segment).segments_idx().end())
      {
        if (verbose)
          std::cerr << "ERROR: Segment[" << i
                    << "] has support line Support_line[" << segment.support_line_idx()
                    << "] which claims it does not support it" << std::endl;
        return false;
      }
    }

    for (std::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
    {
      const Vertex& vertex = m_data.vertex(i);

      if (vertex.segment_idx() == KSR::no_element())
      {
        if (verbose)
          std::cerr << "ERROR: Vertex[" << i
                    << "] is on Segment[-1]" << std::endl;
        return false;
      }
      if (vertex.meta_vertex_idx() == KSR::no_element())
      {
        if (verbose)
          std::cerr << "ERROR: Vertex[" << i
                    << "] has meta vertex Meta_vertex[-1]" << std::endl;
        return false;
      }
      if (m_data.segment_of_vertex(vertex).source_idx() != i
          && m_data.segment_of_vertex(vertex).target_idx() != i)
      {
        if (verbose)
          std::cerr << "ERROR: Vertex[" << i
                    << "] is on Segment[" << vertex.segment_idx()
                    << "] but is neither source nor vertex of it" << std::endl;
        return false;
      }

      if (std::find(m_data.meta_vertex_of_vertex(vertex).vertices_idx().begin(),
                    m_data.meta_vertex_of_vertex(vertex).vertices_idx().end(),
                    i) == m_data.meta_vertex_of_vertex(vertex).vertices_idx().end())
      {
        if (verbose)
          std::cerr << "ERROR: Vertex[" << i
                    << "] has meta vertex Meta_vertex[" << vertex.meta_vertex_idx()
                    << "] which claims it does not contain it" << std::endl;
        return false;
      }
    }

    for (std::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
    {
      const Meta_vertex& meta_vertex = m_data.meta_vertex(i);

      for (KSR::size_t vi : meta_vertex.vertices_idx())
      {
        if (vi == KSR::no_element())
        {
          if (verbose)
            std::cerr << "ERROR: Meta_vertex[" << i
                      << "] contains Vertex[-1]" << std::endl;
          return false;
        }
        if (m_data.vertex(vi).meta_vertex_idx() != i)
        {
          if (verbose)
            std::cerr << "ERROR: Meta_vertex[" << i
                      << "] has vertex Vertex[" << vi
                      << "] which claims to be on by Meta_vertex[" << m_data.vertex(vi).meta_vertex_idx()
                      << "]" << std::endl;
          return false;

        }
      }
    }
    
    return true;
  }

  template <typename OutputIterator>
  OutputIterator output_partition_edges_to_segment_soup (OutputIterator output) const
  {
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
      if (m_data.segment(i).is_active())
        *(output ++) = m_data.segment_2(i);

    return output;
  }

  template <typename VertexOutputIterator, typename FacetOutputIterator>
  void output_partition_cells_to_polygon_soup (VertexOutputIterator vertices,
                                               FacetOutputIterator facets) const
  {
    for (std::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
      *(vertices ++) = m_data.meta_vertex(i).point();
  }

  template <typename MutableFaceGraph>
  bool output_partition_cells_to_face_graph (MutableFaceGraph& mesh) const
  {
    typedef typename boost::graph_traits<MutableFaceGraph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<MutableFaceGraph>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<MutableFaceGraph>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<MutableFaceGraph>::face_descriptor face_descriptor;

    CGAL_static_assertion((CGAL::graph_has_property<MutableFaceGraph, boost::vertex_point_t>::value));
    typedef typename property_map_selector<MutableFaceGraph, CGAL::vertex_point_t>::type VPMap;
    VPMap vpm = get_property_map(boost::vertex_point, mesh);

    std::vector<vertex_descriptor> vdesc;
    vdesc.reserve (m_data.number_of_meta_vertices());

    CGAL_KSR_CERR_1 << "Creating fg vertices" << std::endl;
    for (std::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
    {
      vertex_descriptor vd = add_vertex(mesh);
      put (vpm, vd, m_data.meta_vertex(i).point());
      vdesc.push_back (vd);
    }

    CGAL_KSR_CERR_1 << "Creating fg edges/halfedges" << std::endl;
    std::map<std::pair<KSR::size_t, KSR::size_t>, halfedge_descriptor> hdesc;
    std::set<halfedge_descriptor> is_border_halfedge;
    
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      KSR::size_t source = m_data.source_of_segment(i).meta_vertex_idx();
      KSR::size_t target = m_data.target_of_segment(i).meta_vertex_idx();

      CGAL_assertion (source != target);

      vertex_descriptor v0 = vdesc[source];
      vertex_descriptor v1 = vdesc[target];

      edge_descriptor ed = add_edge(mesh);
      
      halfedge_descriptor hd = halfedge(ed, mesh);
      set_target(hd, v1, mesh);
      halfedge_descriptor opp_hd = opposite(hd, mesh);
      set_target(opp_hd, v0, mesh);
      set_halfedge(v1, hd, mesh);
      set_halfedge(v0, opp_hd, mesh);

      if (m_data.is_bbox_segment(i))
      {
        is_border_halfedge.insert(hd);
        is_border_halfedge.insert(opp_hd);
      }
      
      hdesc.insert (std::make_pair (std::make_pair (source, target), hd));
      hdesc.insert (std::make_pair (std::make_pair (target, source), opp_hd));
    }

    CGAL_KSR_CERR_2 << "* Found " << is_border_halfedge.size() << " border halfedges" << std::endl;
    
    CGAL_KSR_CERR_1 << "Ordering halfedges" << std::endl;
    for (std::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
    {
      const Meta_vertex& meta_vertex = m_data.meta_vertex(i);

      std::vector<KSR::size_t> incident_meta_vertices;
      for (KSR::size_t vertex_idx : meta_vertex.vertices_idx())
      {
        const Vertex& opposite = m_data.opposite_vertex(vertex_idx);
        CGAL_assertion (opposite.meta_vertex_idx() != KSR::no_element());

        if (i == opposite.meta_vertex_idx())
          continue;
        
        incident_meta_vertices.push_back (opposite.meta_vertex_idx());
      }

      std::sort (incident_meta_vertices.begin(), incident_meta_vertices.end(),
                 [&](const KSR::size_t& a, const KSR::size_t& b) -> bool
                 {
                   return (Direction_2 (Segment_2 (meta_vertex.point(), m_data.meta_vertex(a).point()))
                           > Direction_2 (Segment_2 (meta_vertex.point(), m_data.meta_vertex(b).point())));
                 });

      for (std::size_t j = 0; j < incident_meta_vertices.size(); ++ j)
      {
        std::pair<KSR::size_t, KSR::size_t> key0
          = std::make_pair (incident_meta_vertices[j], i);
        std::pair<KSR::size_t, KSR::size_t> key1
          = std::make_pair (incident_meta_vertices[(j+1)%incident_meta_vertices.size()], i);

        CGAL_assertion (hdesc.find(key0) != hdesc.end());
        CGAL_assertion (hdesc.find(key1) != hdesc.end());
        
        halfedge_descriptor h0 = hdesc[key0];
        halfedge_descriptor h1 = hdesc[key1];
        set_next (h0, opposite(h1,mesh),mesh);
      }
    }

    CGAL_KSR_CERR_1 << "Creating faces" << std::endl;
    for (halfedge_descriptor hd : halfedges(mesh))
      set_face (hd, boost::graph_traits<MutableFaceGraph>::null_face(), mesh);

    std::unordered_set<halfedge_descriptor> visited;
    bool found_border_face = false;
    for (halfedge_descriptor halfedge : halfedges(mesh))
    {
      if (!visited.insert(halfedge).second)
        continue;

      // First check if it is border face
      halfedge_descriptor hd = halfedge;
      halfedge_descriptor end = hd;

      bool border = true;

      if (found_border_face)
        border = false;
      else
      {
        do
        {
          // Border face only has border halfedges, as soon as we find one
          // non-border edge, we're done
          if (is_border_halfedge.find(hd)
              == is_border_halfedge.end())
          {
            border = false;
            break;
          }
          hd = next(hd, mesh);
        }
        while (hd != end);

        hd = halfedge;
      }
      
      if (border)
      {
        CGAL_KSR_CERR_2 << "* Found border face" << std::endl;
        found_border_face = true;
        end = hd;
        do
        {
          visited.insert(hd);
          hd = next(hd, mesh);
        }
        while (hd != end);
        continue;
      }
      
      face_descriptor fd = add_face(mesh);
      set_halfedge(fd, hd, mesh);

      end = hd;
      do
      {
        set_face(hd, fd, mesh);
        visited.insert(hd);
        hd = next(hd, mesh);
      }
      while (hd != end);

    }

    return is_valid_face_graph(mesh);
  }

  

private:
  
  void add_bbox_as_segments (const CGAL::Bbox_2& bbox, FT ratio)
  {
    FT xmed = (bbox.xmin() + bbox.xmax()) / 2.;
    FT ymed = (bbox.ymin() + bbox.ymax()) / 2.;
    FT dx = (bbox.xmax() - bbox.xmin()) / 2.;
    FT dy = (bbox.ymax() - bbox.ymin()) / 2.;

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

    m_data.add_meta_vertex (bbox_points[0], 0);
    m_data.add_meta_vertex (bbox_points[1], 1);
    m_data.add_meta_vertex (bbox_points[1], 2);
    m_data.add_meta_vertex (bbox_points[3], 3);
    m_data.add_meta_vertex (bbox_points[3], 4);
    m_data.add_meta_vertex (bbox_points[2], 5);
    m_data.add_meta_vertex (bbox_points[2], 6);
    m_data.add_meta_vertex (bbox_points[0], 7);
  }

  void initialize_vertices_directions (Segment& segment, unsigned int k)
  {
    const Support_line& support_line = m_data.support_line_of_segment (segment);

    Vertex& source = m_data.source_of_segment (segment);
    Vertex& target = m_data.target_of_segment (segment);

    source.remaining_intersections() = k;
    target.remaining_intersections() = k;

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
    std::map<KSR::size_t, std::vector<Point_2> > todo;

    std::size_t nb_inter = 0;
    for (std::size_t i = 4; i < m_data.number_of_segments() - 1; ++ i)
    {
      Segment_2 si_2 = m_data.segment_2(i);

      for (std::size_t j = i+1; j < m_data.number_of_segments(); ++ j)
      {
        Segment_2 sj_2 = m_data.segment_2(j);

        if (!CGAL::do_overlap (si_2.bbox(), sj_2.bbox()))
          continue;

        Point_2 point;
        if (!KSR::intersection_2 (si_2, sj_2, point))
          continue;

        typename std::map<KSR::size_t, std::vector<Point_2> >::iterator
          iter = todo.insert (std::make_pair (KSR::size_t(i), std::vector<Point_2>())).first;
        iter->second.push_back (point);
        
        iter = todo.insert (std::make_pair (KSR::size_t(j), std::vector<Point_2>())).first;
        iter->second.push_back (point);
        ++ nb_inter;
      }
    }

    CGAL_KSR_CERR_2 << "* Found " << nb_inter << " intersection(s) at initialization" << std::endl;

    for (typename std::map<KSR::size_t, std::vector<Point_2> >::iterator
           iter = todo.begin(); iter != todo.end(); ++ iter)
      m_data.cut_segment (iter->first, iter->second);
  }

  bool initialize_queue(FT min_time, FT max_time)
  {
    CGAL_KSR_CERR_1 << "Initializing queue for events in [" << min_time << ";" << max_time << "]" << std::endl;

    m_data.update_positions(min_time);
    
    // First, handle degenerate cases where a collision occur along a
    // same Support_line
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
    {
      const Support_line& support_line = m_data.support_line(i);
      if (support_line.segments_idx().size() < 2)
        continue;

      std::vector<KSR::size_t> vertices_idx;
      vertices_idx.reserve (support_line.segments_idx().size() * 2);
      for (KSR::size_t segment_idx : support_line.segments_idx())
      {
        vertices_idx.push_back (m_data.segment(segment_idx).source_idx());
        vertices_idx.push_back (m_data.segment(segment_idx).target_idx());
      }

      std::sort (vertices_idx.begin(), vertices_idx.end(),
                 [&](const KSR::size_t& a, const KSR::size_t& b) -> bool
                 { return m_data.vertex(a).point(m_data.current_time())
                     < m_data.vertex(b).point(m_data.current_time()); });

      for (std::size_t j = 1; j < vertices_idx.size() - 2; ++ j)
      {
        const Vertex& a = m_data.vertex (vertices_idx[j]);
        const Vertex& b = m_data.vertex (vertices_idx[j+1]);
        if (a.segment_idx() == b.segment_idx())
          continue;
        if (a.is_frozen() && b.is_frozen())
          continue;

        if (a.direction() < 0 || b.direction() > 0)
          continue;

        FT time_to_collision = b.point(m_data.current_time()) - a.point(m_data.current_time());
        if (!a.is_frozen() && ! b.is_frozen())
          time_to_collision /= 2.;

        if (time_to_collision < (max_time-min_time))
        {
          Point_2 point_a = support_line.to_2d(a.point(min_time + time_to_collision));
          Point_2 point_b = support_line.to_2d(b.point(min_time + time_to_collision));
          Point_2 point = CGAL::midpoint (point_a, point_b);
          
          Event ev (Event::PARALLEL, vertices_idx[j], vertices_idx[j+1], point, min_time + time_to_collision);
          CGAL_assertion (is_valid(ev));
          CGAL_KSR_CERR_2 << "* Pushing " << ev << std::endl;
          m_data.queue().push(ev);
        }
      }
    }
    
    // Simulate change of position
    m_data.update_positions(max_time);

    bool still_running = false;

    // Precompute segments and bboxes
    std::vector<Segment_2> segments_2;
    segments_2.reserve (m_data.number_of_segments());
    std::vector<CGAL::Bbox_2> segment_bboxes;
    segment_bboxes.reserve (m_data.number_of_segments());
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      if (m_data.segment(i).is_active())
      {
        segments_2.push_back (m_data.segment_2(i));
        segment_bboxes.push_back (segments_2.back().bbox());
      }
      else
      {
        segments_2.push_back (Segment_2());
        segment_bboxes.push_back (CGAL::Bbox_2());
      }
    }
    std::vector<CGAL::Bbox_2> support_line_bboxes;
    support_line_bboxes.reserve (m_data.number_of_support_lines());
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
      support_line_bboxes.push_back
        (std::accumulate (m_data.support_line(i).segments_idx().begin(),
                          m_data.support_line(i).segments_idx().end(),
                          CGAL::Bbox_2(),
                          [&](const CGAL::Bbox_2& b, const KSR::size_t& segment_idx) -> CGAL::Bbox_2
                          {
                            return b + segment_bboxes[segment_idx];
                          }));

    std::set<std::pair<KSR::size_t, KSR::size_t> > done;
    
    for (std::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
    {
      const Vertex& vertex = m_data.vertex(i);
      if (vertex.is_frozen() || !vertex.is_active())
        continue;

      CGAL_assertion (!m_data.has_meta_vertex(vertex));
      
      still_running = true;

      Segment_2 si (m_data.support_line_of_vertex(vertex).to_2d(vertex.point(min_time)),
                    m_data.point_of_vertex(vertex));
      CGAL::Bbox_2 si_bbox = si.bbox();

      for (std::size_t j = 0; j < m_data.number_of_support_lines(); ++ j)
      {
        if (m_data.segment_of_vertex(vertex).support_line_idx() == j)
          continue;

        const Support_line& support_line = m_data.support_line(j);

        if (!CGAL::do_overlap(si_bbox, support_line_bboxes[j]))
          continue;

        for (KSR::size_t segment_idx : support_line.segments_idx())
        {
          if (!CGAL::do_overlap(si_bbox, segment_bboxes[segment_idx]))
            continue;
          
          Point_2 point;
          if (!KSR::intersection_2 (si, segments_2[segment_idx], point))
            continue;

          Support_line& sli = m_data.support_line_of_vertex(vertex);
//          FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (sli.to_2d(vertex.point(0)), point));
          FT dist = CGAL::abs (vertex.point(0) - sli.to_1d(point));
          FT time = dist / vertex.speed();

          if (time > min_time)
          {
            point = m_data.point_of_vertex(i, time);

            Event ev (Event::BORDER, i, m_data.segment(segment_idx).source_idx(), point, time);

            bool valid = false;
            if (is_valid(ev))
              valid = true;

            if (!valid)
            {
              ev = Event (Event::BORDER, i, m_data.segment(segment_idx).target_idx(), point, time);
              if (is_valid(ev))
                valid = true;
            }
            if (!valid)
            {
              ev = Event (Event::REGULAR, i, segment_idx, point, time);
              if (is_valid(ev))
                valid = true;
            }

            if (valid)
            {
              std::pair<KSR::size_t, KSR::size_t> seg_seg
                = std::make_pair (vertex.segment_idx(), segment_idx);
              if (seg_seg.first > seg_seg.second)
                std::swap (seg_seg.first, seg_seg.second);
              if (!done.insert(seg_seg).second)
                continue;

              CGAL_KSR_CERR_2 << "* Pushing " << ev << std::endl;
              m_data.queue().push(ev);
            }
          }
        }
      }
    }

    m_data.update_positions(min_time);

    return still_running;
  }

  void run()
  {
    CGAL_KSR_CERR_1 << "Unstacking queue" << std::endl;
    
    std::size_t iterations = 0;
    
    Event_queue& queue = m_data.queue();

    static int iter = 0;

    while (!queue.empty())
    {
      Event ev = queue.pop();

      FT current_time = ev.time();

      m_data.update_positions (current_time);

      CGAL_KSR_CERR_2 << "* Applying " << ev << std::endl;

      CGAL_assertion (is_valid(ev));

      if (ev.type() == Event::REGULAR)
        apply_regular_event (ev.vertex_idx(), ev.segment_idx(), ev.intersection());
      else if (ev.type() == Event::BORDER)
        apply_border_event (ev.vertex_idx(), ev.other_vertex_idx(), ev.intersection());
      else // ev.type() == Event::PARALLEL
        apply_parallel_event (ev.vertex_idx(), ev.other_vertex_idx(), ev.intersection());
      
      ++ iterations;

      // std::ofstream out ("dbg.polylines.txt");
      // out.precision(18);
      // output_partition_edges_to_segment_soup
      //   (boost::make_function_output_iterator
      //    ([&](const Segment_2& segment) -> void
      //     {
      //       out << "2 " << segment.source() << " 0 " << segment.target() << " 0" << std::endl;
      //     }));

      // std::cout<<"Press [Enter] to continue . . .";
      // std::cin.get();

      // if (iterations == 6)
      //   break;
    }
  }

  void apply_regular_event (KSR::size_t vertex_idx, KSR::size_t segment_idx, const Point_2& point)
  {
    CGAL_KSR_ASSERT_POINTS_ALMOST_EQUAL (m_data.point_of_vertex (vertex_idx), point);
    
    KSR::size_t new_cut_segment_idx = m_data.cut_segment(segment_idx, point);

    if (!m_data.has_meta_vertex(vertex_idx))
      m_data.add_meta_vertex(point, vertex_idx);

    if (m_data.is_bbox_segment (segment_idx))
      m_data.vertex(vertex_idx).remaining_intersections() = 0;

    redistribute_segment_events (segment_idx, new_cut_segment_idx);

    KSR::size_t new_vertex_idx;
    KSR::size_t new_segment_idx;

    std::tie (new_vertex_idx, new_segment_idx) = freeze_and_propagate_vertex (vertex_idx);
  
    // Transfer events to new moving vertex
    redistribute_vertex_events (vertex_idx, new_vertex_idx);
    redistribute_segment_events (m_data.vertex(vertex_idx).segment_idx(), new_segment_idx);

  }
  
  void apply_border_event (KSR::size_t vertex_idx, KSR::size_t other_vertex_idx, const Point_2& point)
  {
    m_data.connect_vertices(vertex_idx, other_vertex_idx, point);

    for (const std::pair<KSR::size_t, KSR::size_t>& idx
           : { std::make_pair(vertex_idx, other_vertex_idx),
               std::make_pair(other_vertex_idx, vertex_idx) } )
    {
      if (m_data.is_bbox_segment (m_data.vertex(idx.second).segment_idx()))
        m_data.vertex(idx.first).remaining_intersections() = 0;
      KSR::size_t new_vertex_idx;
      KSR::size_t new_segment_idx;

      std::tie (new_vertex_idx, new_segment_idx) = freeze_and_propagate_vertex (idx.first);

      // Transfer events to new moving vertex
      redistribute_vertex_events (idx.first, new_vertex_idx);
      redistribute_segment_events (m_data.vertex(idx.first).segment_idx(), new_segment_idx);
    }
  }

  void apply_parallel_event (KSR::size_t vertex_idx, KSR::size_t other_vertex_idx, const Point_2& point)
  {
    CGAL_assertion_msg (!(m_data.has_meta_vertex(vertex_idx) && m_data.has_meta_vertex(other_vertex_idx)),
                        [&]() -> std::string
                        {
                          return "Meta vertices are located at "
                            + KSR::to_string(m_data.meta_vertex_of_vertex(vertex_idx).point())
                            + " and " + KSR::to_string(m_data.meta_vertex_of_vertex(other_vertex_idx).point());
                        }().c_str());
    
    if (m_data.has_meta_vertex(vertex_idx))
      apply_border_parallel_event (vertex_idx, other_vertex_idx, point);
    else if (m_data.has_meta_vertex(other_vertex_idx))
      apply_border_parallel_event (other_vertex_idx, vertex_idx, point);
    else
    {
      KSR::size_t kept_segment_idx = m_data.vertex(vertex_idx).segment_idx();
      KSR::size_t removed_segment_idx = m_data.vertex(other_vertex_idx).segment_idx();
    
      m_data.merge_segments_of_vertices (vertex_idx, other_vertex_idx);

      transfer_segment_events (removed_segment_idx, kept_segment_idx);
    }

    std::vector<Event> dummy;
    m_data.queue().erase_vertex_events (vertex_idx, dummy);
    m_data.queue().erase_vertex_events (other_vertex_idx, dummy);
  }
  
  void apply_border_parallel_event (KSR::size_t fixed_vertex_idx, KSR::size_t other_vertex_idx, const Point_2& point)
  {
    CGAL_KSR_ASSERT_POINTS_ALMOST_EQUAL (m_data.meta_vertex (fixed_vertex_idx).point(), point);
    
    m_data.vertex(other_vertex_idx).remaining_intersections() = 0;
    m_data.vertex(other_vertex_idx).freeze(m_data.current_time());
    m_data.add_meta_vertex (m_data.meta_vertex(fixed_vertex_idx).point(), other_vertex_idx);
  }

  std::pair<KSR::size_t, KSR::size_t> freeze_and_propagate_vertex (KSR::size_t vertex_idx)
  {
    if (m_data.vertex(vertex_idx).remaining_intersections() != 0)
      m_data.vertex(vertex_idx).remaining_intersections() --;
        
    CGAL_KSR_CERR_3 << "** Remaining intersections = " << m_data.vertex(vertex_idx).remaining_intersections() << std::endl;

    // If there are still intersections to be made
    if (m_data.vertex(vertex_idx).remaining_intersections() != 0)
    {
      // Create a new segment
      KSR::size_t new_moving_vertex_idx = m_data.propagate_segment (vertex_idx);
      m_data.vertex(vertex_idx).freeze(m_data.current_time());
      return std::make_pair(new_moving_vertex_idx, m_data.vertex(new_moving_vertex_idx).segment_idx());
    }
    
    m_data.vertex(vertex_idx).freeze(m_data.current_time());
    return std::make_pair (KSR::no_element(), KSR::no_element());
  }

  void redistribute_vertex_events (KSR::size_t old_vertex, KSR::size_t new_vertex = KSR::no_element())
  {
    CGAL_KSR_CERR_3 << "** Redistribution events of vertex " << old_vertex << std::endl;
    Event_queue& queue = m_data.queue();

    std::vector<Event> events;
    queue.erase_vertex_events (old_vertex, events);

    for (Event& ev : events)
    {
      if (is_valid(ev))
      {
        CGAL_KSR_CERR_4 << "****   - Pushing " << ev << std::endl;
        queue.push (ev);
      }
      else if (new_vertex != KSR::no_element())
      {
        if (ev.vertex_idx() == old_vertex)
          ev.vertex_idx() = new_vertex;
        else
          ev.other_vertex_idx() = new_vertex;
        if (is_valid(ev))
        {
          CGAL_KSR_CERR_4 << "****   - Pushing " << ev << std::endl;
          queue.push (ev);
        }
      }
    }
  }

  void redistribute_segment_events (KSR::size_t old_segment, KSR::size_t new_segment = KSR::no_element())
  {
    CGAL_KSR_CERR_3 << "** Redistribution events of segment " << old_segment << std::endl;
    Event_queue& queue = m_data.queue();

    std::vector<Event> events;
    queue.erase_segment_events (old_segment, events);

    for (Event& ev : events)
    {
      if (is_valid(ev))
      {
        CGAL_KSR_CERR_4 << "****   - Pushing " << ev << std::endl;
        queue.push (ev);
      }
      else if (new_segment != KSR::no_element())
      {
        ev.segment_idx() = new_segment;
        if (is_valid(ev))
        {
          CGAL_KSR_CERR_4 << "****   - Pushing " << ev << std::endl;
          queue.push (ev);
        }
      }
    }
  }

  void transfer_segment_events (KSR::size_t old_segment, KSR::size_t new_segment)
  {
    CGAL_KSR_CERR_3 << "** Transfering segment events " << old_segment << std::endl;
    Event_queue& queue = m_data.queue();

    std::vector<Event> events;
    queue.erase_segment_events (old_segment, events);

    for (Event& ev : events)
    {
      ev.segment_idx() = new_segment;
      CGAL_assertion (is_valid(ev));
      {
        CGAL_KSR_CERR_4 << "****   - Pushing " << ev << std::endl;
        queue.push (ev);
      }
    }
  }

  void transfer_vertex_events_to_segment (KSR::size_t vertex_idx, KSR::size_t segment_idx)
  {
    // TODO ?
  }

  bool is_valid (const Event& ev)
  {
    if (ev.type() == Event::REGULAR)
    {
      Point_2 point = m_data.point_of_vertex(ev.vertex_idx(), ev.time());
      if (point != ev.intersection())
        return false;

      FT point_at_time = m_data.support_line_of_segment(ev.segment_idx()).to_1d (point);
      FT source_at_time = m_data.source_of_segment(ev.segment_idx()).point (ev.time());
      FT target_at_time = m_data.target_of_segment(ev.segment_idx()).point (ev.time());
      return (source_at_time < point_at_time && point_at_time < target_at_time);
    }
    
    if (ev.type() == Event::BORDER)
    {
      Point_2 point = m_data.point_of_vertex(ev.vertex_idx(), ev.time());
      if (point != ev.intersection())
        return false;
      
      FT point_at_time = m_data.support_line_of_vertex(ev.other_vertex_idx()).to_1d (point);
      FT border_at_time = m_data.vertex(ev.other_vertex_idx()).point (ev.time());
      return (point_at_time == border_at_time);
    }
    
    // else ev.type() == Event::PARALLEL
    Point_2 point_a = m_data.point_of_vertex (ev.vertex_idx(), ev.time());
    Point_2 point_b = m_data.point_of_vertex (ev.other_vertex_idx(), ev.time());
    Point_2 point = CGAL::midpoint (point_a, point_b);
    
    return (point == ev.intersection());
  }
};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_2_H
