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

#include <CGAL/box_intersection_d.h>

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

    FT time_step = CGAL::approximate_sqrt(CGAL::squared_distance(Point_2 (bbox.xmin(), bbox.ymin()),
                                                                 Point_2 (bbox.xmax(), bbox.ymax())));
    
    time_step /= 50;
    
    CGAL_KSR_CERR_1 << "Making input segments intersection free" << std::endl;
    make_segments_intersection_free();

    CGAL_assertion(check_integrity(true));


//    m_data.print();

    
    FT min_time = 0;
    while (initialize_queue(min_time, min_time + time_step))
    {
      run();
      min_time += time_step;
    }

    // Prepare output by sorting segments along support lines;
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
    {
      Support_line& support_line = m_data.support_line(i);
      std::sort (support_line.segments_idx().begin(), support_line.segments_idx().end(),
                 [&](const KSR::size_t& a, const KSR::size_t& b) -> bool
                 {
                   return (m_data.source_of_segment(a).point(m_data.current_time())
                           < m_data.source_of_segment(b).point(m_data.current_time()));
                 });
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

      for (KSR::size_t mv : support_line.meta_vertices_idx())
      {
        if (std::find(m_data.meta_vertex(mv).support_lines_idx().begin(),
                      m_data.meta_vertex(mv).support_lines_idx().end(),
                      i) == m_data.meta_vertex(mv).support_lines_idx().end())
        {
          if (verbose)
            std::cerr << "ERROR: Support_line[" << i
                      << "] contains Meta_vertex[" << mv
                      << "] which claims it's not contained by it" << std::endl;
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
      if (m_data.source_of_segment(segment).meta_vertex_idx() != KSR::no_element()
          && m_data.target_of_segment(segment).meta_vertex_idx() != KSR::no_element()
          && m_data.source_of_segment(segment).meta_vertex_idx() == m_data.target_of_segment(segment).meta_vertex_idx())
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
                    << "] which represent the same point "
                    << m_data.point_of_vertex(m_data.source_of_segment(segment)) << std::endl;
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
        std::cerr << m_data.source_of_segment(segment).point(m_data.current_time())
          << " " << m_data.target_of_segment(segment).point(m_data.current_time()) << std::endl;
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
      // if (vertex.meta_vertex_idx() == KSR::no_element())
      // {
      //   if (verbose)
      //     std::cerr << "ERROR: Vertex[" << i
      //               << "] has meta vertex Meta_vertex[-1]" << std::endl;
      //   return false;
      // }
      if (m_data.segment_of_vertex(vertex).source_idx() != i
          && m_data.segment_of_vertex(vertex).target_idx() != i)
      {
        if (verbose)
          std::cerr << "ERROR: Vertex[" << i
                    << "] is on Segment[" << vertex.segment_idx()
                    << "] but is neither source nor vertex of it" << std::endl;
        return false;
      }

    }

    for (std::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
    {
      const Meta_vertex& meta_vertex = m_data.meta_vertex(i);

      for (KSR::size_t sl : meta_vertex.support_lines_idx())
      {
        if (std::find(m_data.support_line(sl).meta_vertices_idx().begin(),
                      m_data.support_line(sl).meta_vertices_idx().end(),
                      i) == m_data.support_line(sl).meta_vertices_idx().end())
        {
          if (verbose)
            std::cerr << "ERROR: Meta_vertex[" << i
                      << "] contains Support_line[" << sl
                      << "] which claims it's not contained by it" << std::endl;
          return false;
        }
      }
    }
    
    return true;
  }

  template <typename OutputIterator>
  OutputIterator output_partition_edges_to_segment_soup (OutputIterator output) const
  {
    std::vector<std::pair<KSR::size_t, KSR::size_t> > meta_edges;
    get_meta_edges (meta_edges);

    for (std::pair<KSR::size_t, KSR::size_t> meta_edge : meta_edges)
      *(output ++) = Segment_2 (m_data.meta_vertex(meta_edge.first).point(),
                                m_data.meta_vertex(meta_edge.second).point());
    
    return output;
  }

  template <typename OutputIterator>
  OutputIterator output_raw_partition_edges_to_segment_soup (OutputIterator output) const
  {
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
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

    std::vector<std::pair<KSR::size_t, KSR::size_t> > meta_edges;
    get_meta_edges (meta_edges);

    CGAL_KSR_CERR_1 << "Creating vertices and edges" << std::endl;
    std::map<KSR::size_t, vertex_descriptor> map_v2v;
    std::map<std::pair<KSR::size_t, KSR::size_t>, halfedge_descriptor> hdesc;
    std::set<halfedge_descriptor> is_border_halfedge;
    for (std::pair<KSR::size_t, KSR::size_t>& meta_edge : meta_edges)
    {
      KSR::size_t source = meta_edge.first;
      KSR::size_t target = meta_edge.second;

      typename std::map<KSR::size_t, vertex_descriptor>::iterator iter;
      bool inserted = false;
      std::tie (iter, inserted) = map_v2v.insert(std::make_pair(source, vertex_descriptor()));
      if (inserted)
      {
        iter->second = add_vertex(mesh);
        put (vpm, iter->second, m_data.meta_vertex(source).point());
      }
      vertex_descriptor v0 = iter->second;
      std::tie (iter, inserted) = map_v2v.insert(std::make_pair(target, vertex_descriptor()));
      if (inserted)
      {
        iter->second = add_vertex(mesh);
        put (vpm, iter->second, m_data.meta_vertex(target).point());
      }
      vertex_descriptor v1 = iter->second;

      edge_descriptor ed = add_edge(mesh);
      
      halfedge_descriptor hd = halfedge(ed, mesh);
      set_target(hd, v1, mesh);
      halfedge_descriptor opp_hd = opposite(hd, mesh);
      set_target(opp_hd, v0, mesh);
      set_halfedge(v1, hd, mesh);
      set_halfedge(v0, opp_hd, mesh);

      if (m_data.is_bbox_meta_edge(source, target))
      {
        is_border_halfedge.insert(hd);
        is_border_halfedge.insert(opp_hd);
      }
      
      hdesc.insert (std::make_pair (std::make_pair (source, target), hd));
      hdesc.insert (std::make_pair (std::make_pair (target, source), opp_hd));
    }

    CGAL_KSR_CERR_2 << "* Found " << is_border_halfedge.size() << " border halfedges" << std::endl;
    
    CGAL_KSR_CERR_1 << "Ordering halfedges" << std::endl;
    for (const std::pair<KSR::size_t, vertex_descriptor> vertex : map_v2v)
    {
      const Meta_vertex& meta_vertex = m_data.meta_vertex(vertex.first);
      
      std::vector<KSR::size_t> incident_meta_vertices;
      for (const std::pair<KSR::size_t, KSR::size_t>& meta_edge : meta_edges)
      {
        if (meta_edge.first == vertex.first)
          incident_meta_vertices.push_back (meta_edge.second);
        else if (meta_edge.second == vertex.first)
          incident_meta_vertices.push_back (meta_edge.first);
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
          = std::make_pair (incident_meta_vertices[j], vertex.first);
        std::pair<KSR::size_t, KSR::size_t> key1
          = std::make_pair (incident_meta_vertices[(j+1)%incident_meta_vertices.size()], vertex.first);

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

    KSR::size_t v0 = m_data.add_meta_vertex (3, 0);
    KSR::size_t v1 = m_data.add_meta_vertex (0, 1);
    KSR::size_t v2 = m_data.add_meta_vertex (1, 2);
    KSR::size_t v3 = m_data.add_meta_vertex (2, 3);

    m_data.attach_vertex_to_meta_vertex (0, 0);
    m_data.attach_vertex_to_meta_vertex (1, 1);
    m_data.attach_vertex_to_meta_vertex (2, 1);
    m_data.attach_vertex_to_meta_vertex (3, 2);
    m_data.attach_vertex_to_meta_vertex (4, 2);
    m_data.attach_vertex_to_meta_vertex (5, 3);
    m_data.attach_vertex_to_meta_vertex (6, 3);
    m_data.attach_vertex_to_meta_vertex (7, 0);

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

  struct Box_with_idx : public CGAL::Box_intersection_d::Box_d<FT,2>
  {
    typedef CGAL::Box_intersection_d::Box_d<FT,2> Base;
    KSR::size_t idx;

    Box_with_idx (const Bbox_2& bbox, KSR::size_t idx)
      : Base(bbox), idx(idx)
    { }
  };

  void make_segments_intersection_free()
  {
    std::vector<std::tuple<Point_2, KSR::size_t, KSR::size_t> > todo;
    std::size_t nb_inter = 0;

    std::vector<Segment_2> segments_2;
    segments_2.reserve (m_data.number_of_segments());
    std::vector<Box_with_idx> boxes;
    boxes.reserve (m_data.number_of_segments());
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      segments_2.push_back (m_data.segment_2(i));
      boxes.push_back (Box_with_idx (segments_2.back().bbox(), i));
    }
    
    CGAL::box_self_intersection_d
      (boxes.begin() + 4, boxes.end(),
       [&](const Box_with_idx& a, const Box_with_idx& b) -> void
       {
         KSR::size_t segment_idx_a = a.idx;
         KSR::size_t segment_idx_b = b.idx;
         
         CGAL_assertion (segment_idx_a != segment_idx_b);
         
         CGAL_assertion (m_data.segment(segment_idx_a).support_line_idx()
                         != m_data.segment(segment_idx_b).support_line_idx());

         Point_2 point;
         if (!KSR::intersection_2 (segments_2[segment_idx_a], segments_2[segment_idx_b], point))
           return;

         todo.push_back (std::make_tuple (point, 
                                          m_data.segment(segment_idx_a).support_line_idx(),
                                          m_data.segment(segment_idx_b).support_line_idx()));
        
         ++ nb_inter;
       });


    CGAL_KSR_CERR_2 << "* Found " << nb_inter << " intersection(s) at initialization" << std::endl;

    std::vector<KSR::size_t> new_meta_vertices;
    
    for (const std::tuple<Point_2, KSR::size_t, KSR::size_t>& t : todo)
      new_meta_vertices.push_back (m_data.add_meta_vertex (get<0>(t), get<1>(t), get<2>(t)));
    
    for (KSR::size_t meta_vertex_idx : new_meta_vertices)
    {
      for (KSR::size_t support_line_idx : m_data.meta_vertex(meta_vertex_idx).support_lines_idx())
      {
        FT position = m_data.position_of_meta_vertex_on_support_line (meta_vertex_idx, support_line_idx);
        for (KSR::size_t segment_idx : m_data.support_line(support_line_idx).segments_idx())
        {
          if (m_data.source_of_segment(segment_idx).point(0) < position
              && position < m_data.target_of_segment(segment_idx).point(0))
          {
            m_data.cut_segment (segment_idx, meta_vertex_idx);
            break;
          }
        }
      }
    }
  }

  bool initialize_queue(FT min_time, FT max_time)
  {
    CGAL_KSR_CERR_1 << "Initializing queue for events in [" << min_time << ";" << max_time << "]" << std::endl;

    m_data.update_positions(max_time);

    bool still_running = false;

    // First, create all new meta vertices at line-line intersections
    // that happened between min_time and max_time
    std::vector<KSR::size_t> new_meta_vertices;

   // Precompute segments and bboxes
    std::vector<Segment_2> segments_2;
    segments_2.reserve (m_data.number_of_segments());
    std::vector<CGAL::Bbox_2> segment_bboxes;
    segment_bboxes.reserve (m_data.number_of_segments());
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      segments_2.push_back (m_data.segment_2(i));
      segment_bboxes.push_back (segments_2.back().bbox());
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
      
    
    for (std::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
    {
      const Vertex& vertex = m_data.vertex(i);
      if (vertex.is_frozen())
        continue;
      
      still_running = true;

      Segment_2 si (m_data.support_line_of_vertex(vertex).to_2d(vertex.point(min_time)),
                    m_data.point_of_vertex(vertex));
      CGAL::Bbox_2 si_bbox = si.bbox();

      for (std::size_t j = 0; j < m_data.number_of_support_lines(); ++ j)
      {
        if (m_data.segment_of_vertex(vertex).support_line_idx() == j)
          continue;

        if (m_data.are_support_lines_connected (m_data.segment_of_vertex(vertex).support_line_idx(),
                                                j))
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
          FT dist = CGAL::approximate_sqrt (CGAL::squared_distance (sli.to_2d(vertex.point(0)), point));
          FT time = dist / vertex.speed();

          if (time > min_time)
          {
            new_meta_vertices.push_back (m_data.add_meta_vertex
                                         (point, j,
                                          m_data.segment_of_vertex(vertex).support_line_idx()));
            break;
          }
        }
      }
    }

    // Make sure structure stays correct
    m_data.update_positions(min_time);
    for (KSR::size_t meta_vertex_idx : new_meta_vertices)
    {
      for (KSR::size_t support_line_idx : m_data.meta_vertex(meta_vertex_idx).support_lines_idx())
      {
        FT position = m_data.position_of_meta_vertex_on_support_line (meta_vertex_idx, support_line_idx);
        for (KSR::size_t segment_idx : m_data.support_line(support_line_idx).segments_idx())
        {
          if (m_data.source_of_segment(segment_idx).point(min_time) < position
              && position < m_data.target_of_segment(segment_idx).point(min_time)
              && !(m_data.source_of_segment(segment_idx).meta_vertex_idx() == meta_vertex_idx
                   || m_data.target_of_segment(segment_idx).meta_vertex_idx() == meta_vertex_idx))
          {
            m_data.cut_segment (segment_idx, meta_vertex_idx);
            break;
          }
        }
      }
    }

    // Second, create all new meta vertices at internal line
    // intersection between two colinear segments
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
    {
      Support_line& support_line = m_data.support_line(i);
      if (support_line.connected_components() < 2)
        continue;

      bool active_vertices = false;
      
      std::vector<KSR::size_t> vertices_idx;
      vertices_idx.reserve (support_line.segments_idx().size() * 2);
      for (KSR::size_t segment_idx : support_line.segments_idx())
      {
        for (KSR::size_t vertex_idx : { m_data.segment(segment_idx).source_idx(),
                                        m_data.segment(segment_idx).target_idx() })
        {
          vertices_idx.push_back (vertex_idx);
          if (!m_data.vertex(vertex_idx).is_frozen())
            active_vertices = true;
        }
      }

      if (!active_vertices)
      {
        support_line.connected_components() = 1;
        continue;
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
        if (a.is_frozen() || b.is_frozen())
          continue;

        if (a.direction() < 0 || b.direction() > 0)
          continue;

        FT time_to_collision = (b.point(m_data.current_time()) - a.point(m_data.current_time())) / 2.;

        if (time_to_collision < (max_time-min_time))
        {
          Point_2 point_a = support_line.to_2d(a.point(min_time + time_to_collision));
          Point_2 point_b = support_line.to_2d(b.point(min_time + time_to_collision));
          Point_2 point = CGAL::midpoint (point_a, point_b);

          KSR::size_t meta_vertex_idx =m_data.add_meta_vertex (point, i);
        }
      }
    }

    // Then compute events along the lines
    
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
    {
      Support_line& support_line = m_data.support_line(i);

      for (KSR::size_t segment_idx : support_line.segments_idx())
      {
        const Segment& segment = m_data.segment(segment_idx);

        for (KSR::size_t vertex_idx : { segment.source_idx() , segment.target_idx() })
        {
          const Vertex& vertex = m_data.vertex(vertex_idx);
          if (vertex.is_frozen())
            continue;

          FT beginning = vertex.point(min_time);
          FT end = vertex.point(max_time);

          auto last_element
            = std::partition (support_line.meta_vertices_idx().begin(),
                              support_line.meta_vertices_idx().end(),
                              [&](const KSR::size_t meta_vertex_idx) -> bool
                              {
                                FT position = m_data.position_of_meta_vertex_on_support_line
                                  (meta_vertex_idx, i);
                                return ((beginning < position & position <= end)
                                        || (end <= position & position < beginning));
                              });

          for (auto it = support_line.meta_vertices_idx().begin(); it != last_element; ++ it)
            m_data.queue().push (Event (vertex_idx, *it,
                                        min_time + CGAL::abs(beginning -
                                                             m_data.position_of_meta_vertex_on_support_line
                                                             (*it, i))));
        }

      }
    }

    if (CGAL_KSR_VERBOSE_LEVEL > 3)
      m_data.queue().print();

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

      apply(ev);
      
      ++ iterations;
    }
  }

  void apply (const Event& ev)
  {
    bool is_meta_vertex_active = m_data.is_meta_vertex_active (ev.meta_vertex_idx());

    CGAL_KSR_CERR_3 << "** Vertex " << ev.vertex_idx() << " is at position "
                    << m_data.vertex(ev.vertex_idx()).point(ev.time()) << std::endl;
    
    // First, attach vertex to meta vertex
    m_data.attach_vertex_to_meta_vertex (ev.vertex_idx(), ev.meta_vertex_idx());

    // Then, check if vertex should be propagated behind
    //  -> if bbox is reached, we don't propagate
    if (m_data.is_bbox_meta_vertex (ev.meta_vertex_idx()))
      m_data.vertex(ev.vertex_idx()).remaining_intersections() = 0;

    //  -> special case for parallel lines, if deadend is reached, we don't propagate
    if (m_data.is_meta_vertex_deadend_of_vertex (ev.meta_vertex_idx(), ev.vertex_idx()))
    {
      CGAL_KSR_CERR_3 << "*** Deadend reached" << std::endl;
      m_data.vertex(ev.vertex_idx()).remaining_intersections() = 0;
    }

    //  -> if the number of K intersections is reached, we don't propagate
    if (is_meta_vertex_active && m_data.vertex(ev.vertex_idx()).remaining_intersections() != 0)
      m_data.vertex(ev.vertex_idx()).remaining_intersections() --;
        
    CGAL_KSR_CERR_3 << "** Remaining intersections = " << m_data.vertex(ev.vertex_idx()).remaining_intersections() << std::endl;

    // If there are still intersections to be made, propagate
    KSR::size_t new_vertex_idx = KSR::no_element();
    if (m_data.vertex(ev.vertex_idx()).remaining_intersections() != 0)
      new_vertex_idx = m_data.propagate_segment (ev.vertex_idx());
    else
      m_data.make_meta_vertex_deadend_of_vertex (ev.meta_vertex_idx(), ev.vertex_idx());
    
    redistribute_vertex_events (ev.vertex_idx(), new_vertex_idx);
    
    m_data.vertex(ev.vertex_idx()).freeze(m_data.current_time());
  }

  void redistribute_vertex_events (KSR::size_t old_vertex, KSR::size_t new_vertex)
  {
    CGAL_KSR_CERR_3 << "** Redistribution events of vertex " << old_vertex << " to " << new_vertex << std::endl;
    Event_queue& queue = m_data.queue();

    std::vector<Event> events;
    queue.erase_vertex_events (old_vertex, events);

    if (new_vertex != KSR::no_element())
      for (Event& ev : events)
      {
        ev.vertex_idx() = new_vertex;
        CGAL_KSR_CERR_4 << "****   - Pushing " << ev << std::endl;
        queue.push (ev);
      }
    else
      for (Event& ev : events)
        if (m_data.is_meta_vertex_deadend_of_vertex (ev.meta_vertex_idx(), ev.vertex_idx()))
        {
          CGAL_KSR_CERR_3 << "*** Remove deadend" << std::endl;
          m_data.make_meta_vertex_no_longer_deadend_of_vertex (ev.meta_vertex_idx(), ev.vertex_idx());
        }
  }

  void get_meta_edges (std::vector<std::pair<KSR::size_t, KSR::size_t> >& meta_edges) const
  {
    for (std::size_t i = 0; i < m_data.number_of_support_lines(); ++ i)
    {
      const Support_line& support_line = m_data.support_line(i);

      CGAL_assertion (support_line.meta_vertices_idx().size() > 1);

      KSR::size_t beginning = KSR::no_element();
      KSR::size_t end = KSR::no_element();
      
      for (KSR::size_t segment_idx : support_line.segments_idx())
      {
        // New segment
        if (beginning == KSR::no_element())
        {
          beginning = m_data.source_of_segment(segment_idx).meta_vertex_idx();
          end = m_data.target_of_segment(segment_idx).meta_vertex_idx();
        }
        else
        {
          // New segment is directly connected and no other line
          // crossed the meta vertex: ignore meta vertex
          if (end == m_data.source_of_segment(segment_idx).meta_vertex_idx()
              && !m_data.is_meta_vertex_intersection (end))
            end = m_data.target_of_segment(segment_idx).meta_vertex_idx();
          // Otherwise, add a vertex and output the segment
          else
          {
            meta_edges.push_back (std::make_pair (beginning, end));
            beginning = m_data.source_of_segment(segment_idx).meta_vertex_idx();
            end = m_data.target_of_segment(segment_idx).meta_vertex_idx();
          }
        }
      }
      meta_edges.push_back (std::make_pair (beginning, end));
    }
  }

};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_2_H
