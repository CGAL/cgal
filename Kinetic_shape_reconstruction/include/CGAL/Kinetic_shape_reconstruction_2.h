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

    CGAL_KSR_CERR << "Adding bbox as segments" << std::endl;
    add_bbox_as_segments (bbox, enlarge_bbox_ratio);

    CGAL_KSR_CERR << "Adding input as segments" << std::endl;
    // Add input as segments
    for (const typename SegmentRange::const_iterator::value_type& vt : segments)
    {
      Segment& segment = m_data.add_segment (get (segment_map, vt));
      initialize_vertices_directions (segment);
    }

    CGAL_KSR_CERR << "Making input segments intersection free" << std::endl;
    make_segments_intersection_free();

    CGAL_KSR_CERR << "Initializing priority queue" << std::endl;
    initialize_queue(k);

    CGAL_KSR_CERR << "Unstacking priority queue" << std::endl;
    run();
  }

  
  template <typename PointRange, typename PointMap, typename VectorMap>
  void reconstruct (const PointRange& points, PointMap point_map, VectorMap normal_map)
  {

  }

  bool check_integrity(bool verbose = false) const
  {
    if (verbose)
      std::cerr << "Checking support lines integrity" << std::endl;

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

    if (verbose)
      std::cerr << "Checking segments integrity" << std::endl;

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

    if (verbose)
      std::cerr << "Checking vertices integrity" << std::endl;
    
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

        std::ofstream test ("metavertex.xyz");
        test << m_data.point_of_vertex(i) << " 0" << std::endl;
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

    if (verbose)
      std::cerr << "Checking meta vertices integrity" << std::endl;
    
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

    CGAL_KSR_CERR << "Creating fg vertices" << std::endl;
    for (std::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
    {
      vertex_descriptor vd = add_vertex(mesh);
      put (vpm, vd, m_data.meta_vertex(i).point());
      vdesc.push_back (vd);
    }

    CGAL_KSR_CERR << "Creating fg edges/halfedges" << std::endl;
    std::map<std::pair<KSR::size_t, KSR::size_t>, halfedge_descriptor> hdesc;
    for (std::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      KSR::size_t source = m_data.source_of_segment(i).meta_vertex_idx();
      KSR::size_t target = m_data.target_of_segment(i).meta_vertex_idx();

      if (source == target) // something fishy here, to check
        continue;

      vertex_descriptor v0 = vdesc[source];
      vertex_descriptor v1 = vdesc[target];

      edge_descriptor ed = add_edge(mesh);
      halfedge_descriptor hd = halfedge(ed, mesh);
      set_target(hd, v1, mesh);
      halfedge_descriptor opp_hd = opposite(hd, mesh);
      set_target(opp_hd, v0, mesh);
      set_halfedge(v1, hd, mesh);
      set_halfedge(v0, opp_hd, mesh);

      hdesc.insert (std::make_pair (std::make_pair (source, target), hd));
      hdesc.insert (std::make_pair (std::make_pair (target, source), opp_hd));
    }
    
    CGAL_KSR_CERR << "Ordering halfedges" << std::endl;
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

    CGAL_KSR_CERR << "Creating faces" << std::endl;
    for (halfedge_descriptor hd : halfedges(mesh))
      set_face (hd, boost::graph_traits<MutableFaceGraph>::null_face(), mesh);

    std::unordered_set<halfedge_descriptor> visited;
    for (halfedge_descriptor hd : halfedges(mesh))
    {
      if (!visited.insert(hd).second)
        continue;

      // First check if it is border face
      halfedge_descriptor end = hd;
      std::size_t nb_border_vdesc = 0;
      do
      {
        if (target (hd,mesh) < 4)
          nb_border_vdesc ++;
        hd = next(hd, mesh);
      }
      while (hd != end);

      if (nb_border_vdesc > 3)
        continue;

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
    for (std::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
      m_data.add_meta_vertex(i);
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

    CGAL_KSR_CERR << "Found " << nb_inter << " intersection(s) at initialization" << std::endl;

    for (typename std::map<KSR::size_t, std::vector<Point_2> >::iterator
           iter = todo.begin(); iter != todo.end(); ++ iter)
      m_data.cut_segment (iter->first, iter->second);
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
    std::size_t iterations = 0;
    while (!m_data.queue_is_empty())
    {
      Event ev = m_data.queue_pop();

      CGAL_KSR_CERR << " * Applying " << ev << std::endl;

      FT current_time = ev.time();

      m_data.update_positions (current_time);

      if (stop_vertex_if_intersection(ev.vertex_idx(), ev.intersection_line_idx()))
      {
        CGAL_KSR_CERR << "  -> Intersection happened" << std::endl;

        m_data.vertex_of_event(ev).remaining_intersections() --;
        if (m_data.is_bbox_segment (ev.intersection_line_idx()))
          m_data.vertex_of_event(ev).remaining_intersections() = 0;
        
        CGAL_KSR_CERR << "  -> Remaining intersections = " << m_data.vertex_of_event(ev).remaining_intersections() << std::endl;
        
        // If there are still intersections to be made
        if (m_data.vertex_of_event(ev).remaining_intersections() != 0)
        {
          // Create a new segment
          Segment& segment = m_data.propagate_segment (ev.vertex_idx());

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
        if (source == point || target == point)
        {
          CGAL_KSR_CERR << "Warning: intersection at point exactly" << std::endl;
        }

        break;
      }
    }

    if (intersected_segment == KSR::no_element()) // No intersection happened
      return false;

    vertex.point() = intersecting.to_1d(point_inter);

    m_data.cut_segment(intersected_segment, point_inter);

    m_data.add_meta_vertex(point_inter, vertex_idx);

    return true;
  }

};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_2_H
