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

#include <CGAL/KSR/Event.h>
#include <CGAL/KSR/Event_queue.h>
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
  typedef typename Data::Intersection_line Intersection_line;
  typedef typename Data::Polygon Polygon;
  typedef typename Data::Segment Segment;
  typedef typename Data::Vertex Vertex;
  
  typedef typename Data::Meta_vertex Meta_vertex;
  
  typedef KSR::Event<Kernel> Event;
  typedef KSR::Event_queue<Kernel> Event_queue;


private:

  Data m_data;

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

    CGAL_KSR_CERR(1) << "Adding bbox as polygons" << std::endl;
    add_bbox_as_polygons (bbox, enlarge_bbox_ratio);

    CGAL_KSR_CERR(1) << "Adding input as polygons" << std::endl;

    KSR::size_t polygon_idx = 0;
    for (const typename PolygonRange::const_iterator::value_type& poly : polygons)
    {
      Polygon& polygon = m_data.add_polygon (get (polygon_map, poly), polygon_idx);
      ++ polygon_idx;
    }

    FT time_step = CGAL::approximate_sqrt(CGAL::squared_distance(Point_3 (bbox.xmin(), bbox.ymin(), bbox.zmin()),
                                                                 Point_3 (bbox.xmax(), bbox.ymax(), bbox.zmax())));
    
    time_step /= 50;
    
    KSR_3::dump (m_data, "init");
    CGAL_KSR_CERR(1) << "Making input polygons intersection free" << std::endl;
    CGAL_assertion(check_integrity(true));
    make_polygons_intersection_free();
    CGAL_assertion(check_integrity(true));
    
    KSR_3::dump (m_data, "iter_0");

    std::size_t iter = 0;
    FT min_time = 0;
    while (initialize_queue(min_time, min_time + time_step))
    {
      run();
      min_time += time_step;

      ++ iter;
//      KSR_3::dump (m_data, "iter_" + std::to_string(iter) + "_");

//      if (iter == 5)
        break;
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
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++ i)
    {
      const Support_plane& support_plane = m_data.support_plane(i);
      for (KSR::size_t p : support_plane.polygons_idx())
      {
        if (p == KSR::no_element())
        {
          if (verbose)
            std::cerr << "ERROR: Support_plane[" << i
                      << "] supports Polygon[-1]" << std::endl;
          return false;
        }
        const Polygon& polygon = m_data.polygon(p);
        if (polygon.support_plane_idx() != i)
        {
          if (verbose)
            std::cerr << "ERROR: Support_plane[" << i
                      << "] supports Polygon[" << p
                      << "] which claims to be supported by Support_plane[" << polygon.support_plane_idx()
                      << "]" << std::endl;
          return false;
        }
      }
      
      for (KSR::size_t l : support_plane.intersection_lines_idx())
      {
        if (l == KSR::no_element())
        {
          if (verbose)
            std::cerr << "ERROR: Support_plane[" << i
                      << "] supports Intersection_line[-1]" << std::endl;
          return false;
        }
        const Intersection_line& intersection_line = m_data.intersection_line(l);
        if (!intersection_line.has_support_plane(i))
        {
          if (verbose)
            std::cerr << "ERROR: Support_plane[" << i
                      << "] supports Intersection_line[" << l
                      << "] which claims it's not contained by it" << std::endl;
          return false;
        }
      }

      for (KSR::size_t mv : support_plane.meta_vertices_idx())
      {
        if (!m_data.meta_vertex(mv).has_support_plane(i))
        {
          if (verbose)
            std::cerr << "ERROR: Support_plane[" << i
                      << "] contains Meta_vertex[" << mv
                      << "] which claims it's not contained by it" << std::endl;
          return false;
        }
      }
    }

    for (KSR::size_t i = 0; i < m_data.number_of_segments(); ++ i)
    {
      const Segment& segment = m_data.segment(i);
      // TODO
    }
    
    for (KSR::size_t i = 0; i < m_data.number_of_polygons(); ++ i)
    {
      const Polygon& polygon = m_data.polygon(i);
      for (KSR::size_t v : polygon.vertices_idx())
      {
        if (v == KSR::no_element())
        {
          if (verbose)
            std::cerr << "ERROR: Polygon[" << i
                      << "] supports Vertex[-1]" << std::endl;
          return false;
        }
        const Vertex& vertex = m_data.vertex(v);
        if (vertex.polygon_idx() != i)
        {
          if (verbose)
            std::cerr << "ERROR: Polygon[" << i
                      << "] supports Vertex[" << v
                      << "] which claims to be supported by Polygon[" << vertex.polygon_idx()
                      << "]" << std::endl;
          return false;
        }
      }

      if (!m_data.support_plane(polygon.support_plane_idx()).has_polygon (i))
      {
        if (verbose)
          std::cerr << "ERROR: Polygon[" << i
                    << "] is supported by Support_plane[" << polygon.support_plane_idx()
                    << "] which claims it does not support it" << std::endl;
        return false;
      }
    }
    
    for (KSR::size_t i = 0; i < m_data.number_of_intersection_lines(); ++ i)
    {
      const Intersection_line& intersection_line = m_data.intersection_line(i);

      for (KSR::size_t p : intersection_line.support_planes_idx())
      {
        if (!m_data.support_plane(p).has_intersection_line(i))
        {
          if (verbose)
            std::cerr << "ERROR: Intersection_line[" << i
                      << "] crossed Support_plane[" << p
                      << "] which claims it does not contain it" << std::endl;
          return false;
        }
      }
      
      for (KSR::size_t s : intersection_line.segments_idx())
      {
        if (m_data.segment(s).intersection_line_idx() != i)
        {
          if (verbose)
            std::cerr << "ERROR: Intersection_line[" << i
                      << "] supports Segment[" << s
                      << "] which claims it's supported by Intersection_line["
                      << m_data.segment(s).intersection_line_idx() << std::endl;
          return false;
        }
      }
    }
    
    for (KSR::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
    {
      const Vertex& vertex = m_data.vertex(i);

      if (!m_data.polygon(vertex.polygon_idx()).has_vertex(i))
      {
        if (verbose)
          std::cerr << "ERROR: Vertex[" << i
                    << "] is supported by Polygon[" << vertex.polygon_idx()
                    << "] which claims it does not contain it" << std::endl;
        return false;
      }

      if (vertex.meta_vertex_idx() == KSR::no_element())
      {
        KSR::size_t meta_vertex_idx = m_data.meta_vertex_exists(m_data.point_of_vertex(vertex));
        if (meta_vertex_idx != KSR::no_element())
        {
          if (verbose)
            std::cerr << "ERROR: Vertex[" << i
                      << "] has no meta vertex but is located on Meta_vertex[" << meta_vertex_idx
                      << "]" << std::endl;
          return false;
        }
      }
    }
    
    for (KSR::size_t i = 0; i < m_data.number_of_meta_vertices(); ++ i)
    {
      const Meta_vertex& meta_vertex = m_data.meta_vertex(i);
      for (KSR::size_t p : meta_vertex.support_planes_idx())
      {
        if (!m_data.support_plane(p).has_meta_vertex(i))
        {
          if (verbose)
            std::cerr << "ERROR: Meta_vertex[" << i
                      << "] has Support_plane[" << p
                      << "] which claims it does not support it" << std::endl;
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
    m_data.update_positions(0.2); 
    for (KSR::size_t i = 0; i < m_data.number_of_vertices(); ++ i)
      *(vertices ++) = m_data.point_of_vertex(i);
    
    std::vector<std::size_t> facet;
    for (KSR::size_t i = 0; i < m_data.number_of_polygons(); ++ i)
    {
#define SKIP_BBOX_FACETS
#ifdef SKIP_BBOX_FACETS
      if (i < 6)
        continue;
#endif
      facet.clear();

      for (KSR::size_t vertex_idx : m_data.polygon(i).vertices_idx())
        facet.push_back (std::size_t(vertex_idx));

      *(facets ++) = facet;
    }
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

    // plane 0       vertex[0]       vertex[1]       vertex[2]       vertex[3]
    facet_points = { bbox_points[0], bbox_points[1], bbox_points[3], bbox_points[2] };
    m_data.add_polygon (facet_points);
    
    // plane 1       vertex[4]       vertex[5]       vertex[6]       vertex[7]
    facet_points = { bbox_points[4], bbox_points[5], bbox_points[7], bbox_points[6] };
    m_data.add_polygon (facet_points);

    // plane 2       vertex[8]       vertex[9]       vertex[10]      vertex[11]
    facet_points = { bbox_points[0], bbox_points[1], bbox_points[5], bbox_points[4] };
    m_data.add_polygon (facet_points);

    // plane 3       vertex[12]      vertex[13]      vertex[14]      vertex[15]
    facet_points = { bbox_points[2], bbox_points[3], bbox_points[7], bbox_points[6] };
    m_data.add_polygon (facet_points);
    
    // plane 4       vertex[16]      vertex[17]      vertex[18]      vertex[19]
    facet_points = { bbox_points[1], bbox_points[5], bbox_points[7], bbox_points[3] };
    m_data.add_polygon (facet_points);
    
    // plane 5       vertex[20]      vertex[21]      vertex[22]      vertex[23]
    facet_points = { bbox_points[0], bbox_points[4], bbox_points[6], bbox_points[2] };
    m_data.add_polygon (facet_points);

    //                                 Point           Planes   Vertices
    m_data.add_meta_vertex_and_attach (bbox_points[0], 0, 2, 5, 0, 8,  20);
    m_data.add_meta_vertex_and_attach (bbox_points[1], 0, 2, 4, 1, 9,  16);
    m_data.add_meta_vertex_and_attach (bbox_points[2], 0, 3, 5, 3, 12, 23);
    m_data.add_meta_vertex_and_attach (bbox_points[3], 0, 3, 4, 2, 13, 19);
    m_data.add_meta_vertex_and_attach (bbox_points[4], 1, 2, 5, 4, 11, 21);
    m_data.add_meta_vertex_and_attach (bbox_points[5], 1, 2, 4, 5, 10, 17);
    m_data.add_meta_vertex_and_attach (bbox_points[6], 1, 3, 5, 7, 15, 22);
    m_data.add_meta_vertex_and_attach (bbox_points[7], 1, 3, 4, 6, 14, 18);

    // Sort meta vertices
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++ i)
    {
      Point_3 centroid;
      KSR::size_t nb = 0;
      for (KSR::size_t meta_vertex_idx : m_data.support_plane(i).meta_vertices_idx())
      {
        centroid = CGAL::barycenter (centroid, nb, m_data.meta_vertex(meta_vertex_idx).point(), 1);
        ++ nb;
      }
      Point_2 centroid_2 = m_data.support_plane(i).to_2d(centroid);
      
      std::sort (m_data.support_plane(i).meta_vertices_idx().begin(),
                 m_data.support_plane(i).meta_vertices_idx().end(),
                 [&] (const KSR::size_t& a, const KSR::size_t& b) -> bool
                 {
                   return (Direction_2 (Segment_2 (centroid_2, m_data.support_plane(i).to_2d(m_data.meta_vertex(a).point())))
                           < Direction_2 (Segment_2 (centroid_2, m_data.support_plane(i).to_2d(m_data.meta_vertex(b).point()))));
                 });
    }

    
    //                            Line                                     Planes
    m_data.add_intersection_line (Line_3 (bbox_points[0], bbox_points[1]), 0, 2);
    m_data.intersection_line(0).meta_vertices_idx().push_back (0);
    m_data.intersection_line(0).meta_vertices_idx().push_back (1);
    m_data.add_segment (0, 0, 1);
    
    m_data.add_intersection_line (Line_3 (bbox_points[1], bbox_points[3]), 0, 4);
    m_data.intersection_line(1).meta_vertices_idx().push_back (1);
    m_data.intersection_line(1).meta_vertices_idx().push_back (3);
    m_data.add_segment (1, 1, 2);
    
    m_data.add_intersection_line (Line_3 (bbox_points[3], bbox_points[2]), 0, 3);
    m_data.intersection_line(2).meta_vertices_idx().push_back (3);
    m_data.intersection_line(2).meta_vertices_idx().push_back (2);
    m_data.add_segment (2, 2, 3);
    
    m_data.add_intersection_line (Line_3 (bbox_points[2], bbox_points[0]), 0, 5);
    m_data.intersection_line(3).meta_vertices_idx().push_back (2);
    m_data.intersection_line(3).meta_vertices_idx().push_back (0);
    m_data.add_segment (3, 3, 0);
    
    m_data.add_intersection_line (Line_3 (bbox_points[4], bbox_points[5]), 1, 2);
    m_data.intersection_line(4).meta_vertices_idx().push_back (4);
    m_data.intersection_line(4).meta_vertices_idx().push_back (5);
    m_data.add_segment (4, 4, 5);
    
    m_data.add_intersection_line (Line_3 (bbox_points[5], bbox_points[7]), 1, 4);
    m_data.intersection_line(5).meta_vertices_idx().push_back (5);
    m_data.intersection_line(5).meta_vertices_idx().push_back (7);
    m_data.add_segment (5, 5, 6);
    
    m_data.add_intersection_line (Line_3 (bbox_points[7], bbox_points[6]), 1, 3);
    m_data.intersection_line(6).meta_vertices_idx().push_back (7);
    m_data.intersection_line(6).meta_vertices_idx().push_back (6);
    m_data.add_segment (6, 6, 7);
    
    m_data.add_intersection_line (Line_3 (bbox_points[6], bbox_points[4]), 1, 5);
    m_data.intersection_line(7).meta_vertices_idx().push_back (6);
    m_data.intersection_line(7).meta_vertices_idx().push_back (4);
    m_data.add_segment (7, 7, 4);
    
    m_data.add_intersection_line (Line_3 (bbox_points[1], bbox_points[5]), 2, 4);
    m_data.intersection_line(8).meta_vertices_idx().push_back (1);
    m_data.intersection_line(8).meta_vertices_idx().push_back (5);
    m_data.add_segment (8, 1, 5);
    
    m_data.add_intersection_line (Line_3 (bbox_points[4], bbox_points[0]), 2, 5);
    m_data.intersection_line(9).meta_vertices_idx().push_back (4);
    m_data.intersection_line(9).meta_vertices_idx().push_back (0);
    m_data.add_segment (9, 4, 0);
    
    m_data.add_intersection_line (Line_3 (bbox_points[3], bbox_points[7]), 3, 4);
    m_data.intersection_line(10).meta_vertices_idx().push_back (3);
    m_data.intersection_line(10).meta_vertices_idx().push_back (7);
    m_data.add_segment (10, 2, 6);
    
    m_data.add_intersection_line (Line_3 (bbox_points[6], bbox_points[2]), 3, 5);
    m_data.intersection_line(11).meta_vertices_idx().push_back (6);
    m_data.intersection_line(11).meta_vertices_idx().push_back (2);
    m_data.add_segment (11, 3, 7);
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

  void make_polygons_intersection_free()
  {
    KSR::vector<std::tuple<Line_3, KSR::size_t, KSR::size_t> > todo;
    KSR::size_t nb_inter = 0;

    KSR::vector<KSR::vector<Point_3> > plane_3;
    plane_3.reserve (m_data.number_of_support_planes());
    KSR::vector<Box_with_idx> boxes;
    boxes.reserve (m_data.number_of_support_planes());
    for (KSR::size_t i = 0; i < m_data.number_of_support_planes(); ++ i)
    {
      plane_3.push_back (m_data.points_of_support_plane(i));
      boxes.push_back (Box_with_idx (CGAL::bbox_3 (plane_3.back().begin(), plane_3.back().end()), i));
    }
    
    CGAL::box_self_intersection_d
      (boxes.begin() + 6, boxes.end(),
       [&](const Box_with_idx& a, const Box_with_idx& b) -> void
       {
         KSR::size_t support_plane_idx_a = a.idx;
         KSR::size_t support_plane_idx_b = b.idx;

         CGAL_assertion (support_plane_idx_a != support_plane_idx_b);
         CGAL_assertion (support_plane_idx_a > 5 && support_plane_idx_b > 5);
         
         Line_3 line;
         if (!KSR::intersection_3 (m_data.support_plane(support_plane_idx_a).plane(),
                                   m_data.support_plane(support_plane_idx_b).plane(),
                                   line))
           return;

         if (KSR::do_intersect (plane_3[support_plane_idx_a], plane_3[support_plane_idx_b]))
         {
           todo.push_back (std::make_tuple (line, support_plane_idx_a, support_plane_idx_b));
           ++ nb_inter;
         }
       });


    CGAL_KSR_CERR(2) << "* Found " << nb_inter << " intersection(s)" << std::endl;

    KSR::Idx_vector new_intersection_lines;

    for (const std::tuple<Line_3, KSR::size_t, KSR::size_t>& t : todo)
      new_intersection_lines.push_back (m_data.add_intersection_line (get<0>(t), get<1>(t), get<2>(t)));

    for (KSR::size_t intersection_line_idx : new_intersection_lines)
    {
      std::cerr << "Intersection_line[" << intersection_line_idx << "]" << std::endl;
      for (KSR::size_t support_plane_idx : m_data.intersection_line(intersection_line_idx).support_planes_idx())
      {
        CGAL_assertion (support_plane_idx > 5);
        KSR::Idx_vector polygons_idx = m_data.support_plane(support_plane_idx).polygons_idx();
        for (KSR::size_t polygon_idx : polygons_idx)
          if (m_data.do_intersect (polygon_idx, m_data.line_on_support_plane (intersection_line_idx, support_plane_idx)))
            m_data.cut_polygon (polygon_idx, intersection_line_idx);
      }
    }
  }

  bool initialize_queue(FT min_time, FT max_time)
  {
    CGAL_KSR_CERR(1) << "Initializing queue for events in [" << min_time << ";" << max_time << "]" << std::endl;

    m_data.update_positions(max_time);

    bool still_running = true;

    return still_running;
  }

  void run()
  {
  }

  void apply (const Event& ev)
  {
  }

  void redistribute_vertex_events (KSR::size_t old_vertex, KSR::size_t new_vertex)
  {
  }

  void get_meta_neighbors (KSR::vector<KSR::vector<KSR::size_t> >& neighbors) const
  {
  }

};



} // namespace CGAL


#endif // CGAL_KINETIC_SHAPE_RECONSTRUCTION_3_H
