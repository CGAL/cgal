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

#ifndef CGAL_KSR_3_DATA_STRUCTURE_H
#define CGAL_KSR_3_DATA_STRUCTURE_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <boost/function_output_iterator.hpp>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/verbosity.h>
#include <CGAL/KSR/debug.h>

#include <CGAL/KSR_3/Support_plane.h>
#include <CGAL/KSR_3/Intersection_line.h>
#include <CGAL/KSR_3/Segment.h>
#include <CGAL/KSR_3/Polygon.h>
#include <CGAL/KSR_3/Vertex.h>

#include <CGAL/KSR_3/Meta_vertex.h>

#include <CGAL/centroid.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Data_structure
{
public:
  
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Direction_2 Direction_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Line_3 Line_3;

  typedef KSR_3::Support_plane<Kernel> Support_plane;
  typedef KSR_3::Intersection_line<Line_3> Intersection_line;
  typedef KSR_3::Segment Segment;
  typedef KSR_3::Polygon Polygon;
  typedef KSR_3::Vertex<Kernel> Vertex;
  
  typedef KSR_3::Meta_vertex<Point_3> Meta_vertex;

  typedef KSR::vector<Support_plane> Support_planes;
  typedef KSR::vector<Intersection_line> Intersection_lines;
  typedef KSR::vector<Segment> Segments;
  typedef KSR::vector<Polygon> Polygons;
  typedef KSR::vector<Vertex> Vertices;
  
  typedef KSR::vector<Meta_vertex> Meta_vertices;

private:

  // Main data structure
  Support_planes m_support_planes;
  Intersection_lines m_intersection_lines;
  Segments m_segments;
  Polygons m_polygons;
  Vertices m_vertices;

  Meta_vertices m_meta_vertices;

  // Helping data structures
  std::map<Point_3, KSR::size_t> m_meta_map;
  
  FT m_current_time;
  
public:

  Data_structure()
    : m_current_time(0)
  { }

  void print() const
  {
    // TODO
  }
  
  const FT& current_time() const { return m_current_time; }

  KSR::size_t number_of_vertices() const { return m_vertices.size(); }
  const Vertex& vertex (KSR::size_t idx) const { return m_vertices[idx]; }
  Vertex& vertex (KSR::size_t idx) { return m_vertices[idx]; }
  KSR::size_t add_vertex (const Vertex& v)
  {
    KSR::size_t out = number_of_vertices();
    m_vertices.push_back(v);
    return out;
  }

  KSR::size_t number_of_polygons() const { return m_polygons.size(); }
  const Polygon& polygon (KSR::size_t idx) const { return m_polygons[idx]; }
  Polygon& polygon (KSR::size_t idx) { return m_polygons[idx]; }
  KSR::size_t add_polygon (const Polygon& p)
  {
    KSR::size_t out = number_of_polygons();
    m_polygons.push_back(p);
    return out;
  }

  KSR::size_t number_of_support_planes() const { return m_support_planes.size(); }
  const Support_plane& support_plane (KSR::size_t idx) const { return m_support_planes[idx]; }
  Support_plane& support_plane (KSR::size_t idx) { return m_support_planes[idx]; }
  KSR::size_t add_support_plane (const Support_plane& new_support_plane)
  {
    KSR::size_t support_plane_idx = KSR::no_element();
    for (KSR::size_t i = 0; i < number_of_support_planes(); ++ i)
      if (new_support_plane == support_plane(i))
      {
        support_plane_idx = i;
        break;
      }

    if (support_plane_idx == KSR::no_element())
    {
      support_plane_idx = number_of_support_planes();
      m_support_planes.push_back (new_support_plane);
    }

    if (number_of_intersection_lines() >= 12) // Intersect planes with bbox... only after the 12 lines of bbox are here!
    {
      std::vector<std::pair<KSR::size_t, Point_3> > intersections;

      Point_3 centroid;
      for (KSR::size_t i = 0; i < 12; ++ i)
      {
        Point_3 point;
        if (!KSR::intersection_3 (support_plane(support_plane_idx).plane(),
                                  intersection_line(i).line(), point))
          continue;

        if (point_is_inside_bbox_section_of_intersection_line (point, i))
        {
          centroid = CGAL::barycenter (centroid, intersections.size(), point, 1);
          intersections.push_back (std::make_pair (i, point));
        }
      }

      Point_2 centroid_2 = support_plane(support_plane_idx).to_2d (centroid);
      std::sort (intersections.begin(), intersections.end(),
                 [&] (const std::pair<KSR::size_t, Point_3>& a,
                      const std::pair<KSR::size_t, Point_3>& b) -> bool
                 {
                   return (Direction_2 (Segment_2 (centroid_2, support_plane(support_plane_idx).to_2d (a.second)))
                           < Direction_2 (Segment_2 (centroid_2, support_plane(support_plane_idx).to_2d (b.second))));
                 });

      for (const std::pair<KSR::size_t, Point_3>& p : intersections)
        add_meta_vertex (p.second, support_plane_idx,
                         intersection_line(p.first).support_planes_idx()[0],
                         intersection_line(p.first).support_planes_idx()[1]);

      for (KSR::size_t i = 0; i < support_plane(support_plane_idx).meta_vertices_idx().size(); ++ i)
        add_intersection_line (support_plane_idx,
                               support_plane(support_plane_idx).meta_vertices_idx()[i],
                               support_plane(support_plane_idx).meta_vertices_idx()
                               [(i+1) % support_plane(support_plane_idx).meta_vertices_idx().size()]);
    }
      
    return support_plane_idx;
  }

  KSR::size_t number_of_intersection_lines() const { return m_intersection_lines.size(); }
  const Intersection_line& intersection_line (KSR::size_t idx) const { return m_intersection_lines[idx]; }
  Intersection_line& intersection_line (KSR::size_t idx) { return m_intersection_lines[idx]; }
  KSR::size_t add_intersection_line (const Intersection_line& l)
  {
    KSR::size_t out = number_of_intersection_lines();
    m_intersection_lines.push_back(l);
    return out;
  }


  KSR::size_t number_of_segments() const { return m_segments.size(); }
  const Segment& segment (KSR::size_t idx) const { return m_segments[idx]; }
  Segment& segment (KSR::size_t idx) { return m_segments[idx]; }
  KSR::size_t add_segment (const Segment& s)
  {
    KSR::size_t out = number_of_segments();
    m_segments.push_back(s);
    return out;
  }

  KSR::size_t number_of_meta_vertices() const { return m_meta_vertices.size(); }
  const Meta_vertex& meta_vertex (KSR::size_t idx) const { return m_meta_vertices[idx]; }
  Meta_vertex& meta_vertex (KSR::size_t idx) { return m_meta_vertices[idx]; }
  KSR::size_t add_meta_vertex (const Meta_vertex& v)
  {
    KSR::size_t out = number_of_meta_vertices();
    m_meta_vertices.push_back(v);
    return out;
  }

  std::string polygon_str (KSR::size_t polygon_idx) const
  {
    std::string out
      = "Polygon[" + std::to_string(polygon_idx)
      + " from " + (polygon(polygon_idx).input_idx() == KSR::no_element() ?
                    "bbox" : std::to_string(polygon(polygon_idx).input_idx()))
      + "](";

    for (KSR::size_t vertex_idx : polygon(polygon_idx).vertices_idx())
      out += " v" + std::to_string(vertex_idx);
    out += " )";
    return out;
  }
  // Vertex/idx -> Point_3
  inline Point_3 point_of_vertex (const Vertex& vertex, FT time) const
  { return support_plane_of_vertex(vertex).to_3d(vertex.point(time)); }
  inline Point_3 point_of_vertex (KSR::size_t vertex_idx, FT time) const
  { return point_of_vertex (m_vertices[vertex_idx], time); }
  inline Point_3 point_of_vertex (const Vertex& vertex) const
  { return point_of_vertex (vertex, m_current_time); }
  inline Point_3 point_of_vertex (KSR::size_t vertex_idx) const
  { return point_of_vertex (vertex_idx, m_current_time); }

  // Vertex/idx -> Vector_3
  inline Vector_3 direction_of_vertex (const Vertex& vertex) const
  { return support_plane_of_vertex(vertex).to_3d(vertex.direction()); }
  inline Vector_3 direction_of_vertex (KSR::size_t vertex_idx) const
  { return direction_of_vertex (vertex(vertex_idx)); }

  // Vertex/ix -> Polygon
  inline const Polygon& polygon_of_vertex (const Vertex& vertex) const
  { return m_polygons[vertex.polygon_idx()]; }
  inline Polygon& polygon_of_vertex (const Vertex& vertex)
  { return m_polygons[vertex.polygon_idx()]; }
  inline const Polygon& polygon_of_vertex (KSR::size_t vertex_idx) const
  { return polygon_of_vertex(vertex(vertex_idx)); }
  inline Polygon& polygon_of_vertex (KSR::size_t vertex_idx)
  { return polygon_of_vertex(vertex(vertex_idx)); }

  // Polygon/idx -> KSR::vector<Point_3>
  inline KSR::vector<Point_3> points_of_support_plane (const Support_plane& support_plane,
                                                       KSR::size_t nb_max = KSR::no_element()) const
  {
    KSR::vector<Point_3> out;
    if (nb_max == KSR::no_element())
      nb_max = support_plane.meta_vertices_idx().size();
    out.reserve (nb_max);
    for (KSR::size_t i = 0; i < nb_max; ++ i)
      out.push_back (meta_vertex(support_plane.meta_vertices_idx()[i]).point());
    return out;
  }
  inline KSR::vector<Point_3> points_of_support_plane (KSR::size_t support_plane_idx, KSR::size_t nb_max = KSR::no_element()) const
  { return points_of_support_plane (support_plane(support_plane_idx), nb_max); }

  // Polygon/idx -> Support_plane
  inline const Support_plane& support_plane_of_polygon (const Polygon& polygon) const
  { return support_plane(polygon.support_plane_idx()); }
  inline Support_plane& support_plane_of_polygon (const Polygon& polygon)
  { return support_plane(polygon.support_plane_idx()); }
  inline const Support_plane& support_plane_of_polygon (KSR::size_t polygon_idx) const
  { return support_plane_of_polygon(polygon(polygon_idx)); }
  inline Support_plane& support_plane_of_polygon (KSR::size_t polygon_idx)
  { return support_plane_of_polygon(polygon(polygon_idx)); }
  
  // Vertex/idx -> Support_plane
  inline const Support_plane& support_plane_of_vertex (const Vertex& vertex) const
  { return support_plane_of_polygon(vertex.polygon_idx()); }
  inline Support_plane& support_plane_of_vertex (const Vertex& vertex)
  { return support_plane_of_polygon(vertex.polygon_idx()); }
  inline const Support_plane& support_plane_of_vertex (KSR::size_t vertex_idx) const
  { return support_plane_of_polygon(vertex(vertex_idx)); }
  inline Support_plane& support_plane_of_vertex (KSR::size_t vertex_idx)
  { return support_plane_of_polygon(vertex(vertex_idx)); }

  // Intersection_line/Support_plane -> Line_2
  inline const Line_2 line_on_support_plane (KSR::size_t intersection_line_idx, KSR::size_t support_plane_idx) const
  { return support_plane(support_plane_idx).to_2d (intersection_line(intersection_line_idx).line()); }

  inline bool is_bbox_polygon (KSR::size_t polygon_idx) const
  { return (polygon(polygon_idx).support_plane_idx() < 6); }
  
  inline bool point_is_inside_bbox_section_of_intersection_line
  (const Point_3& point, KSR::size_t intersection_line_idx) const
  {
    Vector_3 ref (meta_vertex(intersection_line(intersection_line_idx).meta_vertices_idx()[0]).point(),
                  meta_vertex(intersection_line(intersection_line_idx).meta_vertices_idx()[1]).point());
    Vector_3 position (meta_vertex(intersection_line(intersection_line_idx).meta_vertices_idx()[0]).point(),
                       point);

    if (ref * position < 0)
      return false;

    return (position * position) < (ref * ref);
  }
    
  bool has_meta_vertex (const Vertex& vertex) const
  { return vertex.meta_vertex_idx() != KSR::no_element(); }
  bool has_meta_vertex (KSR::size_t vertex_idx) const
  { return has_meta_vertex (m_vertices[vertex_idx]); }
  
  template <typename PointRange>
  Polygon& add_polygon (const PointRange& polygon, KSR::size_t input_idx = KSR::no_element())
  {
    KSR::size_t support_plane_idx = add_support_plane (Support_plane (polygon));

    KSR::size_t polygon_idx = add_polygon (Polygon (input_idx, support_plane_idx));
    support_plane(support_plane_idx).polygons_idx().push_back (polygon_idx);

    // Create ordered polygon

    std::vector<Point_2> points;
    points.reserve (polygon.size());
    for (const Point_3& p : polygon)
      points.push_back (support_plane(support_plane_idx).to_2d(p));
    
    Point_2 centroid = CGAL::centroid (points.begin(), points.end());

    std::sort (points.begin(), points.end(),
               [&](const Point_2& a, const Point_2& b) -> bool
               {
                 return (Direction_2 (Segment_2 (centroid, a))
                         < Direction_2 (Segment_2 (centroid, b)));
               });

    for (const Point_2& p : points)
    {
      KSR::size_t vertex_idx = add_vertex (Vertex (p, polygon_idx));
      m_polygons[polygon_idx].vertices_idx().push_back (vertex_idx);

      if (support_plane_idx > 5) // Bbox doesn't move
      {
        // Initialize direction from center
        m_vertices.back().direction() = KSR::normalize (Vector_2 (centroid, p));
      }
    }

    return m_polygons[polygon_idx];
  }

  KSR::size_t meta_vertex_exists (const Point_3& point) const
  {
    Point_3 p (CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.x()) / CGAL_KSR_SAME_POINT_TOLERANCE),
               CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.y()) / CGAL_KSR_SAME_POINT_TOLERANCE),
               CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.z()) / CGAL_KSR_SAME_POINT_TOLERANCE));
      
    typename std::map<Point_3, KSR::size_t>::const_iterator iter
      = m_meta_map.find(p);
    if (iter != m_meta_map.end())
      return iter->second;
    return KSR::no_element();
  }

  KSR::size_t add_meta_vertex (const Point_3& point,
                               KSR::size_t support_plane_idx_0,
                               KSR::size_t support_plane_idx_1 = KSR::no_element(),
                               KSR::size_t support_plane_idx_2 = KSR::no_element())
  {
    // Avoid several points almost equal
    Point_3 p (CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.x()) / CGAL_KSR_SAME_POINT_TOLERANCE),
               CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.y()) / CGAL_KSR_SAME_POINT_TOLERANCE),
               CGAL_KSR_SAME_POINT_TOLERANCE * std::floor(CGAL::to_double(point.z()) / CGAL_KSR_SAME_POINT_TOLERANCE));
      
    typename std::map<Point_3, KSR::size_t>::iterator iter;
    bool inserted = false;
    std::tie (iter, inserted) = m_meta_map.insert (std::make_pair (p, number_of_meta_vertices()));
    if (inserted)
      add_meta_vertex (Meta_vertex(p));

    KSR::size_t meta_vertex_idx = iter->second;

    CGAL_KSR_CERR(3) << "** Adding meta vertex " << meta_vertex_idx << " between "
                     << support_plane_idx_0 << ", " << support_plane_idx_1 << " and " << support_plane_idx_2
                    << " at point " << p << std::endl;

    for (KSR::size_t support_plane_idx : { support_plane_idx_0, support_plane_idx_1, support_plane_idx_2 })
    {
      if (support_plane_idx != KSR::no_element())
      {
        meta_vertex(meta_vertex_idx).support_planes_idx().insert (support_plane_idx);

        if (std::find(support_plane(support_plane_idx).meta_vertices_idx().begin(),
                      support_plane(support_plane_idx).meta_vertices_idx().end(),
                      meta_vertex_idx) == support_plane(support_plane_idx).meta_vertices_idx().end())
          support_plane(support_plane_idx).meta_vertices_idx().push_back (meta_vertex_idx);
      }
    }

    return meta_vertex_idx;
  }


  void attach_vertex_to_meta_vertex (KSR::size_t vertex_idx, KSR::size_t meta_vertex_idx)
  {
    CGAL_assertion (!has_meta_vertex(vertex_idx));
    CGAL_assertion_msg (meta_vertex(meta_vertex_idx).support_planes_idx().find
                        (polygon_of_vertex(vertex_idx).support_plane_idx())
                        != meta_vertex(meta_vertex_idx).support_planes_idx().end(),
                        "Trying to attach a vertex to a meta vertex not on its support plane");
    vertex(vertex_idx).meta_vertex_idx() = meta_vertex_idx;
  }

  void add_meta_vertex_and_attach (const Point_3& point,
                                   KSR::size_t support_plane_idx_0,
                                   KSR::size_t support_plane_idx_1,
                                   KSR::size_t support_plane_idx_2,
                                   KSR::size_t vertex_idx_0,
                                   KSR::size_t vertex_idx_1,
                                   KSR::size_t vertex_idx_2)
  {
    KSR::size_t meta_vertex_idx = add_meta_vertex
      (point, support_plane_idx_0, support_plane_idx_1, support_plane_idx_2);
    
    attach_vertex_to_meta_vertex (vertex_idx_0, meta_vertex_idx);
    attach_vertex_to_meta_vertex (vertex_idx_1, meta_vertex_idx);
    attach_vertex_to_meta_vertex (vertex_idx_2, meta_vertex_idx);
  }

  KSR::size_t add_intersection_line (KSR::size_t support_plane_idx,
                                     KSR::size_t meta_vertex_idx_0,
                                     KSR::size_t meta_vertex_idx_1)
  {
    KSR::size_t intersection_line_idx
      = add_intersection_line (Intersection_line (Line_3 (meta_vertex(meta_vertex_idx_0).point(),
                                                          meta_vertex(meta_vertex_idx_1).point())));
    
    KSR::size_t common_support_plane_idx = KSR::no_element();
    for (KSR::size_t support_plane_idx_0 : meta_vertex(meta_vertex_idx_0).support_planes_idx())
      for (KSR::size_t support_plane_idx_1 : meta_vertex(meta_vertex_idx_1).support_planes_idx())
        if (support_plane_idx_0 != support_plane_idx
            && support_plane_idx_0 == support_plane_idx_1)
        {
          common_support_plane_idx = support_plane_idx_0;
          break;
        }
    CGAL_assertion (common_support_plane_idx != KSR::no_element());
    CGAL_assertion (common_support_plane_idx < 6);

    intersection_line(intersection_line_idx).meta_vertices_idx().push_back (meta_vertex_idx_0);
    intersection_line(intersection_line_idx).meta_vertices_idx().push_back (meta_vertex_idx_1);

    intersection_line(intersection_line_idx).support_planes_idx().push_back (support_plane_idx);
    intersection_line(intersection_line_idx).support_planes_idx().push_back (common_support_plane_idx);
    support_plane(support_plane_idx).intersection_lines_idx().push_back (intersection_line_idx);
    support_plane(common_support_plane_idx).intersection_lines_idx().push_back (intersection_line_idx);

    KSR::Idx_vector polygons_idx = support_plane(common_support_plane_idx).polygons_idx();
    for (KSR::size_t polygon_idx : polygons_idx)
      if (do_intersect (polygon_idx, line_on_support_plane (intersection_line_idx, common_support_plane_idx)))
        cut_polygon (polygon_idx, intersection_line_idx);

    return intersection_line_idx;
  }

  KSR::size_t add_intersection_line (const Line_3& line, KSR::size_t support_plane_idx_0, KSR::size_t support_plane_idx_1)
  {
    KSR::size_t intersection_line_idx = add_intersection_line (Intersection_line (line));
    if (number_of_intersection_lines() > 12) // After adding bbox, compute intersection with it for all other lines
    {
      Point_3 ref = line.point();
      Vector_3 v = line.to_vector();

      FT pos_min = (std::numeric_limits<double>::max)();
      FT pos_max = -(std::numeric_limits<double>::max)();
      Point_3 source = ref;
      Point_3 target = ref;
      KSR::size_t plane_min = KSR::no_element();
      KSR::size_t plane_max = KSR::no_element();
      for (KSR::size_t i = 0; i < 6; ++ i)
      {
        Point_3 p;
        if (!KSR::intersection_3 (line, points_of_support_plane(i, 4), p))
          continue;

        FT pos = Vector_3 (ref, p) * v;
        if (pos < pos_min)
        {
          source = p;
          pos_min = pos;
          plane_min = i;
        }
        if (pos > pos_max)
        {
          target = p;
          pos_max = pos;
          plane_max = i;
        }
      }

      CGAL_assertion (plane_min != KSR::no_element() && plane_max != KSR:: no_element());

      KSR::size_t vsource = add_meta_vertex (source, plane_min);
      KSR::size_t vtarget = add_meta_vertex (target, plane_max);

      intersection_line(intersection_line_idx).meta_vertices_idx().push_back (vsource);
      intersection_line(intersection_line_idx).meta_vertices_idx().push_back (vtarget);

      for (KSR::size_t support_plane_idx : { support_plane_idx_0, support_plane_idx_1 })
      {
        Line_2 line_2 = support_plane(support_plane_idx).to_2d(line);
        for (KSR::size_t other_intersection_line_idx : support_plane(support_plane_idx).intersection_lines_idx())
        {
          Point_2 point;
          if (!KSR::intersection_2 (line_2,
                                    support_plane(support_plane_idx).to_2d
                                    (intersection_line(other_intersection_line_idx).line()),
                                    point))
            continue;

          Point_3 point_3 = support_plane(support_plane_idx).to_3d(point);

          if (!point_is_inside_bbox_section_of_intersection_line (point_3, intersection_line_idx)
              || !point_is_inside_bbox_section_of_intersection_line (point_3, other_intersection_line_idx))
            continue;

          KSR::size_t meta_vertex_idx = add_meta_vertex (point_3, support_plane_idx);
          intersection_line(intersection_line_idx).meta_vertices_idx().push_back (meta_vertex_idx);
          intersection_line(other_intersection_line_idx).meta_vertices_idx().push_back (meta_vertex_idx);
        }
      }
    }
    
    for (KSR::size_t support_plane_idx : { support_plane_idx_0, support_plane_idx_1 })
    {
      intersection_line(intersection_line_idx).support_planes_idx().push_back (support_plane_idx);
      support_plane(support_plane_idx).intersection_lines_idx().push_back (intersection_line_idx);
    }
    return intersection_line_idx;
  }

  // Add segment on full intersection line, using 2 extrem meta vertices
  KSR::size_t add_segment (KSR::size_t intersection_line_idx, KSR::size_t source_idx, KSR::size_t target_idx)
  {
    return add_segment (Segment (intersection_line_idx, source_idx, target_idx));
  }
  
  void partition (KSR::size_t polygon_idx, const Line_2& line,
                  KSR::Idx_vector& positive_side, KSR::Idx_vector& negative_side) const
  {
    std::vector<bool> has_on_positive_side;
    has_on_positive_side.reserve(polygon(polygon_idx).vertices_idx().size());
    
    for (KSR::size_t vertex_idx : polygon(polygon_idx).vertices_idx())
      has_on_positive_side.push_back (line.has_on_positive_side(vertex(vertex_idx).point(m_current_time)));

    KSR::size_t first_positive = KSR::no_element();
        
    for (std::size_t i = 0; i <= has_on_positive_side.size(); ++ i)
      if (!has_on_positive_side[i] && has_on_positive_side[(i+1) % has_on_positive_side.size()])
      {
        first_positive = KSR::size_t((i+1) % has_on_positive_side.size());
        break;
      }

    KSR::size_t current_position = first_positive;
    do
    {
      if (has_on_positive_side[std::size_t(current_position)])
        positive_side.push_back (polygon(polygon_idx).vertices_idx()[current_position]);
      else
        negative_side.push_back (polygon(polygon_idx).vertices_idx()[current_position]);
      
      current_position = (current_position + 1) % has_on_positive_side.size();
    }
    while (current_position != first_positive);

    CGAL_assertion (!positive_side.empty() && !negative_side.empty());
    CGAL_assertion (positive_side.size() + negative_side.size() >= 3);
  }

  bool do_intersect (KSR::size_t polygon_idx, const Line_2& line) const
  {
    bool positive_side = false, negative_side = false;
    for (KSR::size_t vertex_idx : polygon(polygon_idx).vertices_idx())
    {
      if (line.has_on_positive_side(vertex(vertex_idx).point(m_current_time)))
        positive_side = true;
      else
        negative_side = true;
      if (positive_side && negative_side)
        return true;
    }
    
    return false;
  }

  void cut_polygon (KSR::size_t polygon_idx, KSR::size_t intersection_line_idx)
  {
    CGAL_KSR_CERR(3) << "** Cutting " << polygon_str(polygon_idx) << std::endl;

    Line_2 line_2 = line_on_support_plane (intersection_line_idx, polygon(polygon_idx).support_plane_idx());
    
    KSR::Idx_vector positive_side, negative_side;
    partition (polygon_idx, line_2, positive_side, negative_side);

    Line_2 line_0 (vertex(positive_side.back()).point(m_current_time),
                   vertex(negative_side.front()).point(m_current_time));
    Line_2 line_1 (vertex(negative_side.back()).point(m_current_time),
                   vertex(positive_side.front()).point(m_current_time));
    
    Point_2 inter_0 = KSR::intersection_2<Point_2> (line_0, line_2);
    Point_2 inter_1 = KSR::intersection_2<Point_2> (line_1, line_2);

    // Compute speeds
    Line_2 future_line_0 (vertex(positive_side.back()).point(m_current_time + 1),
                          vertex(negative_side.front()).point(m_current_time + 1));
    Line_2 future_line_1 (vertex(negative_side.back()).point(m_current_time + 1),
                          vertex(positive_side.front()).point(m_current_time + 1));

    Point_2 future_inter_0 = KSR::intersection_2<Point_2> (future_line_0, line_2);
    Point_2 future_inter_1 = KSR::intersection_2<Point_2> (future_line_1, line_2);
    
    Vector_2 direction_0 (inter_0, future_inter_0);
    Vector_2 direction_1 (inter_1, future_inter_1);

    KSR::size_t new_polygon_idx = add_polygon (Polygon(polygon(polygon_idx).input_idx(), polygon(polygon_idx).support_plane_idx()));
    support_plane(polygon(polygon_idx).support_plane_idx()).polygons_idx().push_back (new_polygon_idx);

    for (KSR::size_t vertex_idx : positive_side)
      vertex(vertex_idx).polygon_idx() = polygon_idx;
    for (KSR::size_t vertex_idx : negative_side)
      vertex(vertex_idx).polygon_idx() = new_polygon_idx;

    positive_side.push_back (add_vertex (Vertex (inter_0, polygon_idx)));
    m_vertices.back().intersection_line_idx() = intersection_line_idx;
    m_vertices.back().direction() = direction_0;
    
    positive_side.push_back (add_vertex (Vertex (inter_1, polygon_idx)));
    m_vertices.back().intersection_line_idx() = intersection_line_idx;
    m_vertices.back().direction() = direction_1;

    add_segment (intersection_line_idx, number_of_vertices() - 1, number_of_vertices() - 2);

    negative_side.push_back (add_vertex (Vertex (inter_1, new_polygon_idx)));
    m_vertices.back().intersection_line_idx() = intersection_line_idx;
    m_vertices.back().direction() = direction_1;
    
    negative_side.push_back (add_vertex (Vertex (inter_0, new_polygon_idx)));
    m_vertices.back().intersection_line_idx() = intersection_line_idx;
    m_vertices.back().direction() = direction_0;

    if (direction_0 == CGAL::NULL_VECTOR)
    {
      KSR::size_t meta_vertex_idx = add_meta_vertex (support_plane_of_polygon(polygon_idx).to_3d (inter_0),
                                                     polygon(polygon_idx).support_plane_idx());
      attach_vertex_to_meta_vertex (m_vertices.size() - 4, meta_vertex_idx);
      attach_vertex_to_meta_vertex (m_vertices.size() - 1, meta_vertex_idx);
    }
    
    if (direction_1 == CGAL::NULL_VECTOR)
    {
      KSR::size_t meta_vertex_idx = add_meta_vertex (support_plane_of_polygon(polygon_idx).to_3d (inter_1),
                                                     polygon(polygon_idx).support_plane_idx());
      attach_vertex_to_meta_vertex (m_vertices.size() - 3, meta_vertex_idx);
      attach_vertex_to_meta_vertex (m_vertices.size() - 2, meta_vertex_idx);
    }
    
    add_segment (intersection_line_idx, number_of_vertices() - 1, number_of_vertices() - 2);
    
    polygon(polygon_idx).vertices_idx().swap (positive_side);
    polygon(new_polygon_idx).vertices_idx().swap (negative_side);
    
    CGAL_KSR_CERR(3) << "*** new polygons:";
    for (KSR::size_t i : { polygon_idx, new_polygon_idx })
      CGAL_KSR_CERR(3) << " " << polygon_str(i);
    CGAL_KSR_CERR(3) << std::endl;
  }

  void update_positions (FT time)
  {
    m_current_time = time;
  }
  
};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_DATA_STRUCTURE_H
